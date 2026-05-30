#!/usr/bin/env python3
"""
04_liftover/instruments_to_hg38.py
Lift instrument SNP positions to hg38 for downstream harmonisation.
ARIC_EA and deCODE are already hg38 and pass through unchanged.

Usage:
  python scripts/04_liftover/instruments_to_hg38.py [--cohort ARIC_EA] [--limit N]
"""
import argparse
import logging

import pandas as pd

from scripts.lib.checkpoint import Checkpoint, output_exists
from scripts.lib.config import add_config_arg, load_config
from scripts.lib.liftover import lift_table
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import (
    COHORTS,
    cohort_dir,
    filtered_cis_pqtls_dir,
    filtered_cis_pqtls_hg38_dir,
    instruments_dir,
    instruments_hg38_dir,
    raw_cis_sumstats_dir,
    raw_cis_sumstats_hg38_dir,
)
from scripts.lib.progress import bar
from scripts.lib.sumstats_io import read_norm, write_norm

log = setup_logger("04_liftover")

# Cohorts already in hg38 — pass through with chrom_hg38/pos_hg38 = chrom/pos
# ARIC .glm.linear and deCODE positions are already hg38; do not lift.
HG38_COHORTS = {"deCODE", "ARIC_EA"}

# Cohorts whose cis summary-statistic positions are already in hg38.
CIS_HG38_COHORTS = {"deCODE", "ARIC_EA"}
LIFTED_REQUIRED_COLS = ["seqid", "chrom", "pos", "chrom_hg38", "pos_hg38"]


def _seqid_from_sumstats_path(path) -> str:
    name = path.name
    for suffix in (".tsv.gz", ".tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def lift_cohort(cohort: str, limit: int | None = None) -> int:
    in_dir = instruments_dir(cohort)
    out_dir = instruments_hg38_dir(cohort)
    out_dir.mkdir(parents=True, exist_ok=True)

    tsv_files = sorted(in_dir.glob("*.tsv"))
    if limit:
        tsv_files = tsv_files[:limit]

    cp = Checkpoint(cohort_dir(cohort) / "_state_04.json")
    todo = [f for f in tsv_files if not cp.is_done(f.stem)]

    log.info(f"{cohort}: {len(tsv_files)} instrument files, {len(todo)} to process")
    n_ok = 0
    total_dropped = 0

    for tsv_path in bar(todo, desc=f"{cohort} liftover"):
        seqid = tsv_path.stem
        out_path = out_dir / f"{seqid}.tsv"

        if output_exists(out_path, required_cols=LIFTED_REQUIRED_COLS, min_rows=1):
            cp.mark_done(seqid)
            n_ok += 1
            continue

        df = read_norm(tsv_path)
        if df.empty:
            cp.mark_done(seqid)
            continue

        n_in = len(df)

        if cohort in HG38_COHORTS:
            df = df.copy()
            df["chrom_hg38"] = df["chrom"]
            df["pos_hg38"] = df["pos"]
            df["chrom_hg19"] = None
            df["pos_hg19"] = None
        else:
            df = df.rename(columns={"chrom": "chrom_hg19", "pos": "pos_hg19"})
            # Also add chrom/pos as hg19 for LD-ref usage (plink uses hg19)
            df = lift_table(df, chrom_col="chrom_hg19", pos_col="pos_hg19")
            df["chrom"] = df["chrom_hg19"]
            df["pos"] = df["pos_hg19"]

        n_out = len(df)
        dropped = n_in - n_out
        if dropped and n_in > 0:
            pct = 100 * dropped / n_in
            log.debug(f"{seqid}: {dropped}/{n_in} ({pct:.1f}%) SNPs dropped in liftover")
            total_dropped += dropped

        write_norm(df, out_path)
        cp.mark_done(seqid)
        n_ok += 1

    log.info(f"{cohort}: done. {n_ok} proteins processed. Total SNPs dropped: {total_dropped}")
    return n_ok


def lift_sumstats_file(cohort: str, in_path, out_path) -> int:
    """Lift or pass through one cis summary-statistics file."""
    df = read_norm(in_path)
    if df.empty:
        return 0

    n_in = len(df)
    if cohort not in CIS_HG38_COHORTS:
        df = lift_table(df, chrom_col="chrom", pos_col="pos")
        # lift_table adds chrom_hg38/pos_hg38; overwrite chrom/pos in place.
        df["chrom"] = df.pop("chrom_hg38")
        df["pos"] = df.pop("pos_hg38").astype(int)

    write_norm(df, out_path)
    return n_in - len(df)


def _lift_cis_sumstats_family(
    cohort: str,
    in_dir,
    out_dir,
    state_name: str,
    label: str,
    limit: int | None = None,
) -> int:
    """Lift a family of cis summary-statistics files to hg38."""
    out_dir.mkdir(parents=True, exist_ok=True)

    tsv_files = sorted([*in_dir.glob("*.tsv"), *in_dir.glob("*.tsv.gz")])
    if limit:
        tsv_files = tsv_files[:limit]

    cp = Checkpoint(cohort_dir(cohort) / state_name)
    todo = [f for f in tsv_files if not cp.is_done(_seqid_from_sumstats_path(f))]

    log.info(f"{cohort}: {len(tsv_files)} {label} files, {len(todo)} to lift")
    n_ok = 0
    total_dropped = 0

    for tsv_path in bar(todo, desc=f"{cohort} {label} liftover"):
        seqid = _seqid_from_sumstats_path(tsv_path)
        suffix = ".tsv.gz" if tsv_path.name.endswith(".tsv.gz") else ".tsv"
        out_path = out_dir / f"{seqid}{suffix}"

        if output_exists(out_path, required_cols=["chrom", "pos"], min_rows=1):
            cp.mark_done(seqid)
            n_ok += 1
            continue

        try:
            dropped = lift_sumstats_file(cohort, tsv_path, out_path)
        except pd.errors.EmptyDataError:
            cp.mark_done(seqid)
            continue

        if dropped:
            log.debug(f"{seqid}: {dropped} rows dropped in {label} liftover")
            total_dropped += dropped

        cp.mark_done(seqid)
        n_ok += 1

    log.info(f"{cohort}: {label} liftover done. {n_ok} files. Total rows dropped: {total_dropped}")
    return n_ok


def lift_filtered_cis_pqtls_cohort(cohort: str, limit: int | None = None) -> int:
    return _lift_cis_sumstats_family(
        cohort,
        filtered_cis_pqtls_dir(cohort),
        filtered_cis_pqtls_hg38_dir(cohort),
        "_state_04_filtered_cis.json",
        "filtered_cis_pqtls",
        limit=limit,
    )


def lift_raw_cis_sumstats_cohort(cohort: str, limit: int | None = None) -> int:
    return _lift_cis_sumstats_family(
        cohort,
        raw_cis_sumstats_dir(cohort),
        raw_cis_sumstats_hg38_dir(cohort),
        "_state_04_raw_cis.json",
        "raw_cis_sumstats",
        limit=limit,
    )


def lift_cis_sumstats_cohort(cohort: str, limit: int | None = None) -> int:
    """Compatibility wrapper for the filtered MR-ready cis-pQTL product."""
    return lift_filtered_cis_pqtls_cohort(cohort, limit=limit)


def main() -> None:
    parser = argparse.ArgumentParser(description="Liftover instrument positions to GRCh38")
    parser.add_argument("--cohort", choices=COHORTS + ["all"], default="all")
    parser.add_argument("--limit", type=int, default=None)
    add_config_arg(parser)
    args = parser.parse_args()
    load_config(args.config)  # validate config exists; values not yet used here

    cohorts = COHORTS if args.cohort == "all" else [args.cohort]

    with RunManifest("04_liftover/instruments_to_hg38.py", args=str(args)) as manifest:
        total  = sum(lift_cohort(c, limit=args.limit) for c in cohorts)
        total += sum(lift_filtered_cis_pqtls_cohort(c, limit=args.limit) for c in cohorts)
        total += sum(lift_raw_cis_sumstats_cohort(c, limit=args.limit) for c in cohorts)
        manifest.n_units = total


if __name__ == "__main__":
    main()
