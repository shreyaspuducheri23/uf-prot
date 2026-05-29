#!/usr/bin/env python3
"""
04_liftover/instruments_to_hg38.py
Lift instrument SNP positions to hg38 for downstream harmonisation.
ARIC_EA, UKB_PPP, and deCODE are already hg38 and pass through unchanged.

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
from scripts.lib.paths import COHORTS, cohort_dir, cis_sumstats_dir, cis_sumstats_hg38_dir, instruments_dir, instruments_hg38_dir
from scripts.lib.progress import bar
from scripts.lib.sumstats_io import read_norm, write_norm

log = setup_logger("04_liftover")

# Cohorts already in hg38 — pass through with chrom_hg38/pos_hg38 = chrom/pos
# ARIC .glm.linear and UKB-PPP VCF-format IDs both encode hg38 positions; do not lift.
HG38_COHORTS = {"deCODE", "ARIC_EA", "UKB_PPP"}

# Cohorts whose cis_sumstats positions are already in hg38.
CIS_HG38_COHORTS = {"deCODE", "ARIC_EA", "UKB_PPP"}
LIFTED_REQUIRED_COLS = ["seqid", "chrom", "pos", "chrom_hg38", "pos_hg38"]


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


def lift_cis_sumstats_cohort(cohort: str, limit: int | None = None) -> int:
    """Lift all cis_sumstats files to hg38, writing to cis_sumstats_hg38/.

    hg19 cohorts: positions are lifted via lift_table().
    hg38 cohorts (deCODE): files are passed through unchanged.
    Uses a separate checkpoint (_state_04_cis.json) so the instrument
    checkpoint (_state_04.json) is not disturbed.
    """
    in_dir  = cis_sumstats_dir(cohort)
    out_dir = cis_sumstats_hg38_dir(cohort)
    out_dir.mkdir(parents=True, exist_ok=True)

    tsv_files = sorted(in_dir.glob("*.tsv"))
    if limit:
        tsv_files = tsv_files[:limit]

    cp = Checkpoint(cohort_dir(cohort) / "_state_04_cis.json")
    todo = [f for f in tsv_files if not cp.is_done(f.stem)]

    log.info(f"{cohort}: {len(tsv_files)} cis_sumstats files, {len(todo)} to lift")
    n_ok = 0
    total_dropped = 0

    for tsv_path in bar(todo, desc=f"{cohort} cis_sumstats liftover"):
        seqid    = tsv_path.stem
        out_path = out_dir / f"{seqid}.tsv"

        if output_exists(out_path, required_cols=["chrom", "pos"], min_rows=1):
            cp.mark_done(seqid)
            n_ok += 1
            continue

        df = read_norm(tsv_path)
        if df.empty:
            cp.mark_done(seqid)
            continue

        n_in = len(df)

        if cohort in CIS_HG38_COHORTS:
            pass  # positions already hg38 — write as-is
        else:
            df = lift_table(df, chrom_col="chrom", pos_col="pos")
            # lift_table adds chrom_hg38/pos_hg38; overwrite chrom/pos in place
            df["chrom"] = df.pop("chrom_hg38")
            df["pos"]   = df.pop("pos_hg38").astype(int)

        dropped = n_in - len(df)
        if dropped:
            log.debug(f"{seqid}: {dropped}/{n_in} rows dropped in cis_sumstats liftover")
            total_dropped += dropped

        write_norm(df, out_path)
        cp.mark_done(seqid)
        n_ok += 1

    log.info(f"{cohort}: cis_sumstats liftover done. {n_ok} files. Total rows dropped: {total_dropped}")
    return n_ok


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
        total += sum(lift_cis_sumstats_cohort(c, limit=args.limit) for c in cohorts)
        manifest.n_units = total


if __name__ == "__main__":
    main()
