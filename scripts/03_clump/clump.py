#!/usr/bin/env python3
"""
03_clump/clump.py
LD clumping for all cohorts: 1Mb window, r2<0.001, p<5e-8 vs 1000G EUR LD ref.
Computes F-statistics; outputs per-protein instrument TSVs and a summary.

Usage:
  python scripts/03_clump/clump.py --cohort ARIC_EA [--limit N]
  python scripts/03_clump/clump.py --cohort deCODE
  python scripts/03_clump/clump.py --cohort UKB_PPP
  python scripts/03_clump/clump.py --cohort Fenland
"""
import argparse
import logging

import pandas as pd

from scripts.lib.checkpoint import Checkpoint, output_exists
from scripts.lib.config import add_config_arg, load_config, get_section
from scripts.lib.fstat import add_fstat, WEAK_INSTRUMENT_THRESHOLD
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import (
    COHORTS, cis_sumstats_dir, instruments_dir, cohort_dir
)
from scripts.lib.plink import clump
from scripts.lib.progress import bar
from scripts.lib.sumstats_io import read_norm, write_norm

log = setup_logger("03_clump")

INSTRUMENT_REQUIRED_COLS = ["seqid", "chrom", "pos", "rsid", "pval", "F_stat"]


def clump_cohort(
    cohort: str,
    limit: int | None = None,
    window_kb: int = 1000,
    r2: float = 0.001,
    p1: float = 5e-8,
) -> int:
    in_dir = cis_sumstats_dir(cohort)
    out_dir = instruments_dir(cohort)
    out_dir.mkdir(parents=True, exist_ok=True)

    tsv_files = sorted(in_dir.glob("*.tsv"))
    if limit:
        tsv_files = tsv_files[:limit]

    cp = Checkpoint(cohort_dir(cohort) / "_state_03.json")
    todo = [f for f in tsv_files if not cp.is_done(f.stem)]

    log.info(f"{cohort}: {len(tsv_files)} proteins, {len(todo)} to clump")

    summary_rows = []
    n_ok = 0

    for tsv_path in bar(todo, desc=f"{cohort} clump"):
        seqid = tsv_path.stem
        out_path = out_dir / f"{seqid}.tsv"

        if output_exists(out_path, required_cols=INSTRUMENT_REQUIRED_COLS, min_rows=1):
            cp.mark_done(seqid)
            n_ok += 1
            continue

        try:
            df = read_norm(tsv_path)
        except Exception as exc:
            log.warning(f"{seqid}: read error — {exc}")
            continue

        if df.empty:
            cp.mark_done(seqid)
            continue

        try:
            clumped = clump(df, seqid, window_kb=window_kb, r2=r2, p1=p1)
        except Exception as exc:
            log.warning(f"{seqid}: clumping failed — {exc}")
            continue

        if clumped.empty:
            log.debug(f"{seqid}: 0 instruments after clumping")
            cp.mark_done(seqid)
            continue

        clumped = add_fstat(clumped)
        n_weak = (clumped["F_stat"] < WEAK_INSTRUMENT_THRESHOLD).sum()
        if n_weak:
            log.debug(f"{seqid}: {n_weak} instruments with F<10 (flagged, kept)")

        write_norm(clumped, out_path)
        cp.mark_done(seqid)
        n_ok += 1

        summary_rows.append({
            "seqid": seqid,
            "gene": clumped["gene"].iloc[0] if "gene" in clumped.columns else "",
            "n_instruments": len(clumped),
            "min_pval": clumped["pval"].min(),
            "max_F_stat": clumped["F_stat"].max(),
            "any_weak_F": n_weak > 0,
        })

    # Write summary
    if summary_rows:
        summary_path = cohort_dir(cohort) / "instruments_summary.tsv"
        existing = pd.read_csv(summary_path, sep="\t") if summary_path.exists() else pd.DataFrame()
        updated = pd.concat([existing, pd.DataFrame(summary_rows)], ignore_index=True)
        updated.drop_duplicates("seqid").to_csv(summary_path, sep="\t", index=False)
        log.info(f"Updated instruments_summary.tsv ({len(updated)} proteins)")

    log.info(f"{cohort}: {n_ok} proteins with instruments")
    return n_ok


def main() -> None:
    parser = argparse.ArgumentParser(description="LD clumping for cis-pQTL instruments")
    parser.add_argument("--cohort", choices=COHORTS + ["all"], default="all")
    parser.add_argument("--limit", type=int, default=None)
    add_config_arg(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    clump_cfg = get_section(cfg, "clump")
    fstat_cfg = get_section(cfg, "fstat")
    log.info(f"Clump config: window_kb={clump_cfg['window_kb']}, r2={clump_cfg['r2']}, "
             f"p1={clump_cfg['p1']}, F_threshold={fstat_cfg['weak_threshold']}")

    cohorts = COHORTS if args.cohort == "all" else [args.cohort]

    with RunManifest("03_clump/clump.py", args=str(args)) as manifest:
        total = 0
        for cohort in cohorts:
            total += clump_cohort(
                cohort,
                limit=args.limit,
                window_kb=clump_cfg["window_kb"],
                r2=clump_cfg["r2"],
                p1=clump_cfg["p1"],
            )
        manifest.n_units = total


if __name__ == "__main__":
    main()
