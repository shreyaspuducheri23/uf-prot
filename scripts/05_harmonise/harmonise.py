#!/usr/bin/env python3
"""
05_harmonise/harmonise.py
Join exposure instruments with Kim outcome, search for proxies when SNPs are absent,
then harmonise alleles via TwoSampleMR (via rlib/harmonise.R subprocess).

Usage:
  python scripts/05_harmonise/harmonise.py --cohort ARIC_EA [--limit N]
"""
import argparse
import logging
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from scripts.lib.checkpoint import Checkpoint, output_exists
from scripts.lib.config import add_config_arg, load_config, get_section
from scripts.lib.filters import drop_ambig_palindromes
from scripts.lib.liftover import lift_table
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.outcome import OutcomeLookup, normalize_outcome_row, KIM_N
from scripts.lib.paths import (
    COHORTS, ROOT, cohort_dir, instruments_hg38_dir, harmonised_dir
)
from scripts.lib.plink import find_proxies
from scripts.lib.progress import bar
from scripts.lib.sumstats_io import read_norm, write_norm

log = setup_logger("05_harmonise")

MAF_PROXY_MAX = 0.42


def harmonise_cohort(cohort: str, limit: int | None = None) -> int:
    in_dir = instruments_hg38_dir(cohort)
    out_dir = harmonised_dir(cohort)
    out_dir.mkdir(parents=True, exist_ok=True)

    tsv_files = sorted(in_dir.glob("*.tsv"))
    if limit:
        tsv_files = tsv_files[:limit]

    cp = Checkpoint(cohort_dir(cohort) / "_state_05.json")
    todo = [f for f in tsv_files if not cp.is_done(f.stem)]

    log.info(f"{cohort}: {len(tsv_files)} instrument files, {len(todo)} to harmonise")

    n_ok = 0
    n_proxies_total = 0

    with OutcomeLookup() as outcome:
        for tsv_path in bar(todo, desc=f"{cohort} harmonise"):
            seqid = tsv_path.stem
            out_path = out_dir / f"{seqid}.tsv"

            if output_exists(out_path):
                cp.mark_done(seqid)
                n_ok += 1
                continue

            df = read_norm(tsv_path)
            if df.empty:
                cp.mark_done(seqid)
                continue

            df, n_proxies = _join_outcome(df, outcome)
            n_proxies_total += n_proxies

            if df.empty:
                log.debug(f"{seqid}: 0 SNPs matched in outcome")
                cp.mark_done(seqid)
                continue

            # Call rlib/harmonise.R for canonical allele alignment
            harmonised = _call_harmonise_r(df, seqid)
            if harmonised is None or harmonised.empty:
                cp.mark_done(seqid)
                continue

            write_norm(harmonised, out_path)
            cp.mark_done(seqid)
            n_ok += 1

    log.info(f"{cohort}: {n_ok} proteins harmonised. "
             f"Proxy SNPs used: {n_proxies_total}")
    return n_ok


def _join_outcome(df: pd.DataFrame, outcome: OutcomeLookup) -> tuple[pd.DataFrame, int]:
    """
    For each instrument, look up outcome row by (chrom_hg38, pos_hg38).
    For missing SNPs, attempt proxy search in hg19 coords, then lift proxy to hg38.
    Returns enriched DataFrame with outcome columns appended, and count of proxies used.
    """
    positions = list(zip(df["chrom_hg38"].astype(str), df["pos_hg38"].astype(int)))
    out_df = outcome.fetch_snps(positions)

    # Build lookup: (chrom, pos) → outcome row
    out_lookup: dict[tuple[str, int], dict] = {}
    n_dup_pos = 0
    for _, row in out_df.iterrows():
        key = (str(row["chromosome"]), int(row["base_pair_location"]))
        if key in out_lookup:
            n_dup_pos += 1
        out_lookup[key] = normalize_outcome_row(row)
    if n_dup_pos:
        log.warning(
            f"_join_outcome: {n_dup_pos} duplicate (chrom, pos) positions in outcome — "
            f"last row kept for each duplicate"
        )

    rows = []
    n_proxies = 0
    n_no_rsid = 0
    missing_rsids: list[str] = []
    missing_idx: list[int] = []

    for i, (_, instr) in enumerate(df.iterrows()):
        key = (str(instr["chrom_hg38"]), int(instr["pos_hg38"]))
        if key in out_lookup:
            rows.append({**instr.to_dict(), **out_lookup[key]})
        else:
            rsid = str(instr.get("rsid", "."))
            if rsid != "." and rsid:
                missing_rsids.append(rsid)
                missing_idx.append(i)
            else:
                n_no_rsid += 1

    if n_no_rsid:
        log.warning(
            f"_join_outcome: {n_no_rsid} instruments have no rsid and no position match — "
            f"cannot search for proxies; these are dropped"
        )

    # Proxy search for missing SNPs
    if missing_rsids:
        proxy_map = find_proxies(missing_rsids)  # {target: (proxy_rsid, r2)}
        if proxy_map:
            # Fetch proxies from outcome by rsID
            proxy_rsids = [v[0] for v in proxy_map.values()]
            proxy_out = outcome.fetch_by_rsid(proxy_rsids)
            proxy_lookup: dict[str, dict] = {}
            for _, row in proxy_out.iterrows():
                proxy_lookup[row["rsid"]] = normalize_outcome_row(row)

            for j, i in enumerate(missing_idx):
                target = missing_rsids[j]
                if target not in proxy_map:
                    continue
                proxy_rsid, proxy_r2 = proxy_map[target]
                if proxy_rsid not in proxy_lookup:
                    continue

                instr = df.iloc[i]
                # Check proxy MAF
                proxy_row = proxy_lookup[proxy_rsid]
                eaf = proxy_row.get("EAF_out")
                if eaf is not None:
                    maf = min(eaf, 1 - eaf)
                    if maf > MAF_PROXY_MAX:
                        continue

                merged = {
                    **instr.to_dict(),
                    **proxy_row,
                    "proxy_rsid": proxy_rsid,
                    "proxy_r2": proxy_r2,
                }
                # Update hg38 coords from proxy
                rows.append(merged)
                n_proxies += 1

    if not rows:
        return pd.DataFrame(), n_proxies

    return pd.DataFrame(rows), n_proxies


def _call_harmonise_r(df: pd.DataFrame, seqid: str) -> pd.DataFrame | None:
    """
    Write exposure + outcome portions to temp TSVs and call rlib/harmonise.R.
    Returns harmonised DataFrame or None on failure.
    """
    exp_cols = ["rsid", "chrom", "pos", "EA", "OA", "EAF", "beta", "se", "pval", "N", "seqid", "gene", "uniprot"]
    out_cols = ["rsid", "chrom_hg38", "pos_hg38", "EA_out", "OA_out", "EAF_out",
                "beta_out", "se_out", "pval_out", "N_out"]

    # Rename for harmonise.R expectations
    exp_df = df[[c for c in exp_cols if c in df.columns]].copy()
    out_df = df[[c for c in out_cols if c in df.columns]].copy().rename(columns={
        "rsid": "rsid",
        "chrom_hg38": "chromosome",
        "pos_hg38": "base_pair_location",
        "EA_out": "EA", "OA_out": "OA",
        "EAF_out": "EAF", "beta_out": "beta",
        "se_out": "se", "pval_out": "pval", "N_out": "N",
    })
    out_df["chrom"] = out_df["chromosome"]
    out_df["pos"] = out_df["base_pair_location"]

    with tempfile.TemporaryDirectory(prefix=f"harm_{seqid}_") as tmp:
        tmp = Path(tmp)
        exp_path = tmp / "exposure.tsv"
        out_path = tmp / "outcome.tsv"
        result_path = tmp / "harmonised.tsv"

        exp_df.to_csv(exp_path, sep="\t", index=False)
        out_df.to_csv(out_path, sep="\t", index=False)

        cmd = [
            "Rscript", str(ROOT / "scripts" / "rlib" / "harmonise.R"),
            "--exp", str(exp_path),
            "--out", str(out_path),
            "--result", str(result_path),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            log.warning(f"{seqid}: harmonise.R failed — {result.stderr[:300]}")
            return None

        if not result_path.exists():
            return None

        return pd.read_csv(result_path, sep="\t")


def main() -> None:
    parser = argparse.ArgumentParser(description="Harmonise instruments with Kim outcome")
    parser.add_argument("--cohort", choices=COHORTS + ["all"], default="all")
    parser.add_argument("--limit", type=int, default=None)
    add_config_arg(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    harm_cfg = get_section(cfg, "harmonise")
    global MAF_PROXY_MAX
    MAF_PROXY_MAX = harm_cfg["maf_proxy_max"]

    cohorts = COHORTS if args.cohort == "all" else [args.cohort]

    with RunManifest("05_harmonise/harmonise.py", args=str(args)) as manifest:
        total = sum(harmonise_cohort(c, limit=args.limit) for c in cohorts)
        manifest.n_units = total


if __name__ == "__main__":
    main()
