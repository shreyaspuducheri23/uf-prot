#!/usr/bin/env python3
"""
08_coloc/sharepro.py
Run SharePro colocalization (primary) on candidate proteins.

Builds per-region z-score inputs and LD matrix (plink2 --r square),
then subprocesses into tools/SharePro_coloc/src/sharepro_loc.py.
Requires ±1 Mb regions from 08_coloc/extract_regions.py.

Usage:
  python scripts/08_coloc/sharepro.py [--cohort ARIC_EA] [--limit N]
"""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.lib.checkpoint import Checkpoint
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import COHORTS, COLOC_REGIONS_DIR, SHAREPRO_SCRIPT, cohort_dir
from scripts.lib.plink import r_square_matrix
from scripts.lib.progress import bar

log = setup_logger("08_sharepro")

SHAREPRO_OUT = COLOC_REGIONS_DIR.parent / "sharepro_results.tsv"


def build_zscore_input(df: pd.DataFrame, n: int,
                        beta_col: str = "beta", se_col: str = "se",
                        snp_col: str = "rsid") -> pd.DataFrame:
    """Compute z-scores and return the input table for SharePro."""
    df = df.copy()
    df["z"] = df[beta_col].astype(float) / df[se_col].astype(float)
    df["snp"] = df[snp_col].astype(str)
    df["N"] = n
    return df[["snp", "z", "N"]].dropna()


def _infer_n_exp(exp_df: pd.DataFrame) -> int | None:
    if "N" not in exp_df.columns:
        return None
    n_vals = pd.to_numeric(exp_df["N"], errors="coerce")
    n_vals = n_vals[np.isfinite(n_vals)]
    if n_vals.empty:
        return None
    n_exp = int(np.median(n_vals))
    if n_exp <= 0:
        return None
    return n_exp


def run_sharepro(region_dir: Path, seqid: str, N_out: int) -> tuple[dict | None, str | None]:
    """
    Run SharePro for one protein region.
    Returns (result, None) on success; (None, reason) on failure/skip.
    """
    exp_path = region_dir / "exposure.tsv"
    out_path = region_dir / "outcome.tsv"
    if not exp_path.exists() or not out_path.exists():
        return None, "missing_exposure_or_outcome_file"

    exp_df = pd.read_csv(exp_path, sep="\t", dtype={"chrom": str, "rsid": str})
    out_df = pd.read_csv(out_path, sep="\t", dtype={"chromosome": str, "rsid": str})
    n_exp = _infer_n_exp(exp_df)
    if n_exp is None:
        return None, "invalid_or_missing_N_exp"

    # Merge on rsid to get common SNPs
    out_df = out_df.rename(columns={
        "beta": "beta", "standard_error": "se",
        "effect_allele_frequency": "EAF",
    })
    common_snps = set(exp_df["rsid"]) & set(out_df["rsid"])
    common_snps = {s for s in common_snps if s != "." and pd.notna(s)}

    if len(common_snps) < 5:
        log.debug(f"{seqid}: <5 common SNPs ({len(common_snps)}) — skipping SharePro")
        return None, "insufficient_common_snps"

    exp_sub = exp_df[exp_df["rsid"].isin(common_snps)].drop_duplicates("rsid")
    out_sub = out_df[out_df["rsid"].isin(common_snps)].drop_duplicates("rsid")

    snp_order = sorted(common_snps)
    exp_z = build_zscore_input(exp_sub.set_index("rsid").loc[snp_order].reset_index(),
                                n_exp, snp_col="rsid")
    out_z = build_zscore_input(out_sub.set_index("rsid").loc[snp_order].reset_index(),
                                N_out, snp_col="rsid")
    if exp_z.empty or out_z.empty:
        return None, "empty_sharepro_z_inputs"

    # LD matrix
    try:
        ld_mat = r_square_matrix(snp_order)
    except Exception as exc:
        log.warning(f"{seqid}: LD matrix failed — {exc}")
        return None, "ld_matrix_failed"

    with tempfile.TemporaryDirectory(prefix=f"sharepro_{seqid}_") as tmp:
        tmp = Path(tmp)

        exp_z_path = tmp / "exp_z.txt"
        out_z_path = tmp / "out_z.txt"
        ld_path = tmp / "ld.txt"
        result_path = tmp / "sharepro_result.json"

        exp_z[["snp", "z", "N"]].to_csv(exp_z_path, sep="\t", index=False)
        out_z[["snp", "z", "N"]].to_csv(out_z_path, sep="\t", index=False)
        ld_mat.to_csv(ld_path, sep="\t", header=False, index=False)

        cmd = [
            "python", str(SHAREPRO_SCRIPT),
            "--z1", str(exp_z_path),
            "--z2", str(out_z_path),
            "--ld", str(ld_path),
            "--out", str(result_path),
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            log.warning(f"{seqid}: SharePro failed — {res.stderr[:300]}")
            return None, "sharepro_subprocess_failed"

        if not result_path.exists():
            log.warning(f"{seqid}: SharePro produced no output")
            return None, "sharepro_output_missing"

        with open(result_path) as fh:
            raw = json.load(fh)

        # Extract PP.H4 — SharePro's main colocalization posterior
        pp_h4 = raw.get("PP.H4", raw.get("pp_h4", None))
        return {
            "seqid": seqid,
            "n_snps": len(snp_order),
            "PP_H4": pp_h4,
            "coloc_positive": pp_h4 is not None and pp_h4 >= 0.8,
            "raw": raw,
        }, None


def run_cohort_sharepro(
    cohort: str,
    limit: int | None = None,
    retry_failed: bool = False,
) -> list[dict]:
    region_base = COLOC_REGIONS_DIR / cohort
    if not region_base.exists():
        log.warning(f"{cohort}: no coloc regions directory")
        return []

    candidates = sorted(p for p in region_base.iterdir() if p.is_dir())
    if limit:
        candidates = candidates[:limit]

    cp = Checkpoint(cohort_dir(cohort) / "_state_08_sharepro.json")
    todo = cp.remaining(candidates, key=lambda p: p.name, include_failed=retry_failed)

    log.info(f"{cohort}: {len(candidates)} candidates, {len(todo)} remaining")

    from scripts.lib.outcome import KIM_N

    results = []
    for region_dir in bar(todo, desc=f"{cohort} SharePro"):
        seqid = region_dir.name
        try:
            result, failure_reason = run_sharepro(region_dir, seqid, N_out=KIM_N)
            if result is not None:
                result["cohort"] = cohort
                results.append(result)
                cp.mark_done(seqid)
            else:
                cp.mark_failed(seqid, failure_reason or "sharepro_failed")
        except Exception as exc:
            cp.mark_failed(seqid, f"exception:{exc.__class__.__name__}")
            log.warning(f"{cohort} {seqid}: SharePro exception — {exc}")

    log.info(f"{cohort}: SharePro complete ({len(results)} success, {cp.n_failed} failed)")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="SharePro colocalization")
    parser.add_argument("--cohort", choices=COHORTS + ["all"], default="all")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument(
        "--retry-failed",
        action="store_true",
        help="Include previously failed seqids from checkpoint.",
    )
    args = parser.parse_args()

    cohorts = COHORTS if args.cohort == "all" else [args.cohort]

    with RunManifest("08_coloc/sharepro.py", args=str(args)) as manifest:
        all_results = []
        for cohort in cohorts:
            all_results.extend(
                run_cohort_sharepro(
                    cohort,
                    limit=args.limit,
                    retry_failed=args.retry_failed,
                )
            )

        if all_results:
            cols = ["cohort", "seqid", "n_snps", "PP_H4", "coloc_positive"]
            df = pd.DataFrame(all_results)[cols]
            df.to_csv(SHAREPRO_OUT, sep="\t", index=False)
            log.info(f"SharePro: {len(df)} results written → {SHAREPRO_OUT}")

        manifest.n_units = len(all_results)


if __name__ == "__main__":
    main()
