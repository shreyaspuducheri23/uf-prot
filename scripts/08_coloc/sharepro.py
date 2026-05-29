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


def build_bse_input(df: pd.DataFrame, n: int,
                    beta_col: str = "beta", se_col: str = "se",
                    snp_col: str = "rsid") -> pd.DataFrame:
    """Build SNP/BETA/SE/N table for SharePro input."""
    df = df.copy()
    df["SNP"]  = df[snp_col].astype(str)
    df["BETA"] = df[beta_col].astype(float)
    df["SE"]   = df[se_col].astype(float)
    df["N"]    = n
    return df[["SNP", "BETA", "SE", "N"]].dropna()


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

    # Position-based fallback when rsid matching is insufficient (e.g. UKB_PPP "." rsids)
    if len(common_snps) < 5 and "pos" in exp_df.columns:
        out_pos_map = {
            int(r["base_pair_location"]): str(r["rsid"])
            for _, r in out_df.iterrows()
            if pd.notna(r.get("rsid")) and str(r.get("rsid", "")) not in (".", "")
        }
        for idx, row in exp_df.iterrows():
            try:
                pos = int(row["pos"])
            except (ValueError, KeyError, TypeError):
                continue
            out_rsid = out_pos_map.get(pos)
            if out_rsid and out_rsid not in common_snps:
                exp_df.at[idx, "rsid"] = out_rsid   # relabel so zscore input is consistent
                common_snps.add(out_rsid)
        if len(common_snps) >= 5:
            log.debug(f"{seqid}: position fallback found {len(common_snps)} common SNPs")

    if len(common_snps) < 5:
        log.debug(f"{seqid}: <5 common SNPs ({len(common_snps)}) — skipping SharePro")
        return None, "insufficient_common_snps"

    # LD matrix first — plink may silently drop SNPs absent from the bfile.
    # We derive the final SNP set from ld_mat.index so that the summary stats
    # and LD matrix are always the same size (SharePro asserts this).
    snp_candidates = sorted(common_snps)
    try:
        ld_mat = r_square_matrix(snp_candidates)
    except Exception as exc:
        log.warning(f"{seqid}: LD matrix failed — {exc}")
        return None, "ld_matrix_failed"

    actual_snps = list(ld_mat.index)  # subset plink kept in reference
    if len(actual_snps) < 5:
        log.debug(f"{seqid}: only {len(actual_snps)} SNPs in LD reference — skipping")
        return None, "insufficient_snps_in_ld_reference"

    exp_sub = exp_df[exp_df["rsid"].isin(actual_snps)].drop_duplicates("rsid")
    out_sub = out_df[out_df["rsid"].isin(actual_snps)].drop_duplicates("rsid")

    exp_bse = build_bse_input(exp_sub.set_index("rsid").loc[actual_snps].reset_index(),
                               n_exp, snp_col="rsid")
    out_bse = build_bse_input(out_sub.set_index("rsid").loc[actual_snps].reset_index(),
                               N_out, snp_col="rsid")
    if exp_bse.empty or out_bse.empty:
        return None, "empty_sharepro_inputs"

    with tempfile.TemporaryDirectory(prefix=f"sharepro_{seqid}_") as tmp:
        tmp = Path(tmp)

        exp_bse_path = tmp / "exp_bse.txt"
        out_bse_path = tmp / "out_bse.txt"
        ld_path      = tmp / "ld.txt"
        save_prefix  = tmp / "result"

        exp_bse.to_csv(exp_bse_path, sep="\t", index=False)
        out_bse.to_csv(out_bse_path, sep="\t", index=False)
        ld_mat.to_csv(ld_path, sep="\t", header=False, index=False)

        cmd = [
            "python", str(SHAREPRO_SCRIPT),
            "--z", str(exp_bse_path), str(out_bse_path),
            "--ld", str(ld_path),
            "--save", str(save_prefix),
            "--K", "10",
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            log.warning(f"{seqid}: SharePro failed — {res.stderr[:300]}")
            return None, "sharepro_subprocess_failed"

        result_file = Path(str(save_prefix) + ".sharepro.txt")
        if not result_file.exists():
            log.warning(f"{seqid}: SharePro produced no output")
            return None, "sharepro_output_missing"

        # SharePro writes cs/share/variantProb; share = per-effect-group coloc probability
        sp_df = pd.read_csv(result_file, sep="\t")
        shares = (
            pd.to_numeric(sp_df["share"], errors="coerce").dropna()
            if "share" in sp_df.columns
            else pd.Series([], dtype=float)
        )
        pp_h4 = float(shares.max()) if not shares.empty else 0.0

        return {
            "seqid": seqid,
            "n_snps": len(actual_snps),
            "PP_H4": pp_h4,
            "coloc_positive": pp_h4 >= 0.8,
            "raw": sp_df.to_dict(orient="records"),
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
