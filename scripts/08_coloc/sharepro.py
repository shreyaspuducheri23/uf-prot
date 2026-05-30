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


def _valid_snp_id(value) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip() not in {"", ".", "nan", "None"}


def _variant_key(
    df: pd.DataFrame,
    chrom_col: str,
    pos_col: str,
    ea_col: str,
    oa_col: str,
) -> pd.Series:
    chrom = (
        df[chrom_col]
        .astype("string")
        .str.strip()
        .str.replace(r"^chr", "", case=False, regex=True)
    )
    pos = pd.to_numeric(df[pos_col], errors="coerce").astype("Int64").astype("string")
    ea = df[ea_col].astype("string").str.strip().str.upper()
    oa = df[oa_col].astype("string").str.strip().str.upper()

    valid = (
        chrom.notna()
        & pos.notna()
        & ea.notna()
        & oa.notna()
        & ~chrom.isin(["", "."])
        & ~ea.isin(["", "."])
        & ~oa.isin(["", "."])
    )
    key = chrom + ":" + pos + ":" + ea + ":" + oa
    return key.where(valid)


def _align_sharepro_variants(
    exp_df: pd.DataFrame,
    out_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Return exposure/outcome rows matched by chrom:pos:alleles.

    Outcome rows are aligned to the exposure effect allele. Reverse allele
    matches have outcome beta and EAF flipped. The returned frames include a
    sharepro_snp column with the rsID used for SharePro input and LD lookup.
    """
    required_exp = {"chrom", "pos", "EA", "OA", "beta", "se", "rsid"}
    required_out = {
        "chromosome", "base_pair_location", "effect_allele", "other_allele",
        "beta", "se", "EAF", "rsid",
    }
    if not required_exp.issubset(exp_df.columns) or not required_out.issubset(out_df.columns):
        return pd.DataFrame(), pd.DataFrame()

    exp = exp_df.copy()
    out = out_df.copy()

    exp["match_key"] = _variant_key(exp, "chrom", "pos", "EA", "OA")
    out["key_fwd"] = _variant_key(
        out, "chromosome", "base_pair_location", "effect_allele", "other_allele"
    )
    out["key_rev"] = _variant_key(
        out, "chromosome", "base_pair_location", "other_allele", "effect_allele"
    )

    exp = exp[exp["match_key"].notna()].drop_duplicates("match_key")
    if exp.empty:
        return pd.DataFrame(), pd.DataFrame()

    exp_keys = set(exp["match_key"])
    out["match_key"] = pd.NA
    fwd = out["key_fwd"].isin(exp_keys)
    rev = out["key_rev"].isin(exp_keys) & ~fwd
    out.loc[fwd, "match_key"] = out.loc[fwd, "key_fwd"]
    out.loc[rev, "match_key"] = out.loc[rev, "key_rev"]
    out.loc[rev, "beta"] = -pd.to_numeric(out.loc[rev, "beta"], errors="coerce")
    out.loc[rev, "EAF"] = 1 - pd.to_numeric(out.loc[rev, "EAF"], errors="coerce")
    out = out[out["match_key"].notna()].drop_duplicates("match_key")
    if out.empty:
        return pd.DataFrame(), pd.DataFrame()

    matched_keys = out["match_key"].tolist()
    exp_aligned = exp.set_index("match_key").loc[matched_keys].reset_index()
    out_aligned = out.set_index("match_key").loc[matched_keys].reset_index()

    snp_ids = []
    for exp_rsid, out_rsid in zip(exp_aligned["rsid"], out_aligned["rsid"], strict=True):
        if _valid_snp_id(out_rsid):
            snp_ids.append(str(out_rsid).strip())
        elif _valid_snp_id(exp_rsid):
            snp_ids.append(str(exp_rsid).strip())
        else:
            snp_ids.append(pd.NA)

    exp_aligned["sharepro_snp"] = snp_ids
    out_aligned["sharepro_snp"] = snp_ids
    valid_snp = exp_aligned["sharepro_snp"].notna()
    exp_aligned = exp_aligned[valid_snp].copy()
    out_aligned = out_aligned[valid_snp].copy()

    if exp_aligned.empty:
        return pd.DataFrame(), pd.DataFrame()

    unique_snp = ~exp_aligned["sharepro_snp"].duplicated()
    return exp_aligned[unique_snp].copy(), out_aligned[unique_snp].copy()


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

    out_df = out_df.rename(columns={
        "beta": "beta", "standard_error": "se",
        "effect_allele_frequency": "EAF",
    })

    exp_aligned, out_aligned = _align_sharepro_variants(exp_df, out_df)
    common_snps = set(exp_aligned["sharepro_snp"]) if not exp_aligned.empty else set()

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

    exp_sub = exp_aligned[exp_aligned["sharepro_snp"].isin(actual_snps)].drop_duplicates("sharepro_snp")
    out_sub = out_aligned[out_aligned["sharepro_snp"].isin(actual_snps)].drop_duplicates("sharepro_snp")

    exp_bse = build_bse_input(exp_sub.set_index("sharepro_snp").loc[actual_snps].reset_index(),
                               n_exp, snp_col="sharepro_snp")
    out_bse = build_bse_input(out_sub.set_index("sharepro_snp").loc[actual_snps].reset_index(),
                               N_out, snp_col="sharepro_snp")
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
