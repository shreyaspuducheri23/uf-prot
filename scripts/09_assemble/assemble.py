#!/usr/bin/env python3
"""
09_assemble/assemble.py
Join MR, sensitivity, SharePro, and coloc.abf results into final_results.tsv.
Applies tiering and prints a console summary.

Usage:
  python scripts/09_assemble/assemble.py
"""
import logging
import sys

import pandas as pd

from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import (
    COHORTS, COLOC_DIR, FINAL_RESULTS, MR_ALL_COHORTS, cohort_dir
)

log = setup_logger("09_assemble")

SHAREPRO_PATH = COLOC_DIR / "sharepro_results.tsv"


def load_mr() -> pd.DataFrame:
    if not MR_ALL_COHORTS.exists():
        log.error(f"Missing: {MR_ALL_COHORTS}")
        return pd.DataFrame()
    return pd.read_csv(MR_ALL_COHORTS, sep="\t")


def load_sensitivity() -> pd.DataFrame:
    frames = []
    for cohort in COHORTS:
        path = cohort_dir(cohort) / "sensitivity.tsv"
        if path.exists():
            df = pd.read_csv(path, sep="\t")
            df["cohort"] = cohort
            frames.append(df)
    if not frames:
        log.warning("No sensitivity results found")
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def load_sharepro() -> pd.DataFrame:
    if not SHAREPRO_PATH.exists():
        log.warning(f"SharePro results not found: {SHAREPRO_PATH}")
        return pd.DataFrame()
    return pd.read_csv(SHAREPRO_PATH, sep="\t")


def load_coloc_abf() -> pd.DataFrame:
    frames = []
    for cohort in COHORTS:
        path = COLOC_DIR / f"coloc_abf_{cohort}.tsv"
        if path.exists():
            frames.append(pd.read_csv(path, sep="\t"))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def tier(row: pd.Series) -> str:
    passes_sens = row.get("passes_sensitivity", False)
    if row.get("fdr_pass") and passes_sens and row.get("sharepro_coloc_positive"):
        if row.get("coloc_abf_positive"):
            return "Tier1_replicated"
        return "Tier1"
    if row.get("fdr_pass") and passes_sens:
        return "Tier2"
    if row.get("fdr_pass"):
        return "Tier2_nosens"
    return "Tier3"


def main() -> None:
    with RunManifest("09_assemble/assemble.py") as manifest:
        mr = load_mr()
        if mr.empty:
            log.error("No MR results — run 06_mr first")
            sys.exit(1)

        sens = load_sensitivity()
        sharepro = load_sharepro()
        coloc_abf = load_coloc_abf()

        # Merge sensitivity
        if not sens.empty:
            sens_cols = ["seqid", "cohort", "passes_sensitivity",
                         "Q_pval", "I2", "egger_intercept_pval",
                         "steiger_correct", "direction_consistent"]
            sens_sub = sens[[c for c in sens_cols if c in sens.columns]]
            mr = mr.merge(sens_sub, on=["seqid", "cohort"], how="left")

        # Merge SharePro
        if not sharepro.empty:
            sp_cols = ["seqid", "cohort", "PP_H4", "coloc_positive"]
            sharepro_sub = sharepro[[c for c in sp_cols if c in sharepro.columns]].rename(
                columns={"PP_H4": "sharepro_PP_H4", "coloc_positive": "sharepro_coloc_positive"}
            )
            mr = mr.merge(sharepro_sub, on=["seqid", "cohort"], how="left")

        # Merge coloc.abf
        if not coloc_abf.empty:
            ca_cols = ["seqid", "cohort", "PP_H4", "coloc_positive"]
            coloc_sub = coloc_abf[[c for c in ca_cols if c in coloc_abf.columns]].rename(
                columns={"PP_H4": "coloc_abf_PP_H4", "coloc_positive": "coloc_abf_positive"}
            )
            mr = mr.merge(coloc_sub, on=["seqid", "cohort"], how="left")

        # Proteins without sensitivity results default to False (unknown, not assumed pass).
        if "passes_sensitivity" in mr.columns:
            n_missing = mr["passes_sensitivity"].isna().sum()
            if n_missing:
                log.warning(
                    f"{n_missing} proteins have no sensitivity results — "
                    f"passes_sensitivity defaulted to False"
                )
            mr["passes_sensitivity"] = mr["passes_sensitivity"].fillna(False)

        mr["tier"] = mr.apply(tier, axis=1)

        # Sort by tier, then FDR q
        tier_order = {"Tier1_replicated": 0, "Tier1": 1, "Tier2": 2, "Tier2_nosens": 3, "Tier3": 4}
        mr["_tier_rank"] = mr["tier"].map(tier_order).fillna(99)
        mr = mr.sort_values(["_tier_rank", "fdr_q"]).drop(columns=["_tier_rank"])

        FINAL_RESULTS.parent.mkdir(parents=True, exist_ok=True)
        mr.to_csv(FINAL_RESULTS, sep="\t", index=False)
        log.info(f"Written: {FINAL_RESULTS} ({len(mr)} rows)")

        # Console summary
        print("\n" + "=" * 60)
        print("RESULTS SUMMARY")
        print("=" * 60)
        for t in ["Tier1_replicated", "Tier1", "Tier2", "Tier2_nosens"]:
            subset = mr[mr["tier"] == t]
            if not subset.empty:
                print(f"\n{t}: {len(subset)} protein-cohort associations")
                for _, row in subset.head(10).iterrows():
                    gene = row.get("gene", row["seqid"])
                    or_val = row.get("OR", "NA")
                    p = row.get("pval", "NA")
                    q = row.get("fdr_q", "NA")
                    print(f"  {gene:20s} OR={or_val:.3f}  p={p:.2e}  q={q:.3f}  [{row['cohort']}]")
        print("=" * 60 + "\n")

        manifest.n_units = len(mr)


if __name__ == "__main__":
    main()
