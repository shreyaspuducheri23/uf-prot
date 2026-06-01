#!/usr/bin/env python3
"""
09_assemble/assemble.py
Join MR, sensitivity, SharePro, and coloc.abf results into final_results.tsv.
Prints a per-cohort console summary of FDR-passing associations.

Usage:
  python scripts/09_assemble/assemble.py
"""
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

        mr = mr.sort_values(["cohort", "fdr_q"], na_position="last").reset_index(drop=True)

        FINAL_RESULTS.parent.mkdir(parents=True, exist_ok=True)
        mr.to_csv(FINAL_RESULTS, sep="\t", index=False)
        log.info(f"Written: {FINAL_RESULTS} ({len(mr)} rows)")

        # Console summary
        print("\n" + "=" * 60)
        print("RESULTS SUMMARY")
        print("=" * 60)
        for cohort in sorted(mr["cohort"].dropna().unique()):
            subset = mr[mr["cohort"] == cohort]
            fdr_count = int(subset.get("fdr_pass", pd.Series(False, index=subset.index)).fillna(False).sum())
            print(f"\n{cohort}: {fdr_count} FDR-passing / {len(subset)} total associations")
            top = subset.sort_values("pval", na_position="last").head(5)
            for _, row in top.iterrows():
                gene = row.get("gene", row["seqid"])
                or_val = row.get("OR", float("nan"))
                p = row.get("pval", float("nan"))
                q = row.get("fdr_q", float("nan"))
                print(f"  {gene:20s} OR={or_val:.3f}  p={p:.2e}  q={q:.3f}")
        print("=" * 60 + "\n")

        manifest.n_units = len(mr)


if __name__ == "__main__":
    main()
