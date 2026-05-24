#!/usr/bin/env python3
"""
09_assemble/cross_cohort.py

Cross-cohort gene-level summary.

Reads processed_data/final_results.tsv and produces
processed_data/gene_summary.tsv: one row per gene, sorted by the primary
cohort (UKB_female) p-value.

Analysis structure
------------------
PRIMARY_COHORT ("UKB_female") is the discovery cohort: the outcome GWAS
(Kim et al.) is female-only, so the sex-matched UKB_female pQTL data is the
appropriate primary source.  Gene inclusion requires UKB_female FDR < 0.05.
BH-FDR is applied across UKB_female p-values only.

Three independent cohorts (ARIC_EA, deCODE, Fenland) serve as replication:
they contribute a replication IVW meta-estimate and replication counts, but
do not gate discovery.  UKB_PPP is excluded from replication because it draws
from the same UK Biobank participants as UKB_female (double-counting).

Genes where UKB_female data is absent (protein not in ProteoNexus panel) are
written separately to gene_summary_no_ukb_female.tsv for reference.

Usage:
  uv run python scripts/09_assemble/cross_cohort.py
"""
from __future__ import annotations

import logging
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import COHORTS, FINAL_RESULTS, GENE_SUMMARY, PROCESSED

log = setup_logger("09_assemble_cross_cohort")

PRIMARY_COHORT = "UKB_female"
# UKB_PPP draws from the same UK Biobank participants as UKB_female — including
# it as a replication cohort would double-count the same individuals.
# Only independent cohorts are used for replication.
SENSITIVITY_COHORTS = [c for c in COHORTS if c not in {PRIMARY_COHORT, "UKB_PPP"}]

# Path for genes that had no UKB_female data
GENE_SUMMARY_NO_PRIMARY = PROCESSED / "gene_summary_no_ukb_female.tsv"


# ── Core stats ───────────────────────────────────────────────────────────────

def meta_analysis(betas: List[float], ses: List[float]) -> Dict:
    """Fixed-effects IVW meta-analysis.

    Returns a dict with keys:
      pooled_beta, pooled_se, meta_pval, meta_OR,
      meta_OR_lo95, meta_OR_hi95, I2, Q_pval, direction_consistent
    """
    betas = np.array(betas, dtype=float)
    ses = np.array(ses, dtype=float)

    ws = 1.0 / ses**2
    pooled_beta = np.sum(ws * betas) / np.sum(ws)
    pooled_se = np.sqrt(1.0 / np.sum(ws))
    z = pooled_beta / pooled_se
    meta_pval = float(2.0 * norm.sf(abs(z)))

    n = len(betas)
    Q = float(np.sum(ws * (betas - pooled_beta) ** 2))
    I2 = float(max(0.0, (Q - (n - 1)) / Q) * 100) if Q > 0 else 0.0
    Q_pval = float(chi2.sf(Q, df=n - 1)) if n > 1 else float("nan")

    direction_consistent = bool(len(set(np.sign(betas))) == 1)

    return dict(
        pooled_beta=float(pooled_beta),
        pooled_se=float(pooled_se),
        meta_pval=meta_pval,
        meta_OR=float(np.exp(pooled_beta)),
        meta_OR_lo95=float(np.exp(pooled_beta - 1.96 * pooled_se)),
        meta_OR_hi95=float(np.exp(pooled_beta + 1.96 * pooled_se)),
        I2=I2,
        Q_pval=Q_pval,
        direction_consistent=direction_consistent,
    )


# ── BH-FDR ───────────────────────────────────────────────────────────────────

def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction, returning q-values (same length)."""
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    q = pvals * n / ranks
    q_sorted = q[order]
    for i in range(n - 2, -1, -1):
        q_sorted[i] = min(q_sorted[i], q_sorted[i + 1])
    result = np.empty(n)
    result[order] = q_sorted
    return np.minimum(result, 1.0)


# ── Deduplication ─────────────────────────────────────────────────────────────

def dedup_seqids(gene_df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """For each cohort, keep the seqid with the lowest p-value.

    Returns:
        (deduplicated DataFrame, {cohort: n_seqids_available})
    """
    n_seqids: Dict[str, int] = {}
    kept_rows = []
    for cohort, cdf in gene_df.groupby("cohort"):
        n_seqids[cohort] = len(cdf)
        best = cdf.loc[cdf["pval"].idxmin()]
        kept_rows.append(best)
    return pd.DataFrame(kept_rows), n_seqids


# ── Per-cohort dict builder ──────────────────────────────────────────────────

def _cohort_dict(row: pd.Series) -> dict:
    """Extract per-cohort fields from a deduped result row."""
    return {
        "seqid": row["seqid"],
        "OR": float(row["OR"]),
        "OR_lo95": float(row["OR_lo95"]),
        "OR_hi95": float(row["OR_hi95"]),
        "pval": float(row["pval"]),
        "fdr_q": float(row["fdr_q"]),
        "fdr_pass": bool(row.get("fdr_pass", False)),
        "n_snps": int(row["n_snps"]),
        "beta": float(row["beta"]),
        "se": float(row["se"]),
        "passes_sensitivity": (
            bool(row["passes_sensitivity"])
            if pd.notna(row.get("passes_sensitivity"))
            else float("nan")
        ),
        "steiger_correct": (
            bool(row["steiger_correct"])
            if pd.notna(row.get("steiger_correct"))
            else float("nan")
        ),
        "direction_consistent": (
            bool(row["direction_consistent"])
            if pd.notna(row.get("direction_consistent"))
            else float("nan")
        ),
        "sharepro_PP_H4": (
            float(row["sharepro_PP_H4"])
            if pd.notna(row.get("sharepro_PP_H4"))
            else float("nan")
        ),
        "coloc_abf_PP_H4": (
            float(row["coloc_abf_PP_H4"])
            if pd.notna(row.get("coloc_abf_PP_H4"))
            else float("nan")
        ),
    }


# ── Summary builder ──────────────────────────────────────────────────────────

def build_gene_summary(
    df: pd.DataFrame,
    primary_cohort: str = PRIMARY_COHORT,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build gene-level summary with primary_cohort as the discovery anchor.

    Returns:
        (primary_summary, no_primary_summary)
        primary_summary: genes where primary_cohort passes FDR
        no_primary_summary: genes that passed FDR in other cohorts but have
                            no primary_cohort data (reference only)
    """
    sensitivity_cohorts = [c for c in COHORTS if c != primary_cohort]

    # All genes with any FDR pass
    fdr_genes = set(df.loc[df["fdr_pass"] == True, "gene"])
    df = df[df["gene"].isin(fdr_genes)].copy()

    if df.empty:
        log.warning("No FDR-passing genes found — output will be empty")
        return pd.DataFrame(), pd.DataFrame()

    # Split: genes with / without primary cohort data
    primary_genes = set(
        df.loc[(df["cohort"] == primary_cohort) & (df["fdr_pass"] == True), "gene"]
    )
    no_primary_genes = fdr_genes - primary_genes

    log.info(f"Primary ({primary_cohort}) FDR-passing genes: {len(primary_genes)}")
    log.info(f"Genes with no {primary_cohort} data (reference only): {len(no_primary_genes)}")

    def _build_rows(gene_set: set, require_primary: bool) -> list[dict]:
        rows = []
        for gene in sorted(gene_set):
            gdf = df[df["gene"] == gene].copy()
            gdf_dedup, n_seqids_map = dedup_seqids(gdf)

            # ── Primary cohort data ──────────────────────────────────────────
            primary_row = gdf_dedup[gdf_dedup["cohort"] == primary_cohort]
            if require_primary and primary_row.empty:
                continue

            primary_data: Optional[dict] = None
            if not primary_row.empty:
                primary_data = _cohort_dict(primary_row.iloc[0])
                primary_data["n_seqids"] = n_seqids_map.get(primary_cohort, 1)

            # ── Sensitivity cohort data ──────────────────────────────────────
            sens_rows = gdf_dedup[gdf_dedup["cohort"].isin(sensitivity_cohorts)]
            cohort_data: Dict[str, dict] = {}
            for _, row in sens_rows.iterrows():
                c = row["cohort"]
                cd = _cohort_dict(row)
                cd["n_seqids"] = n_seqids_map.get(c, 1)
                cohort_data[c] = cd

            # Replication meta-analysis (IVW of sensitivity cohorts)
            rep_betas = [cd["beta"] for cd in cohort_data.values()]
            rep_ses   = [cd["se"]   for cd in cohort_data.values()]
            rep_ma    = meta_analysis(rep_betas, rep_ses) if rep_betas else None

            # Pooled meta-analysis (all cohorts, kept for reference)
            all_betas = ([primary_data["beta"]] if primary_data else []) + rep_betas
            all_ses   = ([primary_data["se"]]   if primary_data else []) + rep_ses
            pooled_ma = meta_analysis(all_betas, all_ses) if all_betas else None

            # Replication direction agreement: "x/n" where x = cohorts matching
            # primary direction, n = sensitivity cohorts with data.
            if primary_data and rep_betas:
                primary_sign = np.sign(primary_data["beta"])
                n_agree = sum(1 for b in rep_betas if np.sign(b) == primary_sign)
                direction_consistent_replication = f"{n_agree}/{len(rep_betas)}"
            else:
                direction_consistent_replication = ""

            n_replication_tested  = len(rep_betas)
            n_replication_fdr     = sum(1 for cd in cohort_data.values() if cd["fdr_pass"])
            n_replication_nominal = sum(1 for cd in cohort_data.values() if cd["pval"] < 0.05)

            # ── Colocalization ───────────────────────────────────────────────
            all_cohort_data = ({primary_cohort: primary_data} if primary_data else {})
            all_cohort_data.update(cohort_data)

            sp_pp4  = [cd["sharepro_PP_H4"]  for cd in all_cohort_data.values() if not np.isnan(cd["sharepro_PP_H4"])]
            abf_pp4 = [cd["coloc_abf_PP_H4"] for cd in all_cohort_data.values() if not np.isnan(cd["coloc_abf_PP_H4"])]
            best_sharepro  = max(sp_pp4)  if sp_pp4  else float("nan")
            best_coloc_abf = max(abf_pp4) if abf_pp4 else float("nan")

            n_cohorts_coloc_pos = 0
            coloc_methods_agree = False
            coloc_discordant    = False
            for cd in all_cohort_data.values():
                sp = cd["sharepro_PP_H4"]
                ab = cd["coloc_abf_PP_H4"]
                sp_pos = (not np.isnan(sp)) and sp >= 0.8
                ab_pos = (not np.isnan(ab)) and ab >= 0.8
                sp_neg = (not np.isnan(sp)) and sp < 0.1
                ab_neg = (not np.isnan(ab)) and ab < 0.1
                if sp_pos or ab_pos:
                    n_cohorts_coloc_pos += 1
                if sp_pos and ab_pos:
                    coloc_methods_agree = True
                if (sp_pos and ab_neg) or (ab_pos and sp_neg):
                    coloc_discordant = True

            # ── Assemble output row ──────────────────────────────────────────
            out: dict = {
                "gene": gene,
                "primary_cohort": primary_cohort,

                # Primary (UKB_female) — the discovery estimate
                "primary_OR":       primary_data["OR"]      if primary_data else float("nan"),
                "primary_OR_lo95":  primary_data["OR_lo95"] if primary_data else float("nan"),
                "primary_OR_hi95":  primary_data["OR_hi95"] if primary_data else float("nan"),
                "primary_pval":     primary_data["pval"]    if primary_data else float("nan"),
                "primary_fdr_q":    float("nan"),  # filled after BH pass
                "primary_n_snps":   primary_data["n_snps"]  if primary_data else float("nan"),
                "primary_coloc_sharepro_PP_H4":  primary_data["sharepro_PP_H4"]  if primary_data else float("nan"),
                "primary_coloc_abf_PP_H4":       primary_data["coloc_abf_PP_H4"] if primary_data else float("nan"),

                # Replication (sensitivity cohorts IVW)
                "n_replication_tested":              n_replication_tested,
                "n_replication_fdr":                 n_replication_fdr,
                "n_replication_nominal":             n_replication_nominal,
                "direction_consistent_replication":  direction_consistent_replication,
                "replication_meta_OR":       rep_ma["meta_OR"]       if rep_ma else float("nan"),
                "replication_meta_OR_lo95":  rep_ma["meta_OR_lo95"]  if rep_ma else float("nan"),
                "replication_meta_OR_hi95":  rep_ma["meta_OR_hi95"]  if rep_ma else float("nan"),
                "replication_meta_pval":     rep_ma["meta_pval"]     if rep_ma else float("nan"),
                "replication_I2":            rep_ma["I2"]            if rep_ma else float("nan"),

                # Pooled (all cohorts, reference)
                "pooled_meta_OR":        pooled_ma["meta_OR"]           if pooled_ma else float("nan"),
                "pooled_meta_OR_lo95":   pooled_ma["meta_OR_lo95"]      if pooled_ma else float("nan"),
                "pooled_meta_OR_hi95":   pooled_ma["meta_OR_hi95"]      if pooled_ma else float("nan"),
                "pooled_meta_pval":      pooled_ma["meta_pval"]         if pooled_ma else float("nan"),
                "pooled_I2":             pooled_ma["I2"]                if pooled_ma else float("nan"),
                "pooled_Q_pval":         pooled_ma["Q_pval"]            if pooled_ma else float("nan"),
                "pooled_direction_consistent": pooled_ma["direction_consistent"] if pooled_ma else float("nan"),
                "n_cohorts_tested":      len(all_cohort_data),

                # Colocalization summary (across all cohorts)
                "n_cohorts_coloc_pos":  n_cohorts_coloc_pos,
                "best_sharepro_PP_H4":  best_sharepro,
                "best_coloc_abf_PP_H4": best_coloc_abf,
                "coloc_methods_agree":  coloc_methods_agree,
                "coloc_discordant":     coloc_discordant,
            }

            # Per-cohort columns (primary + all sensitivity cohorts)
            for cohort in COHORTS:
                if cohort == primary_cohort:
                    cd = primary_data or {}
                else:
                    cd = cohort_data.get(cohort, {})
                out[f"seqid_{cohort}"]                  = cd.get("seqid",                float("nan"))
                out[f"n_seqids_{cohort}"]               = cd.get("n_seqids",             float("nan"))
                out[f"OR_{cohort}"]                     = cd.get("OR",                   float("nan"))
                out[f"OR_lo95_{cohort}"]                = cd.get("OR_lo95",              float("nan"))
                out[f"OR_hi95_{cohort}"]                = cd.get("OR_hi95",              float("nan"))
                out[f"pval_{cohort}"]                   = cd.get("pval",                 float("nan"))
                out[f"fdr_q_{cohort}"]                  = cd.get("fdr_q",               float("nan"))
                out[f"n_snps_{cohort}"]                 = cd.get("n_snps",              float("nan"))
                out[f"passes_sensitivity_{cohort}"]     = cd.get("passes_sensitivity",  float("nan"))
                out[f"steiger_correct_{cohort}"]        = cd.get("steiger_correct",     float("nan"))
                out[f"direction_consistent_{cohort}"]   = cd.get("direction_consistent",float("nan"))
                out[f"sharepro_PP_H4_{cohort}"]         = cd.get("sharepro_PP_H4",      float("nan"))
                out[f"coloc_abf_PP_H4_{cohort}"]        = cd.get("coloc_abf_PP_H4",     float("nan"))

            rows.append(out)
        return rows

    # Build primary table (UKB_female FDR-gated)
    primary_rows = _build_rows(primary_genes, require_primary=True)
    primary_df = pd.DataFrame(primary_rows) if primary_rows else pd.DataFrame()

    if not primary_df.empty:
        primary_df["primary_fdr_q"] = _bh_fdr(primary_df["primary_pval"].values)
        primary_df = primary_df.sort_values("primary_pval").reset_index(drop=True)

    # Build reference table (no UKB_female data) using previous all-cohort logic
    no_primary_rows = _build_rows(no_primary_genes, require_primary=False)
    no_primary_df = pd.DataFrame(no_primary_rows) if no_primary_rows else pd.DataFrame()

    if not no_primary_df.empty:
        # Sort by pooled_meta_pval for the reference table
        no_primary_df = no_primary_df.sort_values("pooled_meta_pval").reset_index(drop=True)

    return primary_df, no_primary_df


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    with RunManifest("09_assemble_cross_cohort.py") as manifest:
        if not FINAL_RESULTS.exists():
            log.error(f"Missing: {FINAL_RESULTS}  — run step 9 first")
            sys.exit(1)

        log.info(f"Reading {FINAL_RESULTS}")
        df = pd.read_csv(FINAL_RESULTS, sep="\t")
        log.info(f"  {len(df):,} rows, {df['gene'].nunique():,} genes")

        primary_summary, no_primary_summary = build_gene_summary(df)

        # Write primary results
        n_primary = len(primary_summary)
        GENE_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
        primary_summary.to_csv(GENE_SUMMARY, sep="\t", index=False)
        log.info(f"Wrote {n_primary} primary genes → {GENE_SUMMARY}")

        # Write reference table
        n_no_primary = len(no_primary_summary)
        if n_no_primary > 0:
            no_primary_summary.to_csv(GENE_SUMMARY_NO_PRIMARY, sep="\t", index=False)
            log.info(
                f"Wrote {n_no_primary} genes without {PRIMARY_COHORT} data "
                f"(reference) → {GENE_SUMMARY_NO_PRIMARY}"
            )

        # Console preview
        print(f"\n{'='*72}")
        print(f"Primary discovery ({PRIMARY_COHORT})  —  {n_primary} genes")
        print("=" * 72)
        if not primary_summary.empty:
            preview_cols = [
                "gene", "primary_OR", "primary_pval", "primary_fdr_q",
                "n_replication_fdr", "direction_consistent_replication",
                "best_sharepro_PP_H4", "best_coloc_abf_PP_H4",
            ]
            present = [c for c in preview_cols if c in primary_summary.columns]
            print(primary_summary[present].head(20).to_string(index=False))

        if n_no_primary > 0:
            print(f"\n{'='*72}")
            print(
                f"Reference only (no {PRIMARY_COHORT} data)  —  {n_no_primary} genes\n"
                f"See {GENE_SUMMARY_NO_PRIMARY}"
            )
            print("=" * 72)
            if not no_primary_summary.empty:
                preview_cols = [
                    "gene", "pooled_meta_OR", "pooled_meta_pval",
                    "n_cohorts_tested", "n_replication_fdr",
                ]
                present = [c for c in preview_cols if c in no_primary_summary.columns]
                print(no_primary_summary[present].head(10).to_string(index=False))

        manifest.n_units = n_primary


if __name__ == "__main__":
    main()
