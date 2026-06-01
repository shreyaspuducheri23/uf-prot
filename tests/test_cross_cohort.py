"""Tests for scripts.09_assemble.cross_cohort."""
from __future__ import annotations

import importlib
import math

import numpy as np
import pandas as pd
import pytest

_mod = importlib.import_module("scripts.09_assemble.cross_cohort")
meta_analysis = _mod.meta_analysis
dedup_seqids = _mod.dedup_seqids
build_gene_summary = _mod.build_gene_summary

PRIMARY = "UKB_female"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _make_row(
    gene: str,
    cohort: str,
    seqid: str,
    beta: float,
    se: float,
    pval: float,
    fdr_pass: bool = True,
    n_snps: int = 3,
    fdr_q: float = 0.01,
    passes_sensitivity: object = True,
    steiger_correct: object = True,
    direction_consistent: object = True,
    sharepro_PP_H4: float = float("nan"),
    coloc_abf_PP_H4: float = float("nan"),
) -> dict:
    OR = math.exp(beta)
    return {
        "gene": gene,
        "cohort": cohort,
        "seqid": seqid,
        "beta": beta,
        "se": se,
        "pval": pval,
        "fdr_pass": fdr_pass,
        "fdr_pass_any": fdr_pass,
        "n_snps": n_snps,
        "OR": OR,
        "OR_lo95": math.exp(beta - 1.96 * se),
        "OR_hi95": math.exp(beta + 1.96 * se),
        "fdr_q": fdr_q,
        "method": "Wald ratio" if n_snps == 1 else "Inverse variance weighted",
        "passes_sensitivity": passes_sensitivity,
        "steiger_correct": steiger_correct,
        "direction_consistent": direction_consistent,
        "sharepro_PP_H4": sharepro_PP_H4,
        "coloc_abf_PP_H4": coloc_abf_PP_H4,
        "Q_pval": float("nan"),
        "I2": float("nan"),
        "egger_intercept_pval": float("nan"),
    }


# ── 1. Two-cohort meta-analysis ───────────────────────────────────────────────

def test_meta_analysis_two_cohorts():
    beta1, se1 = -0.5, 0.1
    beta2, se2 = -0.4, 0.2

    res = meta_analysis([beta1, beta2], [se1, se2])

    w1, w2 = 1 / se1**2, 1 / se2**2
    pooled_beta = (w1 * beta1 + w2 * beta2) / (w1 + w2)
    pooled_se = math.sqrt(1 / (w1 + w2))

    assert abs(res["pooled_beta"] - pooled_beta) < 1e-9
    assert abs(res["pooled_se"] - pooled_se) < 1e-9
    assert abs(res["meta_OR"] - math.exp(pooled_beta)) < 1e-9
    assert abs(res["meta_OR_lo95"] - math.exp(pooled_beta - 1.96 * pooled_se)) < 1e-9
    assert abs(res["meta_OR_hi95"] - math.exp(pooled_beta + 1.96 * pooled_se)) < 1e-9

    from scipy.stats import norm
    z = pooled_beta / pooled_se
    expected_p = 2 * norm.sf(abs(z))
    assert abs(res["meta_pval"] - expected_p) < 1e-12

    Q = w1 * (beta1 - pooled_beta) ** 2 + w2 * (beta2 - pooled_beta) ** 2
    expected_I2 = max(0.0, (Q - 1) / Q) * 100
    assert abs(res["I2"] - expected_I2) < 1e-6

    from scipy.stats import chi2
    assert abs(res["Q_pval"] - chi2.sf(Q, df=1)) < 1e-9
    assert res["direction_consistent"] is True


# ── 2. Single-cohort meta-analysis ───────────────────────────────────────────

def test_meta_analysis_single_cohort():
    beta, se = -0.3, 0.07
    res = meta_analysis([beta], [se])

    assert abs(res["pooled_beta"] - beta) < 1e-12
    assert abs(res["pooled_se"] - se) < 1e-12
    assert abs(res["meta_OR"] - math.exp(beta)) < 1e-12
    assert res["I2"] == 0.0
    assert math.isnan(res["Q_pval"]), "Q_pval should be NaN for n=1"


# ── 3. Direction inconsistent ─────────────────────────────────────────────────

def test_direction_inconsistent_flagged():
    res = meta_analysis([0.3, -0.2], [0.05, 0.05])
    assert res["direction_consistent"] is False


# ── 4. Multi-seqid deduplication ──────────────────────────────────────────────

def test_multi_seqid_dedup():
    rows = [
        _make_row("GENE_A", PRIMARY, "seq1", beta=-0.3, se=0.05, pval=1e-9),
        _make_row("GENE_A", PRIMARY, "seq2", beta=-0.2, se=0.06, pval=5e-4),
    ]
    df = pd.DataFrame(rows)

    gene_df = df[df["gene"] == "GENE_A"]
    dedup, n_map = dedup_seqids(gene_df)

    assert len(dedup) == 1
    assert dedup.iloc[0]["seqid"] == "seq1"  # lower pval
    assert n_map[PRIMARY] == 2

    primary_summary, _ = build_gene_summary(df, primary_cohort=PRIMARY)
    assert len(primary_summary) == 1
    assert primary_summary.iloc[0][f"n_seqids_{PRIMARY}"] == 2
    assert primary_summary.iloc[0][f"seqid_{PRIMARY}"] == "seq1"


# ── 5. Only primary-cohort FDR-passing genes in primary table ─────────────────

def test_only_primary_fdr_genes_in_primary_table():
    rows = [
        # Passes in primary cohort
        _make_row("PASSING", PRIMARY, "seq1", beta=-0.3, se=0.05, pval=1e-6, fdr_pass=True),
        # FDR-fails in primary cohort — goes to no_primary or excluded
        _make_row("FAILING", PRIMARY, "seq2", beta=-0.1, se=0.05, pval=0.3, fdr_pass=False),
        # Passes in a sensitivity cohort only — not in primary table
        _make_row("SENS_ONLY", "deCODE", "seq3", beta=-0.4, se=0.05, pval=1e-8, fdr_pass=True),
    ]
    df = pd.DataFrame(rows)
    primary_summary, no_primary_summary = build_gene_summary(df, primary_cohort=PRIMARY)

    primary_genes = set(primary_summary["gene"])
    assert "PASSING" in primary_genes
    assert "FAILING" not in primary_genes
    assert "SENS_ONLY" not in primary_genes

    # SENS_ONLY should appear in the reference table
    no_primary_genes = set(no_primary_summary["gene"]) if not no_primary_summary.empty else set()
    assert "SENS_ONLY" in no_primary_genes


# ── 6. Primary vs replication columns populated correctly ─────────────────────

def test_primary_and_replication_columns():
    """Primary OR comes from UKB_female; replication meta from other cohorts."""
    rows = [
        _make_row("GENE_A", PRIMARY,  "seq_f", beta=-0.4, se=0.05, pval=1e-10, fdr_pass=True),
        _make_row("GENE_A", "deCODE", "seq_d", beta=-0.3, se=0.08, pval=2e-5,  fdr_pass=True),
        _make_row("GENE_A", "ARIC_EA","seq_a", beta=-0.35,se=0.06, pval=5e-8,  fdr_pass=True),
    ]
    df = pd.DataFrame(rows)
    primary_summary, _ = build_gene_summary(df, primary_cohort=PRIMARY)

    assert len(primary_summary) == 1
    row = primary_summary.iloc[0]

    # Primary OR = UKB_female OR
    assert abs(row["primary_OR"] - math.exp(-0.4)) < 1e-6
    assert abs(row["primary_pval"] - 1e-10) < 1e-15

    # Replication meta = IVW of deCODE + ARIC_EA
    rep_ma = meta_analysis([-0.3, -0.35], [0.08, 0.06])
    assert abs(row["replication_meta_OR"] - rep_ma["meta_OR"]) < 1e-6

    # n_replication_fdr: 2 sensitivity cohorts pass FDR
    assert row["n_replication_fdr"] == 2

    # direction_consistent_replication: "x/n" — both sensitivity cohorts agree
    assert row["direction_consistent_replication"] == "2/2"


def test_ukb_ppp_excluded_from_analyses_but_retained():
    """UKB_PPP is retained as a cohort column but excluded from summaries."""
    rows = [
        _make_row("GENE_A", PRIMARY, "seq_f", beta=-0.4, se=0.05, pval=1e-10, fdr_pass=True),
        _make_row("GENE_A", "deCODE", "seq_d", beta=-0.3, se=0.08, pval=2e-5, fdr_pass=True),
        _make_row("GENE_A", "ARIC_EA", "seq_a", beta=-0.35, se=0.06, pval=5e-8, fdr_pass=True),
        _make_row("GENE_A", "UKB_PPP", "seq_u", beta=0.8, se=0.02, pval=1e-30, fdr_pass=True),
    ]
    df = pd.DataFrame(rows)
    primary_summary, _ = build_gene_summary(df, primary_cohort=PRIMARY)
    row = primary_summary.iloc[0]

    rep_ma = meta_analysis([-0.3, -0.35], [0.08, 0.06])
    pooled_ma = meta_analysis([-0.4, -0.3, -0.35], [0.05, 0.08, 0.06])
    assert row["n_replication_tested"] == 2
    assert row["n_replication_fdr"] == 2
    assert row["n_cohorts_tested"] == 3
    assert row["direction_consistent_replication"] == "2/2"
    assert abs(row["replication_meta_OR"] - rep_ma["meta_OR"]) < 1e-6
    assert abs(row["pooled_meta_OR"] - pooled_ma["meta_OR"]) < 1e-6
    assert row["seqid_UKB_PPP"] == "seq_u"
    assert abs(row["OR_UKB_PPP"] - math.exp(0.8)) < 1e-6


# ── 7. Direction flip in replication flagged ─────────────────────────────────

def test_direction_flip_in_replication():
    rows = [
        _make_row("GENE_B", PRIMARY,  "seq_f", beta=0.4,  se=0.05, pval=1e-8, fdr_pass=True),
        _make_row("GENE_B", "deCODE", "seq_d", beta=-0.3, se=0.08, pval=2e-4, fdr_pass=True),
    ]
    df = pd.DataFrame(rows)
    primary_summary, _ = build_gene_summary(df, primary_cohort=PRIMARY)

    row = primary_summary.iloc[0]
    # deCODE disagrees — 0 out of 1
    assert row["direction_consistent_replication"] == "0/1"


# ── 8. Coloc discordant flag ──────────────────────────────────────────────────

def test_coloc_discordant_flag():
    """SharePro ≥ 0.8 but coloc.abf < 0.1 → coloc_discordant=True."""
    rows = [
        _make_row(
            "GENE_X", PRIMARY, "seq1", beta=-0.3, se=0.05, pval=1e-8,
            sharepro_PP_H4=0.85, coloc_abf_PP_H4=0.05,
        ),
    ]
    df = pd.DataFrame(rows)
    primary_summary, _ = build_gene_summary(df, primary_cohort=PRIMARY)

    assert len(primary_summary) == 1
    assert primary_summary.iloc[0]["coloc_discordant"] == True
    assert primary_summary.iloc[0]["coloc_methods_agree"] == False


def test_coloc_methods_agree():
    """Both methods ≥ 0.8 → coloc_methods_agree=True, discordant=False."""
    rows = [
        _make_row(
            "GENE_Y", PRIMARY, "seq1", beta=-0.3, se=0.05, pval=1e-8,
            sharepro_PP_H4=0.90, coloc_abf_PP_H4=0.88,
        ),
    ]
    df = pd.DataFrame(rows)
    primary_summary, _ = build_gene_summary(df, primary_cohort=PRIMARY)

    assert primary_summary.iloc[0]["coloc_methods_agree"] == True
    assert primary_summary.iloc[0]["coloc_discordant"] == False


# ── 9. Sort order is by primary_pval ─────────────────────────────────────────

def test_sorted_by_primary_pval():
    rows = [
        _make_row("GENE_WEAK",   PRIMARY, "seq_w", beta=-0.1, se=0.04, pval=0.04, fdr_pass=True),
        _make_row("GENE_STRONG", PRIMARY, "seq_s", beta=-0.5, se=0.05, pval=1e-15, fdr_pass=True),
        _make_row("GENE_MID",    PRIMARY, "seq_m", beta=-0.3, se=0.05, pval=1e-6,  fdr_pass=True),
    ]
    df = pd.DataFrame(rows)
    primary_summary, _ = build_gene_summary(df, primary_cohort=PRIMARY)

    genes = list(primary_summary["gene"])
    assert genes.index("GENE_STRONG") < genes.index("GENE_MID") < genes.index("GENE_WEAK")


# ── 10. primary_fdr_q applied over primary p-values only ─────────────────────

def test_primary_fdr_q_uses_upstream_per_cohort_q_value():
    rows = [
        _make_row("GENE_A", PRIMARY, "seq1", beta=-0.5, se=0.05, pval=1e-10, fdr_pass=True, fdr_q=0.003),
        _make_row("GENE_B", PRIMARY, "seq2", beta=-0.3, se=0.05, pval=1e-4, fdr_pass=True, fdr_q=0.02),
    ]
    df = pd.DataFrame(rows)
    primary_summary, _ = build_gene_summary(df, primary_cohort=PRIMARY)

    assert "primary_fdr_q" in primary_summary.columns
    assert not primary_summary["primary_fdr_q"].isna().any()
    q_by_gene = dict(zip(primary_summary["gene"], primary_summary["primary_fdr_q"]))
    assert q_by_gene == {"GENE_A": 0.003, "GENE_B": 0.02}
