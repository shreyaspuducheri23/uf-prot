"""Tests for code.lib.fdr"""
import numpy as np
import pandas as pd
import pytest

from scripts.lib.fdr import add_fdr, bh_fdr


class TestBhFdr:
    def test_monotone_after_correction(self):
        pvals = pd.Series([0.04, 0.02, 0.3, 0.001])
        adj = bh_fdr(pvals)
        # Smallest p should have smallest q
        assert adj.iloc[3] <= adj.iloc[1] <= adj.iloc[0] <= adj.iloc[2]

    def test_adjusted_leq_1(self):
        pvals = pd.Series([0.5, 0.8, 0.9, 0.99])
        adj = bh_fdr(pvals)
        assert all(adj <= 1.0)

    def test_nan_preserved(self):
        pvals = pd.Series([0.01, float("nan"), 0.05])
        adj = bh_fdr(pvals)
        assert pd.isna(adj.iloc[1])

    def test_uses_non_nan_count_for_bh_denominator(self):
        pvals = pd.Series([0.01, float("nan"), 0.02])
        adj = bh_fdr(pvals)
        # Non-NaN p-values are [0.01, 0.02], so m=2.
        # q(0.01)=min(0.01*2/1, 0.02*2/2)=0.02; q(0.02)=0.02.
        assert adj.iloc[0] == pytest.approx(0.02)
        assert pd.isna(adj.iloc[1])
        assert adj.iloc[2] == pytest.approx(0.02)

    def test_empty_series(self):
        pvals = pd.Series([], dtype=float)
        adj = bh_fdr(pvals)
        assert len(adj) == 0

    def test_single_value(self):
        pvals = pd.Series([0.03])
        adj = bh_fdr(pvals)
        assert adj.iloc[0] == pytest.approx(0.03)

    def test_all_same_pvals(self):
        pvals = pd.Series([0.05, 0.05, 0.05])
        adj = bh_fdr(pvals)
        # All equal inputs → all equal after correction
        assert len(set(adj.values)) == 1

    def test_known_values(self):
        # BH on [0.01, 0.04, 0.2, 0.4] with n=4
        # ranks: 1,2,3,4; adjusted = p * 4/rank, then monotone ceiling
        pvals = pd.Series([0.01, 0.04, 0.2, 0.4])
        adj = bh_fdr(pvals)
        # 0.01 * 4/1 = 0.04; after monotone: min(0.04, 0.08) = 0.04
        assert adj.iloc[0] == pytest.approx(0.04, abs=1e-10)


    def test_q_equal_to_alpha_is_pass(self):
        # Single p-value: q = p * 1/1 = p. With p=alpha, q==alpha → should PASS (<=, not <).
        pvals = pd.Series([0.05])
        adj = bh_fdr(pvals, alpha=0.05)
        assert adj.iloc[0] == pytest.approx(0.05)

    def test_all_nan_returns_all_nan(self):
        pvals = pd.Series([float("nan"), float("nan")])
        adj = bh_fdr(pvals)
        assert adj.isna().all()

    def test_all_pvals_one(self):
        pvals = pd.Series([1.0, 1.0, 1.0])
        adj = bh_fdr(pvals)
        assert all(adj <= 1.0)

    def test_very_small_pvals_no_underflow(self):
        pvals = pd.Series([1e-300, 1e-299])
        adj = bh_fdr(pvals)
        assert not adj.isna().any()


class TestAddFdr:
    def test_adds_fdr_q_column(self):
        df = pd.DataFrame({"pval": [0.001, 0.05, 0.5]})
        result = add_fdr(df)
        assert "fdr_q" in result.columns
        assert "fdr_pass" in result.columns

    def test_fdr_pass_threshold(self):
        df = pd.DataFrame({"pval": [0.001, 0.5]})
        result = add_fdr(df, alpha=0.05)
        assert result.iloc[0]["fdr_pass"]
        # p=0.5 → adj will be >0.05
        assert not result.iloc[1]["fdr_pass"]

    def test_q_exactly_alpha_is_pass(self):
        # Single p = alpha → q = alpha → fdr_pass must be True (<=, not <)
        df = pd.DataFrame({"pval": [0.05]})
        result = add_fdr(df, alpha=0.05)
        assert result.iloc[0]["fdr_pass"]

    def test_does_not_mutate_input(self):
        df = pd.DataFrame({"pval": [0.01, 0.5]})
        original_cols = list(df.columns)
        add_fdr(df)
        assert list(df.columns) == original_cols
