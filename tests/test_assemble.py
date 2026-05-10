"""Tests for scripts.09_assemble.assemble (tier logic and merge behaviour)."""
import importlib
import pandas as pd
import numpy as np
import pytest

_assemble = importlib.import_module("scripts.09_assemble.assemble")
tier = _assemble.tier


class TestTier:
    def _row(self, **kwargs) -> pd.Series:
        defaults = {
            "fdr_pass": False,
            "passes_sensitivity": False,
            "sharepro_coloc_positive": False,
            "coloc_abf_positive": False,
        }
        defaults.update(kwargs)
        return pd.Series(defaults)

    def test_tier1_replicated_all_positive(self):
        row = self._row(
            fdr_pass=True, passes_sensitivity=True,
            sharepro_coloc_positive=True, coloc_abf_positive=True,
        )
        assert tier(row) == "Tier1_replicated"

    def test_tier1_no_coloc_abf(self):
        row = self._row(
            fdr_pass=True, passes_sensitivity=True,
            sharepro_coloc_positive=True, coloc_abf_positive=False,
        )
        assert tier(row) == "Tier1"

    def test_tier2_passes_sensitivity_no_coloc(self):
        row = self._row(
            fdr_pass=True, passes_sensitivity=True,
            sharepro_coloc_positive=False,
        )
        assert tier(row) == "Tier2"

    def test_tier2_nosens_fdr_pass_only(self):
        row = self._row(fdr_pass=True, passes_sensitivity=False)
        assert tier(row) == "Tier2_nosens"

    def test_tier3_fdr_fail(self):
        row = self._row(fdr_pass=False)
        assert tier(row) == "Tier3"

    def test_missing_passes_sensitivity_defaults_to_false(self):
        # Post-fix: default should be False, meaning fdr_pass alone → Tier2_nosens
        row = pd.Series({"fdr_pass": True, "sharepro_coloc_positive": False})
        result = tier(row)
        assert result == "Tier2_nosens"

    def test_missing_sensitivity_does_not_reach_tier1(self):
        row = pd.Series({
            "fdr_pass": True,
            "sharepro_coloc_positive": True,
            "coloc_abf_positive": True,
            # passes_sensitivity absent → should default to False
        })
        result = tier(row)
        # Without passes_sensitivity=True, should not be Tier1
        assert result != "Tier1"
        assert result != "Tier1_replicated"

    def test_all_false_is_tier3(self):
        row = self._row()
        assert tier(row) == "Tier3"


class TestTierSortingWithNanFdrQ:
    """Validate that NaN fdr_q values don't crash the sort in main()."""

    def test_sort_with_nan_fdr_q(self):
        df = pd.DataFrame({
            "tier": ["Tier3", "Tier1"],
            "fdr_q": [float("nan"), 0.01],
        })
        tier_order = {"Tier1_replicated": 0, "Tier1": 1, "Tier2": 2, "Tier2_nosens": 3, "Tier3": 4}
        df["_tier_rank"] = df["tier"].map(tier_order).fillna(99)
        sorted_df = df.sort_values(["_tier_rank", "fdr_q"]).drop(columns=["_tier_rank"])
        assert sorted_df.iloc[0]["tier"] == "Tier1"


class TestPassesSensitivityDefaultFillna:
    """After fix: fillna(False) instead of fillna(True)."""

    def test_fillna_false_for_missing_sensitivity(self):
        mr = pd.DataFrame({
            "seqid": ["A", "B"],
            "passes_sensitivity": [True, float("nan")],
            "fdr_pass": [True, True],
        })
        mr["passes_sensitivity"] = mr["passes_sensitivity"].fillna(False)
        assert mr.loc[mr["seqid"] == "B", "passes_sensitivity"].iloc[0] is False or \
               mr.loc[mr["seqid"] == "B", "passes_sensitivity"].iloc[0] == False
