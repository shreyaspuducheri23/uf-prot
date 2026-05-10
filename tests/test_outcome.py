"""Tests for code.lib.outcome — requires Kim GWAS tabix index."""
import pytest
import pandas as pd

from scripts.lib.paths import KIM_GWAS


pytestmark = pytest.mark.skipif(
    not KIM_GWAS.exists(),
    reason="Kim GWAS file not present",
)


from scripts.lib.outcome import OutcomeLookup, KIM_N, normalize_outcome_row


class TestOutcomeLookup:
    @pytest.fixture(scope="class")
    def lookup(self):
        with OutcomeLookup() as lkp:
            yield lkp

    def test_fetch_region_returns_dataframe(self, lookup):
        df = lookup.fetch_region("1", 1_000_000, 1_100_000)
        assert isinstance(df, pd.DataFrame)

    def test_fetch_region_columns(self, lookup):
        df = lookup.fetch_region("1", 1_000_000, 1_100_000)
        assert "effect_allele" in df.columns
        assert "beta" in df.columns
        assert "N" in df.columns

    def test_n_constant(self, lookup):
        df = lookup.fetch_region("1", 1_000_000, 1_100_000)
        if not df.empty:
            assert all(df["N"] == KIM_N)

    def test_nonexistent_region_returns_empty(self, lookup):
        df = lookup.fetch_region("999", 1, 100)
        assert df.empty

    def test_fetch_region_positions_in_range(self, lookup):
        start, end = 1_000_000, 1_010_000
        df = lookup.fetch_region("1", start, end)
        if not df.empty:
            assert all(df["base_pair_location"].between(start, end))


class TestNormalizeOutcomeRow:
    def _base_row(self, **kwargs) -> pd.Series:
        defaults = {
            "chromosome": "1",
            "base_pair_location": 1_000_000,
            "rsid": "rs123",
            "effect_allele": "a",
            "other_allele": "t",
            "effect_allele_frequency": 0.35,
            "beta": -0.02,
            "standard_error": 0.005,
            "p_value": 0.0001,
        }
        defaults.update(kwargs)
        return pd.Series(defaults)

    def test_basic_normalization(self):
        result = normalize_outcome_row(self._base_row())
        assert result["EA_out"] == "A"  # uppercased
        assert result["OA_out"] == "T"
        assert result["N_out"] == KIM_N
        assert result["chrom_hg38"] == "1"

    def test_already_uppercase_alleles_preserved(self):
        result = normalize_outcome_row(self._base_row(effect_allele="G", other_allele="C"))
        assert result["EA_out"] == "G"
        assert result["OA_out"] == "C"

    def test_nan_beta_gives_none(self):
        import math
        result = normalize_outcome_row(self._base_row(beta=float("nan")))
        assert result["beta_out"] is None

    def test_nan_eaf_gives_none(self):
        result = normalize_outcome_row(self._base_row(effect_allele_frequency=float("nan")))
        assert result["EAF_out"] is None
