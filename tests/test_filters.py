"""Tests for code.lib.filters"""
import pandas as pd
import pytest

from scripts.lib.filters import (
    cis_window,
    drop_ambig_palindromes,
    exclude_mhc,
    gw_significant,
    maf_above,
)


@pytest.fixture
def base_df():
    return pd.DataFrame({
        "chrom": ["1", "6", "6", "1", "X"],
        "pos":   [100_000, 29_000_000, 10_000, 200_000, 50_000],
        "EA":    ["A", "C", "A", "T", "G"],
        "OA":    ["T", "G", "T", "A", "C"],
        "EAF":   [0.50, 0.30, 0.005, 0.45, 0.20],
        "beta":  [0.1, 0.2, 0.3, 0.4, 0.5],
        "se":    [0.01, 0.02, 0.03, 0.04, 0.05],
        "pval":  [1e-10, 1e-3, 1e-9, 1e-11, 0.04],
    })


class TestMafAbove:
    def test_drops_below_threshold(self, base_df):
        # EAF=0.005 → MAF=0.005 < 0.01
        result = maf_above(base_df, 0.01)
        assert len(result) == 4
        assert all(result["EAF"] != 0.005)

    def test_keeps_high_eaf(self, base_df):
        # EAF=0.95 → MAF=0.05 ≥ 0.01
        df = base_df.copy()
        df["EAF"] = 0.95
        result = maf_above(df, 0.01)
        assert len(result) == 5

    def test_uses_maf_not_eaf(self, base_df):
        # EAF=0.99 → MAF=0.01 — should be kept at threshold 0.01 (≥)
        df = base_df.copy()
        df["EAF"] = 0.99
        result = maf_above(df, 0.01)
        assert len(result) == 5

    def test_empty_input(self):
        df = pd.DataFrame({"EAF": []})
        result = maf_above(df, 0.01)
        assert len(result) == 0

    def test_does_not_mutate_input(self, base_df):
        original_len = len(base_df)
        maf_above(base_df, 0.01)
        assert len(base_df) == original_len


class TestGwSignificant:
    def test_keeps_only_significant(self, base_df):
        result = gw_significant(base_df, p=5e-8)
        assert all(result["pval"] < 5e-8)
        assert len(result) == 3  # rows 0, 2, 3

    def test_none_pass(self, base_df):
        result = gw_significant(base_df, p=1e-15)
        assert len(result) == 0

    def test_all_pass(self, base_df):
        result = gw_significant(base_df, p=1.0)
        assert len(result) == 5


class TestExcludeMhc:
    def test_drops_mhc_hg19(self, base_df):
        # chr6:29M is inside MHC_HG19 (25M–34M)
        result = exclude_mhc(base_df, "hg19")
        assert len(result) == 4
        assert not any((result["chrom"] == "6") & (result["pos"] == 29_000_000))

    def test_keeps_non_mhc_chr6(self, base_df):
        # chr6:10000 is outside MHC
        result = exclude_mhc(base_df, "hg19")
        assert any((result["chrom"] == "6") & (result["pos"] == 10_000))

    def test_mhc_hg38_boundaries(self):
        df = pd.DataFrame({
            "chrom": ["6", "6", "6"],
            "pos":   [28_499_999, 29_000_000, 33_500_001],
            "EA":    ["A", "C", "G"],
            "OA":    ["T", "G", "T"],
        })
        result = exclude_mhc(df, "hg38")
        # First and last should survive; middle is inside MHC_HG38
        assert len(result) == 2

    def test_unknown_build_raises(self, base_df):
        with pytest.raises(ValueError):
            exclude_mhc(base_df, "hg17")


class TestDropAmbigPalindromes:
    def test_drops_at_palindromes_high_maf(self, base_df):
        # Row 0: A/T, EAF=0.5 → MAF=0.5>0.42 → drop
        # Row 3: T/A, EAF=0.45 → MAF=0.45>0.42 → drop
        result = drop_ambig_palindromes(base_df, maf_threshold=0.42)
        assert len(result) == 3

    def test_keeps_cg_palindrome_below_threshold(self):
        df = pd.DataFrame({
            "EA": ["C"], "OA": ["G"], "EAF": [0.3],
        })
        result = drop_ambig_palindromes(df, maf_threshold=0.42)
        assert len(result) == 1

    def test_keeps_non_palindromes(self):
        df = pd.DataFrame({
            "EA": ["A", "G"], "OA": ["C", "T"], "EAF": [0.5, 0.5],
        })
        result = drop_ambig_palindromes(df, maf_threshold=0.42)
        assert len(result) == 2

    def test_case_insensitive(self):
        df = pd.DataFrame({
            "EA": ["a"], "OA": ["T"], "EAF": [0.5],
        })
        result = drop_ambig_palindromes(df, maf_threshold=0.42)
        assert len(result) == 0


class TestCisWindow:
    def test_keeps_within_window(self):
        df = pd.DataFrame({
            "chrom": ["22", "22", "22"],
            "pos":   [25_212_064, 25_512_564, 25_712_565],  # TSS=25_212_564
        })
        result = cis_window(df, tss=25_212_564, gene_chrom="22",
                            build="hg19", kb=500)
        # 25_212_064: |25_212_064 - 25_212_564| = 500 ≤ 500_000 → keep
        # 25_512_564: |25_512_564 - 25_212_564| = 300_000 ≤ 500_000 → keep
        # 25_712_565: |25_712_565 - 25_212_564| = 500_001 > 500_000 → drop
        assert len(result) == 2

    def test_excludes_other_chromosomes(self):
        df = pd.DataFrame({
            "chrom": ["22", "1"],
            "pos":   [25_212_564, 25_212_564],
        })
        result = cis_window(df, tss=25_212_564, gene_chrom="22",
                            build="hg19", kb=500)
        assert len(result) == 1
        assert result["chrom"].iloc[0] == "22"

    def test_boundary_exactly_at_flank_is_kept(self):
        # Variant exactly at TSS + 500_000 should be kept (<=, not <)
        df = pd.DataFrame({"chrom": ["22"], "pos": [25_212_564 + 500_000]})
        result = cis_window(df, tss=25_212_564, gene_chrom="22", build="hg19", kb=500)
        assert len(result) == 1


class TestNaNEdgeCases:
    """Tests for NaN and missing values in filter inputs."""

    def test_maf_above_nan_eaf_dropped(self):
        df = pd.DataFrame({"EAF": [0.3, float("nan")]})
        # NaN EAF → NaN MAF → fails >= threshold → dropped
        result = maf_above(df, 0.01)
        assert len(result) == 1

    def test_drop_ambig_palindromes_nan_eaf_kept(self):
        # NaN EAF → NaN MAF → palindrome check with NaN > 0.42 is False → NOT dropped
        df = pd.DataFrame({"EA": ["A"], "OA": ["T"], "EAF": [float("nan")]})
        result = drop_ambig_palindromes(df, maf_threshold=0.42)
        # palindrome with NaN MAF: NaN > 0.42 is False → kept
        assert len(result) == 1

    def test_maf_above_eaf_exactly_at_threshold_kept(self):
        # EAF=0.01 → MAF=0.01 = threshold; ≥ threshold → kept
        df = pd.DataFrame({"EAF": [0.01]})
        result = maf_above(df, threshold=0.01)
        assert len(result) == 1


class TestAllSixFiltersIntegration:
    """Integration test: all 6 filters applied in sequence produce correct final counts."""

    def test_full_filter_pipeline(self):
        from scripts.lib.filters import (
            cis_window, gw_significant, maf_above, exclude_mhc, drop_ambig_palindromes
        )
        # Deliberately craft rows that each hit a different filter
        df = pd.DataFrame({
            "chrom": ["22", "22", "6",   "22",  "22"],
            "pos":   [
                25_212_064,      # kept: cis, gw_sig, good MAF, not MHC, not palindrome
                30_000_000,      # dropped: outside cis (>25.2M+500k = 25.7M)
                29_000_000,      # dropped: MHC (chr6 inside hg19 MHC)
                25_300_000,      # dropped: p not genome-wide significant
                25_200_000,      # dropped: palindrome with high MAF
            ],
            "EA":    ["A",  "C",  "G",  "A",  "A"],
            "OA":    ["G",  "T",  "A",  "G",  "T"],  # row 4: A/T palindrome
            "EAF":   [0.30, 0.30, 0.30, 0.30, 0.50],  # row 4: MAF=0.5 > 0.42
            "beta":  [0.5,  0.5,  0.5,  0.5,  0.5],
            "se":    [0.05, 0.05, 0.05, 0.05, 0.05],
            "pval":  [1e-9, 1e-9, 1e-9, 0.05, 1e-9],  # row 3: p not gw_sig
            "N":     [7213, 7213, 7213, 7213, 7213],
        })

        tss = 25_212_564
        df2 = cis_window(df, tss, "22", "hg19", kb=500)
        df2 = gw_significant(df2, p=5e-8)
        df2 = maf_above(df2, threshold=0.01)
        df2 = exclude_mhc(df2, "hg19")
        df2 = drop_ambig_palindromes(df2, maf_threshold=0.42)

        assert len(df2) == 1
        assert df2["pos"].iloc[0] == 25_212_064
