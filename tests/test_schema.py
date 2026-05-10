"""Tests for scripts.lib.schema.validate_norm_df."""
import numpy as np
import pandas as pd
import pytest

from scripts.lib.schema import NORM_COLS, validate_norm_df


def _good_df(n: int = 2) -> pd.DataFrame:
    """Return a minimal valid normalized DataFrame."""
    return pd.DataFrame({
        "seqid":   [f"SeqId_{i}" for i in range(n)],
        "gene":    [f"GENE{i}" for i in range(n)],
        "uniprot": [f"Q{i:05d}" for i in range(n)],
        "chrom":   ["1", "22"][:n],
        "pos":     np.array([1_000_000, 2_000_000][:n], dtype="int64"),
        "rsid":    ["rs1", "rs2"][:n],
        "EA":      ["A", "C"][:n],
        "OA":      ["G", "T"][:n],
        "EAF":     [0.3, 0.4][:n],
        "beta":    [0.1, 0.2][:n],
        "se":      [0.01, 0.02][:n],
        "pval":    [1e-9, 1e-8][:n],
        "N":       [7213, 7213][:n],
        "build":   ["hg19", "hg19"][:n],
    })


class TestValidateNormDf:
    def test_valid_dataframe_passes(self):
        validate_norm_df(_good_df())  # should not raise

    def test_missing_column_raises(self):
        df = _good_df()
        df = df.drop(columns=["pval"])
        with pytest.raises(ValueError, match="missing columns"):
            validate_norm_df(df)

    def test_chr_prefix_in_chrom_raises(self):
        df = _good_df()
        df["chrom"] = "chr1"
        with pytest.raises(ValueError, match="invalid chrom values"):
            validate_norm_df(df)

    def test_invalid_chrom_value_raises(self):
        df = _good_df()
        df["chrom"] = "99"
        with pytest.raises(ValueError, match="invalid chrom values"):
            validate_norm_df(df, where="test")

    def test_where_appears_in_error_message(self):
        df = _good_df()
        df = df.drop(columns=["beta"])
        with pytest.raises(ValueError, match="my_caller"):
            validate_norm_df(df, where="my_caller")

    def test_nullable_int64_pos_raises(self):
        df = _good_df()
        df["pos"] = df["pos"].astype("Int64")  # nullable pandas Int64
        with pytest.raises(ValueError, match="nullable dtype"):
            validate_norm_df(df)

    def test_float_pos_raises(self):
        df = _good_df()
        df["pos"] = df["pos"].astype(float)
        with pytest.raises(ValueError, match="int64"):
            validate_norm_df(df)

    def test_lowercase_ea_raises(self):
        df = _good_df()
        df["EA"] = "a"
        with pytest.raises(ValueError, match="lowercase alleles"):
            validate_norm_df(df)

    def test_empty_string_oa_raises(self):
        df = _good_df()
        df["OA"] = ""
        with pytest.raises(ValueError, match="empty strings"):
            validate_norm_df(df)

    def test_x_chromosome_passes(self):
        df = _good_df()
        df["chrom"] = "X"
        validate_norm_df(df)  # should not raise

    def test_single_row_passes(self):
        validate_norm_df(_good_df(n=1))

    def test_all_norm_cols_present(self):
        df = _good_df()
        assert all(c in df.columns for c in NORM_COLS)
