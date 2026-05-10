"""Tests for code.lib.fstat"""
import numpy as np
import pandas as pd
import pytest

from scripts.lib.fstat import WEAK_INSTRUMENT_THRESHOLD, add_fstat, compute_fstat


class TestComputeFstat:
    def test_basic(self):
        df = pd.DataFrame({"beta": [0.1, 0.5], "se": [0.01, 0.1]})
        result = compute_fstat(df)
        assert result.iloc[0] == pytest.approx(100.0)
        assert result.iloc[1] == pytest.approx(25.0)

    def test_zero_se_raises(self):
        df = pd.DataFrame({"beta": [0.1], "se": [0.0]})
        result = compute_fstat(df)
        assert np.isinf(result.iloc[0])

    def test_custom_columns(self):
        df = pd.DataFrame({"b": [2.0], "s": [1.0]})
        result = compute_fstat(df, beta_col="b", se_col="s")
        assert result.iloc[0] == pytest.approx(4.0)


class TestAddFstat:
    def test_adds_column(self):
        df = pd.DataFrame({"beta": [0.3], "se": [0.1]})
        result = add_fstat(df)
        assert "F_stat" in result.columns
        assert result["F_stat"].iloc[0] == pytest.approx(9.0)

    def test_does_not_mutate(self):
        df = pd.DataFrame({"beta": [1.0], "se": [0.5]})
        _ = add_fstat(df)
        assert "F_stat" not in df.columns


    def test_empty_dataframe(self):
        df = pd.DataFrame({"beta": pd.Series([], dtype=float), "se": pd.Series([], dtype=float)})
        result = add_fstat(df)
        assert "F_stat" in result.columns
        assert len(result) == 0


class TestWeakInstrumentThreshold:
    def test_threshold_value(self):
        assert WEAK_INSTRUMENT_THRESHOLD == 10.0


class TestNaNHandling:
    def test_nan_beta_gives_nan_fstat(self):
        df = pd.DataFrame({"beta": [float("nan")], "se": [0.1]})
        result = compute_fstat(df)
        assert pd.isna(result.iloc[0])

    def test_nan_se_gives_nan_fstat(self):
        df = pd.DataFrame({"beta": [0.5], "se": [float("nan")]})
        result = compute_fstat(df)
        assert pd.isna(result.iloc[0])

    def test_negative_se_gives_positive_fstat(self):
        # (beta / -se)^2 = (beta / se)^2 — F-stat is always positive
        df = pd.DataFrame({"beta": [0.5], "se": [-0.1]})
        result = compute_fstat(df)
        assert result.iloc[0] == pytest.approx(25.0)
