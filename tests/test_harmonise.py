"""Tests for scripts.05_harmonise.harmonise (_join_outcome logic)."""
import importlib
import pandas as pd
import numpy as np
import pytest
from unittest.mock import patch, MagicMock

_harmonise_mod = importlib.import_module("scripts.05_harmonise.harmonise")
_join_outcome = _harmonise_mod._join_outcome


def _instrument_df(**kwargs) -> pd.DataFrame:
    defaults = {
        "chrom_hg38": ["22"],
        "pos_hg38": [25_212_564],
        "rsid": ["rs123"],
        "chrom": ["22"],
        "pos": [25_212_564],
        "EA": ["A"], "OA": ["G"],
        "EAF": [0.3], "beta": [0.5], "se": [0.05], "pval": [1e-9],
        "N": [7213], "seqid": ["SeqId_TEST"],
    }
    defaults.update(kwargs)
    return pd.DataFrame(defaults)


def _outcome_df(chrom="22", pos=25_212_564, ea="A", oa="G") -> pd.DataFrame:
    return pd.DataFrame({
        "chromosome": [chrom],
        "base_pair_location": [pos],
        "effect_allele": [ea],
        "other_allele": [oa],
        "beta": [0.2],
        "standard_error": [0.02],
        "effect_allele_frequency": [0.35],
        "p_value": [1e-5],
        "rsid": ["rs123"],
        "rs_id": ["rs123"],
        "hm_coordinate_conversion": [""], "hm_code": [""], "variant_id": [""],
        "N": [434_152],
    })


class TestJoinOutcomePositionMatch:
    def test_direct_position_match(self):
        instr = _instrument_df()
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = _outcome_df()

        result, n_proxies = _join_outcome(instr, mock_outcome)

        assert not result.empty
        assert n_proxies == 0
        assert "EA_out" in result.columns
        assert result["EA_out"].iloc[0] == "A"

    def test_no_match_no_rsid_returns_empty(self):
        instr = _instrument_df(rsid=["."])  # no rsid, no position match
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = pd.DataFrame(columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
        ])

        result, n_proxies = _join_outcome(instr, mock_outcome)

        assert result.empty
        assert n_proxies == 0

    def test_duplicate_positions_in_outcome_logs_warning(self, caplog):
        instr = _instrument_df()
        # Two rows at same position in outcome
        out = pd.concat([_outcome_df(), _outcome_df()], ignore_index=True)
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = out

        with caplog.at_level("WARNING", logger="05_harmonise"):
            result, _ = _join_outcome(instr, mock_outcome)

        assert "duplicate" in caplog.text.lower()


class TestJoinOutcomeProxySearch:
    def test_proxy_used_when_direct_match_missing(self):
        instr = _instrument_df(rsid=["rs999"])
        # No direct position match
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = pd.DataFrame(columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
        ])
        # Proxy returns rs_proxy
        mock_outcome.fetch_by_rsid.return_value = _outcome_df()
        mock_outcome.fetch_by_rsid.return_value["rsid"] = "rs_proxy"

        with patch.object(_harmonise_mod, "find_proxies",
                          return_value={"rs999": ("rs_proxy", 0.95)}):
            result, n_proxies = _join_outcome(instr, mock_outcome)

        assert n_proxies == 1

    def test_no_rsid_instruments_are_dropped_with_warning(self, caplog):
        instr = _instrument_df(rsid=["."])
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = pd.DataFrame(columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
        ])

        with caplog.at_level("WARNING", logger="05_harmonise"):
            result, n_proxies = _join_outcome(instr, mock_outcome)

        assert result.empty
        assert "no rsid" in caplog.text.lower() or "cannot" in caplog.text.lower()

    def test_proxy_high_maf_excluded(self):
        instr = _instrument_df(rsid=["rs999"])
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = pd.DataFrame(columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
        ])
        proxy_row = _outcome_df()
        proxy_row["effect_allele_frequency"] = 0.5
        mock_outcome.fetch_by_rsid.return_value = proxy_row

        with patch.object(_harmonise_mod, "find_proxies",
                          return_value={"rs999": ("rs_proxy", 0.95)}):
            result, n_proxies = _join_outcome(instr, mock_outcome)

        assert n_proxies == 0
        assert result.empty


class TestCallHarmoniseR:
    def test_uses_scripts_rlib_path(self):
        df = _instrument_df()

        def fake_run(cmd, capture_output, text):
            assert cmd[0] == "Rscript"
            assert "scripts/rlib/harmonise.R" in cmd[1]
            result_path = cmd[cmd.index("--result") + 1]
            pd.DataFrame({"seqid": ["SeqId_TEST"], "rsid": ["rs123"]}).to_csv(
                result_path, sep="\t", index=False
            )

            class Result:
                returncode = 0
                stderr = ""

            return Result()

        with patch.object(_harmonise_mod.subprocess, "run", side_effect=fake_run):
            out = _harmonise_mod._call_harmonise_r(df, "SeqId_TEST")

        assert out is not None
        assert not out.empty
