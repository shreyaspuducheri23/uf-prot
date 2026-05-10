"""Tests for scripts.02_cis_pqtl_extract.decode (read_decode_protein logic)."""
import importlib
import pytest
import pandas as pd
from unittest.mock import patch, MagicMock

from scripts.lib.schema import ProteinMeta
from scripts.lib.cis_extract import OUTPUT_COLS

_decode_mod = importlib.import_module("scripts.02_cis_pqtl_extract.decode")
read_decode_protein = _decode_mod.read_decode_protein


@pytest.fixture
def sample_protein():
    return ProteinMeta(
        seqid="10000_28_CRYBB2_CRBB2",
        gene="CRYBB2", uniprot="",
        chrom="22", tss=25_212_564, build="hg38",
        source_cohort="deCODE",
    )


def _fake_rows(n: int = 3) -> list[dict]:
    return [
        {
            "Chrom": "chr22",
            "Pos(hg38)": str(25_212_564 + i * 1000),
            "Name": f"22:25212564+{i}:A:G",
            "rsids": f"rs{100 + i}",
            "effectAllele": "A",
            "otherAllele": "G",
            "Beta": "0.1",
            "Pval": "1e-9",
            "SE": "0.01",
            "N": "35000",
        }
        for i in range(n)
    ]


class TestReadDecodeProtein:
    def test_returns_dataframe_with_expected_cols(self, sample_protein):
        rows = _fake_rows(3)
        eaf_dict = {r["Name"]: 0.30 for r in rows}

        with patch.object(_decode_mod, "_url_map",
                          {sample_protein.seqid: "http://fake"}), \
             patch.object(_decode_mod, "iter_decode_rows",
                          return_value=rows):
            result = read_decode_protein(sample_protein, eaf_dict)

        assert result is not None
        for col in ("chrom", "pos", "EA", "OA", "EAF", "beta", "se", "pval", "N", "rsid"):
            assert col in result.columns, f"missing column: {col}"

    def test_chrom_has_no_chr_prefix(self, sample_protein):
        rows = _fake_rows(1)
        eaf_dict = {r["Name"]: 0.30 for r in rows}

        with patch.object(_decode_mod, "_url_map",
                          {sample_protein.seqid: "http://fake"}), \
             patch.object(_decode_mod, "iter_decode_rows",
                          return_value=rows):
            result = read_decode_protein(sample_protein, eaf_dict)

        assert not result["chrom"].str.startswith("chr").any()

    def test_missing_eaf_rows_dropped_and_logged(self, sample_protein, caplog):
        rows = _fake_rows(3)
        eaf_dict = {}  # none in cache → all dropped

        with patch.object(_decode_mod, "_url_map",
                          {sample_protein.seqid: "http://fake"}), \
             patch.object(_decode_mod, "iter_decode_rows",
                          return_value=rows), \
             caplog.at_level("INFO", logger="02_decode"):
            result = read_decode_protein(sample_protein, eaf_dict)

        assert result is None or (result is not None and len(result) == 0)
        assert "EAF lookup" in caplog.text or result is None

    def test_partial_eaf_cache_some_rows_kept(self, sample_protein):
        rows = _fake_rows(3)
        eaf_dict = {rows[0]["Name"]: 0.30}

        with patch.object(_decode_mod, "_url_map",
                          {sample_protein.seqid: "http://fake"}), \
             patch.object(_decode_mod, "iter_decode_rows",
                          return_value=rows):
            result = read_decode_protein(sample_protein, eaf_dict)

        assert result is not None
        assert len(result) == 1

    def test_missing_url_returns_none(self, sample_protein):
        with patch.object(_decode_mod, "_url_map", {}):
            result = read_decode_protein(sample_protein, {})
        assert result is None

    def test_empty_row_list_returns_none(self, sample_protein):
        with patch.object(_decode_mod, "_url_map",
                          {sample_protein.seqid: "http://fake"}), \
             patch.object(_decode_mod, "iter_decode_rows",
                          return_value=[]):
            result = read_decode_protein(sample_protein, {})
        assert result is None

    def test_n_fills_with_default_when_column_missing(self, sample_protein):
        rows = [
            {
                "Chrom": "chr22", "Pos(hg38)": "25212564",
                "Name": "22:25212564:A:G",
                "rsids": "rs100",
                "effectAllele": "A", "otherAllele": "G",
                "Beta": "0.1", "Pval": "1e-9", "SE": "0.01",
                # No "N" column
            }
        ]
        eaf_dict = {"22:25212564:A:G": 0.30}

        with patch.object(_decode_mod, "_url_map",
                          {sample_protein.seqid: "http://fake"}), \
             patch.object(_decode_mod, "iter_decode_rows",
                          return_value=rows):
            result = read_decode_protein(sample_protein, eaf_dict)

        assert result is not None
        assert result["N"].iloc[0] == 35_000  # default N
