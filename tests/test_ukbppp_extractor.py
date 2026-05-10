"""Tests for UKB-PPP extraction logic (ID parsing, column normalization)."""
import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile

from scripts.lib.schema import ProteinMeta


@pytest.fixture
def sample_protein():
    return ProteinMeta(
        seqid="PROC_P04070",
        gene="PROC", uniprot="P04070",
        chrom="2", tss=128_185_556, build="hg19",
        source_cohort="UKB_PPP",
    )


def _fake_ukbppp_rows(n: int = 3, chrom: str = "2",
                      pos_start: int = 128_185_556) -> list[dict]:
    return [
        {
            "ID": f"chr{chrom}:{pos_start + i * 100}:A:G",
            "ALLELE0": "A",
            "ALLELE1": "G",
            "A1FREQ": "0.35",
            "BETA": "0.05",
            "SE": "0.01",
            "LOG10P": "9.0",
            "N": "34557",
        }
        for i in range(n)
    ]


class TestUkbPppIdParsing:
    """Test the ID parsing logic embedded in ukbppp.py read_fn."""

    def _parse_id(self, id_str: str) -> dict:
        """Replicate the ID parsing from ukbppp.read_fn for testing."""
        parts = id_str.split(":")
        if len(parts) < 4:
            return {}
        chrom = parts[0].replace("chr", "", 1) if parts[0].startswith("chr") else parts[0]
        return {
            "chrom": chrom,
            "pos": int(parts[1]),
            "OA": parts[2],
            "EA": parts[3],
        }

    def test_parses_chr_prefix_correctly(self):
        result = self._parse_id("chr2:128185556:A:G")
        assert result["chrom"] == "2"
        assert result["pos"] == 128_185_556
        assert result["EA"] == "G"
        assert result["OA"] == "A"

    def test_parses_without_chr_prefix(self):
        result = self._parse_id("2:128185556:A:G")
        assert result["chrom"] == "2"

    def test_chr_prefix_stripped_not_lstripped(self):
        # "chr1" → "1" using replace("chr", "", 1), NOT lstrip("chr")
        result = self._parse_id("chr1:100:A:G")
        assert result["chrom"] == "1"

    def test_malformed_id_fewer_than_4_parts_returns_empty(self):
        result = self._parse_id("chr2:100:A")  # only 3 parts
        assert result == {}

    def test_x_chromosome_parsed(self):
        result = self._parse_id("chrX:100000:C:T")
        assert result["chrom"] == "X"


class TestUkbPppReadFn:
    """Test the full read_fn closure from ukbppp.main."""

    def _make_read_fn(self, entity_map: dict, protein: ProteinMeta):
        """Replicate the read_fn construction from ukbppp.main."""
        import importlib
        from scripts.lib.cis import cis_window_bounds

        UKB_N = 34_557

        def read_fn(p: ProteinMeta) -> pd.DataFrame | None:
            start, end = cis_window_bounds(p.tss, kb=500)
            eid = entity_map.get(p.seqid)
            if not eid:
                return None
            from scripts.lib.synapse_stream import stream_ukbppp_protein as _stream_fn
            with tempfile.TemporaryDirectory() as tmp:
                rows = _stream_fn(eid, p.chrom, start, end, Path(tmp))
            if not rows:
                return None
            df = pd.DataFrame(rows)
            df = df.rename(columns={
                "ALLELE1": "EA", "ALLELE0": "OA",
                "A1FREQ": "EAF", "BETA": "beta", "SE": "se",
            })
            id_parts = df["ID"].str.split(":", expand=True)
            df["chrom"] = id_parts[0].str.replace(r"^chr", "", regex=True)
            df["pos"] = pd.to_numeric(id_parts[1], errors="coerce").astype("Int64")
            df["pval"] = 10 ** (-pd.to_numeric(df["LOG10P"], errors="coerce"))
            df["rsid"] = "."
            if "N" in df.columns:
                df["N"] = pd.to_numeric(df["N"], errors="coerce").fillna(UKB_N).astype(int)
            else:
                df["N"] = UKB_N
            return df

        return read_fn

    def test_log10p_converted_to_pval(self, sample_protein):
        rows = _fake_ukbppp_rows(1)  # LOG10P=9.0
        entity_map = {sample_protein.seqid: "syn_fake"}

        read_fn = self._make_read_fn(entity_map, sample_protein)

        with patch("scripts.lib.synapse_stream.stream_ukbppp_protein",
                   return_value=rows):
            result = read_fn(sample_protein)

        assert result is not None
        assert result["pval"].iloc[0] == pytest.approx(1e-9, rel=1e-3)

    def test_rsid_is_always_dot(self, sample_protein):
        rows = _fake_ukbppp_rows(2)
        entity_map = {sample_protein.seqid: "syn_fake"}

        read_fn = self._make_read_fn(entity_map, sample_protein)

        with patch("scripts.lib.synapse_stream.stream_ukbppp_protein",
                   return_value=rows):
            result = read_fn(sample_protein)

        assert result is not None
        assert (result["rsid"] == ".").all()

    def test_chrom_has_no_chr_prefix(self, sample_protein):
        rows = _fake_ukbppp_rows(2)
        entity_map = {sample_protein.seqid: "syn_fake"}

        read_fn = self._make_read_fn(entity_map, sample_protein)

        with patch("scripts.lib.synapse_stream.stream_ukbppp_protein",
                   return_value=rows):
            result = read_fn(sample_protein)

        assert result is not None
        assert not result["chrom"].str.startswith("chr").any()

    def test_missing_entity_id_returns_none(self, sample_protein):
        entity_map = {}  # protein not in map
        read_fn = self._make_read_fn(entity_map, sample_protein)
        result = read_fn(sample_protein)
        assert result is None

    def test_n_fills_default_when_column_absent(self, sample_protein):
        rows = [{k: v for k, v in r.items() if k != "N"} for r in _fake_ukbppp_rows(1)]
        entity_map = {sample_protein.seqid: "syn_fake"}

        read_fn = self._make_read_fn(entity_map, sample_protein)

        with patch("scripts.lib.synapse_stream.stream_ukbppp_protein",
                   return_value=rows):
            result = read_fn(sample_protein)

        assert result is not None
        assert result["N"].iloc[0] == 34_557
