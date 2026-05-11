"""Tests for UKB-PPP extraction logic bound to production transformations."""
import importlib
from pathlib import Path

import pandas as pd
import pytest

from scripts.lib.cis import cis_window_bounds
from scripts.lib.schema import ProteinMeta

_ukb_mod = importlib.import_module("scripts.02_cis_pqtl_extract.ukbppp")


@pytest.fixture
def sample_protein():
    return ProteinMeta(
        seqid="PROC_P04070",
        gene="PROC",
        uniprot="P04070",
        chrom="2",
        tss=128_185_556,
        build="hg19",
        source_cohort="UKB_PPP",
    )


def _fake_ukbppp_rows(n: int = 3, chrom: str = "2", pos_start: int = 128_185_556) -> list[dict]:
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


class TestUkbPppProductionNormalization:
    def test_normalize_rows_parses_id_and_chr_prefix(self):
        rows = _fake_ukbppp_rows(1)
        out = _ukb_mod.normalize_ukbppp_rows(rows)

        assert out is not None
        assert out["chrom"].iloc[0] == "2"
        assert out["pos"].iloc[0] == 128_185_556
        assert out["EA"].iloc[0] == "G"
        assert out["OA"].iloc[0] == "A"

    def test_normalize_rows_converts_log10p_and_applies_rsid_policy(self):
        rows = _fake_ukbppp_rows(1)
        out = _ukb_mod.normalize_ukbppp_rows(rows)

        assert out is not None
        assert out["pval"].iloc[0] == pytest.approx(1e-9, rel=1e-3)
        assert (out["rsid"] == ".").all()

    def test_normalize_rows_n_fallback(self):
        rows = _fake_ukbppp_rows(2)
        rows[0]["N"] = ""
        rows[1].pop("N")

        out = _ukb_mod.normalize_ukbppp_rows(rows)

        assert out is not None
        assert out["N"].iloc[0] == _ukb_mod.UKB_N
        assert out["N"].iloc[1] == _ukb_mod.UKB_N


class TestUkbPppProductionReadFn:
    def test_build_read_fn_uses_streamed_rows_and_window(self, sample_protein):
        calls = []
        expected_start, expected_end = cis_window_bounds(sample_protein.tss, kb=500)

        def fake_stream(entity_id: str, chrom: str, start: int, end: int, tmp: Path):
            calls.append((entity_id, chrom, start, end, tmp))
            return _fake_ukbppp_rows(1, chrom=chrom, pos_start=start)

        read_fn = _ukb_mod.build_read_fn(
            entity_map={sample_protein.seqid: "syn_fake"},
            window_kb=500,
            stream_fn=fake_stream,
        )

        out = read_fn(sample_protein)

        assert out is not None
        assert len(calls) == 1
        eid, chrom, start, end, tmp = calls[0]
        assert eid == "syn_fake"
        assert chrom == sample_protein.chrom
        assert start == expected_start
        assert end == expected_end
        assert tmp.name.startswith("ukb_PROC_P04070_")

    def test_build_read_fn_missing_entity_returns_none(self, sample_protein):
        read_fn = _ukb_mod.build_read_fn(entity_map={}, window_kb=500)
        assert read_fn(sample_protein) is None

    def test_build_read_fn_empty_stream_returns_none(self, sample_protein):
        def fake_stream(_entity_id: str, _chrom: str, _start: int, _end: int, _tmp: Path):
            return []

        read_fn = _ukb_mod.build_read_fn(
            entity_map={sample_protein.seqid: "syn_fake"},
            window_kb=500,
            stream_fn=fake_stream,
        )

        assert read_fn(sample_protein) is None
