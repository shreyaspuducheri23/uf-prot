"""Tests for scripts.lib.cis."""
from unittest.mock import patch

import pandas as pd
import pytest

from scripts.lib.cis import (
    _load_tss_cache,
    _save_tss_cache,
    cis_window_bounds,
    load_aric_tss,
    resolve_tss,
)
from scripts.lib.paths import ARIC_SEQID


class TestCisWindowBounds:
    def test_symmetric_window(self):
        start, end = cis_window_bounds(1_000_000, kb=500)
        assert start == 500_000
        assert end == 1_500_000

    def test_clamped_at_zero(self):
        start, end = cis_window_bounds(100, kb=500)
        assert start == 1
        assert end == 500_100

    def test_1mb_window(self):
        start, end = cis_window_bounds(10_000_000, kb=1000)
        assert end - start == 2_000_000


class TestLoadAricTss:
    @pytest.fixture(scope="class")
    def tss_map(self):
        return load_aric_tss(ARIC_SEQID)

    def test_loads_correct_count(self, tss_map):
        assert len(tss_map) == 4657

    def test_known_entry(self, tss_map):
        assert "SeqId_10000_28" in tss_map
        chrom, tss, uniprot, gene = tss_map["SeqId_10000_28"]
        assert chrom == "22"
        assert tss == 25_212_564
        assert uniprot == "P43320"
        assert gene == "CRYBB2"

    def test_all_have_int_tss(self, tss_map):
        for seqid, (_chrom, tss, _uniprot, _gene) in tss_map.items():
            assert isinstance(tss, int), f"{seqid} has non-int TSS"

    def test_chromosomes_are_strings(self, tss_map):
        for seqid, (chrom, _tss, _uniprot, _gene) in tss_map.items():
            assert isinstance(chrom, str), f"{seqid} chrom is not str"

    def test_invalid_tss_rows_logged_not_silently_skipped(self, tmp_path, caplog):
        tsv_content = (
            "seqid_in_sample\tuniprot_id\tentrezgenesymbol\tchromosome_name\ttranscription_start_site\n"
            "SeqId_GOOD\tQ12345\tGENE1\t22\t25212564\n"
            "SeqId_BAD\tQ99999\tGENE2\t1\tNOT_AN_INT\n"
        )
        path = tmp_path / "seqid.txt"
        path.write_text(tsv_content)
        with caplog.at_level("WARNING", logger="scripts.lib.cis"):
            result = load_aric_tss(path)
        assert "SeqId_GOOD" in result
        assert "SeqId_BAD" not in result
        assert "invalid TSS" in caplog.text.lower() or "skipping" in caplog.text.lower()


class MockResponse:
    def __init__(self, payload=None, status_ok=True):
        self.payload = payload or {}
        self.status_ok = status_ok

    def raise_for_status(self):
        if not self.status_ok:
            raise RuntimeError("404")

    def json(self):
        return self.payload


class TestResolveTss:
    def setup_method(self):
        resolve_tss.cache_clear()

    def _ensembl_payload(self, strand):
        return {
            "seq_region_name": "22",
            "strand": strand,
            "start": 1000,
            "end": 2000,
        }

    def _mock_get(self, ensembl=None, hgnc=None):
        ensembl = ensembl or {}
        hgnc = hgnc or {}

        def side_effect(url, headers=None, timeout=None):
            if "genenames.org" in url:
                if "/fetch/symbol/" in url:
                    symbol = url.rstrip("/").split("/")[-1]
                    docs = hgnc.get(("fetch_symbol", symbol), [])
                    return MockResponse({"response": {"docs": docs}})
                field, symbol = url.rstrip("/").split("/")[-2:]
                docs = [{"symbol": s} for s in hgnc.get((field, symbol), [])]
                return MockResponse({"response": {"docs": docs}})

            symbol = url.split("/lookup/symbol/homo_sapiens/")[1].split("?")[0]
            if symbol not in ensembl:
                return MockResponse(status_ok=False)
            payload = ensembl[symbol]
            if payload is None:
                return MockResponse(status_ok=False)
            return MockResponse(payload)

        return side_effect

    def test_tier_1_strand_plus_1_uses_start(self):
        with patch("scripts.lib.cis.requests.get",
                   side_effect=self._mock_get({"TEST": self._ensembl_payload(1)})):
            result = resolve_tss("TEST", "hg38")
        assert result.resolved
        assert result.tier == 1
        assert result.tss == 1000

    def test_tier_1_strand_minus_1_uses_end(self):
        with patch("scripts.lib.cis.requests.get",
                   side_effect=self._mock_get({"TEST": self._ensembl_payload(-1)})):
            result = resolve_tss("TEST", "hg38")
        assert result.resolved
        assert result.tss == 2000

    def test_tier_1_strand_plus_string_uses_start(self):
        with patch("scripts.lib.cis.requests.get",
                   side_effect=self._mock_get({"TEST": self._ensembl_payload("+")})):
            result = resolve_tss("TEST", "hg38")
        assert result.resolved
        assert result.tss == 1000

    def test_tier_1_bad_strand_falls_through_unresolved(self):
        with patch("scripts.lib.cis.requests.get",
                   side_effect=self._mock_get({"BAD": self._ensembl_payload(0)})):
            result = resolve_tss("BAD", "hg38")
        assert not result.resolved
        assert result.tier == 0

    def test_tier_2_prev_symbol_rescue(self):
        with patch(
            "scripts.lib.cis.requests.get",
            side_effect=self._mock_get(
                {"NEW": self._ensembl_payload(1)},
                {("prev_symbol", "OLD"): ["NEW"]},
            ),
        ):
            result = resolve_tss("OLD", "hg38")
        assert result.resolved
        assert result.tier == 2
        assert result.resolved_symbol == "NEW"
        assert result.requested_symbol == "OLD"
        assert "NEW" in result.attempts

    def test_current_symbol_prev_symbol_rescue_for_legacy_build(self):
        with patch(
            "scripts.lib.cis.requests.get",
            side_effect=self._mock_get(
                {"OLD": self._ensembl_payload(-1)},
                {("fetch_symbol", "NEW"): [{"symbol": "NEW", "prev_symbol": ["OLD"]}]},
            ),
        ):
            result = resolve_tss("NEW", "hg19")
        assert result.resolved
        assert result.tier == 2
        assert result.resolved_symbol == "OLD"
        assert result.tss == 2000

    def test_tier_3_alias_symbol_rescue(self):
        with patch(
            "scripts.lib.cis.requests.get",
            side_effect=self._mock_get(
                {"APPR": self._ensembl_payload(1)},
                {("alias_symbol", "OLD"): ["APPR"]},
            ),
        ):
            result = resolve_tss("OLD", "hg38")
        assert result.resolved
        assert result.tier == 3
        assert result.resolved_symbol == "APPR"

    def test_tier_3_approved_symbol_prev_rescue(self):
        with patch(
            "scripts.lib.cis.requests.get",
            side_effect=self._mock_get(
                {"PREV": self._ensembl_payload(-1)},
                {
                    ("alias_symbol", "OLD"): ["APPR"],
                    ("prev_symbol", "APPR"): ["PREV"],
                },
            ),
        ):
            result = resolve_tss("OLD", "hg38")
        assert result.resolved
        assert result.tier == 3
        assert result.resolved_symbol == "PREV"
        assert result.tss == 2000

    def test_all_tiers_fail_records_attempts(self):
        with patch("scripts.lib.cis.requests.get", side_effect=self._mock_get()):
            result = resolve_tss("MISSING", "hg38")
        assert not result.resolved
        assert result.tier == 0
        assert result.attempts == ("MISSING",)

    def test_override_hg38(self):
        with patch("scripts.lib.cis.requests.get", side_effect=self._mock_get()):
            result = resolve_tss("ALPPL2", "hg38")
        assert result.resolved
        assert result.tier == 4
        assert result.chrom == "2"
        assert result.tss == 232_378_751
        assert "tss_overrides.tsv" in result.source

    def test_override_hg19_guard(self):
        with patch("scripts.lib.cis.requests.get", side_effect=self._mock_get()):
            result = resolve_tss("ALPPL2", "hg19")
        assert not result.resolved


class TestTssCacheRoundTrip:
    def test_legacy_three_column_cache_loads(self, tmp_path):
        cache_path = tmp_path / "_tss.tsv"
        pd.DataFrame([{"gene": "A", "chrom": "1", "tss": 123}]).to_csv(
            cache_path, sep="\t", index=False
        )
        assert _load_tss_cache(cache_path) == {"A": ("1", 123)}

    def test_save_round_trips_provenance_and_integer_null_tier(self, tmp_path):
        cache_path = tmp_path / "_tss.tsv"
        cache = {"A": ("1", 123), "B": ("2", 456)}
        rows = [{
            "gene": "A",
            "chrom": "1",
            "tss": 123,
            "resolved_symbol": "A",
            "tier": 4,
            "source": "manual",
        }]

        _save_tss_cache(cache_path, cache, rows)
        out = pd.read_csv(cache_path, sep="\t", dtype=str)

        assert out.columns.tolist() == ["gene", "chrom", "tss", "resolved_symbol", "tier", "source"]
        assert out.loc[out["gene"] == "A", "tier"].iloc[0] == "4"
        assert out.loc[out["gene"] == "B", "tier"].iloc[0] != "4.0"
