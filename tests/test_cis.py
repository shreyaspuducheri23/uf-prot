"""Tests for code.lib.cis"""
import pytest
import pandas as pd
from unittest.mock import patch

from scripts.lib.cis import cis_window_bounds, load_aric_tss, tss_from_ensembl
from scripts.lib.paths import ARIC_SEQID


class TestCisWindowBounds:
    def test_symmetric_window(self):
        start, end = cis_window_bounds(1_000_000, kb=500)
        assert start == 500_000
        assert end == 1_500_000

    def test_clamped_at_zero(self):
        start, end = cis_window_bounds(100, kb=500)
        assert start == 1  # clamped to 1
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
        for seqid, (chrom, tss, uniprot, gene) in tss_map.items():
            assert isinstance(tss, int), f"{seqid} has non-int TSS"

    def test_chromosomes_are_strings(self, tss_map):
        for seqid, (chrom, tss, uniprot, gene) in tss_map.items():
            assert isinstance(chrom, str), f"{seqid} chrom is not str"

    def test_invalid_tss_rows_logged_not_silently_skipped(self, tmp_path, caplog):
        """Proteins with invalid TSS should produce a WARNING log."""
        import io
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


class TestTssFromEnsemblStrand:
    """Test that strand handling in tss_from_ensembl is correct post-fix."""

    def _mock_response(self, strand):
        return {
            "seq_region_name": "22",
            "strand": strand,
            "start": 1000,
            "end": 2000,
        }

    def test_strand_plus_1_uses_start(self):
        with patch("scripts.lib.cis.requests.get") as mock_get:
            mock_get.return_value.raise_for_status = lambda: None
            mock_get.return_value.json.return_value = self._mock_response(1)
            tss_from_ensembl.cache_clear()
            result = tss_from_ensembl("TESTGENE", "hg38")
        assert result is not None
        _, tss = result
        assert tss == 1000  # start

    def test_strand_minus_1_uses_end(self):
        with patch("scripts.lib.cis.requests.get") as mock_get:
            mock_get.return_value.raise_for_status = lambda: None
            mock_get.return_value.json.return_value = self._mock_response(-1)
            tss_from_ensembl.cache_clear()
            result = tss_from_ensembl("TESTGENE2", "hg38")
        assert result is not None
        _, tss = result
        assert tss == 2000  # end

    def test_strand_plus_string_uses_start(self):
        with patch("scripts.lib.cis.requests.get") as mock_get:
            mock_get.return_value.raise_for_status = lambda: None
            mock_get.return_value.json.return_value = self._mock_response("+")
            tss_from_ensembl.cache_clear()
            result = tss_from_ensembl("TESTGENE3", "hg38")
        assert result is not None
        _, tss = result
        assert tss == 1000

    def test_unexpected_strand_raises(self):
        with patch("scripts.lib.cis.requests.get") as mock_get:
            mock_get.return_value.raise_for_status = lambda: None
            mock_get.return_value.json.return_value = self._mock_response(0)
            tss_from_ensembl.cache_clear()
            result = tss_from_ensembl("TESTGENE4", "hg38")
        # Unexpected strand is caught by the outer try/except → returns None
        assert result is None
