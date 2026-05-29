"""Tests for scripts.02_cis_pqtl_extract.fenland (read_fenland_protein logic)."""
import importlib
import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from pathlib import Path

from scripts.lib.schema import ProteinMeta

_fenland_mod = importlib.import_module("scripts.02_cis_pqtl_extract.fenland")
read_fenland_protein = _fenland_mod.read_fenland_protein


@pytest.fixture
def sample_protein():
    return ProteinMeta(
        seqid="CRYBB2",
        gene="CRYBB2", uniprot="",
        chrom="22", tss=25_212_564, build="hg19",
        source_cohort="Fenland",
    )


def _make_rows(chrom="22", pos=25_212_564, n=2, rsid_col="SNPID") -> list[dict]:
    rows = []
    for i in range(n):
        row = {
            "CHR": chrom,
            "POS": str(pos + i * 100),
            rsid_col: f"rs{200 + i}",
            "EA": "A", "OA": "G",
            "EAF": "0.35",
            "BETA": "0.05",
            "SE": "0.01",
            "P": "1e-10",
            "N": "10708",
        }
        rows.append(row)
    return rows


class TestReadFenlandProtein:
    def test_returns_dataframe_with_standard_columns(self, sample_protein, tmp_path):
        entity_map = {"CRYBB2": [("syn_fake_1", "CRYBB2.txt.gz")]}
        rows = _make_rows()

        with patch.object(_fenland_mod, "stream_fenland_protein",
                          return_value=rows):
            result = read_fenland_protein(
                sample_protein, entity_map,
                cis_start=25_000_000, cis_end=25_700_000,
            )

        assert result is not None
        for col in ("chrom", "pos", "rsid", "EA", "OA", "EAF", "beta", "se", "pval", "N"):
            assert col in result.columns

    def test_snpid_column_renamed_to_rsid(self, sample_protein):
        entity_map = {"CRYBB2": [("syn1", "f.txt.gz")]}
        rows = _make_rows(rsid_col="SNPID")

        with patch.object(_fenland_mod, "stream_fenland_protein",
                          return_value=rows):
            result = read_fenland_protein(sample_protein, entity_map, 25_000_000, 25_700_000)

        assert "rsid" in result.columns

    def test_rsid_column_is_preserved(self, sample_protein):
        entity_map = {"CRYBB2": [("syn1", "f.txt.gz")]}
        rows = _make_rows(rsid_col="rsid")

        with patch.object(_fenland_mod, "stream_fenland_protein",
                          return_value=rows):
            result = read_fenland_protein(sample_protein, entity_map, 25_000_000, 25_700_000)

        assert "rsid" in result.columns

    def test_snp_column_as_rsid_fallback(self, sample_protein):
        entity_map = {"CRYBB2": [("syn1", "f.txt.gz")]}
        rows = _make_rows(rsid_col="SNP")

        with patch.object(_fenland_mod, "stream_fenland_protein",
                          return_value=rows):
            result = read_fenland_protein(sample_protein, entity_map, 25_000_000, 25_700_000)

        assert "rsid" in result.columns

    def test_missing_rsid_column_defaults_to_dot(self, sample_protein):
        entity_map = {"CRYBB2": [("syn1", "f.txt.gz")]}
        rows = [{"CHR": "22", "POS": "25212564", "EA": "A", "OA": "G",
                 "EAF": "0.3", "BETA": "0.1", "SE": "0.01", "P": "1e-9", "N": "10708"}]

        with patch.object(_fenland_mod, "stream_fenland_protein",
                          return_value=rows):
            result = read_fenland_protein(sample_protein, entity_map, 25_000_000, 25_700_000)

        assert "rsid" in result.columns
        assert (result["rsid"] == ".").all()

    def test_chrom_has_no_chr_prefix(self, sample_protein):
        entity_map = {"CRYBB2": [("syn1", "f.txt.gz")]}
        rows = [{"CHR": "chr22", "POS": "25212564", "rsid": "rs1",
                 "EA": "A", "OA": "G", "EAF": "0.3",
                 "BETA": "0.1", "SE": "0.01", "P": "1e-9", "N": "10708"}]

        with patch.object(_fenland_mod, "stream_fenland_protein",
                          return_value=rows):
            result = read_fenland_protein(sample_protein, entity_map, 25_000_000, 25_700_000)

        assert not result["chrom"].str.startswith("chr").any()

    def test_multi_file_aggregation(self, sample_protein):
        entity_map = {"CRYBB2": [("syn1", "f1.txt.gz"), ("syn2", "f2.txt.gz")]}
        rows_per_file = _make_rows(n=2)

        with patch.object(_fenland_mod, "stream_fenland_protein",
                          return_value=rows_per_file):
            result = read_fenland_protein(sample_protein, entity_map, 25_000_000, 25_700_000)

        assert result is not None
        assert len(result) == 4  # 2 rows × 2 files

    def test_no_files_returns_none(self, sample_protein):
        entity_map = {"CRYBB2": []}
        result = read_fenland_protein(sample_protein, entity_map, 25_000_000, 25_700_000)
        assert result is None

    def test_all_empty_files_returns_none(self, sample_protein):
        entity_map = {"CRYBB2": [("syn1", "f.txt.gz")]}

        with patch.object(_fenland_mod, "stream_fenland_protein",
                          return_value=[]):
            result = read_fenland_protein(sample_protein, entity_map, 25_000_000, 25_700_000)

        assert result is None

    def test_n_filled_with_default_when_missing(self, sample_protein):
        entity_map = {"CRYBB2": [("syn1", "f.txt.gz")]}
        rows = [{"CHR": "22", "POS": "25212564", "rsid": "rs1",
                 "EA": "A", "OA": "G", "EAF": "0.3",
                 "BETA": "0.1", "SE": "0.01", "P": "1e-9"}]  # no "N"

        with patch.object(_fenland_mod, "stream_fenland_protein",
                          return_value=rows):
            result = read_fenland_protein(sample_protein, entity_map, 25_000_000, 25_700_000)

        assert result["N"].iloc[0] == 10_708

    def test_source_n_is_preserved_when_present(self, sample_protein):
        entity_map = {"CRYBB2": [("syn1", "f.txt.gz")]}
        rows = _make_rows(n=1)
        rows[0]["N"] = "9999"

        with patch.object(_fenland_mod, "stream_fenland_protein",
                          return_value=rows):
            result = read_fenland_protein(sample_protein, entity_map, 25_000_000, 25_700_000)

        assert result["N"].iloc[0] == 9999


# ── Real-data integration test ────────────────────────────────────────────────
# Downloads one Fenland file from Synapse and exercises the full parsing path.
# Cached in tests/_fenland_cache/ after first run; skipped if Synapse unavailable.

_CACHE_DIR = Path(__file__).parent / "_fenland_cache"
_ABL1_ENTITY   = "syn52965334"   # ABL1_3341_33.txt.gz
_ABL1_FILENAME = "ABL1_3341_33.txt.gz"
_ABL1_GENE     = "ABL1"
_ABL1_CHROM    = "9"
_ABL1_TSS      = 133_589_333
_CIS_WINDOW_KB = 500


def _synapse_available() -> bool:
    try:
        import synapseclient
        syn = synapseclient.Synapse()
        syn.login(silent=True)
        return True
    except Exception:
        return False


@pytest.mark.slow
@pytest.mark.skipif(not _synapse_available(), reason="Synapse credentials not available")
class TestStreamFenlandRealDownload:
    """Parses a real Fenland file to catch format regressions (gzip vs plaintext, column names)."""

    @pytest.fixture(scope="class")
    def cached_file(self):
        import shutil
        _CACHE_DIR.mkdir(exist_ok=True)
        cached = _CACHE_DIR / _ABL1_FILENAME
        if not cached.exists():
            from scripts.lib.synapse_stream import download_entity
            import tempfile
            with tempfile.TemporaryDirectory() as tmp:
                path = download_entity(_ABL1_ENTITY, Path(tmp))
                shutil.copy(path, cached)
        return cached

    def _stream(self, cached_file: Path) -> list[dict]:
        import shutil, tempfile
        from scripts.lib.synapse_stream import stream_fenland_protein
        from scripts.lib.cis import cis_window_bounds
        from unittest.mock import patch
        start, end = cis_window_bounds(_ABL1_TSS, kb=_CIS_WINDOW_KB)
        with tempfile.TemporaryDirectory() as tmp:
            file_copy = Path(tmp) / cached_file.name
            shutil.copy(cached_file, file_copy)
            with patch("scripts.lib.synapse_stream.download_entity", return_value=file_copy):
                return stream_fenland_protein(_ABL1_ENTITY, _ABL1_CHROM, start, end, Path(tmp))

    def test_returns_nonempty_rows_in_cis_window(self, cached_file):
        rows = self._stream(cached_file)
        assert len(rows) > 0, "Expected cis-window rows for ABL1 but got none"

    def test_expected_columns_present(self, cached_file):
        rows = self._stream(cached_file)
        assert rows
        # Actual METAL meta-analysis format used by Fenland
        expected = {"chr", "pos", "rsid", "Allele1", "Allele2", "Freq1",
                    "Effect", "StdErr", "Pvalue", "TotalSampleSize"}
        actual = set(rows[0].keys())
        missing = expected - actual
        assert not missing, f"Columns missing from real Fenland file: {missing}. Got: {actual}"
