"""Tests for UKB-PPP extraction logic bound to production transformations."""
import importlib
import tempfile
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


# ── Real-data integration test ────────────────────────────────────────────────
# Downloads one UKB-PPP tar from Synapse and verifies the full parsing path.
# Cached in tests/_ukbppp_cache/ after first run; skipped if Synapse login fails.

_CACHE_DIR = Path(__file__).parent / "_ukbppp_cache"
# A1BG (chr19:58864865) — first entity in syn51365303
_A1BG_ENTITY = "syn52363617"
_A1BG_CHROM = "19"
_A1BG_TSS = 58_864_865
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
class TestStreamUkbpppRealDownload:
    """Parses an actual UKB-PPP tar to catch format regressions early."""

    @pytest.fixture(scope="class")
    def cached_tar(self):
        _CACHE_DIR.mkdir(exist_ok=True)
        cached = _CACHE_DIR / "A1BG_P04217_OID30771_v1_Inflammation_II.tar"
        if not cached.exists():
            from scripts.lib.synapse_stream import download_entity
            import shutil
            with tempfile.TemporaryDirectory() as tmp:
                path = download_entity(_A1BG_ENTITY, Path(tmp))
                shutil.copy(path, cached)
        return cached

    def _stream(self, cached_tar: Path) -> list[dict]:
        """Copy cached tar into a temp dir (stream_ukbppp_protein deletes after use)."""
        import shutil
        from scripts.lib.synapse_stream import stream_ukbppp_protein
        from scripts.lib.cis import cis_window_bounds
        from unittest.mock import patch
        start, end = cis_window_bounds(_A1BG_TSS, kb=_CIS_WINDOW_KB)
        with tempfile.TemporaryDirectory() as tmp:
            tar_copy = Path(tmp) / cached_tar.name
            shutil.copy(cached_tar, tar_copy)
            with patch("scripts.lib.synapse_stream.download_entity", return_value=tar_copy):
                return stream_ukbppp_protein("syn_fake", _A1BG_CHROM, start, end, Path(tmp))

    def test_returns_nonempty_rows_in_cis_window(self, cached_tar):
        rows = self._stream(cached_tar)
        assert len(rows) > 0, "Expected cis-window rows for A1BG but got none"

    def test_id_field_parseable_as_chrom_pos(self, cached_tar):
        from scripts.lib.cis import cis_window_bounds
        start, end = cis_window_bounds(_A1BG_TSS, kb=_CIS_WINDOW_KB)
        rows = self._stream(cached_tar)
        for row in rows[:5]:
            parts = row["ID"].split(":")
            assert len(parts) >= 2, f"Unparseable ID: {row['ID']}"
            chrom = parts[0].lstrip("chr")
            assert chrom == _A1BG_CHROM
            pos = int(parts[1])
            assert start <= pos <= end
