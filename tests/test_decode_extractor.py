"""Tests for scripts.02_cis_pqtl_extract.decode (read_decode_protein logic)."""
import importlib
import json
import gzip
import io
import pytest
import pandas as pd
from unittest.mock import patch, MagicMock, call

from scripts.lib.schema import ProteinMeta
from scripts.lib.cis_extract import OUTPUT_COLS

_decode_mod = importlib.import_module("scripts.02_cis_pqtl_extract.decode")
read_decode_protein = _decode_mod.read_decode_protein

# Columns the rename dict in read_decode_protein maps FROM (must exist in real files).
_EXPECTED_SOURCE_COLS = {"Chrom", "Pos", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval", "SE", "ImpMAF"}


@pytest.fixture
def sample_protein():
    return ProteinMeta(
        seqid="10000_28_CRYBB2_CRBB2",
        gene="CRYBB2", uniprot="",
        chrom="22", tss=25_212_564, build="hg38",
        source_cohort="deCODE",
    )


def _fake_df(n: int = 3, impmaf: str = "0.12") -> pd.DataFrame:
    rows = [
        {
            "Chrom": "chr22",
            "Pos": str(25_212_564 + i * 1000),
            "Name": f"22:25212564+{i}:A:G",
            "rsids": f"rs{100 + i}",
            "effectAllele": "A",
            "otherAllele": "G",
            "Beta": "0.1",
            "Pval": "1e-9",
            "SE": "0.01",
            "N": "35000",
            "ImpMAF": impmaf,
        }
        for i in range(n)
    ]
    return pd.DataFrame(rows)


class TestBuildProteinList:
    def test_uses_protein_names_from_s3_key_index(self, tmp_path):
        cohort = "deCODE_test"
        cohort_path = tmp_path / cohort
        cohort_path.mkdir()
        (cohort_path / "_tss_hg38.tsv").write_text(
            "gene\tchrom\ttss\tresolved_symbol\ttier\tsource\n"
            "CRYBB2\t22\t25212564\tCRYBB2\t1\tcache\n"
            "RAF1\t3\t12645699\tRAF1\t1\tcache\n"
        )
        protein_names = [
            "10000_28_CRYBB2_CRBB2",
            "10001_7_RAF1_c_Raf",
            "malformed",
        ]

        with patch.object(_decode_mod, "COHORT", cohort), \
             patch.object(_decode_mod, "cohort_dir", lambda c: tmp_path / c), \
             patch.object(_decode_mod, "resolve_tss") as mock_resolve:
            proteins = _decode_mod.build_protein_list(protein_names, build="hg38")

        assert [p.seqid for p in proteins] == [
            "10000_28_CRYBB2_CRBB2",
            "10001_7_RAF1_c_Raf",
        ]
        assert [p.gene for p in proteins] == ["CRYBB2", "RAF1"]
        assert [p.chrom for p in proteins] == ["22", "3"]
        mock_resolve.assert_not_called()


class TestReadDecodeProtein:
    def _call(self, protein, df, n_default=_decode_mod._DEFAULT_N):
        with patch.object(_decode_mod, "_s3_key_map", {protein.seqid: "fake/key"}), \
             patch.object(_decode_mod, "_load_decode_raw_df", return_value=df):
            return read_decode_protein(protein, n_default=n_default)

    def test_returns_dataframe_with_expected_cols(self, sample_protein):
        result = self._call(sample_protein, _fake_df(3))
        assert result is not None
        for col in ("chrom", "pos", "EA", "OA", "EAF", "beta", "se", "pval", "N", "rsid"):
            assert col in result.columns, f"missing column: {col}"

    def test_eaf_comes_from_impmaf(self, sample_protein):
        result = self._call(sample_protein, _fake_df(1, impmaf="0.23"))
        assert result is not None
        assert abs(result["EAF"].iloc[0] - 0.23) < 1e-9

    def test_eaf_is_numeric(self, sample_protein):
        result = self._call(sample_protein, _fake_df(3))
        assert result is not None
        assert pd.api.types.is_float_dtype(result["EAF"])

    def test_row_with_nan_impmaf_dropped(self, sample_protein):
        df = _fake_df(3)
        df.loc[1, "ImpMAF"] = None  # middle row has no ImpMAF
        result = self._call(sample_protein, df)
        assert result is not None
        assert len(result) == 2

    def test_all_nan_impmaf_returns_none(self, sample_protein):
        df = _fake_df(3, impmaf="not_a_number")
        result = self._call(sample_protein, df)
        assert result is None

    def test_chrom_has_no_chr_prefix(self, sample_protein):
        result = self._call(sample_protein, _fake_df(1))
        assert not result["chrom"].str.startswith("chr").any()

    def test_missing_key_in_s3_map_returns_none(self, sample_protein):
        with patch.object(_decode_mod, "_s3_key_map", {}):
            result = read_decode_protein(sample_protein)
        assert result is None

    def test_none_from_load_returns_none(self, sample_protein):
        with patch.object(_decode_mod, "_s3_key_map", {sample_protein.seqid: "fake/key"}), \
             patch.object(_decode_mod, "_load_decode_raw_df", return_value=None):
            result = read_decode_protein(sample_protein)
        assert result is None

    def test_n_fills_with_default_when_column_missing(self, sample_protein):
        rows = [{
            "Chrom": "chr22", "Pos": "25212564",
            "Name": "22:25212564:A:G", "rsids": "rs100",
            "effectAllele": "A", "otherAllele": "G",
            "Beta": "0.1", "Pval": "1e-9", "SE": "0.01",
            "ImpMAF": "0.15",
            # No "N" column
        }]
        result = self._call(sample_protein, pd.DataFrame(rows))
        assert result is not None
        assert result["N"].iloc[0] == 35_559

    def test_source_n_is_preserved_when_present(self, sample_protein):
        df = _fake_df(1)
        df.loc[0, "N"] = "35938"
        result = self._call(sample_protein, df)
        assert result is not None
        assert result["N"].iloc[0] == 35_938


# ---------------------------------------------------------------------------
# _build_s3_key_index
# ---------------------------------------------------------------------------

class TestBuildS3KeyIndex:
    """Tests for _build_s3_key_index — the one-time S3 listing that maps
    core_name → full S3 key, with a JSON disk cache."""

    def _invoke(self, prefix, pattern, tmp_path, mock_s3_index=None):
        """Call _build_s3_key_index with COHORT pointing at tmp_path."""
        cohort = "deCODE_test"
        (tmp_path / cohort).mkdir(parents=True, exist_ok=True)
        with patch.object(_decode_mod, "COHORT", cohort), \
             patch.object(_decode_mod, "cohort_dir", lambda c: tmp_path / c), \
             patch.object(_decode_mod, "_get_s3_client", return_value=mock_s3_index):
            return _decode_mod._build_s3_key_index(prefix, pattern)

    def _make_paginator(self, keys: list[str]) -> MagicMock:
        """Mock paginator that returns all keys in a single page."""
        page = {"Contents": [{"Key": k} for k in keys]}
        pager = MagicMock()
        pager.paginate.return_value = [page]
        s3 = MagicMock()
        s3.get_paginator.return_value = pager
        return s3

    def test_loads_from_cache_without_hitting_s3(self, tmp_path):
        cache_dir = tmp_path / "deCODE_test"
        cache_dir.mkdir()
        cached = {"10000_28_CRYBB2_CRBB2": "final_somascan_raw/Proteomics_PC0_10000_28_CRYBB2_CRBB2_07082019.txt.gz"}
        (cache_dir / "_s3_key_index.json").write_text(json.dumps(cached))

        mock_s3 = MagicMock()
        result = self._invoke("final_somascan_raw", r"Proteomics_PC0_(.+)_\d{8}\.txt\.gz", tmp_path, mock_s3)

        assert result == cached
        mock_s3.get_paginator.assert_not_called()

    def test_builds_index_from_s3_listing(self, tmp_path):
        keys = [
            "final_somascan_raw/Proteomics_PC0_10000_28_CRYBB2_CRBB2_07082019.txt.gz",
            "final_somascan_raw/Proteomics_PC0_10001_7_RAF1_c_Raf_07082019.txt.gz",
        ]
        s3 = self._make_paginator(keys)
        result = self._invoke("final_somascan_raw", r"Proteomics_PC0_(.+)_\d{8}\.txt\.gz", tmp_path, s3)

        assert result["10000_28_CRYBB2_CRBB2"] == keys[0]
        assert result["10001_7_RAF1_c_Raf"] == keys[1]

    def test_non_matching_keys_excluded(self, tmp_path):
        keys = [
            "final_somascan_raw/Proteomics_PC0_GENE_07082019.txt.gz",
            "final_somascan_raw/README.txt",
            "final_somascan_raw/Proteomics_SMP_GENE_07082019.txt.gz",  # wrong prefix
        ]
        s3 = self._make_paginator(keys)
        result = self._invoke("final_somascan_raw", r"Proteomics_PC0_(.+)_\d{8}\.txt\.gz", tmp_path, s3)

        assert "GENE" in result
        assert len(result) == 1  # README and SMP file excluded

    def test_cache_written_after_s3_build(self, tmp_path):
        keys = ["final_somascan_raw/Proteomics_PC0_MYGENE_20230101.txt.gz"]
        s3 = self._make_paginator(keys)
        self._invoke("final_somascan_raw", r"Proteomics_PC0_(.+)_\d{8}\.txt\.gz", tmp_path, s3)

        cache_file = tmp_path / "deCODE_test" / "_s3_key_index.json"
        assert cache_file.exists()
        on_disk = json.loads(cache_file.read_text())
        assert "MYGENE" in on_disk

    def test_empty_bucket_prefix_returns_empty_dict(self, tmp_path):
        s3 = self._make_paginator([])
        result = self._invoke("final_somascan_raw", r"Proteomics_PC0_(.+)_\d{8}\.txt\.gz", tmp_path, s3)
        assert result == {}

    def test_paginator_called_with_correct_bucket_and_prefix(self, tmp_path):
        s3 = self._make_paginator([])
        self._invoke("my_prefix", r"Proteomics_PC0_(.+)_\d{8}\.txt\.gz", tmp_path, s3)
        s3.get_paginator.assert_called_once_with("list_objects_v2")
        s3.get_paginator.return_value.paginate.assert_called_once_with(
            Bucket=_decode_mod._S3_BUCKET, Prefix="my_prefix/"
        )

    def test_smp_pattern_extracts_correct_core(self, tmp_path):
        keys = ["final_somascan_smp/Proteomics_SMP_PC0_10000_28_CRYBB2_CRBB2_07082019.txt.gz"]
        s3 = self._make_paginator(keys)
        result = self._invoke("final_somascan_smp", r"Proteomics_SMP_PC0_(.+)_\d{8}\.txt\.gz", tmp_path, s3)
        assert "10000_28_CRYBB2_CRBB2" in result

    def test_multiple_pages_all_indexed(self, tmp_path):
        page1 = {"Contents": [{"Key": "pfx/Proteomics_PC0_GENE_A_20230101.txt.gz"}]}
        page2 = {"Contents": [{"Key": "pfx/Proteomics_PC0_GENE_B_20230101.txt.gz"}]}
        pager = MagicMock()
        pager.paginate.return_value = [page1, page2]
        s3 = MagicMock()
        s3.get_paginator.return_value = pager
        result = self._invoke("pfx", r"Proteomics_PC0_(.+)_\d{8}\.txt\.gz", tmp_path, s3)
        assert "GENE_A" in result and "GENE_B" in result


# ---------------------------------------------------------------------------
# _get_s3_client caching
# ---------------------------------------------------------------------------

class TestGetS3Client:
    def test_returns_boto3_client(self):
        import scripts.lib.decode_stream as ds
        original = ds._S3_CLIENT
        try:
            ds._S3_CLIENT = None
            with patch("boto3.client") as mock_boto:
                mock_boto.return_value = MagicMock()
                from scripts.lib.decode_stream import _get_s3_client
                client = _get_s3_client("https://ep", "key", "secret")
                mock_boto.assert_called_once_with(
                    "s3",
                    endpoint_url="https://ep",
                    aws_access_key_id="key",
                    aws_secret_access_key="secret",
                    region_name="us-east-1",
                )
                assert client is mock_boto.return_value
        finally:
            ds._S3_CLIENT = original

    def test_second_call_returns_cached_instance(self):
        import scripts.lib.decode_stream as ds
        original = ds._S3_CLIENT
        try:
            ds._S3_CLIENT = None
            with patch("boto3.client") as mock_boto:
                mock_boto.return_value = MagicMock()
                from scripts.lib.decode_stream import _get_s3_client
                c1 = _get_s3_client("https://ep", "key", "secret")
                c2 = _get_s3_client("https://ep", "key", "secret")
                assert c1 is c2
                assert mock_boto.call_count == 1
        finally:
            ds._S3_CLIENT = original


@pytest.mark.network
def test_decode_file_columns_match_rename_dict(real_s3):
    """Contract test: stream the header from a real deCODE S3 file and assert
    every column the rename dict depends on is actually present."""
    body = real_s3.get_object(Bucket=_S3_BUCKET, Key=_S3_RAW_KEY)["Body"]
    with gzip.open(io.BufferedReader(body), "rt") as fh:
        actual_cols = set(fh.readline().strip().split("\t"))

    missing = _EXPECTED_SOURCE_COLS - actual_cols
    assert not missing, (
        f"deCODE file format changed — columns missing from real file: {missing}. "
        f"Actual columns: {actual_cols}"
    )


# ---------------------------------------------------------------------------
# Real S3 integration tests
# ---------------------------------------------------------------------------

_S3_ENDPOINT   = _decode_mod._S3_ENDPOINT
_S3_BUCKET     = _decode_mod._S3_BUCKET
_S3_ACCESS_KEY = _decode_mod._S3_ACCESS_KEY
_S3_SECRET_KEY = _decode_mod._S3_SECRET_KEY

# Well-known protein on chr1 (near start of file → fast early abort)
_S3_TEST_PROTEIN  = "10015_119_KCNAB2_KCAB2"
_S3_TEST_CHROM    = "1"
_S3_TEST_TSS      = 5_990_927
_S3_TEST_WINDOW   = 500_000
_S3_RAW_KEY       = "final_somascan_raw/Proteomics_PC0_10015_119_KCNAB2_KCAB2_07082019.txt.gz"
_S3_SMP_KEY       = "final_somascan_smp/Proteomics_SMP_PC0_10015_119_KCNAB2_KCAB2_10032022.txt.gz"
_STREAM_USECOLS   = frozenset(["Chrom", "Pos", "Name", "Beta", "Pval", "SE", "N"])


@pytest.fixture(scope="module")
def real_s3():
    from scripts.lib.decode_stream import _get_s3_client
    import scripts.lib.decode_stream as ds
    orig = ds._S3_CLIENT
    ds._S3_CLIENT = None
    client = _get_s3_client(_S3_ENDPOINT, _S3_ACCESS_KEY, _S3_SECRET_KEY)
    yield client
    ds._S3_CLIENT = orig


@pytest.mark.network
class TestStreamS3CisRowsRealData:
    """Integration tests against the real deCODE S3 bucket."""

    def test_raw_returns_nonempty_rows(self, real_s3):
        from scripts.lib.decode_stream import stream_s3_cis_rows
        rows = list(stream_s3_cis_rows(
            real_s3, _S3_BUCKET, _S3_RAW_KEY,
            _S3_TEST_CHROM, _S3_TEST_TSS, _S3_TEST_WINDOW, _STREAM_USECOLS,
        ))
        assert len(rows) > 100, f"Expected >100 cis rows, got {len(rows)}"

    def test_raw_all_rows_on_target_chrom(self, real_s3):
        from scripts.lib.decode_stream import stream_s3_cis_rows
        rows = list(stream_s3_cis_rows(
            real_s3, _S3_BUCKET, _S3_RAW_KEY,
            _S3_TEST_CHROM, _S3_TEST_TSS, _S3_TEST_WINDOW, _STREAM_USECOLS,
        ))
        chroms = {r["Chrom"].lstrip("chr") for r in rows}
        assert chroms == {_S3_TEST_CHROM}, f"Unexpected chroms: {chroms}"

    def test_raw_all_rows_within_window(self, real_s3):
        from scripts.lib.decode_stream import stream_s3_cis_rows
        rows = list(stream_s3_cis_rows(
            real_s3, _S3_BUCKET, _S3_RAW_KEY,
            _S3_TEST_CHROM, _S3_TEST_TSS, _S3_TEST_WINDOW, _STREAM_USECOLS,
        ))
        low  = _S3_TEST_TSS - _S3_TEST_WINDOW
        high = _S3_TEST_TSS + _S3_TEST_WINDOW
        for r in rows:
            pos = int(r["Pos"])
            assert low <= pos <= high, f"Row at pos {pos} outside [{low}, {high}]"

    def test_raw_rows_have_expected_columns(self, real_s3):
        from scripts.lib.decode_stream import stream_s3_cis_rows
        rows = list(stream_s3_cis_rows(
            real_s3, _S3_BUCKET, _S3_RAW_KEY,
            _S3_TEST_CHROM, _S3_TEST_TSS, _S3_TEST_WINDOW, _STREAM_USECOLS,
        ))
        assert rows
        assert set(rows[0].keys()) == _STREAM_USECOLS

    def test_raw_numeric_fields_parseable(self, real_s3):
        from scripts.lib.decode_stream import stream_s3_cis_rows
        import math
        rows = list(stream_s3_cis_rows(
            real_s3, _S3_BUCKET, _S3_RAW_KEY,
            _S3_TEST_CHROM, _S3_TEST_TSS, _S3_TEST_WINDOW, _STREAM_USECOLS,
        ))
        for r in rows[:20]:
            assert not math.isnan(float(r["Beta"]))
            assert 0 < float(r["Pval"]) <= 1
            assert float(r["SE"]) > 0
            assert int(r["N"]) > 0

    def test_smp_returns_nonempty_rows(self, real_s3):
        from scripts.lib.decode_stream import stream_s3_cis_rows
        rows = list(stream_s3_cis_rows(
            real_s3, _S3_BUCKET, _S3_SMP_KEY,
            _S3_TEST_CHROM, _S3_TEST_TSS, _S3_TEST_WINDOW, _STREAM_USECOLS,
        ))
        assert len(rows) > 100

    def test_raw_and_smp_share_same_variants(self, real_s3):
        """Both normalizations cover the same genomic window."""
        from scripts.lib.decode_stream import stream_s3_cis_rows
        raw = list(stream_s3_cis_rows(
            real_s3, _S3_BUCKET, _S3_RAW_KEY,
            _S3_TEST_CHROM, _S3_TEST_TSS, _S3_TEST_WINDOW, _STREAM_USECOLS,
        ))
        smp = list(stream_s3_cis_rows(
            real_s3, _S3_BUCKET, _S3_SMP_KEY,
            _S3_TEST_CHROM, _S3_TEST_TSS, _S3_TEST_WINDOW, _STREAM_USECOLS,
        ))
        raw_names = {r["Name"] for r in raw}
        smp_names = {r["Name"] for r in smp}
        # Overlap should be very high (same variants, different normalization)
        overlap = len(raw_names & smp_names) / max(len(raw_names), len(smp_names))
        assert overlap > 0.95, f"raw/smp variant overlap only {overlap:.1%}"

    def test_raw_and_smp_betas_differ(self, real_s3):
        """SMP normalization changes Beta values — confirms we're reading different files."""
        from scripts.lib.decode_stream import stream_s3_cis_rows
        raw = {r["Name"]: float(r["Beta"]) for r in stream_s3_cis_rows(
            real_s3, _S3_BUCKET, _S3_RAW_KEY,
            _S3_TEST_CHROM, _S3_TEST_TSS, _S3_TEST_WINDOW, _STREAM_USECOLS,
        )}
        smp = {r["Name"]: float(r["Beta"]) for r in stream_s3_cis_rows(
            real_s3, _S3_BUCKET, _S3_SMP_KEY,
            _S3_TEST_CHROM, _S3_TEST_TSS, _S3_TEST_WINDOW, _STREAM_USECOLS,
        )}
        common = set(raw) & set(smp)
        n_differ = sum(1 for k in common if abs(raw[k] - smp[k]) > 1e-6)
        assert n_differ > len(common) * 0.5, "Expected majority of Betas to differ between raw and SMP"

    def test_s3_key_index_lists_expected_protein(self, real_s3, tmp_path):
        """Live S3 listing contains the test protein in both raw and SMP prefixes."""
        import re
        for prefix, pat in [
            ("final_somascan_raw", r"Proteomics_PC0_(.+)_\d{8}\.txt\.gz"),
            ("final_somascan_smp", r"Proteomics_SMP_PC0_(.+)_\d{8}\.txt\.gz"),
        ]:
            regex = re.compile(pat)
            found = False
            for page in real_s3.get_paginator("list_objects_v2").paginate(
                Bucket=_S3_BUCKET, Prefix=f"{prefix}/"
            ):
                for obj in page.get("Contents", []):
                    m = regex.match(obj["Key"].split("/")[-1])
                    if m and m.group(1) == _S3_TEST_PROTEIN:
                        found = True
                        break
                if found:
                    break
            assert found, f"{_S3_TEST_PROTEIN} not found in s3://{_S3_BUCKET}/{prefix}/"

    def test_live_s3_key_index_builds_protein_list_without_bulk_urls(self, real_s3, tmp_path):
        """Live integration for the production path: S3 index keys → ProteinMeta list."""
        cohort = "deCODE_live_test"
        cohort_path = tmp_path / cohort
        cohort_path.mkdir()

        with patch.object(_decode_mod, "COHORT", cohort), \
             patch.object(_decode_mod, "cohort_dir", lambda c: tmp_path / c), \
             patch.object(_decode_mod, "_get_s3_client", return_value=real_s3):
            key_index = _decode_mod._build_s3_key_index(
                "final_somascan_raw",
                r"Proteomics_PC0_(.+)_\d{8}\.txt\.gz",
            )

            genes = sorted({name.split("_")[2] for name in key_index if len(name.split("_")) >= 3})
            (cohort_path / "_tss_hg38.tsv").write_text(
                "gene\tchrom\ttss\tresolved_symbol\ttier\tsource\n"
                + "".join(f"{gene}\t1\t1000000\t{gene}\t1\ttest\n" for gene in genes)
            )

            with patch.object(_decode_mod, "resolve_tss") as mock_resolve:
                proteins = _decode_mod.build_protein_list(key_index.keys(), build="hg38")

        assert _S3_TEST_PROTEIN in key_index
        assert len(proteins) == len(key_index)
        assert {p.seqid for p in proteins} == set(key_index)
        assert next(p for p in proteins if p.seqid == _S3_TEST_PROTEIN).gene == "KCNAB2"
        mock_resolve.assert_not_called()
