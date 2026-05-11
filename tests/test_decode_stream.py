"""Tests for code.lib.decode_stream"""
import gzip
import io
import pytest
from unittest.mock import patch, MagicMock

from scripts.lib.decode_stream import parse_bulk_urls, iter_decode_rows, stream_s3_cis_rows
from scripts.lib.paths import DECODE_URLS


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

class _RawBytes(io.RawIOBase):
    """Minimal RawIOBase wrapper for bytes — compatible with io.BufferedReader."""
    def __init__(self, data: bytes):
        self._buf = memoryview(data)
        self._pos = 0

    def readinto(self, b):
        avail = len(self._buf) - self._pos
        n = min(len(b), avail)
        b[:n] = self._buf[self._pos: self._pos + n]
        self._pos += n
        return n

    def readable(self):
        return True


def _make_gz_tsv(rows: list[dict], header: list[str]) -> bytes:
    """Serialize rows to a gzip-compressed tab-separated byte string."""
    lines = "\t".join(header) + "\n"
    for row in rows:
        lines += "\t".join(str(row.get(c, "")) for c in header) + "\n"
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as fh:
        fh.write(lines.encode())
    return buf.getvalue()


def _mock_s3(rows: list[dict], header: list[str]) -> MagicMock:
    """Return a mock boto3 S3 client whose get_object returns the given rows."""
    gz = _make_gz_tsv(rows, header)
    client = MagicMock()
    client.get_object.return_value = {"Body": _RawBytes(gz)}
    return client


_HEADER = ["Chrom", "Pos", "Name", "Beta", "Pval", "SE", "N"]
_USECOLS = frozenset(_HEADER)


def _row(chrom, pos, name="v", beta="0.1", pval="1e-10", se="0.01", n="35000"):
    return {"Chrom": chrom, "Pos": str(pos), "Name": name,
            "Beta": beta, "Pval": pval, "SE": se, "N": n}


class TestParseBulkUrls:
    @pytest.fixture(scope="class")
    def urls(self):
        if not DECODE_URLS.exists():
            pytest.fail("deCODE bulk_urls.txt not present")
        return parse_bulk_urls(DECODE_URLS)

    def test_returns_list_of_tuples(self, urls):
        assert isinstance(urls, list)
        assert len(urls) > 0
        assert all(isinstance(t, tuple) and len(t) == 2 for t in urls)

    def test_only_digit_prefix_files(self, urls):
        for name, _ in urls:
            assert name[0].isdigit(), f"Non-digit prefix: {name}"

    def test_no_md5_files(self, urls):
        for name, _ in urls:
            assert not name.endswith(".md5sum")

    def test_gene_parseable_from_name(self, urls):
        for name, _ in urls[:10]:
            parts = name.split("_")
            assert len(parts) >= 3, f"Name doesn't have ≥3 parts: {name}"


class TestIterDecodeRows:
    def test_yields_dicts(self):
        data = b"ColA\tColB\nval1\tval2\nval3\tval4\n"
        gz_bytes = io.BytesIO()
        with gzip.GzipFile(fileobj=gz_bytes, mode="wb") as f:
            f.write(data)
        gz_content = gz_bytes.getvalue()

        with patch("scripts.lib.decode_stream.download_bytes", return_value=gz_content):
            rows = list(iter_decode_rows("http://fake.url"))

        assert rows == [{"ColA": "val1", "ColB": "val2"},
                        {"ColA": "val3", "ColB": "val4"}]

    def test_empty_file(self):
        data = b"ColA\tColB\n"
        gz_bytes = io.BytesIO()
        with gzip.GzipFile(fileobj=gz_bytes, mode="wb") as f:
            f.write(data)

        with patch("scripts.lib.decode_stream.download_bytes", return_value=gz_bytes.getvalue()):
            rows = list(iter_decode_rows("http://fake.url"))

        assert rows == []


# ---------------------------------------------------------------------------
# stream_s3_cis_rows
# ---------------------------------------------------------------------------

class TestStreamS3CisRows:
    # target: chr1, tss=1_000_000, window=200_000  →  low=800_000  high=1_200_000

    def _call(self, rows, target_chrom="1", tss=1_000_000, window=200_000, usecols=_USECOLS):
        s3 = _mock_s3(rows, _HEADER)
        return list(stream_s3_cis_rows(s3, "bucket", "key", target_chrom, tss, window, usecols))

    def test_yields_rows_inside_window(self):
        rows = [
            _row("chr1", 800_000, "in_low"),
            _row("chr1", 1_000_000, "in_mid"),
            _row("chr1", 1_200_000, "in_high"),
        ]
        result = self._call(rows)
        assert {r["Name"] for r in result} == {"in_low", "in_mid", "in_high"}

    def test_skips_rows_below_low_bound(self):
        rows = [
            _row("chr1", 500_000, "before_window"),
            _row("chr1", 1_000_000, "in_window"),
        ]
        result = self._call(rows)
        assert [r["Name"] for r in result] == ["in_window"]

    def test_aborts_after_high_bound(self):
        rows = [
            _row("chr1", 1_000_000, "in_window"),
            _row("chr1", 1_300_000, "past_window"),  # pos > high → break
            _row("chr1", 1_050_000, "never_reached"),  # sorted order violated but after break
        ]
        result = self._call(rows)
        assert [r["Name"] for r in result] == ["in_window"]

    def test_aborts_when_past_target_chrom(self):
        rows = [
            _row("chr1", 1_000_000, "in_window"),
            _row("chr2", 500_000, "chr2_row"),  # seen_target=True → break
        ]
        result = self._call(rows)
        assert [r["Name"] for r in result] == ["in_window"]

    def test_skips_pre_target_chrom_rows(self):
        rows = [
            _row("chr1", 1_000_000, "chr1_row"),  # pre-target, skip
            _row("chr2", 1_000_000, "in_window"),
        ]
        result = self._call(rows, target_chrom="2")
        assert [r["Name"] for r in result] == ["in_window"]

    def test_strips_chr_prefix_in_comparison(self):
        """chr22 in file is matched when target_chrom='22'."""
        rows = [_row("chr22", 1_000_000, "hit")]
        result = self._call(rows, target_chrom="22")
        assert len(result) == 1

    def test_only_usecols_in_output(self):
        rows = [_row("chr1", 1_000_000, "v")]
        usecols = frozenset(["Name", "Beta"])
        result = self._call(rows, usecols=usecols)
        assert result[0].keys() == {"Name", "Beta"}

    def test_empty_file_returns_empty(self):
        assert self._call([]) == []

    def test_no_rows_on_target_chrom_returns_empty(self):
        rows = [_row("chr2", 1_000_000, "chr2_only")]
        assert self._call(rows, target_chrom="1") == []

    def test_s3_get_object_called_with_correct_params(self):
        rows = [_row("chr1", 1_000_000)]
        s3 = _mock_s3(rows, _HEADER)
        list(stream_s3_cis_rows(s3, "my-bucket", "prefix/key.txt.gz", "1", 1_000_000, 200_000, _USECOLS))
        s3.get_object.assert_called_once_with(Bucket="my-bucket", Key="prefix/key.txt.gz")

    def test_row_at_exact_low_bound_is_included(self):
        rows = [_row("chr1", 800_000, "exact_low")]   # pos == tss - window
        result = self._call(rows)
        assert len(result) == 1 and result[0]["Name"] == "exact_low"

    def test_row_at_exact_high_bound_is_included(self):
        rows = [_row("chr1", 1_200_000, "exact_high")]  # pos == tss + window
        result = self._call(rows)
        assert len(result) == 1 and result[0]["Name"] == "exact_high"

    def test_row_one_past_high_bound_excluded(self):
        rows = [
            _row("chr1", 1_200_000, "at_high"),
            _row("chr1", 1_200_001, "past_high"),
        ]
        result = self._call(rows)
        assert [r["Name"] for r in result] == ["at_high"]

    def test_malformed_row_short_columns_skipped(self):
        """Rows with too few columns don't crash; valid rows still yielded."""
        header = _HEADER
        # Build raw content manually to inject a short row
        lines = "\t".join(header) + "\n"
        lines += "chr1\t1000000\n"  # only 2 fields — short row
        lines += "\t".join([str(v) for v in ["chr1", "1050000", "good_row", "0.1", "1e-9", "0.01", "35000"]]) + "\n"
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb") as fh:
            fh.write(lines.encode())
        s3 = MagicMock()
        s3.get_object.return_value = {"Body": _RawBytes(buf.getvalue())}
        result = list(stream_s3_cis_rows(s3, "b", "k", "1", 1_000_000, 200_000, _USECOLS))
        names = [r["Name"] for r in result]
        assert "good_row" in names

    def test_invalid_pos_row_skipped(self):
        """Non-integer Pos values don't crash; other rows still yielded."""
        rows = [
            _row("chr1", "not_a_number", "bad_pos"),
            _row("chr1", 1_000_000, "good_pos"),
        ]
        result = self._call(rows)
        names = [r["Name"] for r in result]
        assert "good_pos" in names
        assert "bad_pos" not in names

    def test_no_chr_prefix_in_file_matches_target(self):
        """Files where Chrom is '1' not 'chr1' still match target_chrom='1'."""
        rows = [_row("1", 1_000_000, "no_prefix")]
        result = self._call(rows, target_chrom="1")
        assert len(result) == 1

    def test_zero_window_only_exact_tss(self):
        rows = [
            _row("chr1", 999_999, "before_tss"),
            _row("chr1", 1_000_000, "at_tss"),
            _row("chr1", 1_000_001, "after_tss"),
        ]
        result = self._call(rows, window=0)
        assert [r["Name"] for r in result] == ["at_tss"]

    def test_multiple_rows_in_window_all_yielded(self):
        rows = [_row("chr1", 900_000 + i * 50_000, f"v{i}") for i in range(7)]
        # low=800000 high=1200000; positions: 900000,950000,1000000,1050000,1100000,1150000,1200000
        result = self._call(rows)
        assert len(result) == 7

    def test_correct_values_in_yielded_row(self):
        rows = [_row("chr1", 1_000_000, name="rs123", beta="-0.5", pval="3e-12", se="0.025", n="40000")]
        result = self._call(rows)
        assert len(result) == 1
        r = result[0]
        assert r["Name"] == "rs123"
        assert r["Beta"] == "-0.5"
        assert r["Pval"] == "3e-12"
        assert r["N"] == "40000"

    def test_multiple_chroms_before_target_all_skipped(self):
        rows = [
            _row("chr1", 1_000_000, "chr1_row"),
            _row("chr2", 1_000_000, "chr2_row"),
            _row("chr3", 1_000_000, "chr3_target"),
            _row("chr4", 1_000_000, "chr4_row"),
        ]
        result = self._call(rows, target_chrom="3")
        assert [r["Name"] for r in result] == ["chr3_target"]

    def test_extra_columns_in_file_not_in_usecols_excluded(self):
        header_extra = _HEADER + ["ExtraCol1", "ExtraCol2"]
        rows_extra = [dict(_row("chr1", 1_000_000, "v"), ExtraCol1="x", ExtraCol2="y")]
        gz = _make_gz_tsv(rows_extra, header_extra)
        s3 = MagicMock()
        s3.get_object.return_value = {"Body": _RawBytes(gz)}
        result = list(stream_s3_cis_rows(s3, "b", "k", "1", 1_000_000, 200_000, _USECOLS))
        assert len(result) == 1
        assert "ExtraCol1" not in result[0]
        assert "ExtraCol2" not in result[0]

    def test_usecols_subset_returns_only_requested(self):
        rows = [_row("chr1", 1_000_000, "v")]
        result = self._call(rows, usecols=frozenset(["Name", "Pval", "N"]))
        assert result[0].keys() == {"Name", "Pval", "N"}
