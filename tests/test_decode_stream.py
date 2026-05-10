"""Tests for code.lib.decode_stream"""
import gzip
import io
import pytest
from unittest.mock import patch, MagicMock

from scripts.lib.decode_stream import parse_bulk_urls, iter_decode_rows
from scripts.lib.paths import DECODE_URLS


class TestParseBulkUrls:
    @pytest.fixture(scope="class")
    def urls(self):
        if not DECODE_URLS.exists():
            pytest.skip("deCODE bulk_urls.txt not present")
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
