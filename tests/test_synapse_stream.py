"""Tests for code.lib.synapse_stream"""
import gzip
import io
import tarfile
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

from scripts.lib.synapse_stream import iter_gz_rows, iter_tar_gz_rows, stream_ukbppp_protein


def _make_gz(content: bytes) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write(content)
    return buf.getvalue()


def _make_tar(members: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tf:
        for name, data in members.items():
            info = tarfile.TarInfo(name=name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
    return buf.getvalue()


class TestIterGzRows:
    def test_yields_header_keyed_dicts(self):
        data = b"A\tB\nval1\tval2\n"
        rows = list(iter_gz_rows(_make_gz(data)))
        assert rows == [{"A": "val1", "B": "val2"}]

    def test_skips_mismatched_rows(self):
        data = b"A\tB\nval1\tval2\nextra_col1\textra_col2\textra_col3\n"
        rows = list(iter_gz_rows(_make_gz(data)))
        assert len(rows) == 1

    def test_empty_gz(self):
        data = b"A\tB\n"
        rows = list(iter_gz_rows(_make_gz(data)))
        assert rows == []


class TestIterTarGzRows:
    def test_reads_gz_members(self, tmp_path):
        gz_content = _make_gz(b"COL1\tCOL2\nrow1a\trow1b\n")
        tar_bytes = _make_tar({"chr1.gz": gz_content})
        tar_path = tmp_path / "test.tar"
        tar_path.write_bytes(tar_bytes)

        rows = list(iter_tar_gz_rows(tar_path))
        assert len(rows) == 1
        member, row = rows[0]
        assert member == "chr1.gz"
        assert row == {"COL1": "row1a", "COL2": "row1b"}

    def test_skips_non_gz_members(self, tmp_path):
        tar_bytes = _make_tar({"readme.txt": b"hello"})
        tar_path = tmp_path / "test.tar"
        tar_path.write_bytes(tar_bytes)

        rows = list(iter_tar_gz_rows(tar_path))
        assert rows == []

    def test_member_filter(self, tmp_path):
        gz1 = _make_gz(b"ID\n1:100:A:T\n")
        gz2 = _make_gz(b"ID\n1:200:C:G\n")
        tar_bytes = _make_tar({"chr1.gz": gz1, "chr2.gz": gz2})
        tar_path = tmp_path / "test.tar"
        tar_path.write_bytes(tar_bytes)

        rows = list(iter_tar_gz_rows(tar_path, member_filter=lambda n: "chr1" in n))
        assert len(rows) == 1


class TestStreamUkbpppProtein:
    def test_filters_to_cis_window(self, tmp_path):
        # UKB ID format: chr<c>:<pos>:ref:alt; real files are space-delimited
        content = b"ID BETA\n1:100000:A:T 0.1\n1:200000:A:T 0.2\n2:100000:G:C 0.3\n"
        gz_content = _make_gz(content)
        tar_bytes = _make_tar({"chr1.gz": gz_content, "chr2.gz": _make_gz(
            b"ID BETA\n2:100000:G:C 0.3\n"
        )})

        # Download mock returns our tar
        with patch("scripts.lib.synapse_stream.download_entity") as mock_dl:
            tar_path = tmp_path / "fake.tar"
            tar_path.write_bytes(tar_bytes)
            mock_dl.return_value = tar_path

            rows = stream_ukbppp_protein(
                "syn123",
                cis_chrom="1", cis_start=50_000, cis_end=150_000,
                tmp_dir=tmp_path,
            )

        assert len(rows) == 1
        assert rows[0]["ID"] == "1:100000:A:T"
