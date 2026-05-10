"""Tests for code.lib.sumstats_io"""
import pandas as pd
import pytest

from scripts.lib.sumstats_io import read_norm, write_norm


class TestWriteReadNorm:
    def test_roundtrip(self, tmp_path):
        df = pd.DataFrame({
            "seqid": ["SeqId_1"], "gene": ["BRCA2"],
            "chrom": ["13"], "pos": [32_914_437],
            "rsid": ["rs1234"], "EA": ["A"], "OA": ["G"],
            "EAF": [0.3], "beta": [0.05], "se": [0.01],
            "pval": [1e-8], "N": [7213], "build": ["hg19"],
        })
        out = tmp_path / "test.tsv"
        write_norm(df, out)
        loaded = read_norm(out)
        assert list(loaded.columns) == list(df.columns)
        assert loaded["chrom"].iloc[0] == "13"
        assert loaded["rsid"].iloc[0] == "rs1234"

    def test_creates_parent_dir(self, tmp_path):
        out = tmp_path / "subdir" / "test.tsv"
        df = pd.DataFrame({"beta": [0.1]})
        write_norm(df, out)
        assert out.exists()

    def test_chrom_preserved_as_str(self, tmp_path):
        df = pd.DataFrame({"chrom": ["1", "22", "X"]})
        out = tmp_path / "chrom.tsv"
        write_norm(df, out)
        loaded = read_norm(out)
        assert loaded["chrom"].tolist() == ["1", "22", "X"]
