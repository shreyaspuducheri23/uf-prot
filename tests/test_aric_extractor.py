"""Tests for code/02_cis_pqtl_extract/aric.py"""
import importlib
import sys
import pytest
import pandas as pd
from pathlib import Path

from scripts.lib.cis import load_aric_tss
from scripts.lib.schema import ProteinMeta
from scripts.lib.paths import ARIC_EA_DIR, ARIC_SEQID

aric = importlib.import_module("scripts.02_cis_pqtl_extract.aric")


class TestReadAricProtein:
    @pytest.fixture(scope="class")
    def sample_protein(self):
        tss_map = load_aric_tss(ARIC_SEQID)
        seqid = "SeqId_10000_28"
        chrom, tss, uniprot, gene = tss_map[seqid]
        return ProteinMeta(seqid=seqid, gene=gene, uniprot=uniprot,
                           chrom=chrom, tss=tss, build="hg38",
                           source_cohort="ARIC_EA")

    def test_returns_dataframe_or_none(self, sample_protein):
        result = aric.read_aric_protein(sample_protein)
        assert result is None or isinstance(result, pd.DataFrame)

    def test_required_columns_present(self, sample_protein):
        result = aric.read_aric_protein(sample_protein)
        if result is None:
            pytest.fail("File not found for sample protein")
        for col in ["chrom", "pos", "EA", "OA", "EAF", "beta", "se", "pval", "N"]:
            assert col in result.columns, f"Missing column: {col}"

    def test_rsid_dot_for_non_rs(self, sample_protein):
        result = aric.read_aric_protein(sample_protein)
        if result is None:
            pytest.fail("File not found")
        non_rs = result[~result["rsid"].str.startswith("rs")]
        assert all(non_rs["rsid"] == ".")

    def test_chrom_is_string(self, sample_protein):
        result = aric.read_aric_protein(sample_protein)
        if result is None:
            pytest.fail("File not found")
        assert pd.api.types.is_string_dtype(result["chrom"])

    def test_missing_seqid_returns_none(self):
        protein = ProteinMeta(
            seqid="SeqId_NONEXISTENT", gene="FAKE", uniprot="",
            chrom="1", tss=1_000_000, build="hg38", source_cohort="ARIC_EA",
        )
        result = aric.read_aric_protein(protein)
        assert result is None


class TestLoadAricProteins:
    def test_loads_4657_proteins(self):
        proteins = aric.load_aric_proteins()
        assert len(proteins) == 4657

    def test_all_are_protein_meta(self):
        proteins = aric.load_aric_proteins()
        for p in proteins[:10]:
            assert isinstance(p, ProteinMeta)
            assert p.build == "hg38"
            assert p.source_cohort == "ARIC_EA"
