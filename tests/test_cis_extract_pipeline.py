"""Tests for code.lib.cis_extract (cohort-agnostic pipeline)."""
import importlib
import pandas as pd
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from scripts.lib.schema import ProteinMeta
from scripts.lib.cis_extract import _apply_filters, _normalize, run_extraction, OUTPUT_COLS


@pytest.fixture
def sample_protein():
    return ProteinMeta(
        seqid="SeqId_TEST", gene="TESTGENE", uniprot="Q12345",
        chrom="22", tss=25_212_564, build="hg19", source_cohort="ARIC_EA",
    )


@pytest.fixture
def raw_cis_df():
    """Simulated raw sumstats with variants spanning the cis window and beyond."""
    return pd.DataFrame({
        "chrom": ["22", "22", "22", "22"],
        "pos":   [25_212_064, 25_712_000, 26_000_000, 1_000],
        "EA":    ["A", "C", "G", "T"],
        "OA":    ["G", "T", "A", "C"],
        "EAF":   [0.30, 0.25, 0.45, 0.10],
        "beta":  [0.5,  0.3,  0.2,  0.1],
        "se":    [0.05, 0.03, 0.02, 0.01],
        "pval":  [1e-9, 5e-10, 0.05, 1e-12],
        "N":     [7213, 7213, 7213, 7213],
    })


class TestApplyFilters:
    def test_cis_window_applied(self, sample_protein, raw_cis_df):
        result = _apply_filters(raw_cis_df, sample_protein)
        # pos=26_000_000: |26M - 25.2M| = 787_436 > 500_000 → outside cis → drop
        # pos=1_000: on chr22 but 25M away → drop
        for pos in result["pos"]:
            assert abs(pos - sample_protein.tss) <= 500_000

    def test_gw_significance_applied(self, sample_protein, raw_cis_df):
        result = _apply_filters(raw_cis_df, sample_protein)
        # p=0.05 is not genome-wide significant
        assert all(result["pval"] < 5e-8)

    def test_empty_on_no_variants(self, sample_protein):
        df = pd.DataFrame({
            "chrom": ["22"], "pos": [30_000_000],
            "EA": ["A"], "OA": ["G"], "EAF": [0.3],
            "beta": [0.1], "se": [0.01], "pval": [1e-10], "N": [7213],
        })
        result = _apply_filters(df, sample_protein)
        assert result.empty


class TestNormalize:
    def test_adds_metadata_columns(self, sample_protein):
        df = pd.DataFrame({
            "chrom": ["22"], "pos": [25_212_564],
            "rsid": ["rs123"], "EA": ["A"], "OA": ["G"],
            "EAF": [0.3], "beta": [0.5], "se": [0.05], "pval": [1e-9], "N": [7213],
        })
        result = _normalize(df, sample_protein)
        assert result["seqid"].iloc[0] == "SeqId_TEST"
        assert result["gene"].iloc[0] == "TESTGENE"
        assert result["uniprot"].iloc[0] == "Q12345"
        assert result["build"].iloc[0] == "hg19"

    def test_output_col_order(self, sample_protein):
        from scripts.lib.cis_extract import OUTPUT_COLS
        df = pd.DataFrame({c: ["x"] for c in OUTPUT_COLS})
        result = _normalize(df, sample_protein)
        assert list(result.columns) == OUTPUT_COLS


class TestApplyFiltersAllSix:
    """Verify that all 6 filters in _apply_filters are actually applied."""

    def test_maf_filter_applied(self, sample_protein):
        df = pd.DataFrame({
            "chrom": ["22"], "pos": [25_212_064],
            "EA": ["A"], "OA": ["G"], "EAF": [0.001],  # MAF too low → dropped
            "beta": [0.5], "se": [0.05], "pval": [1e-9], "N": [7213],
        })
        result = _apply_filters(df, sample_protein)
        assert result.empty

    def test_mhc_filter_applied(self):
        protein_chr6 = ProteinMeta(
            seqid="SeqId_MHC", gene="GENE_MHC", uniprot="Q99999",
            chrom="6", tss=29_000_000, build="hg19", source_cohort="ARIC_EA",
        )
        df = pd.DataFrame({
            "chrom": ["6"], "pos": [29_000_000],  # inside MHC hg19
            "EA": ["A"], "OA": ["G"], "EAF": [0.3],
            "beta": [0.5], "se": [0.05], "pval": [1e-9], "N": [7213],
        })
        result = _apply_filters(df, protein_chr6)
        assert result.empty

    def test_palindrome_filter_applied(self, sample_protein):
        df = pd.DataFrame({
            "chrom": ["22"], "pos": [25_212_064],
            "EA": ["A"], "OA": ["T"], "EAF": [0.50],  # palindrome, MAF=0.5 > 0.42 → dropped
            "beta": [0.5], "se": [0.05], "pval": [1e-9], "N": [7213],
        })
        result = _apply_filters(df, sample_protein)
        assert result.empty

    def test_variant_with_zero_survivors_checkpointed(self, tmp_path, sample_protein):
        """Proteins with 0 variants after filtering should still be checkpointed."""
        cohort = "ARIC_EA"
        out_dir = tmp_path / "cis_sumstats"
        out_dir.mkdir(parents=True)
        cp_dir = tmp_path / "cohort_state"
        cp_dir.mkdir(parents=True)

        def read_fn(p):
            return pd.DataFrame({
                "chrom": ["22"], "pos": [30_000_000],  # outside cis window → all dropped
                "EA": ["A"], "OA": ["G"], "EAF": [0.3],
                "beta": [0.5], "se": [0.05], "pval": [1e-9], "N": [7213],
            })

        with patch("scripts.lib.cis_extract.cis_sumstats_dir", return_value=out_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=cp_dir):
            n_ok = run_extraction(cohort, [sample_protein], read_fn)

        assert n_ok == 0  # no cis-pQTL proteins
        # Output file should NOT exist (protein was empty after filters)
        assert not (out_dir / f"{sample_protein.seqid}.tsv").exists()


class TestRunExtraction:
    def test_checkpointing_skips_done(self, tmp_path, sample_protein):
        """Proteins already in checkpoint should not call read_fn."""
        cohort = "ARIC_EA"
        out_dir = tmp_path / "processed_data" / cohort / "cis_sumstats"
        out_dir.mkdir(parents=True)

        # Pre-create a minimally valid output file to simulate completed protein
        pd.DataFrame(
            [{col: "x" for col in OUTPUT_COLS}]
        ).to_csv(out_dir / f"{sample_protein.seqid}.tsv", sep="\t", index=False)

        call_count = {"n": 0}
        def read_fn(p):
            call_count["n"] += 1
            return None

        with patch("scripts.lib.cis_extract.cis_sumstats_dir", return_value=out_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=tmp_path / "processed_data" / cohort):
            run_extraction(cohort, [sample_protein], read_fn)

        # read_fn should not have been called since output exists
        assert call_count["n"] == 0

    def test_parallel_workers_produce_same_result(self, tmp_path, sample_protein):
        """Workers=2 should produce the same output files as workers=1."""
        cohort = "ARIC_EA"

        def make_out_dir(suffix):
            d = tmp_path / f"cis_sumstats_{suffix}"
            d.mkdir(parents=True)
            return d

        def make_read_fn():
            return lambda p: pd.DataFrame({
                "chrom": ["22"], "pos": [25_212_564],
                "rsid": ["rs123"],
                "EA": ["A"], "OA": ["G"], "EAF": [0.3],
                "beta": [0.5], "se": [0.05], "pval": [1e-9], "N": [7213],
            })

        for workers in (1, 2):
            out_dir = make_out_dir(workers)
            with patch("scripts.lib.cis_extract.cis_sumstats_dir", return_value=out_dir), \
                 patch("scripts.lib.cis_extract.cohort_dir",
                       return_value=tmp_path / f"state_{workers}"):
                run_extraction(cohort, [sample_protein], make_read_fn(), workers=workers)
            files = list(out_dir.glob("*.tsv"))
            assert len(files) == 1, f"workers={workers}: expected 1 output file"

    def test_parallel_writes_protein_index(self, tmp_path, sample_protein):
        cohort = "ARIC_EA"
        out_dir = tmp_path / "cis_sumstats"
        state_dir = tmp_path / "processed_data" / cohort
        out_dir.mkdir(parents=True)
        state_dir.mkdir(parents=True, exist_ok=True)

        def read_fn(_p):
            return pd.DataFrame({
                "chrom": ["22"],
                "pos": [25_212_564],
                "rsid": ["rs123"],
                "EA": ["A"],
                "OA": ["G"],
                "EAF": [0.3],
                "beta": [0.5],
                "se": [0.05],
                "pval": [1e-9],
                "N": [7213],
            })

        with patch("scripts.lib.cis_extract.cis_sumstats_dir", return_value=out_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=state_dir):
            run_extraction(cohort, [sample_protein], read_fn, workers=2)

        index_path = state_dir / "protein_index.tsv"
        assert index_path.exists()
        index_df = pd.read_csv(index_path, sep="\t")
        assert list(index_df.columns) == ["seqid", "gene", "uniprot", "chrom", "tss", "build"]
        assert len(index_df) == 1
        assert index_df.loc[0, "seqid"] == sample_protein.seqid
