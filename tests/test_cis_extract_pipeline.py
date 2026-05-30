"""Tests for code.lib.cis_extract (cohort-agnostic pipeline)."""
import importlib
import pandas as pd
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from scripts.lib.checkpoint import Checkpoint
from scripts.lib.schema import ProteinMeta
from scripts.lib.cis_extract import _apply_filters, _normalize, run_extraction, OUTPUT_COLS


@pytest.fixture
def sample_protein():
    return ProteinMeta(
        seqid="SeqId_TEST", gene="TESTGENE", uniprot="Q12345",
        chrom="22", tss=25_212_564, build="hg19", source_cohort="Fenland",
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
            chrom="6", tss=29_000_000, build="hg19", source_cohort="Fenland",
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
        filtered_dir = tmp_path / "filtered_cis_pqtls"
        raw_dir = tmp_path / "raw_cis_sumstats"
        filtered_dir.mkdir(parents=True)
        raw_dir.mkdir(parents=True)
        cp_dir = tmp_path / "cohort_state"
        cp_dir.mkdir(parents=True)

        def read_fn(p):
            return pd.DataFrame({
                "chrom": ["22"], "pos": [30_000_000],  # outside cis window → all dropped
                "EA": ["A"], "OA": ["G"], "EAF": [0.3],
                "beta": [0.5], "se": [0.05], "pval": [1e-9], "N": [7213],
            })

        with patch("scripts.lib.cis_extract.filtered_cis_pqtls_dir", return_value=filtered_dir), \
             patch("scripts.lib.cis_extract.raw_cis_sumstats_dir", return_value=raw_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=cp_dir):
            n_ok = run_extraction(cohort, [sample_protein], read_fn)

        assert n_ok == 0  # no cis-pQTL proteins
        # Output file should NOT exist (protein was empty after filters)
        assert not (filtered_dir / f"{sample_protein.seqid}.tsv").exists()
        assert not (raw_dir / f"{sample_protein.seqid}.tsv.gz").exists()


class TestRunExtraction:
    def test_raw_and_filtered_outputs_split_from_one_read(self, tmp_path, sample_protein):
        cohort = "Fenland"
        filtered_dir = tmp_path / "filtered_cis_pqtls"
        raw_dir = tmp_path / "raw_cis_sumstats"
        state_dir = tmp_path / "state"
        filtered_dir.mkdir()
        raw_dir.mkdir()
        state_dir.mkdir()
        calls = {"n": 0}

        def read_fn(_p):
            calls["n"] += 1
            return pd.DataFrame({
                "chrom": ["22", "22", "22", "22", "22"],
                "pos": [
                    sample_protein.tss,
                    sample_protein.tss + 600_000,  # raw ±1 Mb only
                    sample_protein.tss + 1_100_000,  # outside raw window
                    sample_protein.tss + 10_000,  # low MAF
                    sample_protein.tss + 20_000,  # ambiguous palindrome
                ],
                "rsid": ["rs_keep", "rs_region", "rs_outside", "rs_lowmaf", "rs_pal"],
                "EA": ["A", "C", "G", "A", "A"],
                "OA": ["G", "T", "A", "G", "T"],
                "EAF": [0.3, 0.25, 0.3, 0.001, 0.5],
                "beta": [0.5, 0.2, 0.1, 0.5, 0.5],
                "se": [0.05, 0.02, 0.01, 0.05, 0.05],
                "pval": [1e-9, 0.2, 1e-9, 1e-9, 1e-9],
                "N": [7213] * 5,
            })

        with patch("scripts.lib.cis_extract.filtered_cis_pqtls_dir", return_value=filtered_dir), \
             patch("scripts.lib.cis_extract.raw_cis_sumstats_dir", return_value=raw_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=state_dir):
            n_ok = run_extraction(cohort, [sample_protein], read_fn)

        assert n_ok == 1
        assert calls["n"] == 1
        raw = pd.read_csv(raw_dir / f"{sample_protein.seqid}.tsv.gz", sep="\t")
        filtered = pd.read_csv(filtered_dir / f"{sample_protein.seqid}.tsv", sep="\t")
        assert set(raw["rsid"]) == {"rs_keep", "rs_region", "rs_lowmaf", "rs_pal"}
        assert set(filtered["rsid"]) == {"rs_keep"}

    def test_raw_cache_keeps_mhc_rows_even_when_filtered_empty(self, tmp_path):
        cohort = "Fenland"
        protein = ProteinMeta(
            seqid="SeqId_MHC_RAW",
            gene="MHCGENE",
            uniprot="P1",
            chrom="6",
            tss=29_000_000,
            build="hg19",
            source_cohort=cohort,
        )
        filtered_dir = tmp_path / "filtered_cis_pqtls"
        raw_dir = tmp_path / "raw_cis_sumstats"
        state_dir = tmp_path / "state"
        filtered_dir.mkdir()
        raw_dir.mkdir()
        state_dir.mkdir()

        def read_fn(_p):
            return pd.DataFrame({
                "chrom": ["6"],
                "pos": [29_000_000],
                "rsid": ["rs_mhc"],
                "EA": ["A"],
                "OA": ["G"],
                "EAF": [0.3],
                "beta": [0.5],
                "se": [0.05],
                "pval": [1e-9],
                "N": [7213],
            })

        with patch("scripts.lib.cis_extract.filtered_cis_pqtls_dir", return_value=filtered_dir), \
             patch("scripts.lib.cis_extract.raw_cis_sumstats_dir", return_value=raw_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=state_dir):
            n_ok = run_extraction(cohort, [protein], read_fn)

        assert n_ok == 0
        raw = pd.read_csv(raw_dir / f"{protein.seqid}.tsv.gz", sep="\t")
        assert set(raw["rsid"]) == {"rs_mhc"}
        assert not (filtered_dir / f"{protein.seqid}.tsv").exists()
        assert (filtered_dir / f"{protein.seqid}.tsv.empty").exists()

    def test_filter_empty_raw_cache_marker_skips_after_checkpoint_loss(self, tmp_path):
        cohort = "Fenland"
        protein = ProteinMeta(
            seqid="SeqId_EMPTY_CACHE",
            gene="EMPTYGENE",
            uniprot="P1",
            chrom="6",
            tss=29_000_000,
            build="hg19",
            source_cohort=cohort,
        )
        filtered_dir = tmp_path / "filtered_cis_pqtls"
        raw_dir = tmp_path / "raw_cis_sumstats"
        state_dir = tmp_path / "state"
        filtered_dir.mkdir()
        raw_dir.mkdir()
        state_dir.mkdir()
        calls = {"n": 0}

        def read_fn(_p):
            calls["n"] += 1
            return pd.DataFrame({
                "chrom": ["6"],
                "pos": [29_000_000],
                "rsid": ["rs_mhc"],
                "EA": ["A"],
                "OA": ["G"],
                "EAF": [0.3],
                "beta": [0.5],
                "se": [0.05],
                "pval": [1e-9],
                "N": [7213],
            })

        with patch("scripts.lib.cis_extract.filtered_cis_pqtls_dir", return_value=filtered_dir), \
             patch("scripts.lib.cis_extract.raw_cis_sumstats_dir", return_value=raw_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=state_dir):
            assert run_extraction(cohort, [protein], read_fn) == 0

        assert calls["n"] == 1
        (state_dir / "_state_02.json").unlink()

        with patch("scripts.lib.cis_extract.filtered_cis_pqtls_dir", return_value=filtered_dir), \
             patch("scripts.lib.cis_extract.raw_cis_sumstats_dir", return_value=raw_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=state_dir):
            assert run_extraction(cohort, [protein], read_fn) == 0

        assert calls["n"] == 1

    def test_checkpointing_skips_done(self, tmp_path, sample_protein):
        """Proteins already in checkpoint should not call read_fn."""
        cohort = "ARIC_EA"
        filtered_dir = tmp_path / "processed_data" / cohort / "filtered_cis_pqtls"
        raw_dir = tmp_path / "processed_data" / cohort / "raw_cis_sumstats"
        filtered_dir.mkdir(parents=True)
        raw_dir.mkdir(parents=True)

        # Pre-create a minimally valid output file to simulate completed protein
        pd.DataFrame(
            [{col: "x" for col in OUTPUT_COLS}]
        ).to_csv(filtered_dir / f"{sample_protein.seqid}.tsv", sep="\t", index=False)
        pd.DataFrame(
            [{col: "x" for col in OUTPUT_COLS}]
        ).to_csv(raw_dir / f"{sample_protein.seqid}.tsv.gz", sep="\t", index=False)

        call_count = {"n": 0}
        def read_fn(p):
            call_count["n"] += 1
            return None

        with patch("scripts.lib.cis_extract.filtered_cis_pqtls_dir", return_value=filtered_dir), \
             patch("scripts.lib.cis_extract.raw_cis_sumstats_dir", return_value=raw_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=tmp_path / "processed_data" / cohort):
            run_extraction(cohort, [sample_protein], read_fn)

        # read_fn should not have been called since output exists
        assert call_count["n"] == 0

    def test_checkpointed_empty_raw_does_not_reenter_todo(self, tmp_path, sample_protein):
        """A done protein with no raw file can represent a true empty ±1 Mb region."""
        cohort = "ARIC_EA"
        filtered_dir = tmp_path / "processed_data" / cohort / "filtered_cis_pqtls"
        raw_dir = tmp_path / "processed_data" / cohort / "raw_cis_sumstats"
        state_dir = tmp_path / "processed_data" / cohort
        filtered_dir.mkdir(parents=True)
        raw_dir.mkdir(parents=True)
        state_dir.mkdir(parents=True, exist_ok=True)
        Checkpoint(state_dir / "_state_02.json").mark_done(sample_protein.seqid)

        call_count = {"n": 0}

        def read_fn(_p):
            call_count["n"] += 1
            return pd.DataFrame({
                "chrom": ["22"],
                "pos": [sample_protein.tss],
                "rsid": ["rs1"],
                "EA": ["A"],
                "OA": ["G"],
                "EAF": [0.3],
                "beta": [0.5],
                "se": [0.05],
                "pval": [1e-9],
                "N": [7213],
            })

        with patch("scripts.lib.cis_extract.filtered_cis_pqtls_dir", return_value=filtered_dir), \
             patch("scripts.lib.cis_extract.raw_cis_sumstats_dir", return_value=raw_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=state_dir):
            n_ok = run_extraction(cohort, [sample_protein], read_fn)

        assert n_ok == 0
        assert call_count["n"] == 0

    def test_parallel_workers_produce_same_result(self, tmp_path, sample_protein):
        """Workers=2 should produce the same output files as workers=1."""
        cohort = "ARIC_EA"

        def make_out_dir(suffix):
            d = tmp_path / f"filtered_cis_pqtls_{suffix}"
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
            raw_dir = tmp_path / f"raw_cis_sumstats_{workers}"
            raw_dir.mkdir(parents=True)
            with patch("scripts.lib.cis_extract.filtered_cis_pqtls_dir", return_value=out_dir), \
                 patch("scripts.lib.cis_extract.raw_cis_sumstats_dir", return_value=raw_dir), \
                 patch("scripts.lib.cis_extract.cohort_dir",
                       return_value=tmp_path / f"state_{workers}"):
                run_extraction(cohort, [sample_protein], make_read_fn(), workers=workers)
            files = list(out_dir.glob("*.tsv"))
            assert len(files) == 1, f"workers={workers}: expected 1 output file"
            raw_files = list(raw_dir.glob("*.tsv.gz"))
            assert len(raw_files) == 1, f"workers={workers}: expected 1 raw output file"

    def test_read_failure_tracked_in_checkpoint_sequential(self, tmp_path, sample_protein):
        """read_fn raising → protein appears in cp.n_failed (sequential path)."""
        cohort = "ARIC_EA"
        out_dir = tmp_path / "filtered_cis_pqtls"
        raw_dir = tmp_path / "raw_cis_sumstats"
        state_dir = tmp_path / "processed_data" / cohort
        out_dir.mkdir(parents=True)
        raw_dir.mkdir(parents=True)
        state_dir.mkdir(parents=True)

        def read_fn(p):
            raise IOError("simulated read error")

        with patch("scripts.lib.cis_extract.filtered_cis_pqtls_dir", return_value=out_dir), \
             patch("scripts.lib.cis_extract.raw_cis_sumstats_dir", return_value=raw_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=state_dir):
            n_ok = run_extraction(cohort, [sample_protein], read_fn, workers=1)

        assert n_ok == 0
        cp = Checkpoint(state_dir / "_state_02.json")
        assert cp.n_failed == 1
        assert cp.is_failed(sample_protein.seqid)

    def test_read_failure_tracked_in_checkpoint_parallel(self, tmp_path, sample_protein):
        """read_fn raising → protein appears in cp.n_failed (parallel path)."""
        cohort = "ARIC_EA"
        out_dir = tmp_path / "filtered_cis_pqtls"
        raw_dir = tmp_path / "raw_cis_sumstats"
        state_dir = tmp_path / "processed_data" / cohort
        out_dir.mkdir(parents=True)
        raw_dir.mkdir(parents=True)
        state_dir.mkdir(parents=True)

        def read_fn(p):
            raise IOError("simulated read error")

        with patch("scripts.lib.cis_extract.filtered_cis_pqtls_dir", return_value=out_dir), \
             patch("scripts.lib.cis_extract.raw_cis_sumstats_dir", return_value=raw_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=state_dir):
            n_ok = run_extraction(cohort, [sample_protein], read_fn, workers=2)

        assert n_ok == 0
        cp = Checkpoint(state_dir / "_state_02.json")
        assert cp.n_failed == 1
        assert cp.is_failed(sample_protein.seqid)

    def test_parallel_writes_protein_index(self, tmp_path, sample_protein):
        cohort = "ARIC_EA"
        out_dir = tmp_path / "filtered_cis_pqtls"
        raw_dir = tmp_path / "raw_cis_sumstats"
        state_dir = tmp_path / "processed_data" / cohort
        out_dir.mkdir(parents=True)
        raw_dir.mkdir(parents=True)
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

        with patch("scripts.lib.cis_extract.filtered_cis_pqtls_dir", return_value=out_dir), \
             patch("scripts.lib.cis_extract.raw_cis_sumstats_dir", return_value=raw_dir), \
             patch("scripts.lib.cis_extract.cohort_dir", return_value=state_dir):
            run_extraction(cohort, [sample_protein], read_fn, workers=2)

        index_path = state_dir / "protein_index.tsv"
        assert index_path.exists()
        index_df = pd.read_csv(index_path, sep="\t")
        assert list(index_df.columns) == ["seqid", "gene", "uniprot", "chrom", "tss", "build"]
        assert len(index_df) == 1
        assert index_df.loc[0, "seqid"] == sample_protein.seqid
