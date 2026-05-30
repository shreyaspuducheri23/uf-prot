"""Tests for scripts.08_coloc.extract_regions."""
import importlib
import json
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

_mod = importlib.import_module("scripts.08_coloc.extract_regions")


def _write_candidates(base: Path, seqid: str, build: str = "hg19") -> None:
    base.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"seqid": [seqid], "fdr_pass": [True]}).to_csv(
        base / "mr_results.tsv", sep="\t", index=False
    )
    pd.DataFrame({"seqid": [seqid], "chrom": ["1"], "tss": [100_000], "build": [build]}).to_csv(
        base / "protein_index.tsv", sep="\t", index=False
    )


def _write_raw_cis_hg38(raw_dir: Path, seqid: str, df: pd.DataFrame) -> None:
    raw_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(raw_dir / f"{seqid}.tsv.gz", sep="\t", index=False)


class FakeOutcome:
    def __init__(self, out_df: pd.DataFrame):
        self.out_df = out_df
        self.calls = []

    def fetch_region(self, chrom: str, start: int, end: int) -> pd.DataFrame:
        self.calls.append((chrom, start, end))
        return self.out_df.copy()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return None


def test_hg19_region_lifted_before_outcome_query(tmp_path):
    cohort = "ARIC_EA"
    seqid = "SeqId_1"
    cohort_base = tmp_path / cohort
    _write_candidates(cohort_base, seqid, build="hg19")

    exp_df = pd.DataFrame(
        {
            "chrom": ["1"],
            "pos": [100_000],
            "rsid": ["rs1"],
            "EA": ["A"],
            "OA": ["G"],
            "EAF": [0.2],
            "beta": [0.1],
            "se": [0.01],
            "pval": [1e-9],
            "N": [1000],
            "seqid": [seqid],
            "gene": ["G1"],
            "uniprot": ["P1"],
            "build": ["hg19"],
        }
    )
    out_df = pd.DataFrame(
        {
            "chromosome": ["1"],
            "base_pair_location": [200_000],
            "rsid": ["rs1"],
            "beta": [0.2],
            "standard_error": [0.02],
            "effect_allele_frequency": [0.25],
            "p_value": [1e-6],
            "effect_allele": ["A"],
            "other_allele": ["G"],
        }
    )
    fake_outcome = FakeOutcome(out_df)
    raw_dir = tmp_path / "raw_cis_sumstats_hg38"
    _write_raw_cis_hg38(raw_dir, seqid, exp_df)

    def fake_lift_position(chrom: str, pos: int):
        if pos == 1_100_000:
            return chrom, 2_100_000
        return chrom, pos + 1_000_000

    with patch.object(_mod, "cohort_dir", return_value=cohort_base), \
         patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "raw_cis_sumstats_hg38_dir", return_value=raw_dir), \
         patch.object(_mod, "OutcomeLookup", return_value=fake_outcome), \
         patch.object(_mod, "lift_position", side_effect=fake_lift_position):
        n_ok = _mod.extract_cohort_regions(cohort)

    assert n_ok == 1
    assert fake_outcome.calls == [("1", 1_000_001, 2_100_000)]


def test_empty_outcome_region_marked_failed_not_done(tmp_path):
    cohort = "ARIC_EA"
    seqid = "SeqId_fail"
    cohort_base = tmp_path / cohort
    _write_candidates(cohort_base, seqid, build="hg38")

    exp_df = pd.DataFrame(
        {
            "chrom": ["1"],
            "pos": [100_000],
            "rsid": ["rs1"],
            "EA": ["A"],
            "OA": ["G"],
            "EAF": [0.2],
            "beta": [0.1],
            "se": [0.01],
            "pval": [1e-9],
            "N": [1000],
            "seqid": [seqid],
            "gene": ["G1"],
            "uniprot": ["P1"],
            "build": ["hg38"],
        }
    )
    fake_outcome = FakeOutcome(pd.DataFrame())
    raw_dir = tmp_path / "raw_cis_sumstats_hg38"
    _write_raw_cis_hg38(raw_dir, seqid, exp_df)

    with patch.object(_mod, "cohort_dir", return_value=cohort_base), \
         patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "raw_cis_sumstats_hg38_dir", return_value=raw_dir), \
         patch.object(_mod, "OutcomeLookup", return_value=fake_outcome):
        n_ok = _mod.extract_cohort_regions(cohort)

    assert n_ok == 0
    state = json.loads((cohort_base / "_state_08_regions.json").read_text())
    assert state["done"] == []
    assert state["status"][seqid]["state"] == "failed"


def test_reads_raw_hg38_cache_and_retains_variant_absent_from_filtered_file(tmp_path):
    cohort = "Fenland"
    seqid = "SeqId_raw"
    cohort_base = tmp_path / cohort
    _write_candidates(cohort_base, seqid, build="hg38")

    raw_df = pd.DataFrame({
        "chrom": ["1", "1"],
        "pos": [100_000, 150_000],
        "rsid": ["rs_filtered", "rs_raw_only"],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
        "EAF": [0.2, 0.3],
        "beta": [0.1, 0.05],
        "se": [0.01, 0.02],
        "pval": [1e-9, 0.4],
        "N": [1000, 1000],
        "seqid": [seqid, seqid],
        "gene": ["G1", "G1"],
        "uniprot": ["P1", "P1"],
        "build": ["hg38", "hg38"],
    })
    raw_dir = tmp_path / "raw_cis_sumstats_hg38"
    _write_raw_cis_hg38(raw_dir, seqid, raw_df)
    out_df = pd.DataFrame({"chromosome": ["1"], "base_pair_location": [100_000]})
    fake_outcome = FakeOutcome(out_df)

    with patch.object(_mod, "cohort_dir", return_value=cohort_base), \
         patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "raw_cis_sumstats_hg38_dir", return_value=raw_dir), \
         patch.object(_mod, "_recover_raw_cis_hg38", side_effect=AssertionError("should not recover")), \
         patch.object(_mod, "OutcomeLookup", return_value=fake_outcome):
        n_ok = _mod.extract_cohort_regions(cohort)

    assert n_ok == 1
    written = pd.read_csv(tmp_path / "regions" / cohort / seqid / "exposure.tsv", sep="\t")
    assert set(written["rsid"]) == {"rs_filtered", "rs_raw_only"}


def test_missing_raw_cache_uses_one_protein_recovery(tmp_path):
    cohort = "UKB_PPP"
    seqid = "SeqId_recover"
    cohort_base = tmp_path / cohort
    _write_candidates(cohort_base, seqid, build="hg19")

    recovered_df = pd.DataFrame({
        "chrom": ["1"],
        "pos": [100_000],
        "rsid": ["rs_recovered"],
        "EA": ["A"],
        "OA": ["G"],
        "EAF": [0.2],
        "beta": [0.1],
        "se": [0.01],
        "pval": [0.2],
        "N": [1000],
        "seqid": [seqid],
        "gene": ["G1"],
        "uniprot": ["P1"],
        "build": ["hg38"],
    })
    out_df = pd.DataFrame({"chromosome": ["1"], "base_pair_location": [100_000]})
    fake_outcome = FakeOutcome(out_df)

    with patch.object(_mod, "cohort_dir", return_value=cohort_base), \
         patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "raw_cis_sumstats_hg38_dir", return_value=tmp_path / "missing_raw"), \
         patch.object(_mod, "_recover_raw_cis_hg38", return_value=recovered_df) as recover, \
         patch.object(_mod, "OutcomeLookup", return_value=fake_outcome), \
         patch.object(_mod, "lift_position", side_effect=lambda chrom, pos: (chrom, pos)):
        n_ok = _mod.extract_cohort_regions(cohort)

    assert n_ok == 1
    recover.assert_called_once()
    written = pd.read_csv(tmp_path / "regions" / cohort / seqid / "exposure.tsv", sep="\t")
    assert set(written["rsid"]) == {"rs_recovered"}


def test_recover_raw_cis_hg38_removes_partial_liftover_output(tmp_path):
    cohort = "UKB_PPP"
    seqid = "SeqId_partial"
    raw_dir = tmp_path / "raw_cis_sumstats_hg38"
    native_path = tmp_path / "native.tsv.gz"
    raw = pd.DataFrame({
        "chrom": ["1"],
        "pos": [100_000],
        "rsid": ["rs1"],
        "EA": ["A"],
        "OA": ["G"],
        "EAF": [0.2],
        "beta": [0.1],
        "se": [0.01],
        "pval": [0.2],
        "N": [1000],
    })

    class FakeLiftover:
        @staticmethod
        def lift_sumstats_file(_cohort, _native_path, out_path):
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text("chrom\tpos\n1\t100\n")
            raise OSError("simulated liftover write failure")

    with patch.object(_mod, "raw_cis_sumstats_hg38_dir", return_value=raw_dir), \
         patch.object(_mod, "_extract_raw_native_region", return_value=raw), \
         patch.object(_mod, "write_raw_cis_cache", return_value=native_path), \
         patch.object(_mod.importlib, "import_module", return_value=FakeLiftover):
        with pytest.raises(OSError, match="simulated liftover"):
            _mod._recover_raw_cis_hg38(cohort, seqid, "1", 100_000, "hg19")

    assert not (raw_dir / f"{seqid}.tsv.gz").exists()
    assert not list(raw_dir.glob("*.tmp.tsv.gz"))


def test_extract_regions_materializes_same_rsid_different_positions_without_matching(tmp_path):
    cohort = "Fenland"
    seqid = "SeqId_same_rsid_diff_pos"
    cohort_base = tmp_path / cohort
    _write_candidates(cohort_base, seqid, build="hg38")

    exp_df = pd.DataFrame({
        "chrom": ["1"],
        "pos": [100_000],
        "rsid": ["rs_shared"],
        "EA": ["A"],
        "OA": ["G"],
        "EAF": [0.2],
        "beta": [0.1],
        "se": [0.01],
        "pval": [1e-9],
        "N": [1000],
        "seqid": [seqid],
        "gene": ["G1"],
        "uniprot": ["P1"],
        "build": ["hg38"],
    })
    raw_dir = tmp_path / "raw_cis_sumstats_hg38"
    _write_raw_cis_hg38(raw_dir, seqid, exp_df)
    out_df = pd.DataFrame({
        "chromosome": ["1"],
        "base_pair_location": [101_000],
        "rsid": ["rs_shared"],
        "effect_allele": ["A"],
        "other_allele": ["G"],
        "beta": [0.2],
        "standard_error": [0.02],
        "effect_allele_frequency": [0.3],
        "p_value": [1e-6],
    })
    fake_outcome = FakeOutcome(out_df)

    with patch.object(_mod, "cohort_dir", return_value=cohort_base), \
         patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "raw_cis_sumstats_hg38_dir", return_value=raw_dir), \
         patch.object(_mod, "OutcomeLookup", return_value=fake_outcome):
        assert _mod.extract_cohort_regions(cohort) == 1

    region_dir = tmp_path / "regions" / cohort / seqid
    exposure = pd.read_csv(region_dir / "exposure.tsv", sep="\t")
    outcome = pd.read_csv(region_dir / "outcome.tsv", sep="\t")
    assert exposure["rsid"].iloc[0] == "rs_shared"
    assert outcome["rsid"].iloc[0] == "rs_shared"
    assert exposure["pos"].iloc[0] == 100_000
    assert outcome["base_pair_location"].iloc[0] == 101_000


def test_extract_regions_preserves_recovered_exposure_columns_for_downstream_alignment(tmp_path):
    cohort = "UKB_PPP"
    seqid = "SeqId_lowercase_recovered"
    cohort_base = tmp_path / cohort
    _write_candidates(cohort_base, seqid, build="hg38")
    recovered_df = pd.DataFrame({
        "chrom": ["1"],
        "pos": [100_000],
        "rsid": ["rs_lower"],
        "EA": ["a"],
        "OA": ["g"],
        "EAF": [0.2],
        "beta": [0.1],
        "se": [0.01],
        "pval": [0.2],
        "N": [1000],
        "seqid": [seqid],
        "gene": ["G1"],
        "uniprot": ["P1"],
        "build": ["hg38"],
    })
    out_df = pd.DataFrame({
        "chromosome": ["1"],
        "base_pair_location": [100_000],
        "rsid": ["rs_lower"],
        "effect_allele": ["A"],
        "other_allele": ["G"],
        "beta": [0.2],
        "standard_error": [0.02],
        "effect_allele_frequency": [0.3],
        "p_value": [1e-6],
    })

    with patch.object(_mod, "cohort_dir", return_value=cohort_base), \
         patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "raw_cis_sumstats_hg38_dir", return_value=tmp_path / "missing_raw"), \
         patch.object(_mod, "_recover_raw_cis_hg38", return_value=recovered_df), \
         patch.object(_mod, "OutcomeLookup", return_value=FakeOutcome(out_df)):
        assert _mod.extract_cohort_regions(cohort) == 1

    exposure = pd.read_csv(tmp_path / "regions" / cohort / seqid / "exposure.tsv", sep="\t")
    required = {"chrom", "pos", "rsid", "EA", "OA", "beta", "se", "EAF", "N"}
    assert required.issubset(exposure.columns)
    assert exposure["EA"].iloc[0] == "a"
    assert exposure["OA"].iloc[0] == "g"


def test_extract_regions_preserves_duplicate_outcome_positions_for_coloc_alignment(tmp_path):
    cohort = "Fenland"
    seqid = "SeqId_duplicate_outcome"
    cohort_base = tmp_path / cohort
    _write_candidates(cohort_base, seqid, build="hg38")

    exp_df = pd.DataFrame({
        "chrom": ["1"],
        "pos": [100_000],
        "rsid": ["rs_dup"],
        "EA": ["A"],
        "OA": ["G"],
        "EAF": [0.2],
        "beta": [0.1],
        "se": [0.01],
        "pval": [1e-9],
        "N": [1000],
        "seqid": [seqid],
        "gene": ["G1"],
        "uniprot": ["P1"],
        "build": ["hg38"],
    })
    raw_dir = tmp_path / "raw_cis_sumstats_hg38"
    _write_raw_cis_hg38(raw_dir, seqid, exp_df)
    out_df = pd.DataFrame({
        "chromosome": ["1", "1"],
        "base_pair_location": [100_000, 100_000],
        "rsid": ["rs_dup_a", "rs_dup_b"],
        "effect_allele": ["A", "C"],
        "other_allele": ["G", "T"],
        "beta": [0.2, 0.9],
        "standard_error": [0.02, 0.03],
        "effect_allele_frequency": [0.3, 0.4],
        "p_value": [1e-6, 2e-6],
    })

    with patch.object(_mod, "cohort_dir", return_value=cohort_base), \
         patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "raw_cis_sumstats_hg38_dir", return_value=raw_dir), \
         patch.object(_mod, "OutcomeLookup", return_value=FakeOutcome(out_df)):
        assert _mod.extract_cohort_regions(cohort) == 1

    outcome = pd.read_csv(tmp_path / "regions" / cohort / seqid / "outcome.tsv", sep="\t")
    assert len(outcome) == 2
    assert outcome["base_pair_location"].tolist() == [100_000, 100_000]
    assert set(zip(outcome["effect_allele"], outcome["other_allele"])) == {("A", "G"), ("C", "T")}
