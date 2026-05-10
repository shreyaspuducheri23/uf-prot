"""Tests for scripts.08_coloc.extract_regions."""
import importlib
import json
from pathlib import Path
from unittest.mock import patch

import pandas as pd

_mod = importlib.import_module("scripts.08_coloc.extract_regions")


def _write_candidates(base: Path, seqid: str, build: str = "hg19") -> None:
    base.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"seqid": [seqid], "fdr_pass": [True]}).to_csv(
        base / "mr_results.tsv", sep="\t", index=False
    )
    pd.DataFrame({"seqid": [seqid], "chrom": ["1"], "tss": [100_000], "build": [build]}).to_csv(
        base / "protein_index.tsv", sep="\t", index=False
    )


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

    def fake_lift_position(chrom: str, pos: int):
        if pos == 1_100_000:
            return chrom, 2_100_000
        return chrom, pos + 1_000_000

    with patch.object(_mod, "cohort_dir", return_value=cohort_base), \
         patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "extract_aric_region", return_value=exp_df), \
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

    with patch.object(_mod, "cohort_dir", return_value=cohort_base), \
         patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "extract_aric_region", return_value=exp_df), \
         patch.object(_mod, "OutcomeLookup", return_value=fake_outcome):
        n_ok = _mod.extract_cohort_regions(cohort)

    assert n_ok == 0
    state = json.loads((cohort_base / "_state_08_regions.json").read_text())
    assert state["done"] == []
    assert state["status"][seqid]["state"] == "failed"
