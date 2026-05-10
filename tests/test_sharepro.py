"""Tests for scripts.08_coloc.sharepro."""
import importlib
import json
from pathlib import Path
from unittest.mock import patch

import pandas as pd

_mod = importlib.import_module("scripts.08_coloc.sharepro")


def _write_region(region_dir: Path, with_n: bool) -> None:
    region_dir.mkdir(parents=True, exist_ok=True)
    exp = pd.DataFrame(
        {
            "rsid": [f"rs{i}" for i in range(1, 7)],
            "beta": [0.2] * 6,
            "se": [0.05] * 6,
        }
    )
    if with_n:
        exp["N"] = [1000] * 6
    exp.to_csv(region_dir / "exposure.tsv", sep="\t", index=False)

    out = pd.DataFrame(
        {
            "rsid": [f"rs{i}" for i in range(1, 7)],
            "beta": [0.1] * 6,
            "standard_error": [0.02] * 6,
            "effect_allele_frequency": [0.3] * 6,
        }
    )
    out.to_csv(region_dir / "outcome.tsv", sep="\t", index=False)


def test_run_sharepro_rejects_missing_n_exp(tmp_path):
    region_dir = tmp_path / "SeqId_A"
    _write_region(region_dir, with_n=False)
    result, reason = _mod.run_sharepro(region_dir, "SeqId_A", N_out=1000)
    assert result is None
    assert reason == "invalid_or_missing_N_exp"


def test_failed_checkpoint_skipped_unless_retry_failed(tmp_path):
    cohort = "ARIC_EA"
    seqid = "SeqId_B"
    region_base = tmp_path / "regions" / cohort / seqid
    _write_region(region_base, with_n=False)

    cohort_state = tmp_path / cohort
    cohort_state.mkdir(parents=True, exist_ok=True)

    with patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "cohort_dir", return_value=cohort_state):
        _mod.run_cohort_sharepro(cohort)

    state_path = cohort_state / "_state_08_sharepro.json"
    state = json.loads(state_path.read_text())
    assert state["status"][seqid]["state"] == "failed"

    with patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "cohort_dir", return_value=cohort_state), \
         patch.object(_mod, "run_sharepro", side_effect=AssertionError("should not rerun failed")):
        results = _mod.run_cohort_sharepro(cohort, retry_failed=False)
    assert results == []

    with patch.object(_mod, "COLOC_REGIONS_DIR", tmp_path / "regions"), \
         patch.object(_mod, "cohort_dir", return_value=cohort_state), \
         patch.object(
             _mod,
             "run_sharepro",
             return_value=(
                 {
                     "seqid": seqid,
                     "n_snps": 6,
                     "PP_H4": 0.91,
                     "coloc_positive": True,
                     "raw": {},
                 },
                 None,
             ),
         ):
        retried = _mod.run_cohort_sharepro(cohort, retry_failed=True)
    assert len(retried) == 1
