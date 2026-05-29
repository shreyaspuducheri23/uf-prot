"""Tests for scripts.08_coloc.sharepro."""
import importlib
import json
import subprocess
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

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


# ---------------------------------------------------------------------------
# Interface regression tests
# ---------------------------------------------------------------------------

def _make_region(tmp_path: Path, snps: list[str]) -> Path:
    """Minimal region dir with matching exposure/outcome rsids."""
    region = tmp_path / "region"
    region.mkdir(exist_ok=True)
    exp = pd.DataFrame({
        "rsid": snps,
        "beta": [0.1] * len(snps),
        "se":   [0.05] * len(snps),
        "N":    [1000] * len(snps),
        "pos":  list(range(100, 100 + len(snps))),
    })
    out = pd.DataFrame({
        "rsid":                    snps,
        "beta":                    [0.1] * len(snps),
        "standard_error":          [0.05] * len(snps),
        "effect_allele_frequency": [0.3] * len(snps),
        "base_pair_location":      list(range(100, 100 + len(snps))),
    })
    exp.to_csv(region / "exposure.tsv", sep="\t", index=False)
    out.to_csv(region / "outcome.tsv",  sep="\t", index=False)
    return region


def _identity_ld(snp_list):
    n = len(snp_list)
    return pd.DataFrame(
        [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)],
        index=snp_list, columns=snp_list,
    )


def test_sharepro_script_path_exists():
    """The SharePro entry-point must exist on disk (catches wrong path)."""
    from scripts.lib.paths import SHAREPRO_SCRIPT
    assert SHAREPRO_SCRIPT.exists(), (
        f"SharePro script not found: {SHAREPRO_SCRIPT}\n"
        "Expected: tools/SharePro_coloc/src/SharePro/sharepro_coloc.py"
    )


def test_build_bse_input_columns():
    """Input builder must produce SNP/BETA/SE/N — not snp/z/N."""
    df = pd.DataFrame({
        "rsid": ["rs1", "rs2"],
        "beta": [0.1, -0.2],
        "se":   [0.05, 0.04],
        "N":    [1000, 1000],
    })
    result = _mod.build_bse_input(df, n=1000)
    assert list(result.columns) == ["SNP", "BETA", "SE", "N"], (
        f"Expected SNP/BETA/SE/N, got {list(result.columns)}"
    )
    assert "z" not in result.columns, "z-score column must not be present"


def test_run_sharepro_uses_z_and_save_flags(tmp_path, monkeypatch):
    """run_sharepro must use --z and --save flags (not --z1/--z2/--out)."""
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(list(cmd))
        save_prefix = cmd[cmd.index("--save") + 1]
        Path(save_prefix + ".sharepro.txt").write_text("cs\tshare\tvariantProb\n")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    snps = [f"rs{i}" for i in range(10)]
    monkeypatch.setattr(_mod.subprocess, "run", fake_run)
    monkeypatch.setattr(_mod, "r_square_matrix", lambda sl: _identity_ld(sl))

    _mod.run_sharepro(_make_region(tmp_path, snps), "test_seqid", N_out=5000)

    assert calls, "subprocess.run was never called"
    cmd = calls[0]
    assert "--z" in cmd,    f"Expected --z flag, got: {cmd}"
    assert "--save" in cmd, f"Expected --save flag, got: {cmd}"
    assert "--z1" not in cmd,  "Deprecated --z1 must not appear"
    assert "--z2" not in cmd,  "Deprecated --z2 must not appear"
    assert "--out" not in cmd, "Deprecated --out must not appear"


def test_run_sharepro_parses_sharepro_txt_not_json(tmp_path, monkeypatch):
    """run_sharepro must read .sharepro.txt 'share' column, not a JSON PP.H4 key."""
    def fake_run(cmd, **kwargs):
        save_prefix = cmd[cmd.index("--save") + 1]
        Path(save_prefix + ".sharepro.txt").write_text(
            "cs\tshare\tvariantProb\n"
            "rs1/rs2\t0.92\t0.6/0.4\n"
        )
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    snps = [f"rs{i}" for i in range(10)]
    monkeypatch.setattr(_mod.subprocess, "run", fake_run)
    monkeypatch.setattr(_mod, "r_square_matrix", lambda sl: _identity_ld(sl))

    result, err = _mod.run_sharepro(_make_region(tmp_path, snps), "test_seqid", N_out=5000)

    assert result is not None, f"Expected success, got failure: {err}"
    assert abs(result["PP_H4"] - 0.92) < 1e-6, f"Expected PP_H4=0.92, got {result['PP_H4']}"
    assert result["coloc_positive"] is True


def test_run_sharepro_empty_output_gives_pp_h4_zero(tmp_path, monkeypatch):
    """Header-only .sharepro.txt (no effect groups found) → PP_H4 = 0, not an error."""
    def fake_run(cmd, **kwargs):
        save_prefix = cmd[cmd.index("--save") + 1]
        Path(save_prefix + ".sharepro.txt").write_text("cs\tshare\tvariantProb\n")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    snps = [f"rs{i}" for i in range(10)]
    monkeypatch.setattr(_mod.subprocess, "run", fake_run)
    monkeypatch.setattr(_mod, "r_square_matrix", lambda sl: _identity_ld(sl))

    result, err = _mod.run_sharepro(_make_region(tmp_path, snps), "test_seqid", N_out=5000)

    assert result is not None, f"Expected success dict, got failure: {err}"
    assert result["PP_H4"] == 0.0
    assert result["coloc_positive"] is False


def test_run_sharepro_ld_snp_subset_aligns_with_sumstats(tmp_path, monkeypatch):
    """When plink drops SNPs from the LD matrix, sumstats must be restricted to the
    same subset — not passed at full size (which causes SharePro AssertionError)."""
    snps = [f"rs{i}" for i in range(10)]
    # Simulate plink dropping the last 3 SNPs (only 7 survive in LD)
    surviving = snps[:7]

    captured_row_counts = []

    def fake_run(cmd, **kwargs):
        # Read the exp_bse file *now*, while the tmpdir is still alive
        exp_bse_path = cmd[cmd.index("--z") + 1]
        n_rows = len(pd.read_csv(exp_bse_path, sep="\t"))
        captured_row_counts.append(n_rows)
        save_prefix = cmd[cmd.index("--save") + 1]
        Path(save_prefix + ".sharepro.txt").write_text("cs\tshare\tvariantProb\n")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(_mod.subprocess, "run", fake_run)
    # LD stub returns only 7×7 matrix
    monkeypatch.setattr(_mod, "r_square_matrix", lambda sl: _identity_ld(surviving))

    result, err = _mod.run_sharepro(_make_region(tmp_path, snps), "test_seqid", N_out=5000)

    assert captured_row_counts, "subprocess.run was never called"
    assert captured_row_counts[0] == 7, (
        f"exp_bse must have 7 rows (matching LD), got {captured_row_counts[0]}"
    )
