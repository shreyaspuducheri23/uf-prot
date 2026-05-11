"""Fail-fast checks for mandatory local assets/tooling required by this test suite."""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from scripts.lib.paths import (
    ARIC_EA_DIR,
    ARIC_SEQID,
    CHAIN_HG19_TO_HG38,
    DECODE_URLS,
    KIM_GWAS,
    LD_REF_PREFIX,
)

PLINK2_CMD = "plink2"
PLINK2_FALLBACK = Path("/Users/spuduch/Research/MR_IA/plink2_mac_arm64_20260228/plink2")


def _assert_exists(path: Path, label: str) -> None:
    if not path.exists():
        pytest.fail(f"Missing required {label}: {path}")

def _is_runnable(exec_path: Path) -> bool:
    try:
        probe = subprocess.run([str(exec_path), "--help"], capture_output=True, text=True)
    except OSError:
        return False
    return probe.returncode in {0, 1}

def _cmd_is_runnable(cmd: str) -> bool:
    try:
        probe = subprocess.run([cmd, "--help"], capture_output=True, text=True)
    except OSError:
        return False
    return probe.returncode in {0, 1}


def test_required_local_assets_present() -> None:
    _assert_exists(ARIC_SEQID, "ARIC seqid index")
    _assert_exists(ARIC_EA_DIR, "ARIC EA directory")
    _assert_exists(DECODE_URLS, "deCODE bulk URL list")
    _assert_exists(CHAIN_HG19_TO_HG38, "hg19->hg38 liftover chain")
    _assert_exists(KIM_GWAS, "Kim GWAS bgz file")
    _assert_exists(Path(f"{KIM_GWAS}.tbi"), "Kim GWAS tabix index")
    _assert_exists(Path(f"{LD_REF_PREFIX}.bed"), "LD reference .bed")
    _assert_exists(Path(f"{LD_REF_PREFIX}.bim"), "LD reference .bim")
    _assert_exists(Path(f"{LD_REF_PREFIX}.fam"), "LD reference .fam")
    _assert_exists(Path(f"{LD_REF_PREFIX}.snplist"), "LD reference snplist")


def test_required_executables_available() -> None:
    if not _cmd_is_runnable(PLINK2_CMD):
        _assert_exists(PLINK2_FALLBACK, "PLINK2 fallback binary")
        if not _is_runnable(PLINK2_FALLBACK):
            pytest.fail(f"PLINK2 is not runnable on PATH and fallback is also not runnable: {PLINK2_FALLBACK}")

    r_probe = subprocess.run(["Rscript", "--version"], capture_output=True, text=True)
    if r_probe.returncode != 0:
        pytest.fail("Rscript is required but not runnable")

    two_sample_mr_probe = subprocess.run(
        ["Rscript", "-e", 'quit(save="no", status=ifelse(requireNamespace("TwoSampleMR", quietly=TRUE),0,42))'],
        capture_output=True,
        text=True,
    )
    if two_sample_mr_probe.returncode == 42:
        pytest.fail("R package TwoSampleMR is required but not installed")
    if two_sample_mr_probe.returncode != 0:
        pytest.fail("Failed while probing TwoSampleMR availability via Rscript")

    gwas_backend_probe = subprocess.run(
        [
            "Rscript",
            "-e",
            'quit(save="no", status=ifelse(requireNamespace("GWASBrewer", quietly=TRUE),0,42))',
        ],
        capture_output=True,
        text=True,
    )
    if gwas_backend_probe.returncode == 42:
        pytest.fail("R package GWASBrewer is required for oracle simulation tests")
    if gwas_backend_probe.returncode != 0:
        pytest.fail("Failed while probing GWASBrewer availability via Rscript")


def test_repo_python_virtualenv_present() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    py = repo_root / ".venv" / "bin" / "python"
    _assert_exists(py, "repo virtualenv Python (.venv/bin/python)")
