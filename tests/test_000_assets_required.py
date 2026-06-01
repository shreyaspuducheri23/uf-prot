"""Fail-fast checks for mandatory local assets/tooling required by this test suite."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from scripts.lib.paths import (
    ARIC_EA_DIR,
    ARIC_SEQID,
    CHAIN_HG19_TO_HG38,
    KIM_GWAS,
    LD_REF_PREFIX,
)


def _assert_exists(path: Path, label: str) -> None:
    if not path.exists():
        pytest.fail(f"Missing required {label}: {path}")


def test_required_local_assets_present() -> None:
    _assert_exists(ARIC_SEQID, "ARIC seqid index")
    _assert_exists(ARIC_EA_DIR, "ARIC EA directory")
    _assert_exists(CHAIN_HG19_TO_HG38, "hg19->hg38 liftover chain")
    _assert_exists(KIM_GWAS, "Kim GWAS bgz file")
    _assert_exists(Path(f"{KIM_GWAS}.tbi"), "Kim GWAS tabix index")
    _assert_exists(Path(f"{LD_REF_PREFIX}.bed"), "LD reference .bed")
    _assert_exists(Path(f"{LD_REF_PREFIX}.bim"), "LD reference .bim")
    _assert_exists(Path(f"{LD_REF_PREFIX}.fam"), "LD reference .fam")
    _assert_exists(Path(f"{LD_REF_PREFIX}.snplist"), "LD reference snplist")


def test_plink2_on_path() -> None:
    if shutil.which("plink2") is None:
        pytest.fail("plink2 is not on $PATH — clumping and proxy steps will silently produce 0 instruments")


def test_required_executables_available() -> None:

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
