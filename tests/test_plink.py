"""Tests for scripts.lib.plink helpers used by clumping/proxy workflows."""
import csv
import subprocess
from pathlib import Path

import importlib
import pytest

_plink = importlib.import_module("scripts.lib.plink")


class _FakeCompleted:
    def __init__(self, stdout: str = "", stderr: str = ""):
        self.returncode = 0
        self.stdout = stdout
        self.stderr = stderr


def test_in_phase_allele_map_parses_plink_output():
    def fake_run(_cmd, cwd=None):
        return _FakeCompleted(
            stderr=(
                "--ld rs_target rs_proxy:\n\n"
                "rs_target alleles:\n"
                "  MAJOR = G\n"
                "  MINOR = A\n"
                "rs_proxy alleles:\n"
                "  MAJOR = T\n"
                "  MINOR = C\n"
                "  Major alleles are in phase with each other.\n"
            )
        )

    original = _plink._run
    _plink._run = fake_run
    try:
        mapping = _plink.in_phase_allele_map("rs_target", "rs_proxy")
    finally:
        _plink._run = original

    assert mapping == {"G": "T", "A": "C"}


def test_in_phase_allele_map_returns_none_when_unparseable():
    def fake_run(_cmd, cwd=None):
        return _FakeCompleted(stderr="No valid LD statistics computed\n")

    original = _plink._run
    _plink._run = fake_run
    try:
        mapping = _plink.in_phase_allele_map("rs_target", "rs_proxy")
    finally:
        _plink._run = original

    assert mapping is None


def test_find_proxies_picks_highest_r2_per_target():
    def fake_run(cmd, cwd=None):
        out_prefix = Path(cmd[cmd.index("--out") + 1])
        vcor_file = out_prefix.with_suffix(".vcor")
        vcor_file.write_text(
            "#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n"
            "1\t100\trsT\t1\t200\trsA\t0.81\n"
            "1\t100\trsT\t1\t300\trsB\t0.95\n"
            "1\t100\trsT\t1\t400\trsC\t0.95\n"
            "1\t100\trsT\t1\t500\trsT\t1.00\n"
        )
        return _FakeCompleted()

    original = _plink._run
    _plink._run = fake_run
    try:
        proxies = _plink.find_proxies(["rsT"])
    finally:
        _plink._run = original

    assert proxies == {"rsT": ("rsB", 0.95)}


def test_find_proxies_matches_ground_truth_for_rs10757278():
    target = "rs10757278"
    ground_truth_path = Path(__file__).resolve().parent / "proxy_test.txt"
    assert ground_truth_path.exists(), f"Missing ground-truth file: {ground_truth_path}"

    # Prefer PATH plink2; fallback to known local binary in this workstation.
    plink_exec = _plink._PLINK2
    try:
        probe = subprocess.run([plink_exec, "--help"], capture_output=True, text=True)
    except OSError:
        probe = None
    if probe is None or probe.returncode not in {0, 1}:
        fallback = Path("/Users/spuduch/Research/MR_IA/plink2_mac_arm64_20260228/plink2")
        if not fallback.exists():
            pytest.fail(f"plink2 unavailable on PATH and fallback missing: {fallback}")
        plink_exec = str(fallback)

    expected: dict[str, float] = {}
    with ground_truth_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rsid = str(row["RS_Number"]).strip()
            if rsid == target:
                continue
            try:
                r2 = float(row["R2"])
            except (TypeError, ValueError):
                continue
            if r2 >= 0.8:
                expected[rsid] = r2
    assert expected, f"No >=0.8 proxy candidates in {ground_truth_path}"

    original_exec = _plink._PLINK2
    _plink._PLINK2 = plink_exec
    try:
        proxy_map = _plink.find_proxies([target], r2_threshold=0.8)
    finally:
        _plink._PLINK2 = original_exec

    assert target in proxy_map, f"No proxy returned for {target}"
    proxy_rsid, proxy_r2 = proxy_map[target]

    best_r2 = max(expected.values())
    best_candidates = sorted(rsid for rsid, r2 in expected.items() if r2 == best_r2)
    expected_proxy = best_candidates[0]

    assert proxy_rsid == expected_proxy
    assert proxy_r2 == pytest.approx(best_r2)
