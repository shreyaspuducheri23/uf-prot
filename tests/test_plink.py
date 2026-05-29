"""Tests for scripts.lib.plink helpers used by clumping/proxy workflows."""
import csv
import textwrap
from pathlib import Path
from unittest.mock import patch

import importlib
import pandas as pd
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

    proxy_map = _plink.find_proxies([target], r2_threshold=0.8)

    assert target in proxy_map, f"No proxy returned for {target}"
    proxy_rsid, proxy_r2 = proxy_map[target]

    best_r2 = max(expected.values())
    best_candidates = sorted(rsid for rsid, r2 in expected.items() if r2 == best_r2)
    expected_proxy = best_candidates[0]

    assert proxy_rsid == expected_proxy
    assert proxy_r2 == pytest.approx(best_r2)


def test_bim_pos_to_rsid_returns_mapping(tmp_path):
    bim = tmp_path / "ref.bim"
    bim.write_text(
        "1\trs111\t0\t100000\tA\tG\n"
        "1\trs222\t0\t200000\tC\tT\n"
        "2\trs333\t0\t500000\tA\tC\n"
    )
    bfile = tmp_path / "ref"
    _plink._bim_pos_to_rsid.cache_clear()
    mapping = _plink._bim_pos_to_rsid(bfile)
    assert mapping[("1", 100000)] == "rs111"
    assert mapping[("2", 500000)] == "rs333"
    assert ("1", 999999) not in mapping


def test_clump_annotates_missing_rsids_from_bim(tmp_path):
    bim = tmp_path / "ref.bim"
    bim.write_text(
        "22\trs999\t0\t25212564\tA\tG\n"
        "22\trs888\t0\t25222564\tC\tT\n"
    )
    bfile = tmp_path / "ref"

    sumstats = pd.DataFrame({
        "seqid":  ["SeqId_TEST"] * 2,
        "chrom":  ["22", "22"],
        "pos":    [25_212_564, 25_222_564],
        "rsid":   [".", "."],
        "pval":   [1e-10, 1e-9],
        "EA":     ["A", "C"],
        "OA":     ["G", "T"],
        "beta":   [0.5, 0.4],
        "se":     [0.05, 0.05],
        "EAF":    [0.3, 0.4],
        "N":      [30000, 30000],
    })

    clumps_content = "SNP P\nrs999 1e-10\n"

    def fake_run(cmd, cwd=None):
        out_prefix = Path(cmd[cmd.index("--out") + 1])
        out_prefix.with_suffix(".clumps").write_text(clumps_content)
        return _FakeCompleted()

    _plink._bim_pos_to_rsid.cache_clear()
    original_run = _plink._run
    _plink._run = fake_run
    try:
        result = _plink.clump(sumstats, "SeqId_TEST", bfile=bfile)
    finally:
        _plink._run = original_run
        _plink._bim_pos_to_rsid.cache_clear()

    assert len(result) == 1
    assert result["rsid"].iloc[0] == "rs999"


def test_clump_drops_variants_not_in_bim(tmp_path):
    bim = tmp_path / "ref.bim"
    bim.write_text("22\trs999\t0\t25212564\tA\tG\n")
    bfile = tmp_path / "ref"

    sumstats = pd.DataFrame({
        "seqid":  ["SeqId_TEST"],
        "chrom":  ["22"],
        "pos":    [99999999],  # not in bim
        "rsid":   ["."],
        "pval":   [1e-10],
        "EA":     ["A"],
        "OA":     ["G"],
        "beta":   [0.5],
        "se":     [0.05],
        "EAF":    [0.3],
        "N":      [30000],
    })

    _plink._bim_pos_to_rsid.cache_clear()
    try:
        result = _plink.clump(sumstats, "SeqId_TEST", bfile=bfile)
    finally:
        _plink._bim_pos_to_rsid.cache_clear()

    assert result.empty


# ---------------------------------------------------------------------------
# r_square_matrix regression tests
# ---------------------------------------------------------------------------

def test_r_square_matrix_uses_r_unphased(monkeypatch):
    """r_square_matrix must call --r-unphased, not --r (removed in plink2)."""
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        raise FileNotFoundError("no matrix")  # stop after the plink call

    monkeypatch.setattr(_plink, "_run", fake_run)
    with pytest.raises((FileNotFoundError, Exception)):
        _plink.r_square_matrix(["rs1260326"])

    assert calls, "plink was not called"
    assert "--r-unphased" in calls[0], \
        f"Expected --r-unphased in plink cmd, got: {calls[0]}"
    assert "--r" not in [a for a in calls[0] if a == "--r"], \
        "Deprecated --r flag must not be used"


def test_r_square_matrix_reads_vars_file_for_snp_order(tmp_path, monkeypatch):
    """r_square_matrix must use the .vars file for index/columns, not the input list."""
    matrix_content = "1.0\t0.3\n0.3\t1.0\n"
    vars_content = "rs780094\nrs1260326\n"  # different order from input

    def fake_run(cmd, **kwargs):
        out_prefix = cmd[cmd.index("--out") + 1]
        mf = Path(out_prefix + ".unphased.vcor1")
        mf.write_text(matrix_content)
        Path(str(mf) + ".vars").write_text(vars_content)

    monkeypatch.setattr(_plink, "_run", fake_run)
    monkeypatch.setattr(_plink, "LD_REF_PREFIX", tmp_path / "dummy")

    result = _plink.r_square_matrix(["rs1260326", "rs780094"])
    assert list(result.index)   == ["rs780094", "rs1260326"]
    assert list(result.columns) == ["rs780094", "rs1260326"]


def test_r_square_matrix_uses_maf_filter(monkeypatch):
    """r_square_matrix must pass --maf 0.01 to plink to exclude monomorphic variants."""
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        raise FileNotFoundError("no matrix")

    monkeypatch.setattr(_plink, "_run", fake_run)
    with pytest.raises((FileNotFoundError, Exception)):
        _plink.r_square_matrix(["rs1260326"])

    assert calls, "plink was not called"
    assert "--maf" in calls[0], f"Expected --maf flag in plink cmd: {calls[0]}"
    maf_val = calls[0][calls[0].index("--maf") + 1]
    assert float(maf_val) >= 0.01, f"MAF threshold should be >= 0.01, got {maf_val}"


def test_r_square_matrix_drops_residual_nan_snps(tmp_path, monkeypatch):
    """Residual NaN after --maf (rare clustered variants) must be dropped, not raised.

    Real pattern: NaN is confined to the problem-SNP rows/columns; clean SNPs
    have no NaN in their rows.  Dropping the flagged rows/columns yields a clean matrix.
    """
    # rs_a and rs_b have NaN with each other (intra-cluster rare variants)
    # but normal correlations with rs_clean — matching the real observed pattern.
    matrix_content = (
        "1.0\t0.3\t0.2\t0.1\n"   # rs_clean1 row — all OK
        "0.3\t1.0\t0.4\t0.2\n"   # rs_clean2 row — all OK
        "0.2\t0.4\t1.0\tnan\n"   # rs_a row — NaN with rs_b only
        "0.1\t0.2\tnan\t1.0\n"   # rs_b row — NaN with rs_a only
    )
    vars_content = "rs_clean1\nrs_clean2\nrs_a\nrs_b\n"

    def fake_run(cmd, **kwargs):
        out_prefix = cmd[cmd.index("--out") + 1]
        mf = Path(out_prefix + ".unphased.vcor1")
        mf.write_text(matrix_content)
        Path(str(mf) + ".vars").write_text(vars_content)

    monkeypatch.setattr(_plink, "_run", fake_run)
    monkeypatch.setattr(_plink, "LD_REF_PREFIX", tmp_path / "dummy")

    result = _plink.r_square_matrix(["rs_clean1", "rs_clean2", "rs_a", "rs_b"])
    # Problem SNPs dropped; clean SNPs survive
    assert "rs_a" not in result.index, "rs_a must be dropped (has NaN)"
    assert "rs_b" not in result.index, "rs_b must be dropped (has NaN)"
    assert list(result.index) == ["rs_clean1", "rs_clean2"]
    assert not result.isna().any().any(), "No NaN should remain"
