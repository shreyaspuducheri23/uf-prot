"""Tests for scripts.05_harmonise.harmonise (_join_outcome logic)."""
import importlib
import pandas as pd
import numpy as np
import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
import subprocess

_harmonise_mod = importlib.import_module("scripts.05_harmonise.harmonise")
_join_outcome = _harmonise_mod._join_outcome


def _instrument_df(**kwargs) -> pd.DataFrame:
    defaults = {
        "chrom_hg38": ["22"],
        "pos_hg38": [25_212_564],
        "rsid": ["rs123"],
        "chrom": ["22"],
        "pos": [25_212_564],
        "EA": ["A"], "OA": ["G"],
        "EAF": [0.3], "beta": [0.5], "se": [0.05], "pval": [1e-9],
        "N": [7213], "seqid": ["SeqId_TEST"],
    }
    defaults.update(kwargs)
    return pd.DataFrame(defaults)


def _outcome_df(chrom="22", pos=25_212_564, ea="A", oa="G") -> pd.DataFrame:
    return pd.DataFrame({
        "chromosome": [chrom],
        "base_pair_location": [pos],
        "effect_allele": [ea],
        "other_allele": [oa],
        "beta": [0.2],
        "standard_error": [0.02],
        "effect_allele_frequency": [0.35],
        "p_value": [1e-5],
        "rsid": ["rs123"],
        "rs_id": ["rs123"],
        "hm_coordinate_conversion": [""], "hm_code": [""], "variant_id": [""],
        "N": [434_152],
    })


class TestJoinOutcomePositionMatch:
    def test_direct_position_match(self):
        instr = _instrument_df()
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = _outcome_df()

        result, n_proxies = _join_outcome(instr, mock_outcome)

        assert not result.empty
        assert n_proxies == 0
        assert "EA_out" in result.columns
        assert result["EA_out"].iloc[0] == "A"

    def test_no_match_no_rsid_returns_empty(self):
        instr = _instrument_df(rsid=["."])  # no rsid, no position match
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = pd.DataFrame(columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
        ])

        result, n_proxies = _join_outcome(instr, mock_outcome)

        assert result.empty
        assert n_proxies == 0

    def test_duplicate_positions_in_outcome_logs_warning(self, caplog):
        instr = _instrument_df()
        # Two rows at same position in outcome
        out = pd.concat([_outcome_df(), _outcome_df()], ignore_index=True)
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = out

        with caplog.at_level("WARNING", logger="05_harmonise"):
            result, _ = _join_outcome(instr, mock_outcome)

        assert "duplicate" in caplog.text.lower()


class TestJoinOutcomeProxySearch:
    def test_proxy_used_when_direct_match_missing(self):
        instr = _instrument_df(rsid=["rs999"])
        # No direct position match
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = pd.DataFrame(columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
        ])
        # Proxy returns rs_proxy
        mock_outcome.fetch_by_rsid.return_value = _outcome_df()
        mock_outcome.fetch_by_rsid.return_value["rsid"] = "rs_proxy"

        with patch.object(_harmonise_mod, "find_proxies",
                          return_value={"rs999": ("rs_proxy", 0.95)}), \
             patch.object(_harmonise_mod, "in_phase_allele_map",
                          return_value={"A": "A", "G": "G"}):
            result, n_proxies = _join_outcome(instr, mock_outcome)

        assert n_proxies == 1
        assert not result.empty
        # Target SNP identity is preserved for harmonisation.
        assert result["rsid"].iloc[0] == "rs999"
        # Proxy provenance is retained separately.
        assert result["proxy_rsid"].iloc[0] == "rs_proxy"
        assert bool(result["proxy_used"].iloc[0]) is True

    def test_no_rsid_instruments_are_dropped_with_warning(self, caplog):
        instr = _instrument_df(rsid=["."])
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = pd.DataFrame(columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
        ])

        with caplog.at_level("WARNING", logger="05_harmonise"):
            result, n_proxies = _join_outcome(instr, mock_outcome)

        assert result.empty
        assert "no rsid" in caplog.text.lower() or "cannot" in caplog.text.lower()

    def test_proxy_high_maf_excluded(self):
        instr = _instrument_df(rsid=["rs999"])
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = pd.DataFrame(columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
        ])
        proxy_row = _outcome_df()
        proxy_row["effect_allele_frequency"] = 0.5
        mock_outcome.fetch_by_rsid.return_value = proxy_row

        with patch.object(_harmonise_mod, "find_proxies",
                          return_value={"rs999": ("rs_proxy", 0.95)}), \
             patch.object(_harmonise_mod, "in_phase_allele_map",
                          return_value={"A": "A", "G": "G"}):
            result, n_proxies = _join_outcome(instr, mock_outcome)

        assert n_proxies == 0
        assert result.empty

    def test_proxy_flip_applied_when_proxy_effect_allele_maps_to_target_oa(self):
        instr = _instrument_df(rsid=["rs999"], EA=["A"], OA=["G"])
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = pd.DataFrame(columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
        ])
        proxy_row = _outcome_df()
        proxy_row["rsid"] = "rs_proxy"
        proxy_row["effect_allele"] = "G"  # maps to target OA in mocked phase map
        proxy_row["beta"] = 0.2
        proxy_row["effect_allele_frequency"] = 0.35
        mock_outcome.fetch_by_rsid.return_value = proxy_row

        with patch.object(_harmonise_mod, "find_proxies",
                          return_value={"rs999": ("rs_proxy", 0.95)}), \
             patch.object(_harmonise_mod, "in_phase_allele_map",
                          return_value={"A": "A", "G": "G"}):
            result, n_proxies = _join_outcome(instr, mock_outcome)

        assert n_proxies == 1
        assert result["beta_out"].iloc[0] == pytest.approx(-0.2)
        assert result["EAF_out"].iloc[0] == pytest.approx(0.65)
        assert bool(result["proxy_flip"].iloc[0]) is True

    def test_proxy_dropped_when_phase_alignment_unavailable(self):
        instr = _instrument_df(rsid=["rs999"])
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = pd.DataFrame(columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
        ])
        proxy_row = _outcome_df()
        proxy_row["rsid"] = "rs_proxy"
        mock_outcome.fetch_by_rsid.return_value = proxy_row

        with patch.object(_harmonise_mod, "find_proxies",
                          return_value={"rs999": ("rs_proxy", 0.95)}), \
             patch.object(_harmonise_mod, "in_phase_allele_map", return_value=None):
            result, n_proxies = _join_outcome(instr, mock_outcome)

        assert n_proxies == 0
        assert result.empty


class TestCallHarmoniseR:
    def test_uses_scripts_rlib_path(self):
        df = _instrument_df()

        def fake_run(cmd, capture_output, text):
            assert cmd[0] == "Rscript"
            assert "scripts/rlib/harmonise.R" in cmd[1]
            result_path = cmd[cmd.index("--result") + 1]
            pd.DataFrame({"seqid": ["SeqId_TEST"], "rsid": ["rs123"]}).to_csv(
                result_path, sep="\t", index=False
            )

            class Result:
                returncode = 0
                stderr = ""

            return Result()

        with patch.object(_harmonise_mod.subprocess, "run", side_effect=fake_run):
            out = _harmonise_mod._call_harmonise_r(df, "SeqId_TEST")

        assert out is not None
        assert not out.empty

    def test_proxy_metadata_does_not_break_r_inputs(self):
        df = _instrument_df(
            rsid=["rs_target"],
            chrom_hg38=["22"],
            pos_hg38=[25_212_564],
            EA_out=["A"],
            OA_out=["G"],
            EAF_out=[0.31],
            beta_out=[0.12],
            se_out=[0.02],
            pval_out=[1e-6],
            N_out=[434_152],
            proxy_rsid=["rs_proxy"],
            proxy_r2=[0.91],
            proxy_used=[True],
            proxy_flip=[False],
            outcome_rsid=["rs_target"],
            outcome_chrom_hg38=["22"],
            outcome_pos_hg38=[25_212_570],
        )

        def fake_run(cmd, capture_output, text):
            exp_path = cmd[cmd.index("--exp") + 1]
            out_path = cmd[cmd.index("--out") + 1]
            result_path = cmd[cmd.index("--result") + 1]

            exp_in = pd.read_csv(exp_path, sep="\t")
            out_in = pd.read_csv(out_path, sep="\t")

            # Target SNP identity should be preserved for harmonisation.
            assert exp_in["rsid"].iloc[0] == "rs_target"
            assert out_in["rsid"].iloc[0] == "rs_target"
            # Proxy metadata should not be required/forwarded to R inputs.
            assert "proxy_rsid" not in exp_in.columns
            assert "proxy_rsid" not in out_in.columns

            pd.DataFrame({"seqid": ["SeqId_TEST"], "rsid": ["rs_target"]}).to_csv(
                result_path, sep="\t", index=False
            )

            class Result:
                returncode = 0
                stderr = ""

            return Result()

        with patch.object(_harmonise_mod.subprocess, "run", side_effect=fake_run):
            out = _harmonise_mod._call_harmonise_r(df, "SeqId_TEST")

        assert out is not None
        assert not out.empty

    def test_proxy_flipped_values_flow_into_r_outcome_input(self):
        instr = _instrument_df(rsid=["rs999"], EA=["A"], OA=["G"])
        mock_outcome = MagicMock()
        mock_outcome.fetch_snps.return_value = pd.DataFrame(columns=[
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
        ])
        proxy_row = _outcome_df()
        proxy_row["rsid"] = "rs_proxy"
        proxy_row["effect_allele"] = "G"  # maps to target OA in mocked phase map
        proxy_row["beta"] = 0.2
        proxy_row["effect_allele_frequency"] = 0.35
        mock_outcome.fetch_by_rsid.return_value = proxy_row

        with patch.object(_harmonise_mod, "find_proxies",
                          return_value={"rs999": ("rs_proxy", 0.95)}), \
             patch.object(_harmonise_mod, "in_phase_allele_map",
                          return_value={"A": "A", "G": "G"}):
            joined, _ = _join_outcome(instr, mock_outcome)

        assert not joined.empty

        def fake_run(cmd, capture_output, text):
            out_path = cmd[cmd.index("--out") + 1]
            result_path = cmd[cmd.index("--result") + 1]
            out_in = pd.read_csv(out_path, sep="\t")
            # Flip should have been applied before calling R.
            assert out_in["beta"].iloc[0] == pytest.approx(-0.2)
            assert out_in["EAF"].iloc[0] == pytest.approx(0.65)

            pd.DataFrame({"seqid": ["SeqId_TEST"], "rsid": ["rs999"]}).to_csv(
                result_path, sep="\t", index=False
            )

            class Result:
                returncode = 0
                stderr = ""

            return Result()

        with patch.object(_harmonise_mod.subprocess, "run", side_effect=fake_run):
            out = _harmonise_mod._call_harmonise_r(joined, "SeqId_TEST")

        assert out is not None
        assert not out.empty


def test_harmonise_cohort_proxy_path_blackbox_real_plink(tmp_path):
    seqid = "SeqId_TEST"
    cohort = "ARIC_EA"
    plink_mod = importlib.import_module("scripts.lib.plink")
    paths_mod = importlib.import_module("scripts.lib.paths")

    # Real PLINK2 + LD reference guardrails.
    plink_exec = str(plink_mod._PLINK2)
    try:
        probe = subprocess.run([plink_exec, "--help"], capture_output=True, text=True)
    except OSError:
        probe = None
    if probe is None or probe.returncode not in {0, 1}:
        fallback = Path("/Users/spuduch/Research/MR_IA/plink2_mac_arm64_20260228/plink2")
        if not fallback.exists():
            pytest.fail(f"PLINK2 is not runnable via PATH and fallback is missing: {fallback}")
        plink_mod._PLINK2 = str(fallback)
    ld_prefix = Path(paths_mod.LD_REF_PREFIX)
    if not Path(f"{ld_prefix}.bed").exists():
        pytest.fail(f"LD reference bed file missing: {ld_prefix}.bed")
    snplist = Path(f"{ld_prefix}.snplist")
    if not snplist.exists():
        pytest.fail(f"LD reference snplist missing: {snplist}")

    rsids: list[str] = []
    with snplist.open() as fh:
        for line in fh:
            rs = line.strip()
            if rs.startswith("rs"):
                rsids.append(rs)
            if len(rsids) >= 300:
                break
    if len(rsids) < 2:
        pytest.fail("Not enough rsIDs in LD reference snplist")

    proxy_map = plink_mod.find_proxies(rsids)
    if not proxy_map:
        pytest.fail("No proxy SNPs found in sampled LD-reference rsIDs")

    chosen_target = None
    chosen_proxy = None
    chosen_phase_map = None
    for target, (proxy, _r2) in list(proxy_map.items())[:40]:
        phase_map = plink_mod.in_phase_allele_map(target, proxy)
        if not phase_map:
            continue
        keys = [a for a in phase_map.keys() if len(a) == 1 and a in {"A", "C", "G", "T"}]
        vals = [phase_map[a] for a in keys if phase_map[a] in {"A", "C", "G", "T"}]
        if len(keys) == 2 and len(set(keys)) == 2 and len(vals) == 2:
            chosen_target = target
            chosen_proxy = proxy
            chosen_phase_map = phase_map
            break
    if not chosen_target or not chosen_proxy or not chosen_phase_map:
        pytest.fail("Could not find target/proxy pair with parseable single-base phase map")

    # Build an instrument row whose alleles are the two target alleles in phase map.
    target_alleles = list(chosen_phase_map.keys())[:2]
    target_ea = target_alleles[0]
    target_oa = target_alleles[1]

    in_dir = tmp_path / "instruments_hg38"
    out_dir = tmp_path / "harmonised"
    state_dir = tmp_path / "state"
    in_dir.mkdir(parents=True, exist_ok=True)
    state_dir.mkdir(parents=True, exist_ok=True)

    instr = pd.DataFrame({
        "seqid": [seqid],
        "gene": ["TESTGENE"],
        "uniprot": ["Q12345"],
        "chrom": ["1"],
        "pos": [100],
        "chrom_hg38": ["1"],
        "pos_hg38": [100],
        "rsid": [chosen_target],
        "EA": [target_ea],
        "OA": [target_oa],
        "EAF": [0.30],
        "beta": [0.50],
        "se": [0.05],
        "pval": [1e-9],
        "N": [7213],
        "build": ["hg19"],
    })
    instr.to_csv(in_dir / f"{seqid}.tsv", sep="\t", index=False)

    class FakeOutcome:
        def fetch_snps(self, positions):
            # Force proxy branch: no direct positional matches.
            return pd.DataFrame(columns=[
                "chromosome", "base_pair_location", "effect_allele", "other_allele",
                "beta", "standard_error", "effect_allele_frequency", "p_value",
                "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
            ])

        def fetch_by_rsid(self, rsids):
            assert rsids == [chosen_proxy]
            proxy_effect_allele = chosen_phase_map[target_oa]
            return pd.DataFrame({
                "chromosome": ["1"],
                "base_pair_location": [150],
                "effect_allele": [proxy_effect_allele],  # maps to target OA -> flip expected
                "other_allele": ["A" if proxy_effect_allele != "A" else "C"],
                "beta": [0.3],
                "standard_error": [0.02],
                "effect_allele_frequency": [0.2],
                "p_value": [1e-6],
                "rsid": [chosen_proxy],
                "rs_id": [chosen_proxy],
                "hm_coordinate_conversion": [""],
                "hm_code": [""],
                "variant_id": [""],
            })

        def __enter__(self):
            return self

        def __exit__(self, *_):
            return None

    real_subprocess_run = _harmonise_mod.subprocess.run

    def fake_run(cmd, *args, **kwargs):
        # Let PLINK subprocesses execute normally against real binary.
        if "--exp" not in cmd:
            return real_subprocess_run(cmd, *args, **kwargs)

        exp_path = Path(cmd[cmd.index("--exp") + 1])
        out_path = Path(cmd[cmd.index("--out") + 1])
        result_path = Path(cmd[cmd.index("--result") + 1])

        exp_in = pd.read_csv(exp_path, sep="\t")
        out_in = pd.read_csv(out_path, sep="\t")

        # Proxy was detected and transformed to target-allele orientation.
        assert exp_in["rsid"].iloc[0] == chosen_target
        assert out_in["rsid"].iloc[0] == chosen_target
        assert out_in["beta"].iloc[0] == pytest.approx(-0.3)
        assert out_in["EAF"].iloc[0] == pytest.approx(0.8)

        pd.DataFrame({"seqid": [seqid], "rsid": [chosen_target]}).to_csv(
            result_path, sep="\t", index=False
        )

        class Result:
            returncode = 0
            stderr = ""
            stdout = ""

        return Result()

    with patch.object(_harmonise_mod, "instruments_hg38_dir", return_value=in_dir), \
         patch.object(_harmonise_mod, "harmonised_dir", return_value=out_dir), \
         patch.object(_harmonise_mod, "cohort_dir", return_value=state_dir), \
         patch.object(_harmonise_mod, "OutcomeLookup", return_value=FakeOutcome()), \
         patch.object(_harmonise_mod.subprocess, "run", side_effect=fake_run):
        n_ok = _harmonise_mod.harmonise_cohort(cohort)

    assert n_ok == 1
    out_file = out_dir / f"{seqid}.tsv"
    assert out_file.exists()
