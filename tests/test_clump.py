"""Tests for scripts.03_clump.clump (clump_cohort logic), with mocked PLINK."""
import importlib
import pandas as pd
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

_clump_mod = importlib.import_module("scripts.03_clump.clump")
clump_cohort = _clump_mod.clump_cohort


def _write_cis_tsv(path: Path, seqid: str = "SeqId_TEST", n_variants: int = 3) -> None:
    df = pd.DataFrame({
        "seqid":   [seqid] * n_variants,
        "gene":    ["TESTGENE"] * n_variants,
        "uniprot": ["Q12345"] * n_variants,
        "chrom":   ["22"] * n_variants,
        "pos":     [25_212_564 + i * 10_000 for i in range(n_variants)],
        "rsid":    [f"rs{100 + i}" for i in range(n_variants)],
        "EA":      ["A"] * n_variants,
        "OA":      ["G"] * n_variants,
        "EAF":     [0.3] * n_variants,
        "beta":    [0.5] * n_variants,
        "se":      [0.05] * n_variants,
        "pval":    [1e-9] * n_variants,
        "N":       [7213] * n_variants,
        "build":   ["hg19"] * n_variants,
    })
    df.to_csv(path, sep="\t", index=False)


def _make_clump_output(rsids: list[str]) -> pd.DataFrame:
    """Fake PLINK clump output with one lead SNP per rsid."""
    return pd.DataFrame({"SNP": rsids, "P": [1e-9] * len(rsids)})


class TestClumpCohort:
    def test_single_protein_produces_instrument_tsv(self, tmp_path):
        seqid = "SeqId_TEST"
        in_dir = tmp_path / "filtered_cis_pqtls"
        out_dir = tmp_path / "instruments"
        state_dir = tmp_path / "state"
        in_dir.mkdir()
        state_dir.mkdir()

        _write_cis_tsv(in_dir / f"{seqid}.tsv", seqid=seqid)

        def fake_clump(df, seqid_arg, **kwargs):
            return df.iloc[:1].copy()

        with patch.object(_clump_mod, "filtered_cis_pqtls_dir", return_value=in_dir), \
             patch.object(_clump_mod, "instruments_dir", return_value=out_dir), \
             patch.object(_clump_mod, "cohort_dir", return_value=state_dir), \
             patch.object(_clump_mod, "clump", side_effect=fake_clump):
            n = clump_cohort("ARIC_EA")

        assert n == 1
        out_files = list(out_dir.glob("*.tsv"))
        assert len(out_files) == 1
        result = pd.read_csv(out_files[0], sep="\t")
        assert "F_stat" in result.columns

    def test_zero_cis_pqtls_protein_skipped_gracefully(self, tmp_path):
        seqid = "SeqId_EMPTY"
        in_dir = tmp_path / "filtered_cis_pqtls"
        out_dir = tmp_path / "instruments"
        state_dir = tmp_path / "state"
        in_dir.mkdir()
        state_dir.mkdir()

        df = pd.DataFrame({
            "seqid": [seqid], "gene": ["G"], "uniprot": ["Q1"], "chrom": ["22"],
            "pos": [25_212_564], "rsid": ["."], "EA": ["A"], "OA": ["G"],
            "EAF": [0.3], "beta": [0.5], "se": [0.05], "pval": [1e-9],
            "N": [7213], "build": ["hg19"],
        })
        (in_dir / f"{seqid}.tsv").write_text(df.to_csv(sep="\t", index=False))

        with patch.object(_clump_mod, "filtered_cis_pqtls_dir", return_value=in_dir), \
             patch.object(_clump_mod, "instruments_dir", return_value=out_dir), \
             patch.object(_clump_mod, "cohort_dir", return_value=state_dir):
            n = clump_cohort("ARIC_EA")

        assert n == 0
        assert not list(out_dir.glob("*.tsv"))

    def test_weak_instrument_flagged_not_dropped(self, tmp_path):
        seqid = "SeqId_WEAK"
        in_dir = tmp_path / "filtered_cis_pqtls"
        out_dir = tmp_path / "instruments"
        state_dir = tmp_path / "state"
        in_dir.mkdir()
        state_dir.mkdir()

        df = pd.DataFrame({
            "seqid": [seqid], "gene": ["G"], "uniprot": ["Q1"], "chrom": ["22"],
            "pos": [25_212_564], "rsid": ["rs100"], "EA": ["A"], "OA": ["G"],
            "EAF": [0.3],
            "beta": [0.1],  # F = (0.1/0.05)^2 = 4 < 10 → weak
            "se": [0.05],
            "pval": [1e-9], "N": [7213], "build": ["hg19"],
        })
        (in_dir / f"{seqid}.tsv").write_text(df.to_csv(sep="\t", index=False))

        def fake_clump(df_in, seqid_arg, **kwargs):
            return df_in.copy()

        with patch.object(_clump_mod, "filtered_cis_pqtls_dir", return_value=in_dir), \
             patch.object(_clump_mod, "instruments_dir", return_value=out_dir), \
             patch.object(_clump_mod, "cohort_dir", return_value=state_dir), \
             patch.object(_clump_mod, "clump", side_effect=fake_clump):
            n = clump_cohort("ARIC_EA")

        assert n == 1
        result = pd.read_csv(out_dir / f"{seqid}.tsv", sep="\t")
        assert result["F_stat"].iloc[0] == pytest.approx(4.0)
        assert len(result) == 1

    def test_clump_config_passed_to_plink_wrapper(self, tmp_path):
        seqid = "SeqId_CFG"
        in_dir = tmp_path / "filtered_cis_pqtls"
        out_dir = tmp_path / "instruments"
        state_dir = tmp_path / "state"
        in_dir.mkdir()
        state_dir.mkdir()
        _write_cis_tsv(in_dir / f"{seqid}.tsv", seqid=seqid)

        seen = {}

        def fake_clump(df, seqid_arg, **kwargs):
            seen.update(kwargs)
            return df.iloc[:1].copy()

        with patch.object(_clump_mod, "filtered_cis_pqtls_dir", return_value=in_dir), \
             patch.object(_clump_mod, "instruments_dir", return_value=out_dir), \
             patch.object(_clump_mod, "cohort_dir", return_value=state_dir), \
             patch.object(_clump_mod, "clump", side_effect=fake_clump):
            n = clump_cohort("ARIC_EA", window_kb=2500, r2=0.02, p1=1e-6)

        assert n == 1
        assert seen["window_kb"] == 2500
        assert seen["r2"] == 0.02
        assert seen["p1"] == 1e-6
