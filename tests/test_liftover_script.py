"""Tests for scripts.04_liftover.instruments_to_hg38 (lift_cohort logic)."""
import importlib
import pandas as pd
import numpy as np
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

_liftover_script = importlib.import_module("scripts.04_liftover.instruments_to_hg38")
lift_cohort = _liftover_script.lift_cohort


def _make_instrument_tsv(tmp_path: Path, seqid: str = "SeqId_TEST") -> Path:
    out = tmp_path / f"{seqid}.tsv"
    df = pd.DataFrame({
        "seqid": [seqid, seqid],
        "gene": ["GENE", "GENE"],
        "uniprot": ["Q12345", "Q12345"],
        "chrom": ["22", "22"],
        "pos": [25_212_564, 25_300_000],
        "rsid": ["rs1", "rs2"],
        "EA": ["A", "C"],
        "OA": ["G", "T"],
        "EAF": [0.3, 0.4],
        "beta": [0.1, 0.2],
        "se": [0.01, 0.02],
        "pval": [1e-9, 1e-10],
        "N": [7213, 7213],
        "build": ["hg19", "hg19"],
    })
    df.to_csv(out, sep="\t", index=False)
    return out


class TestLiftCohortHg38Passthrough:
    def test_hg38_cohort_columns_passthrough(self, tmp_path):
        """deCODE (hg38) should add chrom_hg38/pos_hg38 == chrom/pos without liftover."""
        in_dir = tmp_path / "instruments"
        out_dir = tmp_path / "instruments_hg38"
        state_dir = tmp_path / "state"
        in_dir.mkdir()
        state_dir.mkdir()

        tsv = _make_instrument_tsv(in_dir)
        df_in = pd.read_csv(tsv, sep="\t")
        df_in["build"] = "hg38"
        df_in.to_csv(tsv, sep="\t", index=False)

        with patch.object(_liftover_script, "instruments_dir", return_value=in_dir), \
             patch.object(_liftover_script, "instruments_hg38_dir", return_value=out_dir), \
             patch.object(_liftover_script, "cohort_dir", return_value=state_dir):
            n = lift_cohort("deCODE")

        assert n == 1
        out_files = list(out_dir.glob("*.tsv"))
        assert len(out_files) == 1
        df_out = pd.read_csv(out_files[0], sep="\t")
        assert "chrom_hg38" in df_out.columns
        assert "pos_hg38" in df_out.columns
        assert (df_out["chrom_hg38"] == df_out["chrom"]).all()
        assert (df_out["pos_hg38"] == df_out["pos"]).all()


class TestLiftCohortZeroDivisionGuard:
    def test_n_in_zero_does_not_crash(self, tmp_path):
        """Empty instrument file should not cause ZeroDivisionError in drop-pct logging."""
        in_dir = tmp_path / "instruments"
        out_dir = tmp_path / "instruments_hg38"
        state_dir = tmp_path / "state"
        in_dir.mkdir()
        state_dir.mkdir()

        # Write an empty (header-only) TSV
        from scripts.lib.schema import NORM_COLS
        tsv = in_dir / "SeqId_EMPTY.tsv"
        pd.DataFrame(columns=NORM_COLS).to_csv(tsv, sep="\t", index=False)

        with patch.object(_liftover_script, "instruments_dir", return_value=in_dir), \
             patch.object(_liftover_script, "instruments_hg38_dir", return_value=out_dir), \
             patch.object(_liftover_script, "cohort_dir", return_value=state_dir):
            n = lift_cohort("ARIC_EA")

        assert n == 0


class TestLiftCohortHg19Liftover:
    def test_hg19_cohort_calls_lift_table(self, tmp_path):
        """Fenland (genuinely hg19) should call lift_table and produce hg38 columns."""
        in_dir = tmp_path / "instruments"
        out_dir = tmp_path / "instruments_hg38"
        state_dir = tmp_path / "state"
        in_dir.mkdir()
        state_dir.mkdir()

        _make_instrument_tsv(in_dir)

        # Mock lift_table to return input with added hg38 columns
        def fake_lift_table(df, chrom_col, pos_col, **kwargs):
            df = df.copy()
            df["chrom_hg38"] = df[chrom_col]
            df["pos_hg38"] = df[pos_col]
            return df

        with patch.object(_liftover_script, "instruments_dir", return_value=in_dir), \
             patch.object(_liftover_script, "instruments_hg38_dir", return_value=out_dir), \
             patch.object(_liftover_script, "cohort_dir", return_value=state_dir), \
             patch.object(_liftover_script, "lift_table", side_effect=fake_lift_table):
            n = lift_cohort("Fenland")

        assert n == 1
        out_files = list(out_dir.glob("*.tsv"))
        assert len(out_files) == 1
        df_out = pd.read_csv(out_files[0], sep="\t")
        assert "chrom_hg38" in df_out.columns
        assert "pos_hg38" in df_out.columns


def test_hg38_cohorts_includes_aric_and_ukbppp():
    """ARIC_EA and UKB_PPP must be in HG38_COHORTS — their positions are already hg38."""
    assert "ARIC_EA" in _liftover_script.HG38_COHORTS, \
        "ARIC_EA positions are hg38; must not be lifted"
    assert "UKB_PPP" in _liftover_script.HG38_COHORTS, \
        "UKB_PPP positions are hg38; must not be lifted"


def test_aric_ukbppp_instruments_pass_through_without_shift(tmp_path):
    """lift_cohort for ARIC_EA/UKB_PPP must not shift positions (pos == pos_hg38)."""
    for cohort in ("ARIC_EA", "UKB_PPP"):
        in_dir = tmp_path / cohort / "instruments"
        out_dir = tmp_path / cohort / "instruments_hg38"
        state_dir = tmp_path / cohort
        in_dir.mkdir(parents=True)

        df = pd.DataFrame({
            "seqid": ["S"], "chrom": ["2"], "pos": [119_390_517],
            "rsid": ["rs1"], "EA": ["A"], "OA": ["G"],
            "EAF": [0.3], "beta": [0.1], "se": [0.01],
            "pval": [1e-9], "N": [7000], "build": ["hg38"],
        })
        df.to_csv(in_dir / "S.tsv", sep="\t", index=False)

        with patch.object(_liftover_script, "instruments_dir", return_value=in_dir), \
             patch.object(_liftover_script, "instruments_hg38_dir", return_value=out_dir), \
             patch.object(_liftover_script, "cohort_dir", return_value=state_dir):
            lift_cohort(cohort)

        result = pd.read_csv(next(out_dir.glob("*.tsv")), sep="\t")
        assert (result["pos"] == result["pos_hg38"]).all(), \
            f"{cohort}: pos and pos_hg38 must be equal (no liftover applied)"


def test_fenland_instruments_do_shift_positions(tmp_path):
    """lift_cohort for Fenland (genuinely hg19) must shift positions."""
    in_dir  = tmp_path / "instruments"
    out_dir = tmp_path / "instruments_hg38"
    state_dir = tmp_path
    in_dir.mkdir()

    df = pd.DataFrame({
        "seqid": ["S"], "chrom": ["4"], "pos": [100_092_382],
        "rsid": ["rs1"], "EA": ["A"], "OA": ["G"],
        "EAF": [0.3], "beta": [0.1], "se": [0.01],
        "pval": [1e-9], "N": [10000], "build": ["hg19"],
    })
    df.to_csv(in_dir / "S.tsv", sep="\t", index=False)

    def fake_lift_table(df, chrom_col, pos_col, **kwargs):
        """Simulate a position shift (hg19→hg38 adds ~200 kb for chr4)."""
        df = df.copy()
        df["chrom_hg38"] = df[chrom_col]
        df["pos_hg38"] = df[pos_col] + 200_000
        return df

    with patch.object(_liftover_script, "instruments_dir", return_value=in_dir), \
         patch.object(_liftover_script, "instruments_hg38_dir", return_value=out_dir), \
         patch.object(_liftover_script, "cohort_dir", return_value=state_dir), \
         patch.object(_liftover_script, "lift_table", side_effect=fake_lift_table):
        lift_cohort("Fenland")

    result = pd.read_csv(next(out_dir.glob("*.tsv")), sep="\t")
    assert (result["pos"] != result["pos_hg38"]).any(), \
        "Fenland: hg19 positions must be shifted by liftover"
