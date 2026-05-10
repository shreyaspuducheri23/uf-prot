"""
End-to-end smoke test: 2 synthetic proteins through extraction → clump → liftover → harmonise → assemble.

All external dependencies (PLINK, R, Synapse, tabix) are mocked.
This test validates that the pipeline wiring is correct and produces
a final results table with expected structure.
"""
import importlib
import json
import pandas as pd
import numpy as np
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from scripts.lib.schema import ProteinMeta, NORM_COLS
from scripts.lib.cis_extract import run_extraction
from scripts.lib.fdr import add_fdr

# Import numbered-directory modules via importlib (identifiers can't start with digits)
_clump_mod = importlib.import_module("scripts.03_clump.clump")
_liftover_mod = importlib.import_module("scripts.04_liftover.instruments_to_hg38")
_harmonise_mod = importlib.import_module("scripts.05_harmonise.harmonise")
_assemble_mod = importlib.import_module("scripts.09_assemble.assemble")


# ─── Fixture proteins ────────────────────────────────────────────────────────

PROTEINS = [
    ProteinMeta(
        seqid="SeqId_P1", gene="GENE1", uniprot="Q11111",
        chrom="22", tss=25_212_564, build="hg19", source_cohort="ARIC_EA",
    ),
    ProteinMeta(
        seqid="SeqId_P2", gene="GENE2", uniprot="Q22222",
        chrom="1", tss=10_000_000, build="hg19", source_cohort="ARIC_EA",
    ),
]


def _good_cis_df(protein: ProteinMeta) -> pd.DataFrame:
    """Minimal valid cis sumstats for one protein."""
    tss = protein.tss
    return pd.DataFrame({
        "chrom": [protein.chrom],
        "pos":   [tss],
        "rsid":  ["rs12345"],
        "EA":    ["A"],
        "OA":    ["G"],
        "EAF":   [0.30],
        "beta":  [0.50],
        "se":    [0.05],
        "pval":  [1e-9],
        "N":     [7213],
    })


# ─── Step 02: extraction ─────────────────────────────────────────────────────

@pytest.fixture
def extracted_dir(tmp_path):
    """Run extraction for both proteins, return the output directory."""
    out_dir = tmp_path / "cis_sumstats"
    out_dir.mkdir(parents=True)

    with patch("scripts.lib.cis_extract.cis_sumstats_dir", return_value=out_dir), \
         patch("scripts.lib.cis_extract.cohort_dir", return_value=tmp_path / "state"):
        run_extraction(
            "ARIC_EA", PROTEINS,
            read_fn=_good_cis_df,
        )
    return out_dir


# ─── Step 03: clump (mocked PLINK) ───────────────────────────────────────────

@pytest.fixture
def instruments_dir(tmp_path, extracted_dir):
    """Clump extracted proteins; return instruments directory."""
    clump_cohort = _clump_mod.clump_cohort

    out_dir = tmp_path / "instruments"
    state_dir = tmp_path / "clump_state"
    state_dir.mkdir(parents=True)

    def fake_clump(df, seqid, **kwargs):
        return df.copy()  # trivially: all variants are lead SNPs

    with patch.object(_clump_mod, "cis_sumstats_dir", return_value=extracted_dir), \
         patch.object(_clump_mod, "instruments_dir", return_value=out_dir), \
         patch.object(_clump_mod, "cohort_dir", return_value=state_dir), \
         patch.object(_clump_mod, "clump", side_effect=fake_clump):
        clump_cohort("ARIC_EA")

    return out_dir


# ─── Step 04: liftover (mocked) ──────────────────────────────────────────────

@pytest.fixture
def instruments_hg38_dir(tmp_path, instruments_dir):
    """Lift instrument files; return hg38 directory."""
    lift_cohort = _liftover_mod.lift_cohort

    out_dir = tmp_path / "instruments_hg38"
    state_dir = tmp_path / "liftover_state"
    state_dir.mkdir(parents=True)

    def fake_lift_table(df, chrom_col="chrom", pos_col="pos", **kwargs):
        df = df.copy()
        df["chrom_hg38"] = df[chrom_col]
        df["pos_hg38"] = df[pos_col]
        return df

    with patch.object(_liftover_mod, "instruments_dir", return_value=instruments_dir), \
         patch.object(_liftover_mod, "instruments_hg38_dir", return_value=out_dir), \
         patch.object(_liftover_mod, "cohort_dir", return_value=state_dir), \
         patch.object(_liftover_mod, "lift_table", side_effect=fake_lift_table):
        lift_cohort("ARIC_EA")

    return out_dir


# ─── Step 05: harmonise (mocked outcome + mocked R) ─────────────────────────

def _fake_outcome_row(chrom, pos):
    return {
        "chrom_hg38": str(chrom),
        "pos_hg38": int(pos),
        "rsid": "rs12345",
        "EA_out": "A", "OA_out": "G",
        "EAF_out": 0.32,
        "beta_out": 0.18, "se_out": 0.02,
        "pval_out": 1e-6,
        "N_out": 434_152,
    }


@pytest.fixture
def harmonised_dir(tmp_path, instruments_hg38_dir):
    harmonise_cohort = _harmonise_mod.harmonise_cohort

    out_dir = tmp_path / "harmonised"
    state_dir = tmp_path / "harm_state"
    state_dir.mkdir(parents=True)

    class FakeOutcome:
        def fetch_snps(self, positions):
            rows = []
            for chrom, pos in positions:
                rows.append({
                    "chromosome": str(chrom), "base_pair_location": pos,
                    "effect_allele": "A", "other_allele": "G",
                    "beta": "0.18", "standard_error": "0.02",
                    "effect_allele_frequency": "0.32",
                    "p_value": "1e-6",
                    "rsid": "rs12345", "rs_id": "rs12345",
                    "hm_coordinate_conversion": "", "hm_code": "", "variant_id": "",
                })
            cols = [
                "chromosome", "base_pair_location", "effect_allele", "other_allele",
                "beta", "standard_error", "effect_allele_frequency", "p_value",
                "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
            ]
            df = pd.DataFrame(rows, columns=cols)
            df["base_pair_location"] = df["base_pair_location"].astype(int)
            df["N"] = 434_152
            return df

        def fetch_by_rsid(self, rsids):
            return pd.DataFrame()

        def __enter__(self):
            return self

        def __exit__(self, *_):
            pass

    def fake_harmonise_r(df, seqid):
        # Return a minimal harmonised DataFrame
        return pd.DataFrame({
            "seqid": [seqid], "rsid": ["rs12345"],
            "EA": ["A"], "OA": ["G"],
            "beta_exp": [0.5], "se_exp": [0.05], "pval_exp": [1e-9], "N_exp": [7213],
            "EA_out": ["A"], "OA_out": ["G"],
            "beta_out": [0.18], "se_out": [0.02], "pval_out": [1e-6], "N_out": [434_152],
        })

    with patch.object(_harmonise_mod, "instruments_hg38_dir", return_value=instruments_hg38_dir), \
         patch.object(_harmonise_mod, "harmonised_dir", return_value=out_dir), \
         patch.object(_harmonise_mod, "cohort_dir", return_value=state_dir), \
         patch.object(_harmonise_mod, "OutcomeLookup", return_value=FakeOutcome()), \
         patch.object(_harmonise_mod, "find_proxies", return_value={}), \
         patch.object(_harmonise_mod, "_call_harmonise_r", side_effect=fake_harmonise_r):
        harmonise_cohort("ARIC_EA")

    return out_dir


# ─── Step 09: assemble ───────────────────────────────────────────────────────

class TestEndToEndSmoke:
    def test_pipeline_produces_results_for_both_proteins(
        self, tmp_path, harmonised_dir
    ):
        """Verify that both proteins have output files after harmonisation."""
        files = list(harmonised_dir.glob("*.tsv"))
        seqids = {f.stem for f in files}
        assert "SeqId_P1" in seqids
        assert "SeqId_P2" in seqids

    def test_fdr_tiering_produces_expected_columns(self, harmonised_dir):
        """Read harmonised outputs, add FDR, apply tier — expect tier column."""
        tier = _assemble_mod.tier

        frames = []
        for tsv in harmonised_dir.glob("*.tsv"):
            df = pd.read_csv(tsv, sep="\t")
            df["seqid"] = tsv.stem
            frames.append(df)

        mr = pd.concat(frames, ignore_index=True)
        mr["pval"] = mr.get("pval_exp", pd.Series([1e-9] * len(mr)))
        mr = add_fdr(mr, pval_col="pval", alpha=0.05)
        mr["passes_sensitivity"] = True
        mr["sharepro_coloc_positive"] = True
        mr["coloc_abf_positive"] = True
        mr["tier"] = mr.apply(tier, axis=1)

        assert "tier" in mr.columns
        assert len(mr) == 2  # one row per protein
        assert all(mr["tier"] == "Tier1_replicated")

    def test_extraction_output_has_required_norm_cols(self, extracted_dir):
        """Extracted TSVs must have all NORM_COLS."""
        for tsv in extracted_dir.glob("*.tsv"):
            df = pd.read_csv(tsv, sep="\t")
            missing = [c for c in NORM_COLS if c not in df.columns]
            assert not missing, f"{tsv.name}: missing columns {missing}"
