"""Deterministic simulation benchmark suite for target-recovery stress testing."""
import importlib
import math
from dataclasses import dataclass
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from scripts.lib.cis_extract import run_extraction
from scripts.lib.fdr import add_fdr
from scripts.lib.schema import ProteinMeta

_clump_mod = importlib.import_module("scripts.03_clump.clump")
_liftover_mod = importlib.import_module("scripts.04_liftover.instruments_to_hg38")
_harm_mod = importlib.import_module("scripts.05_harmonise.harmonise")
_assemble_mod = importlib.import_module("scripts.09_assemble.assemble")


def _norm_sf(z: float) -> float:
    return 0.5 * math.erfc(z / math.sqrt(2))


@dataclass(frozen=True)
class Scenario:
    name: str
    mediation_scale: float
    pleiotropy_beta: float
    coloc_for_causal: bool


def _build_proteins() -> tuple[list[ProteinMeta], dict[str, str]]:
    proteins: list[ProteinMeta] = []
    labels: dict[str, str] = {}
    for i in range(12):
        seqid = f"SeqId_{i:02d}"
        proteins.append(
            ProteinMeta(
                seqid=seqid,
                gene=f"G{i:02d}",
                uniprot=f"P{i:05d}",
                chrom="1",
                tss=1_000_000 + i * 50_000,
                build="hg19",
                source_cohort="ARIC_EA",
            )
        )
        if i < 4:
            labels[seqid] = "causal"
        elif i < 6:
            labels[seqid] = "pleiotropic"
        elif i < 8:
            labels[seqid] = "weak"
        else:
            labels[seqid] = "null"
    return proteins, labels


def _run_scenario(tmp_path: Path, scenario: Scenario, seed: int = 42) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    proteins, labels = _build_proteins()
    cohort = "ARIC_EA"
    root = tmp_path / scenario.name
    cohort_dir = root / cohort
    effect_by_pos: dict[int, float] = {}

    def read_fn(protein: ProteinMeta) -> pd.DataFrame:
        rows = []
        label = labels[protein.seqid]
        for j in range(8):
            pos = protein.tss - 350_000 + j * 100_000
            if label == "weak":
                beta_exp = float(rng.normal(0.06, 0.015))
                se_exp = float(rng.uniform(0.04, 0.06))
            else:
                beta_exp = float(rng.normal(0.25, 0.05))
                se_exp = float(rng.uniform(0.03, 0.06))
            pval_exp = 1e-10 if j < 4 else 1e-7
            rows.append(
                {
                    "chrom": protein.chrom,
                    "pos": pos,
                    "rsid": f"rs{protein.seqid}_{j}",
                    "EA": "A",
                    "OA": "G",
                    "EAF": 0.2 + 0.05 * (j % 3),
                    "beta": beta_exp,
                    "se": se_exp,
                    "pval": pval_exp,
                    "N": 7000,
                }
            )
            if label in {"causal", "weak"}:
                effect = scenario.mediation_scale * beta_exp + float(rng.normal(0.0, 0.02))
            elif label == "pleiotropic":
                effect = scenario.pleiotropy_beta + float(rng.normal(0.0, 0.02))
            else:
                effect = float(rng.normal(0.0, 0.02))
            effect_by_pos[pos] = effect
        return pd.DataFrame(rows)

    cis_dir = cohort_dir / "cis_sumstats"
    with patch("scripts.lib.cis_extract.cis_sumstats_dir", return_value=cis_dir), \
         patch("scripts.lib.cis_extract.cohort_dir", return_value=cohort_dir):
        run_extraction(cohort, proteins, read_fn)

    def fake_clump(df: pd.DataFrame, seqid: str, **kwargs) -> pd.DataFrame:
        return df.sort_values("pval").head(3).copy()

    inst_dir = cohort_dir / "instruments"
    with patch.object(_clump_mod, "cis_sumstats_dir", return_value=cis_dir), \
         patch.object(_clump_mod, "instruments_dir", return_value=inst_dir), \
         patch.object(_clump_mod, "cohort_dir", return_value=cohort_dir), \
         patch.object(_clump_mod, "clump", side_effect=fake_clump):
        _clump_mod.clump_cohort(cohort, window_kb=1000, r2=0.001, p1=5e-8)

    hg38_dir = cohort_dir / "instruments_hg38"

    def fake_lift_table(df: pd.DataFrame, chrom_col: str = "chrom", pos_col: str = "pos", **kwargs):
        lifted = df.copy()
        lifted["chrom_hg38"] = lifted[chrom_col]
        lifted["pos_hg38"] = lifted[pos_col]
        return lifted

    with patch.object(_liftover_mod, "instruments_dir", return_value=inst_dir), \
         patch.object(_liftover_mod, "instruments_hg38_dir", return_value=hg38_dir), \
         patch.object(_liftover_mod, "cohort_dir", return_value=cohort_dir), \
         patch.object(_liftover_mod, "lift_table", side_effect=fake_lift_table):
        _liftover_mod.lift_cohort(cohort)

    harm_dir = cohort_dir / "harmonised"

    class FakeOutcome:
        def fetch_snps(self, positions):
            rows = []
            for chrom, pos in positions:
                rows.append(
                    {
                        "chromosome": str(chrom),
                        "base_pair_location": int(pos),
                        "effect_allele": "A",
                        "other_allele": "G",
                        "beta": effect_by_pos[int(pos)],
                        "standard_error": 0.03,
                        "effect_allele_frequency": 0.28,
                        "p_value": 1e-6,
                        "rsid": f"rs{int(pos)}",
                        "rs_id": f"rs{int(pos)}",
                        "hm_coordinate_conversion": "",
                        "hm_code": "",
                        "variant_id": "",
                        "N": 434_152,
                    }
                )
            return pd.DataFrame(rows)

        def fetch_by_rsid(self, rsids):
            return pd.DataFrame()

        def __enter__(self):
            return self

        def __exit__(self, *_):
            return None

    def fake_harmonise_r(df: pd.DataFrame, seqid: str) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "seqid": seqid,
                "rsid": df["rsid"],
                "EA": df["EA"],
                "OA": df["OA"],
                "beta_exp": df["beta"],
                "se_exp": df["se"],
                "pval_exp": df["pval"],
                "N_exp": df["N"],
                "EA_out": df["EA_out"],
                "OA_out": df["OA_out"],
                "beta_out": df["beta_out"],
                "se_out": df["se_out"],
                "pval_out": df["pval_out"],
                "N_out": df["N_out"],
            }
        )

    with patch.object(_harm_mod, "instruments_hg38_dir", return_value=hg38_dir), \
         patch.object(_harm_mod, "harmonised_dir", return_value=harm_dir), \
         patch.object(_harm_mod, "cohort_dir", return_value=cohort_dir), \
         patch.object(_harm_mod, "OutcomeLookup", return_value=FakeOutcome()), \
         patch.object(_harm_mod, "find_proxies", return_value={}), \
         patch.object(_harm_mod, "_call_harmonise_r", side_effect=fake_harmonise_r):
        _harm_mod.harmonise_cohort(cohort)

    rows = []
    for path in sorted(harm_dir.glob("*.tsv")):
        df = pd.read_csv(path, sep="\t")
        mean_beta = float(df["beta_out"].mean())
        z = abs(mean_beta) / 0.03 * math.sqrt(len(df))
        pval = 2 * _norm_sf(z)
        mean_f = float(((df["beta_exp"] / df["se_exp"]) ** 2).mean())
        label = labels[path.stem]
        passes_sens = label in {"causal", "null"}
        rows.append(
            {
                "seqid": path.stem,
                "class": label,
                "cohort": cohort,
                "pval": pval,
                "mean_F": mean_f,
                "passes_sensitivity": passes_sens,
                "sharepro_coloc_positive": scenario.coloc_for_causal and label == "causal",
                "coloc_abf_positive": scenario.coloc_for_causal and label == "causal",
            }
        )

    mr = add_fdr(pd.DataFrame(rows), pval_col="pval", alpha=0.05)
    mr["tier"] = mr.apply(_assemble_mod.tier, axis=1)
    return mr


def test_shared_signal_recovers_reasonable_targets(tmp_path):
    mr = _run_scenario(
        tmp_path,
        Scenario(
            name="shared_signal",
            mediation_scale=0.9,
            pleiotropy_beta=0.20,
            coloc_for_causal=True,
        ),
    )
    causal = mr[mr["class"] == "causal"]
    nulls = mr[mr["class"] == "null"]
    assert (causal["tier"] == "Tier1_replicated").mean() >= 0.75
    assert (nulls["tier"] == "Tier3").mean() >= 0.75


def test_coloc_mismatch_downgrades_causal_targets(tmp_path):
    mr = _run_scenario(
        tmp_path,
        Scenario(
            name="coloc_mismatch",
            mediation_scale=0.9,
            pleiotropy_beta=0.20,
            coloc_for_causal=False,
        ),
    )
    causal = mr[mr["class"] == "causal"]
    assert (causal["tier"].isin(["Tier1", "Tier1_replicated"])).sum() == 0
    assert (causal["tier"] == "Tier2").mean() >= 0.75


def test_pleiotropy_not_promoted_to_primary_targets(tmp_path):
    mr = _run_scenario(
        tmp_path,
        Scenario(
            name="pleiotropy_guardrail",
            mediation_scale=0.9,
            pleiotropy_beta=0.25,
            coloc_for_causal=True,
        ),
    )
    pleio = mr[mr["class"] == "pleiotropic"]
    assert (pleio["fdr_pass"]).mean() >= 0.5
    assert (pleio["tier"].isin(["Tier1", "Tier1_replicated", "Tier2"])).sum() == 0


def test_weak_instrument_targets_not_prioritized(tmp_path):
    mr = _run_scenario(
        tmp_path,
        Scenario(
            name="weak_instrument_guardrail",
            mediation_scale=0.7,
            pleiotropy_beta=0.20,
            coloc_for_causal=True,
        ),
    )
    weak = mr[mr["class"] == "weak"]
    assert (weak["mean_F"] < 10).mean() >= 0.5
    assert (weak["tier"].isin(["Tier1", "Tier1_replicated", "Tier2"])).sum() == 0
