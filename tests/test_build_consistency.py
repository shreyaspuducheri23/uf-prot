"""Tests that build metadata is internally consistent across the pipeline."""
import importlib
import pytest
import pandas as pd
from pathlib import Path


# ---------------------------------------------------------------------------
# Unit: BUILD constants in step 2 scripts must match native coordinates
# ---------------------------------------------------------------------------

def test_aric_build_constant_is_hg38():
    mod = importlib.import_module("scripts.02_cis_pqtl_extract.aric")
    assert mod.BUILD == "hg38", (
        "ARIC_EA .glm.linear files are in hg38 — BUILD must reflect that. "
        "Wrong value causes instruments to be incorrectly lifted in step 4."
    )


def test_ukbppp_build_constant_is_hg19():
    mod = importlib.import_module("scripts.02_cis_pqtl_extract.ukbppp")
    assert mod.BUILD == "hg19", (
        "UKB-PPP positions are native hg19/GRCh37 and must be lifted in step 4."
    )


def test_extractor_build_constants_match_config(pipeline_cfg):
    extractors = {
        "ARIC_EA": "scripts.02_cis_pqtl_extract.aric",
        "deCODE": "scripts.02_cis_pqtl_extract.decode",
        "UKB_PPP": "scripts.02_cis_pqtl_extract.ukbppp",
        "Fenland": "scripts.02_cis_pqtl_extract.fenland",
        "UKB_female": "scripts.02_cis_pqtl_extract.ukb_female",
    }
    for cohort, module_name in extractors.items():
        mod = importlib.import_module(module_name)
        assert mod.BUILD == pipeline_cfg["cohorts"][cohort]["build"]


def test_hg38_liftover_passthrough_sets_match_extractor_builds(pipeline_cfg):
    liftover = importlib.import_module("scripts.04_liftover.instruments_to_hg38")
    expected_hg38 = {
        cohort for cohort, cfg in pipeline_cfg["cohorts"].items()
        if cfg["build"] == "hg38"
    }
    assert liftover.HG38_COHORTS == expected_hg38
    assert liftover.CIS_HG38_COHORTS == expected_hg38


# ---------------------------------------------------------------------------
# Integration: lead SNP must be within ±window_kb of TSS in the same build
# ---------------------------------------------------------------------------

PROCESSED = Path("processed_data")
WINDOW_KB = 500   # step 2 extraction window
TOLERANCE = 50_000  # allow for TSS imprecision / transcript model differences


@pytest.mark.skipif(
    not (PROCESSED / "ARIC_EA" / "protein_index.tsv").exists(),
    reason="processed_data not present",
)
@pytest.mark.parametrize("cohort", ["ARIC_EA", "UKB_PPP", "Fenland", "deCODE"])
def test_lead_snp_within_tss_window(cohort):
    """For each cohort, sample proteins and verify the lead instrument is within
    ±window_kb of the TSS. Fails if build metadata causes coordinate mismatch."""
    idx_path = PROCESSED / cohort / "protein_index.tsv"
    inst_dir = PROCESSED / cohort / "instruments"
    if not idx_path.exists() or not inst_dir.exists():
        pytest.skip(f"{cohort}: missing protein_index or instruments dir")

    idx = pd.read_csv(idx_path, sep="\t", dtype=str)
    sample = idx.sample(min(10, len(idx)), random_state=42)
    failures = []

    for _, row in sample.iterrows():
        seqid = row["seqid"]
        tss   = int(row["tss"])
        inst_path = inst_dir / f"{seqid}.tsv"
        if not inst_path.exists():
            continue
        inst = pd.read_csv(inst_path, sep="\t", dtype={"chrom": str})
        if inst.empty:
            continue
        lead_pos = int(inst.iloc[0]["pos"])
        dist = abs(lead_pos - tss)
        if dist > WINDOW_KB * 1000 + TOLERANCE:
            failures.append(
                f"{seqid}: lead pos={lead_pos:,}, tss={tss:,}, dist={dist:,} "
                f"(>{WINDOW_KB * 1000 + TOLERANCE:,}) — possible build mismatch"
            )

    assert not failures, (
        f"{cohort}: {len(failures)} proteins have lead SNP far from TSS:\n"
        + "\n".join(failures)
    )


@pytest.mark.skipif(
    not (PROCESSED / "ARIC_EA" / "instruments_hg38").exists(),
    reason="processed_data not present",
)
@pytest.mark.parametrize("cohort", ["ARIC_EA", "deCODE"])
def test_instruments_hg38_pos_equals_pos_hg38_for_hg38_cohorts(cohort):
    """For hg38 cohorts, liftover must be a no-op: pos == pos_hg38."""
    _mod = importlib.import_module("scripts.04_liftover.instruments_to_hg38")
    HG38_COHORTS = _mod.HG38_COHORTS
    assert cohort in HG38_COHORTS, f"{cohort} should be in HG38_COHORTS"

    hg38_dir = PROCESSED / cohort / "instruments_hg38"
    if not hg38_dir.exists():
        pytest.skip(f"{cohort}: instruments_hg38 dir missing")

    files = list(hg38_dir.glob("*.tsv"))[:10]  # sample 10
    for f in files:
        df = pd.read_csv(f, sep="\t", dtype={"chrom": str, "chrom_hg38": str})
        if df.empty:
            continue
        bad = df[df["pos"].astype(int) != df["pos_hg38"].astype(int)]
        assert bad.empty, (
            f"{cohort} {f.stem}: {len(bad)} rows have pos != pos_hg38 "
            f"(liftover incorrectly applied to already-hg38 positions)"
        )
