"""GWASBrewer oracle scenario matrix for pipeline behavior validation."""

from __future__ import annotations

import pandas as pd
import pytest

from tests.gwas_oracle_harness import (
    compute_oracle_mr_and_tiers,
    generate_oracle_dataset,
    run_full_oracle_02_to_09,
    run_oracle_pipeline,
)


def _assert_non_empty(summary: pd.DataFrame) -> None:
    assert not summary.empty, "Oracle scenario produced empty summary"


def test_oracle_proxy_branch_executes_and_is_detected(tmp_path):
    bundle = generate_oracle_dataset(tmp_path, scenario="proxy_branch", seed=21, n_proteins=12, snps_per_protein=8)
    artifacts = run_oracle_pipeline(tmp_path, bundle, force_proxy_branch=True, use_real_harmonise_r=False)
    summary = compute_oracle_mr_and_tiers(artifacts)
    _assert_non_empty(summary)

    expected_proxy = summary[summary["expected_proxy"] == True]
    assert not expected_proxy.empty, "Manifest expected proxy proteins but none were present"
    proxy_seen = expected_proxy["proxy_used_any"].mean()
    assert proxy_seen >= 0.6, f"Proxy path under-executed: {proxy_seen:.3f}"

    # Downstream branch effect: harmonised outputs must explicitly carry proxy provenance.
    harm_dir = artifacts["harmonised_dir"]
    proxy_flag_seen = False
    for path in harm_dir.glob("*.tsv"):
        df = pd.read_csv(path, sep="\t")
        if "proxy_used" in df.columns and df["proxy_used"].astype(bool).any():
            proxy_flag_seen = True
            break
    assert proxy_flag_seen, "force_proxy_branch run produced no harmonised proxy_used rows"


@pytest.mark.gwas_oracle_slow
def test_oracle_coloc_mismatch_downgrades_causal_targets(tmp_path):
    bundle = generate_oracle_dataset(tmp_path, scenario="coloc_mismatch", seed=31, n_proteins=12, snps_per_protein=8)
    artifacts = run_oracle_pipeline(tmp_path, bundle, force_proxy_branch=False, use_real_harmonise_r=False)
    summary = compute_oracle_mr_and_tiers(artifacts)
    _assert_non_empty(summary)

    causal = summary[summary["class"] == "causal"]
    assert (causal["tier"].isin(["Tier1", "Tier1_replicated"])).sum() == 0
    assert (causal["tier"] == "Tier2").mean() >= 0.6


@pytest.mark.gwas_oracle_slow
def test_oracle_pleiotropy_guardrail(tmp_path):
    bundle = generate_oracle_dataset(tmp_path, scenario="pleiotropy_guardrail", seed=41, n_proteins=12, snps_per_protein=8)
    artifacts = run_oracle_pipeline(tmp_path, bundle, force_proxy_branch=False, use_real_harmonise_r=False)
    summary = compute_oracle_mr_and_tiers(artifacts)
    _assert_non_empty(summary)

    pleio = summary[summary["class"] == "pleiotropic"]
    assert not pleio.empty
    promoted = pleio["tier"].isin(["Tier1", "Tier1_replicated", "Tier2"]).mean()
    assert promoted <= 0.25, f"Pleiotropic proteins promoted too often: {promoted:.3f}"


@pytest.mark.gwas_oracle_slow
def test_oracle_weak_and_null_control(tmp_path):
    bundle = generate_oracle_dataset(tmp_path, scenario="weak_instrument_guardrail", seed=51, n_proteins=12, snps_per_protein=8)
    artifacts = run_oracle_pipeline(tmp_path, bundle, force_proxy_branch=False, use_real_harmonise_r=False)
    summary = compute_oracle_mr_and_tiers(artifacts)
    _assert_non_empty(summary)

    manifest = bundle.manifest.copy()
    weak_expected = set(manifest.loc[manifest["class"] == "weak", "seqid"].astype(str))
    null_expected = set(manifest.loc[manifest["class"] == "null", "seqid"].astype(str))
    seen = set(summary["seqid"].astype(str))

    # Weak/null proteins may be filtered out before MR (expected); if present, they must remain deprioritized.
    weak = summary[summary["class"] == "weak"]
    if not weak.empty:
        weak_top = weak["tier"].isin(["Tier1", "Tier1_replicated", "Tier2"]).mean()
        assert weak_top <= 0.25, f"Weak instruments promoted too often: {weak_top:.3f}"
    else:
        assert weak_expected.isdisjoint(seen), "Weak class missing from summary but unexpectedly present as non-weak rows"

    nulls = summary[summary["class"] == "null"]
    if not nulls.empty:
        null_fp = nulls["fdr_pass"].mean()
        assert null_fp <= 0.25, f"Null false-positive rate too high: {null_fp:.3f}"
    else:
        assert null_expected.isdisjoint(seen), "Null class missing from summary but unexpectedly present as non-null rows"


@pytest.mark.gwas_oracle_slow
def test_oracle_full_02_to_09_integration(tmp_path):
    bundle = generate_oracle_dataset(tmp_path, scenario="full_pipeline", seed=61, n_proteins=12, snps_per_protein=8)
    final_path = run_full_oracle_02_to_09(tmp_path, bundle)

    final = pd.read_csv(final_path, sep="\t")
    assert not final.empty, "final_results.tsv is empty"
    assert "tier" in final.columns

    # Integrated stability check: tiered output should contain at least one non-Tier3 protein.
    assert (final["tier"] != "Tier3").any(), "No prioritized proteins found in full synthetic run"

    mr_path = final_path.parent / "ARIC_EA" / "mr_results.tsv"
    assert mr_path.exists(), "mr_results.tsv missing after full 02→09 run"
    mr = pd.read_csv(mr_path, sep="\t")
    assert "gene" in mr.columns, "mr_results.tsv missing gene column"
    assert mr["gene"].notna().all(), "mr_results.tsv contains missing gene metadata"


@pytest.mark.gwas_oracle_slow
def test_oracle_close_loci_duplicate_position_behavior(tmp_path, caplog):
    bundle = generate_oracle_dataset(
        tmp_path,
        scenario="shared_signal",
        seed=71,
        n_proteins=12,
        snps_per_protein=8,
        locus_spacing_bp=50_000,
    )
    with caplog.at_level("WARNING", logger="05_harmonise"):
        artifacts = run_oracle_pipeline(tmp_path, bundle, force_proxy_branch=False, use_real_harmonise_r=False)
    summary = compute_oracle_mr_and_tiers(artifacts)
    _assert_non_empty(summary)
    assert "duplicate" in caplog.text.lower()
