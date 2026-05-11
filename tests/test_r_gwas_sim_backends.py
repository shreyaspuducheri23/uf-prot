"""Fast GWASBrewer oracle smoke tests (always-on)."""

from __future__ import annotations

import pandas as pd

from tests.gwas_oracle_harness import (
    compute_oracle_mr_and_tiers,
    generate_oracle_dataset,
    run_oracle_pipeline,
)


def test_gwasbrewer_oracle_smoke_is_deterministic(tmp_path):
    b1 = generate_oracle_dataset(tmp_path, scenario="shared_signal", seed=123, n_proteins=12, snps_per_protein=8)
    b2 = generate_oracle_dataset(tmp_path, scenario="shared_signal", seed=123, n_proteins=12, snps_per_protein=8)

    pd.testing.assert_frame_equal(
        b1.manifest.sort_values("seqid").reset_index(drop=True),
        b2.manifest.sort_values("seqid").reset_index(drop=True),
    )
    pd.testing.assert_frame_equal(
        b1.exposure.sort_values(["seqid", "rsid"]).reset_index(drop=True),
        b2.exposure.sort_values(["seqid", "rsid"]).reset_index(drop=True),
    )


def test_gwasbrewer_oracle_smoke_tier_and_mr_calibration(tmp_path):
    bundle = generate_oracle_dataset(tmp_path, scenario="shared_signal", seed=7, n_proteins=12, snps_per_protein=8)
    artifacts = run_oracle_pipeline(tmp_path, bundle, force_proxy_branch=False, use_real_harmonise_r=False)
    summary = compute_oracle_mr_and_tiers(artifacts)

    causal = summary[summary["class"] == "causal"]
    assert not causal.empty, "Expected causal proteins in oracle manifest"

    top_tiers = {"Tier1_replicated", "Tier1", "Tier2"}
    causal_recall = causal["tier"].isin(top_tiers).mean()
    assert causal_recall >= 0.75, f"Causal top-tier recall too low: {causal_recall:.3f}"

    nulls = summary[summary["class"] == "null"]
    if not nulls.empty:
        null_top_rate = nulls["tier"].isin(top_tiers).mean()
        assert null_top_rate <= 0.25, f"Null promotion too high: {null_top_rate:.3f}"

    sign_ok = (causal["mr_beta"] * causal["expected_effect_sign"] > 0).mean()
    assert sign_ok >= 0.75, f"Causal sign agreement too low: {sign_ok:.3f}"

    cal_ok = (
        (causal["mr_beta"] - causal["target_mr_beta"]).abs()
        <= causal["beta_tolerance"]
    ).mean()
    assert cal_ok >= 0.75, f"Causal MR calibration too loose: {cal_ok:.3f}"


def test_gwasbrewer_oracle_smoke_proxy_mode_requires_proxy_rows(tmp_path):
    bundle = generate_oracle_dataset(tmp_path, scenario="proxy_branch", seed=21, n_proteins=12, snps_per_protein=8)
    artifacts = run_oracle_pipeline(tmp_path, bundle, force_proxy_branch=True, use_real_harmonise_r=False)
    summary = compute_oracle_mr_and_tiers(artifacts)

    expected_proxy = summary[summary["expected_proxy"] == True]
    assert not expected_proxy.empty, "Expected proxy-marked proteins in proxy_branch scenario"
    assert expected_proxy["proxy_used_any"].all(), "Proxy branch expected but proxy_used flag missing"
