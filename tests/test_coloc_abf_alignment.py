"""Synthetic tests for coloc.abf allele/key alignment."""
import shutil
import subprocess

import pandas as pd
import pytest


pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")


def _run_align(tmp_path, exp: pd.DataFrame, out: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    exp_path = tmp_path / "exposure.tsv"
    out_path = tmp_path / "outcome.tsv"
    aligned_exp_path = tmp_path / "aligned_exposure.tsv"
    aligned_out_path = tmp_path / "aligned_outcome.tsv"
    exp.to_csv(exp_path, sep="\t", index=False)
    out.to_csv(out_path, sep="\t", index=False)

    code = f"""
source("scripts/rlib/coloc_align.R")
exp_df <- read.delim("{exp_path}", stringsAsFactors = FALSE, check.names = FALSE)
out_df <- read.delim("{out_path}", stringsAsFactors = FALSE, check.names = FALSE)
aligned <- align_coloc_region(exp_df, out_df)
exp_out <- aligned$exp
out_out <- aligned$out
exp_out$snp_key <- aligned$snp
out_out$snp_key <- aligned$snp
write.table(exp_out, "{aligned_exp_path}", sep = "\\t", row.names = FALSE, quote = FALSE)
write.table(out_out, "{aligned_out_path}", sep = "\\t", row.names = FALSE, quote = FALSE)
"""
    res = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    return (
        pd.read_csv(aligned_exp_path, sep="\t"),
        pd.read_csv(aligned_out_path, sep="\t"),
    )


def _exp(snps=None, alleles=None, positions=None) -> pd.DataFrame:
    snps = snps or [f"rs{i}" for i in range(1, 7)]
    alleles = alleles or [("A", "G")] * len(snps)
    positions = positions or list(range(100, 100 + len(snps)))
    return pd.DataFrame({
        "chrom": ["1"] * len(snps),
        "pos": positions,
        "rsid": snps,
        "EA": [ea for ea, _oa in alleles],
        "OA": [oa for _ea, oa in alleles],
        "beta": [0.2] * len(snps),
        "se": [0.05] * len(snps),
        "EAF": [0.4] * len(snps),
        "N": [1000] * len(snps),
    })


def _out(snps=None, alleles=None, positions=None, betas=None, eafs=None) -> pd.DataFrame:
    snps = snps or [f"rs{i}" for i in range(1, 7)]
    alleles = alleles or [("A", "G")] * len(snps)
    positions = positions or list(range(100, 100 + len(snps)))
    betas = betas or [0.7] * len(snps)
    eafs = eafs or [0.8] * len(snps)
    return pd.DataFrame({
        "chromosome": ["1"] * len(snps),
        "base_pair_location": positions,
        "rsid": snps,
        "effect_allele": [ea for ea, _oa in alleles],
        "other_allele": [oa for _ea, oa in alleles],
        "beta": betas,
        "standard_error": [0.03] * len(snps),
        "effect_allele_frequency": eafs,
    })


def test_coloc_abf_forward_matches_preserve_outcome_effects(tmp_path):
    aligned_exp, aligned_out = _run_align(tmp_path, _exp(), _out())

    assert aligned_exp["snp_key"].tolist() == aligned_out["snp_key"].tolist()
    assert aligned_out["beta"].tolist() == pytest.approx([0.7] * 6)
    assert aligned_out["effect_allele_frequency"].tolist() == pytest.approx([0.8] * 6)


def test_coloc_abf_reverse_matches_flip_outcome_effects(tmp_path):
    aligned_exp, aligned_out = _run_align(
        tmp_path,
        _exp(alleles=[("A", "G")] * 6),
        _out(alleles=[("G", "A")] * 6),
    )

    assert len(aligned_exp) == 6
    assert aligned_out["beta"].tolist() == pytest.approx([-0.7] * 6)
    assert aligned_out["effect_allele_frequency"].tolist() == pytest.approx([0.2] * 6)


def test_coloc_abf_excludes_incompatible_same_rsid_variants(tmp_path):
    aligned_exp, aligned_out = _run_align(
        tmp_path,
        _exp(alleles=[("A", "G")] * 6),
        _out(alleles=[("C", "T")] * 6),
    )

    assert aligned_exp.empty
    assert aligned_out.empty


def test_coloc_abf_matches_by_position_and_alleles_not_rsid_only(tmp_path):
    aligned_exp, aligned_out = _run_align(
        tmp_path,
        _exp(snps=[f"exp_rs{i}" for i in range(1, 7)]),
        _out(snps=[f"out_rs{i}" for i in range(1, 7)]),
    )

    assert len(aligned_exp) == 6
    assert aligned_exp["rsid"].tolist() != aligned_out["rsid"].tolist()
    assert aligned_exp["snp_key"].tolist() == aligned_out["snp_key"].tolist()


def test_coloc_abf_duplicate_keys_keep_first_and_align_order(tmp_path):
    exp = _exp(
        snps=["rs_dup_exp", "rs_dup_exp_late", "rs_unique1", "rs_unique2", "rs_unique3", "rs_unique4"],
        positions=[100, 100, 101, 102, 103, 104],
    )
    exp.loc[1, "beta"] = 9.9
    out = _out(
        snps=["rs_dup_out", "rs_dup_out_late", "rs_unique1", "rs_unique2", "rs_unique3", "rs_unique4"],
        positions=[100, 100, 101, 102, 103, 104],
        betas=[0.7, 9.9, 0.71, 0.72, 0.73, 0.74],
    )

    aligned_exp, aligned_out = _run_align(tmp_path, exp, out)

    assert aligned_exp["rsid"].tolist()[0] == "rs_dup_exp"
    assert aligned_out["rsid"].tolist()[0] == "rs_dup_out"
    assert aligned_exp["beta"].tolist()[0] == pytest.approx(0.2)
    assert aligned_out["beta"].tolist()[0] == pytest.approx(0.7)
    assert aligned_exp["snp_key"].tolist() == aligned_out["snp_key"].tolist()
