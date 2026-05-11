"""Shared GWASBrewer oracle fixtures for deterministic pipeline validation."""
from __future__ import annotations

import importlib
import math
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pandas as pd

from scripts.lib.cis_extract import run_extraction
from scripts.lib.fdr import add_fdr
from scripts.lib.schema import ProteinMeta

_clump_mod = importlib.import_module("scripts.03_clump.clump")
_liftover_mod = importlib.import_module("scripts.04_liftover.instruments_to_hg38")
_harm_mod = importlib.import_module("scripts.05_harmonise.harmonise")
_assemble_mod = importlib.import_module("scripts.09_assemble.assemble")


@dataclass(frozen=True)
class OracleBundle:
    scenario: str
    seed: int
    root: Path
    exposure: pd.DataFrame
    outcome: pd.DataFrame
    manifest: pd.DataFrame


def _norm_sf(z: float) -> float:
    return 0.5 * math.erfc(z / math.sqrt(2))


def generate_oracle_dataset(
    tmp_path: Path,
    *,
    scenario: str,
    seed: int = 42,
    n_proteins: int = 12,
    snps_per_protein: int = 8,
    locus_spacing_bp: int = 2_000_000,
) -> OracleBundle:
    """Generate deterministic exposure/outcome GWAS and oracle manifest via GWASBrewer."""
    root = tmp_path / f"oracle_{scenario}_{seed}"
    root.mkdir(parents=True, exist_ok=True)

    script_path = Path(__file__).parent / "r" / "gwasbrewer_oracle_generate.R"
    cmd = [
        "Rscript",
        str(script_path),
        str(root),
        scenario,
        str(seed),
        str(n_proteins),
        str(snps_per_protein),
        str(locus_spacing_bp),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"GWASBrewer generator failed\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")

    exposure = pd.read_csv(root / "exposure_gwas.tsv", sep="\t", keep_default_na=False)
    outcome = pd.read_csv(root / "outcome_gwas.tsv", sep="\t", keep_default_na=False)
    manifest = pd.read_csv(root / "oracle_manifest.tsv", sep="\t", keep_default_na=False)

    if exposure.empty or outcome.empty or manifest.empty:
        raise ValueError("Generator produced empty outputs")

    return OracleBundle(
        scenario=scenario,
        seed=seed,
        root=root,
        exposure=exposure,
        outcome=outcome,
        manifest=manifest,
    )


def _build_proteins(manifest: pd.DataFrame) -> list[ProteinMeta]:
    proteins: list[ProteinMeta] = []
    for row in manifest.itertuples(index=False):
        proteins.append(
            ProteinMeta(
                seqid=str(row.seqid),
                gene=str(row.gene),
                uniprot=str(row.uniprot),
                chrom=str(row.chrom),
                tss=int(row.tss),
                build="hg19",
                source_cohort="ARIC_EA",
            )
        )
    return proteins


def run_oracle_pipeline(
    tmp_path: Path,
    bundle: OracleBundle,
    *,
    force_proxy_branch: bool = False,
    use_real_harmonise_r: bool = False,
) -> dict[str, Any]:
    """Run deterministic 02→05 pipeline stages using oracle-generated synthetic data."""
    cohort = "ARIC_EA"
    run_root = tmp_path / f"run_{bundle.scenario}_{bundle.seed}_{'proxy' if force_proxy_branch else 'direct'}"
    cohort_root = run_root / "processed_data" / cohort

    proteins = _build_proteins(bundle.manifest)
    exp_by_seqid = {k: g.copy() for k, g in bundle.exposure.groupby("seqid", sort=False)}

    def read_fn(protein: ProteinMeta) -> pd.DataFrame:
        df = exp_by_seqid.get(protein.seqid)
        if df is None:
            return pd.DataFrame()
        cols = ["chrom", "pos", "rsid", "EA", "OA", "EAF", "beta", "se", "pval", "N"]
        out = df[cols].copy()
        out["chrom"] = out["chrom"].astype(str)
        out["pos"] = out["pos"].astype(int)
        return out

    cis_cfg = {
        "window_kb": 500,
        "pval_gw": 5e-8,
        "maf_min": 0.01,
        "palindrome_maf_max": 0.42,
    }

    cis_dir = cohort_root / "cis_sumstats"
    with patch("scripts.lib.cis_extract.cis_sumstats_dir", return_value=cis_dir), patch(
        "scripts.lib.cis_extract.cohort_dir", return_value=cohort_root
    ):
        run_extraction(cohort, proteins, read_fn, cfg=cis_cfg)

    def fake_clump(df: pd.DataFrame, seqid: str, **kwargs) -> pd.DataFrame:
        return df.sort_values("pval", kind="stable").head(3).copy()

    inst_dir = cohort_root / "instruments"
    with patch.object(_clump_mod, "cis_sumstats_dir", return_value=cis_dir), patch.object(
        _clump_mod, "instruments_dir", return_value=inst_dir
    ), patch.object(_clump_mod, "cohort_dir", return_value=cohort_root), patch.object(
        _clump_mod, "clump", side_effect=fake_clump
    ):
        _clump_mod.clump_cohort(cohort, window_kb=1000, r2=0.001, p1=5e-8)

    hg38_dir = cohort_root / "instruments_hg38"

    def fake_lift_table(df: pd.DataFrame, chrom_col: str = "chrom", pos_col: str = "pos", **kwargs):
        lifted = df.copy()
        lifted["chrom_hg38"] = lifted[chrom_col].astype(str)
        lifted["pos_hg38"] = lifted[pos_col].astype(int)
        return lifted

    with patch.object(_liftover_mod, "instruments_dir", return_value=inst_dir), patch.object(
        _liftover_mod, "instruments_hg38_dir", return_value=hg38_dir
    ), patch.object(_liftover_mod, "cohort_dir", return_value=cohort_root), patch.object(
        _liftover_mod, "lift_table", side_effect=fake_lift_table
    ):
        _liftover_mod.lift_cohort(cohort)

    harm_dir = cohort_root / "harmonised"

    direct_lookup: dict[tuple[str, int], list[dict[str, Any]]] = {}
    by_rsid: dict[str, dict[str, Any]] = {}
    for row in bundle.outcome.to_dict(orient="records"):
        key = (str(row["chromosome"]), int(row["base_pair_location"]))
        direct_lookup.setdefault(key, []).append(row)
        by_rsid[str(row["rsid"])] = row

    proxy_targets: set[str] = set()
    proxy_map: dict[str, tuple[str, float]] = {}
    proxy_rows: dict[str, dict[str, Any]] = {}
    if force_proxy_branch:
        proxy_seqids = set(
            bundle.manifest.loc[bundle.manifest["expected_proxy"] == True, "seqid"].astype(str).tolist()
        )
        for row in bundle.exposure.to_dict(orient="records"):
            rsid = str(row["rsid"])
            if str(row["seqid"]) not in proxy_seqids:
                continue
            proxy_rsid = f"{rsid}_proxy"
            proxy_targets.add(rsid)
            proxy_map[rsid] = (proxy_rsid, 0.92)
            src = by_rsid[rsid]
            proxy_rows[proxy_rsid] = {
                **src,
                "rsid": proxy_rsid,
                "rs_id": proxy_rsid,
                "base_pair_location": int(src["base_pair_location"]) + 1,
            }

    class FakeOutcome:
        def fetch_snps(self, positions):
            rows = []
            for chrom, pos in positions:
                matched = direct_lookup.get((str(chrom), int(pos)), [])
                if not matched:
                    continue
                for row in matched:
                    if force_proxy_branch and str(row["rsid"]) in proxy_targets:
                        continue
                    rows.append(row)
            return pd.DataFrame(rows)

        def fetch_by_rsid(self, rsids):
            rows = [proxy_rows[r] for r in rsids if r in proxy_rows]
            return pd.DataFrame(rows)

        def __enter__(self):
            return self

        def __exit__(self, *_):
            return None

    def fake_harmonise_r(df: pd.DataFrame, seqid: str) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "seqid": seqid,
                "rsid": df["rsid"],
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
                "proxy_used": df.get("proxy_used", False),
            }
        )

    patchers = [
        patch.object(_harm_mod, "instruments_hg38_dir", return_value=hg38_dir),
        patch.object(_harm_mod, "harmonised_dir", return_value=harm_dir),
        patch.object(_harm_mod, "cohort_dir", return_value=cohort_root),
        patch.object(_harm_mod, "OutcomeLookup", return_value=FakeOutcome()),
        patch.object(_harm_mod, "find_proxies", return_value=proxy_map if force_proxy_branch else {}),
        patch.object(_harm_mod, "in_phase_allele_map", return_value={"A": "A", "G": "G"}),
    ]
    if not use_real_harmonise_r:
        patchers.append(patch.object(_harm_mod, "_call_harmonise_r", side_effect=fake_harmonise_r))

    with patchers[0], patchers[1], patchers[2], patchers[3], patchers[4], patchers[5]:
        if use_real_harmonise_r:
            _harm_mod.harmonise_cohort(cohort)
        else:
            with patchers[6]:
                _harm_mod.harmonise_cohort(cohort)

    return {
        "run_root": run_root,
        "cohort_root": cohort_root,
        "harmonised_dir": harm_dir,
        "manifest": bundle.manifest.copy(),
    }


def compute_oracle_mr_and_tiers(run_artifacts: dict[str, Any]) -> pd.DataFrame:
    """Compute deterministic MR-like summaries + tiers for oracle assertions."""
    harm_dir: Path = run_artifacts["harmonised_dir"]
    manifest = run_artifacts["manifest"].copy()
    cohort = "ARIC_EA"

    rows: list[dict[str, Any]] = []
    for path in sorted(harm_dir.glob("*.tsv")):
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue

        beta_exp = pd.to_numeric(df["beta_exp"], errors="coerce")
        beta_out = pd.to_numeric(df["beta_out"], errors="coerce")
        se_out = pd.to_numeric(df["se_out"], errors="coerce")

        valid = pd.DataFrame(
            {
                "bx": pd.to_numeric(beta_exp, errors="coerce"),
                "by": pd.to_numeric(beta_out, errors="coerce"),
                "se_y": pd.to_numeric(se_out, errors="coerce"),
            }
        ).dropna()
        valid = valid[valid["bx"] != 0]
        valid = valid[valid["se_y"] > 0]
        if valid.empty:
            mr_beta = float("nan")
            pval = 1.0
        else:
            w = 1.0 / (valid["se_y"] ** 2)
            denom = float((w * (valid["bx"] ** 2)).sum())
            if denom <= 0:
                mr_beta = float("nan")
                pval = 1.0
            else:
                mr_beta = float((w * valid["bx"] * valid["by"]).sum() / denom)
                se_ivw = math.sqrt(1.0 / denom)
                z = abs(mr_beta) / max(se_ivw, 1e-12)
                pval = 2 * _norm_sf(z)

        class_row = manifest[manifest["seqid"] == path.stem]
        if class_row.empty:
            continue
        class_label = str(class_row.iloc[0]["class"])
        coloc_expected = bool(class_row.iloc[0]["expected_coloc_positive"])

        rows.append(
            {
                "seqid": path.stem,
                "class": class_label,
                "cohort": cohort,
                "pval": pval,
                "mr_beta": mr_beta,
                "mean_F": float(((beta_exp / pd.to_numeric(df["se_exp"], errors="coerce")) ** 2).mean()),
                "passes_sensitivity": class_label in {"causal", "null"},
                "sharepro_coloc_positive": coloc_expected and class_label == "causal",
                "coloc_abf_positive": coloc_expected and class_label == "causal",
                "proxy_used_any": bool(df.get("proxy_used", pd.Series(False)).astype(bool).any()),
            }
        )

    summary = pd.DataFrame(rows)
    if summary.empty:
        raise ValueError("No harmonised outputs were produced for oracle scenario")

    summary = add_fdr(summary, pval_col="pval", alpha=0.05)
    summary["tier"] = summary.apply(_assemble_mod.tier, axis=1)
    return summary.merge(manifest, on="seqid", how="left", suffixes=("", "_oracle"))


def run_full_oracle_02_to_09(tmp_path: Path, bundle: OracleBundle) -> Path:
    """Run an integrated 02→09 synthetic flow (02-05 synthetic + real 06/07/09 scripts)."""
    artifacts = run_oracle_pipeline(
        tmp_path,
        bundle,
        force_proxy_branch=False,
        use_real_harmonise_r=True,
    )
    run_root: Path = artifacts["run_root"]

    repo_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env.update(
        {
            "LEIO_ROOT": str(repo_root),
            "LEIO_PROCESSED_DIR": str(run_root / "processed_data"),
            "LEIO_LOGS_DIR": str(run_root / "logs"),
        }
    )
    cmds = [
        ["Rscript", str(repo_root / "scripts" / "06_mr" / "run_mr.R"), "--cohort", "ARIC_EA"],
        ["Rscript", str(repo_root / "scripts" / "07_sensitivity" / "run_sensitivity.R"), "--cohort", "ARIC_EA"],
        ["python", str(repo_root / "scripts" / "09_assemble" / "assemble.py")],
    ]

    for cmd in cmds:
        proc = subprocess.run(cmd, cwd=run_root, env=env, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"Command failed: {' '.join(cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
            )

    final_path = run_root / "processed_data" / "final_results.tsv"
    if not final_path.exists() or final_path.stat().st_size == 0:
        raise ValueError("final_results.tsv missing or empty after 09 assemble")
    return final_path
