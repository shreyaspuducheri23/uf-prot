"""Real end-to-end smoke test using actual pipeline CLIs and external tools.

This test intentionally avoids monkeypatching PLINK, liftover, outcome lookup,
or R harmonisation. It builds a tiny synthetic ARIC input and Kim outcome file,
then runs steps 02→05 as black-box commands.
"""

import importlib
import json
import os
import subprocess
from pathlib import Path

import pandas as pd
import pysam
import pytest

from scripts.lib.fdr import add_fdr
from scripts.lib.schema import NORM_COLS

_assemble_mod = importlib.import_module("scripts.09_assemble.assemble")


PROTEINS = [
    {
        "seqid": "SeqId_P1",
        "gene": "GENE1",
        "uniprot": "Q11111",
        "chrom": "2",
        "tss": 136_608_646,
    },
    {
        "seqid": "SeqId_P2",
        "gene": "GENE2",
        "uniprot": "Q22222",
        "chrom": "13",
        "tss": 32_914_437,
    },
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _run(cmd: list[str], env: dict[str, str], cwd: Path) -> None:
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=cwd)
    assert proc.returncode == 0, (
        f"Command failed: {' '.join(cmd)}\n"
        f"STDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}"
    )


def _load_rsids(ld_ref_dir: Path, n: int = 2) -> list[str]:
    snplist = ld_ref_dir / "data_maf0.01_rs.snplist"
    rsids: list[str] = []
    with snplist.open() as fh:
        for line in fh:
            rsid = line.strip()
            if rsid.startswith("rs"):
                rsids.append(rsid)
            if len(rsids) == n:
                break
    if len(rsids) < n:
        raise RuntimeError(f"Expected at least {n} rsIDs in {snplist}")
    return rsids


def _write_config(path: Path) -> None:
    cfg = {
        "_meta": {"version": "1.0"},
        "cis_extract": {
            "window_kb": 500,
            "pval_gw": 5e-8,
            "maf_min": 0.01,
            "palindrome_maf_max": 0.42,
        },
        "clump": {"window_kb": 1000, "r2": 0.001, "p1": 5e-8},
        "fstat": {"weak_threshold": 10.0},
        "harmonise": {"maf_proxy_max": 0.42, "proxy_r2_min": 0.8},
        "outcome": {"kim_N": 434152},
        "mhc": {"hg19": [25000000, 34000000], "hg38": [28500000, 33500000]},
        "cohorts": {
            "ARIC_EA": {"N": None, "build": "hg38"},
            "Fenland": {"N": 10708, "build": "hg19"},
            "deCODE": {"N_default": 35000, "build": "hg38"},
            "UKB_PPP": {"N": None, "build": "hg19"},
            "UKB_female": {"N": None, "build": "hg19"},
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cfg))


def _write_synthetic_aric(raw_dir: Path, rsids: list[str]) -> None:
    aric_dir = raw_dir / "ARIC"
    ea_dir = aric_dir / "EA"
    ea_dir.mkdir(parents=True, exist_ok=True)

    seqid_df = pd.DataFrame(
        {
            "seqid_in_sample": [p["seqid"] for p in PROTEINS],
            "uniprot_id": [p["uniprot"] for p in PROTEINS],
            "entrezgenesymbol": [p["gene"] for p in PROTEINS],
            "chromosome_name": [p["chrom"] for p in PROTEINS],
            "transcription_start_site": [p["tss"] for p in PROTEINS],
        }
    )
    seqid_df.to_csv(aric_dir / "seqid.txt", sep="\t", index=False)

    for protein, rsid in zip(PROTEINS, rsids):
        # One valid GW-significant cis SNP and one non-significant SNP.
        df = pd.DataFrame(
            {
                "#CHROM": [int(protein["chrom"]), int(protein["chrom"])],
                "POS": [protein["tss"], protein["tss"] + 1_000],
                "ID": [rsid, "rs999999999"],
                "A1": ["A", "A"],
                "REF": ["G", "G"],
                "A1_FREQ": [0.25, 0.25],
                "BETA": [0.5, 0.1],
                "SE": [0.05, 0.05],
                "P": [1e-10, 1e-4],
                "OBS_CT": [7213, 7213],
                "TEST": ["ADD", "ADD"],
            }
        )
        df.to_csv(ea_dir / f"{protein['seqid']}.PHENO1.glm.linear", sep="\t", index=False)


def _write_synthetic_kim(raw_dir: Path, instruments_hg38_dir: Path) -> Path:
    kim_dir = raw_dir / "kim_fibroid_gwas"
    kim_dir.mkdir(parents=True, exist_ok=True)

    rows: list[list[object]] = []
    for tsv in sorted(instruments_hg38_dir.glob("*.tsv")):
        instr = pd.read_csv(tsv, sep="\t")
        for _, row in instr.iterrows():
            rows.append(
                [
                    str(row["chrom_hg38"]).lstrip("chr"),
                    int(row["pos_hg38"]),
                    "A",  # effect_allele
                    "G",  # other_allele
                    0.18,  # beta
                    0.02,  # standard_error
                    0.32,  # effect_allele_frequency
                    1e-6,  # p_value
                    str(row["rsid"]),
                    str(row["rsid"]),
                    "",  # hm_coordinate_conversion
                    "",  # hm_code
                    "",  # variant_id
                ]
            )

    df = pd.DataFrame(
        rows,
        columns=[
            "chromosome",
            "base_pair_location",
            "effect_allele",
            "other_allele",
            "beta",
            "standard_error",
            "effect_allele_frequency",
            "p_value",
            "rsid",
            "rs_id",
            "hm_coordinate_conversion",
            "hm_code",
            "variant_id",
        ],
    )
    df = df.sort_values(["chromosome", "base_pair_location"], kind="stable")

    plain_tsv = kim_dir / "GCST90461958.h.tsv"
    df.to_csv(plain_tsv, sep="\t", index=False, header=False)

    gz_path = Path(
        pysam.tabix_index(
            str(plain_tsv),
            seq_col=0,
            start_col=1,
            end_col=1,
            force=True,
            zerobased=False,
        )
    )
    assert gz_path.exists()
    assert (gz_path.parent / f"{gz_path.name}.tbi").exists()
    return gz_path


@pytest.fixture(scope="module")
def real_run(tmp_path_factory):
    repo_root = _repo_root()
    py = repo_root / ".venv" / "bin" / "python"

    ld_ref_dir = repo_root / "data" / "ld_ref" / "ld_files"
    plink2_fallback = Path("/Users/spuduch/Research/MR_IA/plink2_mac_arm64_20260228/plink2")
    ld_prefix = ld_ref_dir / "data_maf0.01_rs"
    chain = repo_root / "data" / "ref" / "hg19ToHg38.over.chain.gz"

    if not py.exists():
        pytest.fail("Missing .venv/bin/python")
    if not Path(f"{ld_prefix}.bed").exists():
        pytest.fail("LD reference not found at README paths under data/ld_ref/ld_files")
    if not chain.exists():
        pytest.fail("Liftover chain file missing: data/ref/hg19ToHg38.over.chain.gz")
    try:
        plink_probe = subprocess.run(["plink2", "--help"], capture_output=True, text=True)
    except OSError:
        plink_probe = None
    if plink_probe is None or plink_probe.returncode not in {0, 1}:
        if not plink2_fallback.exists():
            pytest.fail(
                "plink2 is not available on PATH and fallback binary not found "
                "at /Users/spuduch/Research/MR_IA/plink2_mac_arm64_20260228/plink2"
            )

    r_check = subprocess.run(
        ["Rscript", "-e", 'quit(save="no", status=ifelse(requireNamespace("TwoSampleMR", quietly=TRUE),0,42))'],
        capture_output=True,
        text=True,
    )
    if r_check.returncode == 42:
        pytest.fail("TwoSampleMR is not installed in this environment")
    if r_check.returncode != 0:
        pytest.fail("Rscript unavailable for harmonisation")

    run_root = tmp_path_factory.mktemp("e2e_real")
    raw_dir = run_root / "raw"
    processed_dir = run_root / "processed"
    logs_dir = run_root / "logs"
    config_path = run_root / "config" / "pipeline.json"

    rsids = _load_rsids(ld_ref_dir, n=2)
    _write_synthetic_aric(raw_dir, rsids)
    _write_config(config_path)

    env = os.environ.copy()
    env.update(
        {
            "LEIO_ROOT": str(repo_root),
            "LEIO_RAW_DIR": str(raw_dir),
            "LEIO_PROCESSED_DIR": str(processed_dir),
            "LEIO_LOGS_DIR": str(logs_dir),
            "LEIO_LD_REF_DIR": str(ld_ref_dir),
            "LEIO_REF_DIR": str(repo_root / "data" / "ref"),
        }
    )
    if plink_probe is None or plink_probe.returncode not in {0, 1}:
        env["PATH"] = f"{plink2_fallback.parent}:{env.get('PATH', '')}"

    _run(
        [
            str(py),
            "scripts/02_cis_pqtl_extract/aric.py",
            "--limit",
            "2",
            "--config",
            str(config_path),
        ],
        env=env,
        cwd=repo_root,
    )
    _run(
        [
            str(py),
            "scripts/03_clump/clump.py",
            "--cohort",
            "ARIC_EA",
            "--config",
            str(config_path),
        ],
        env=env,
        cwd=repo_root,
    )
    _run(
        [
            str(py),
            "scripts/04_liftover/instruments_to_hg38.py",
            "--cohort",
            "ARIC_EA",
            "--config",
            str(config_path),
        ],
        env=env,
        cwd=repo_root,
    )

    instruments_hg38_dir = processed_dir / "ARIC_EA" / "instruments_hg38"
    _write_synthetic_kim(raw_dir, instruments_hg38_dir)

    _run(
        [
            str(py),
            "scripts/05_harmonise/harmonise.py",
            "--cohort",
            "ARIC_EA",
            "--config",
            str(config_path),
        ],
        env=env,
        cwd=repo_root,
    )

    return {
        "extracted_dir": processed_dir / "ARIC_EA" / "cis_sumstats",
        "harmonised_dir": processed_dir / "ARIC_EA" / "harmonised",
    }


class TestEndToEndSmoke:
    def test_pipeline_produces_results_for_both_proteins(self, real_run):
        files = list(real_run["harmonised_dir"].glob("*.tsv"))
        seqids = {f.stem for f in files}
        assert "SeqId_P1" in seqids
        assert "SeqId_P2" in seqids

    def test_fdr_tiering_produces_expected_columns(self, real_run):
        tier = _assemble_mod.tier

        frames = []
        for tsv in real_run["harmonised_dir"].glob("*.tsv"):
            df = pd.read_csv(tsv, sep="\t")
            df["seqid"] = tsv.stem
            frames.append(df)

        mr = pd.concat(frames, ignore_index=True)
        mr["pval"] = mr.get("pval.exposure", pd.Series([1e-9] * len(mr)))
        mr = add_fdr(mr, pval_col="pval", alpha=0.05)
        mr["passes_sensitivity"] = True
        mr["sharepro_coloc_positive"] = True
        mr["coloc_abf_positive"] = True
        mr["tier"] = mr.apply(tier, axis=1)

        assert "tier" in mr.columns
        assert len(mr) == 2
        assert all(mr["tier"] == "Tier1_replicated")

    def test_extraction_output_has_required_norm_cols(self, real_run):
        for tsv in real_run["extracted_dir"].glob("*.tsv"):
            df = pd.read_csv(tsv, sep="\t")
            missing = [c for c in NORM_COLS if c not in df.columns]
            assert not missing, f"{tsv.name}: missing columns {missing}"
