"""Black-box CLI wiring test using a fake PLINK binary."""
import json
import os
import subprocess
import textwrap
from pathlib import Path

import pandas as pd


def _write_fake_plink2(path: Path) -> None:
    script = textwrap.dedent(
        """\
        #!/usr/bin/env python3
        import csv
        import sys

        argv = sys.argv
        if "--clump" not in argv:
            sys.exit(0)

        def arg(flag):
            return argv[argv.index(flag) + 1]

        assoc = arg("--clump")
        out_prefix = arg("--out")
        p1 = float(arg("--clump-p1"))

        kept = []
        with open(assoc, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\\t")
            for row in reader:
                if float(row["P"]) <= p1:
                    kept.append(row)

        if kept:
            with open(out_prefix + ".clumps", "w") as out:
                out.write("ID P\\n")
                for row in kept:
                    out.write(f"{row['SNP']} {row['P']}\\n")
        """
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(script)
    path.chmod(0o755)


def _write_synthetic_aric(raw_dir: Path) -> None:
    aric_dir = raw_dir / "ARIC"
    ea_dir = aric_dir / "EA"
    ea_dir.mkdir(parents=True, exist_ok=True)

    seqid = "SeqId_TEST"
    seqid_index = pd.DataFrame(
        {
            "seqid_in_sample": [seqid],
            "uniprot_id": ["P00001"],
            "entrezgenesymbol": ["GTEST"],
            "chromosome_name": ["1"],
            "transcription_start_site": [1_000_000],
        }
    )
    seqid_index.to_csv(aric_dir / "seqid.txt", sep="\t", index=False)

    # Two variants should survive cis+GW filters, one should fail p-value, one should fail cis window.
    df = pd.DataFrame(
        {
            "#CHROM": [1, 1, 1, 1],
            "POS": [1_000_000, 1_000_100, 1_000_200, 2_000_000],
            "ID": ["rs1", "rs2", "rs3", "rs4"],
            "A1": ["A", "A", "A", "A"],
            "REF": ["G", "G", "G", "G"],
            "A1_FREQ": [0.2, 0.25, 0.3, 0.2],
            "BETA": [0.5, 0.4, 0.3, 0.6],
            "SE": [0.05, 0.05, 0.05, 0.05],
            "P": [1e-12, 1e-9, 1e-6, 1e-12],
            "OBS_CT": [7213, 7213, 7213, 7213],
            "TEST": ["ADD", "ADD", "ADD", "ADD"],
        }
    )
    df.to_csv(ea_dir / f"{seqid}.PHENO1.glm.linear", sep="\t", index=False)


def _write_config(path: Path, clump_p1: float) -> None:
    cfg = {
        "_meta": {"version": "1.0"},
        "cis_extract": {
            "window_kb": 500,
            "pval_gw": 5e-8,
            "maf_min": 0.01,
            "palindrome_maf_max": 0.42,
        },
        "clump": {"window_kb": 1000, "r2": 0.001, "p1": clump_p1},
        "fstat": {"weak_threshold": 10.0},
        "harmonise": {"maf_proxy_max": 0.42, "proxy_r2_min": 0.8},
        "outcome": {"kim_N": 434152},
        "mhc": {"hg19": [25000000, 34000000], "hg38": [28500000, 33500000]},
        "cohorts": {
            "ARIC_EA": {"N": None, "build": "hg38"},
            "Fenland": {"N": 10708, "build": "hg19"},
            "deCODE": {"N_default": 35000, "build": "hg38"},
            "UKB_PPP": {"N": None, "build": "hg19"},
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cfg))


def _run_blackbox(tmp_path: Path, clump_p1: float) -> int:
    run_root = tmp_path / f"blackbox_{str(clump_p1).replace('.', '_')}"
    raw_dir = run_root / "raw"
    processed_dir = run_root / "processed"
    logs_dir = run_root / "logs"
    ld_ref_dir = run_root / "ld_ref" / "ld_files"
    config_path = run_root / "config" / "pipeline.json"
    plink2_path = run_root / "bin" / "plink2"

    _write_synthetic_aric(raw_dir)
    _write_config(config_path, clump_p1)
    _write_fake_plink2(plink2_path)
    ld_ref_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env.update(
        {
            "LEIO_RAW_DIR": str(raw_dir),
            "LEIO_PROCESSED_DIR": str(processed_dir),
            "LEIO_LOGS_DIR": str(logs_dir),
            "LEIO_LD_REF_DIR": str(ld_ref_dir),
            "PATH": f"{plink2_path.parent}:{os.environ.get('PATH', '')}",
        }
    )

    cmd_extract = [
        "uv",
        "run",
        "python",
        "scripts/02_cis_pqtl_extract/aric.py",
        "--limit",
        "1",
        "--config",
        str(config_path),
    ]
    cmd_clump = [
        "uv",
        "run",
        "python",
        "scripts/03_clump/clump.py",
        "--cohort",
        "ARIC_EA",
        "--config",
        str(config_path),
    ]
    res_extract = subprocess.run(cmd_extract, capture_output=True, text=True, env=env)
    assert res_extract.returncode == 0, res_extract.stderr + "\n" + res_extract.stdout
    res_clump = subprocess.run(cmd_clump, capture_output=True, text=True, env=env)
    assert res_clump.returncode == 0, res_clump.stderr + "\n" + res_clump.stdout

    cis_path = processed_dir / "ARIC_EA" / "cis_sumstats" / "SeqId_TEST.tsv"
    assert cis_path.exists()
    cis_df = pd.read_csv(cis_path, sep="\t")
    assert set(cis_df["rsid"]) == {"rs1", "rs2"}

    inst_path = processed_dir / "ARIC_EA" / "instruments" / "SeqId_TEST.tsv"
    assert inst_path.exists()
    inst_df = pd.read_csv(inst_path, sep="\t")
    return len(inst_df)


def test_fake_plink_blackbox_cli_respects_clump_p1_config_wiring_only(tmp_path):
    """This test validates config plumbing into clump CLI, not real PLINK LD semantics."""
    strict_n = _run_blackbox(tmp_path, clump_p1=1e-10)
    relaxed_n = _run_blackbox(tmp_path, clump_p1=1e-8)
    assert strict_n == 1
    assert relaxed_n == 2
