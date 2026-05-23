#!/usr/bin/env python3
"""
02_cis_pqtl_extract/aric.py
Extract cis-pQTLs from ARIC EA PLINK2 .glm.linear files.

Usage:
  python scripts/02_cis_pqtl_extract/aric.py [--limit N] [--workers N]
"""
import argparse
import glob
import logging
from pathlib import Path

import pandas as pd

from scripts.lib.cis import load_aric_tss
from scripts.lib.config import add_config_arg, load_config, get_section
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import ARIC_EA_DIR, ARIC_SEQID
from scripts.lib.schema import ProteinMeta
from scripts.lib.cis_extract import run_extraction

log = setup_logger("02_aric")

COHORT = "ARIC_EA"
BUILD = "hg38"  # seqid.txt TSS and .glm.linear positions are both hg38
ARIC_N = 7_213  # approximate OBS_CT


def load_aric_proteins() -> list[ProteinMeta]:
    tss_map = load_aric_tss(ARIC_SEQID)  # {seqid: (chrom, tss, uniprot, gene)}
    proteins = []
    for seqid, (chrom, tss, uniprot, gene) in tss_map.items():
        proteins.append(ProteinMeta(
            seqid=seqid, gene=gene, uniprot=uniprot,
            chrom=chrom, tss=tss, build=BUILD, source_cohort=COHORT,
        ))
    log.info(f"ARIC EA: {len(proteins)} proteins from seqid.txt")
    return proteins


def read_aric_protein(protein: ProteinMeta) -> pd.DataFrame | None:
    """Read one ARIC EA .glm.linear file and normalize columns."""
    pattern = str(ARIC_EA_DIR / f"{protein.seqid}.PHENO1.glm.linear")
    matches = glob.glob(pattern)
    if not matches:
        log.warning(f"{protein.seqid}: no PLINK2 file found at {pattern}")
        return None
    if len(matches) > 1:
        log.warning(f"{protein.seqid}: multiple files matched {pattern} — using first: {matches[0]}")

    # PLINK2 header starts with #CHROM — read without comment skipping,
    # then strip the '#' from the column name.
    df = pd.read_csv(matches[0], sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    df = df[df["TEST"] == "ADD"].copy()
    df = df.rename(columns={
        "CHROM": "chrom",
        "POS": "pos",
        "ID": "rsid",
        "A1": "EA",
        "REF": "OA",
        "A1_FREQ": "EAF",
        "BETA": "beta",
        "SE": "se",
        "P": "pval",
        "OBS_CT": "N",
    })
    # REF is OA (A1 is effect allele)
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = df["pos"].astype(int)
    df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
    df["N"] = df["N"].fillna(ARIC_N).astype(int)

    # rsid: use ID column if starts with 'rs', else '.'
    df["rsid"] = df["rsid"].where(df["rsid"].str.startswith("rs"), ".")

    return df


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract ARIC EA cis-pQTLs")
    parser.add_argument("--limit", type=int, default=None,
                        help="Process at most N proteins (for testing)")
    add_config_arg(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    cis_cfg = get_section(cfg, "cis_extract")

    with RunManifest("02_cis_pqtl_extract/aric.py") as manifest:
        proteins = load_aric_proteins()
        n = run_extraction(COHORT, proteins, read_aric_protein, limit=args.limit, cfg=cis_cfg)
        manifest.n_units = n


if __name__ == "__main__":
    main()
