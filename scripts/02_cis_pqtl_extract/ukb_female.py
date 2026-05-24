#!/usr/bin/env python3
"""
02_cis_pqtl_extract/ukb_female.py
Phase 2: Extract cis-pQTLs from ProteoNexus UKB-female pre-filtered TSVs.

Reads from processed_data/UKB_female/cis_raw/{GENE}.tsv (written by
protonexus_unpack.py) via plain pd.read_csv() — no tar, no gzip.
Applies the standard filter pipeline (p < 5e-8, MAF, MHC, palindromes)
and writes to processed_data/UKB_female/cis_sumstats/{GENE}.tsv.

Usage:
  python scripts/02_cis_pqtl_extract/ukb_female.py [--workers N] [--limit N]
"""
import argparse
import logging
from pathlib import Path
from typing import Callable

import pandas as pd

from scripts.lib.config import add_config_arg, load_config, get_section
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import UKB_FEMALE_CIS_RAW, cohort_dir
from scripts.lib.schema import ProteinMeta
from scripts.lib.cis_extract import run_extraction

log = setup_logger("02e_ukb_female")

COHORT = "UKB_female"
BUILD = "hg19"


def _load_tss_cache(cache_path: Path) -> dict[str, tuple[str, int]]:
    """Load {gene_upper: (chrom, tss)} from the TSS cache written by Phase 1."""
    if not cache_path.exists():
        return {}
    try:
        df = pd.read_csv(cache_path, sep="\t", dtype=str)
        result: dict[str, tuple[str, int]] = {}
        for _, row in df.iterrows():
            try:
                result[row["gene"].upper()] = (row["chrom"], int(row["tss"]))
            except (ValueError, KeyError):
                pass
        return result
    except Exception as exc:
        log.warning(f"TSS cache read error: {exc}")
        return {}


def build_protein_list() -> list[ProteinMeta]:
    """
    Scan cis_raw/ for existing TSVs written by protonexus_unpack.py.
    Returns a ProteinMeta for each TSV found.
    """
    tss_cache_path = cohort_dir(COHORT) / "_tss_hg19.tsv"
    tss_cache = _load_tss_cache(tss_cache_path)

    proteins: list[ProteinMeta] = []
    missing_tss: list[str] = []

    for tsv_path in sorted(UKB_FEMALE_CIS_RAW.glob("*.tsv")):
        gene = tsv_path.stem  # e.g. "GMPR2"
        if gene not in tss_cache:
            log.debug(f"No TSS cached for {gene} — skipping")
            missing_tss.append(gene)
            continue

        chrom, tss = tss_cache[gene]
        proteins.append(ProteinMeta(
            seqid=gene,
            gene=gene,
            uniprot="",
            chrom=str(chrom),
            tss=tss,
            build=BUILD,
            source_cohort=COHORT,
        ))

    log.info(
        f"UKB_female: {len(proteins)} proteins from cis_raw "
        f"({len(missing_tss)} skipped — no TSS)"
    )
    return proteins


def normalize_protonexus_rows(rows: list[dict]) -> pd.DataFrame | None:
    """
    Map GEMMA column names to the shared pipeline schema.

    GEMMA → pipeline:
      chr      → chrom  (strip "chr" prefix; cast to str)
      ps       → pos    (int)
      rs       → rsid
      allele1  → EA     (effect allele)
      allele0  → OA     (other allele)
      af       → EAF
      beta     → beta
      se       → se
      p_wald   → pval
      n_obs    → N      (per-SNP integer)
    """
    if not rows:
        return None
    df = pd.DataFrame(rows)
    if df.empty:
        return None

    df = df.rename(columns={
        "rs":      "rsid",
        "allele1": "EA",
        "allele0": "OA",
        "af":      "EAF",
        "beta":    "beta",
        "se":      "se",
        "p_wald":  "pval",
    })

    # chrom: strip "chr" prefix and cast to str
    df["chrom"] = df["chr"].astype(str).str.replace(r"^chr", "", regex=True)

    # pos: integer
    df["pos"] = pd.to_numeric(df["ps"], errors="coerce").astype("int64")

    # N: per-SNP from n_obs (always present in GEMMA output)
    df["N"] = pd.to_numeric(df["n_obs"], errors="coerce").astype("Int64").astype(int)

    # Uppercase alleles
    df["EA"] = df["EA"].astype(str).str.upper()
    df["OA"] = df["OA"].astype(str).str.upper()

    # Drop unused GEMMA columns
    for col in ("chr", "ps", "n_mis", "pip_susie", "fwer", "n_obs"):
        if col in df.columns:
            df = df.drop(columns=[col])

    return df


def build_read_fn() -> Callable[[ProteinMeta], pd.DataFrame | None]:
    """
    Return a read function that reads the pre-filtered cis TSV and
    normalizes it.  No tar, no gzip — plain pd.read_csv().
    """
    def read_fn(protein: ProteinMeta) -> pd.DataFrame | None:
        path = UKB_FEMALE_CIS_RAW / f"{protein.seqid}.tsv"
        if not path.exists():
            log.debug(f"{protein.seqid}: cis_raw TSV not found — {path}")
            return None
        try:
            df = pd.read_csv(path, sep="\t", dtype={"chr": str, "rs": str})
        except Exception as exc:
            log.warning(f"{protein.seqid}: failed to read {path}: {exc}")
            return None
        return normalize_protonexus_rows(df.to_dict("records"))

    return read_fn


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract UKB-female cis-pQTLs from ProteoNexus cis_raw TSVs"
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int, default=None)
    add_config_arg(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    cis_cfg = get_section(cfg, "cis_extract")

    with RunManifest("02_cis_pqtl_extract/ukb_female.py") as manifest:
        proteins = build_protein_list()
        read_fn = build_read_fn()
        n = run_extraction(
            COHORT,
            proteins[:args.limit] if args.limit else proteins,
            read_fn,
            workers=args.workers,
            cfg=cis_cfg,
        )
        manifest.n_units = n


if __name__ == "__main__":
    main()
