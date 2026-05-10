#!/usr/bin/env python3
"""
02_cis_pqtl_extract/ukbppp.py
Extract cis-pQTLs from UKB-PPP via Synapse streaming (syn51365303).

Each protein is stored as a ~550 MB tar. Downloads one tar at a time,
extracts per-chromosome .gz files in memory, filters to cis positions, saves, deletes.

Usage:
  python scripts/02_cis_pqtl_extract/ukbppp.py [--workers N] [--limit N]
"""
import argparse
import logging
import tempfile
from pathlib import Path

import pandas as pd

from scripts.lib.cis import tss_from_ensembl
from scripts.lib.config import add_config_arg, load_config, get_section
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import cohort_dir
from scripts.lib.progress import bar, counter
from scripts.lib.schema import ProteinMeta
from scripts.lib.synapse_stream import list_folder, stream_ukbppp_protein
from scripts.lib.cis_extract import run_extraction

log = setup_logger("02_ukbppp")

COHORT = "UKB_PPP"
BUILD = "hg19"
SYNAPSE_FOLDER = "syn51365303"
UKB_N = 34_557


def load_ukbppp_manifest() -> list[tuple[str, str, str]]:
    """
    List Synapse folder, infer gene from entity name.
    Returns [(entity_id, protein_name, gene_symbol), ...].
    Entity names follow UKB-PPP conventions (e.g. 'PROC_P04070.tar').
    """
    entities = list_folder(SYNAPSE_FOLDER)
    manifest = []
    for e in entities:
        name = e["name"]  # e.g. 'PROC_P04070.tar' or '<GENE>_<UNIPROT>.tar'
        stem = name.removesuffix(".tar")
        parts = stem.split("_")
        gene = parts[0] if parts else stem
        manifest.append((e["id"], stem, gene))
    log.info(f"UKB-PPP Synapse folder: {len(manifest)} entities")
    return manifest


def build_protein_list(manifest: list[tuple[str, str, str]]) -> tuple[list[ProteinMeta], dict[str, str]]:
    """Returns (proteins, entity_id_map {seqid: entity_id})."""
    entity_map: dict[str, str] = {}
    proteins = []

    tss_cache_path = cohort_dir(COHORT) / "_tss_hg19.tsv"
    tss_cache: dict[str, tuple[str, int]] = {}
    if tss_cache_path.exists():
        df = pd.read_csv(tss_cache_path, sep="\t", dtype=str)
        for _, row in df.iterrows():
            try:
                tss_cache[row["gene"]] = (row["chrom"], int(row["tss"]))
            except (ValueError, KeyError):
                pass

    new_rows = []
    for entity_id, protein_name, gene in bar(manifest, desc="UKB-PPP TSS lookup"):
        if gene not in tss_cache:
            result = tss_from_ensembl(gene, BUILD)
            if result:
                tss_cache[gene] = result
                new_rows.append({"gene": gene, "chrom": result[0], "tss": result[1]})
            else:
                log.debug(f"TSS not found for UKB-PPP gene {gene}")
                continue

        chrom, tss = tss_cache[gene]
        protein = ProteinMeta(
            seqid=protein_name, gene=gene, uniprot=protein_name.split("_")[-1] if "_" in protein_name else "",
            chrom=str(chrom), tss=tss, build=BUILD, source_cohort=COHORT,
        )
        proteins.append(protein)
        entity_map[protein_name] = entity_id

    if new_rows:
        existing = pd.read_csv(tss_cache_path, sep="\t") if tss_cache_path.exists() else pd.DataFrame()
        updated = pd.concat([existing, pd.DataFrame(new_rows)], ignore_index=True)
        updated.drop_duplicates("gene").to_csv(tss_cache_path, sep="\t", index=False)

    return proteins, entity_map


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract UKB-PPP cis-pQTLs")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int, default=None)
    add_config_arg(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    cis_cfg = get_section(cfg, "cis_extract")
    window_kb = cis_cfg["window_kb"]

    with RunManifest("02_cis_pqtl_extract/ukbppp.py") as manifest:
        manifest_list = load_ukbppp_manifest()
        proteins, entity_map = build_protein_list(manifest_list)

        # cis positions keyed by (chrom, pos) in hg19 for each protein
        # We pass the protein object and let stream_ukbppp_protein do position-based filtering
        # via the ±500kb window; the stream function accepts a position set.
        # For UKB-PPP we pre-build per-protein cis windows here.

        from scripts.lib.cis import cis_window_bounds
        from scripts.lib.filters import exclude_mhc

        def read_fn(protein: ProteinMeta) -> pd.DataFrame | None:
            start, end = cis_window_bounds(protein.tss, kb=window_kb)
            eid = entity_map.get(protein.seqid)
            if not eid:
                return None
            with tempfile.TemporaryDirectory(prefix=f"ukb_{protein.seqid}_") as tmp:
                rows = stream_ukbppp_protein(eid, protein.chrom, start, end, Path(tmp))
            if not rows:
                return None
            df = pd.DataFrame(rows)
            df = df.rename(columns={
                "ALLELE1": "EA", "ALLELE0": "OA",
                "A1FREQ": "EAF", "BETA": "beta", "SE": "se",
            })
            id_parts = df["ID"].str.split(":", expand=True)
            df["chrom"] = id_parts[0].str.replace(r"^chr", "", regex=True)
            df["pos"] = pd.to_numeric(id_parts[1], errors="coerce").astype("Int64")
            df["pval"] = 10 ** (-pd.to_numeric(df["LOG10P"], errors="coerce"))
            df["rsid"] = "."
            if "N" in df.columns:
                df["N"] = pd.to_numeric(df["N"], errors="coerce").fillna(UKB_N).astype(int)
            else:
                df["N"] = UKB_N
            return df

        n = run_extraction(COHORT, proteins[:args.limit] if args.limit else proteins,
                           read_fn, workers=args.workers, cfg=cis_cfg)
        manifest.n_units = n


if __name__ == "__main__":
    main()
