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
from typing import Callable

import pandas as pd

from scripts.lib.cis import _append_unresolved, _load_tss_cache, _save_tss_cache, resolve_tss
from scripts.lib.config import (
    add_config_arg, load_config, get_section, get_cohort_build, get_cohort_sample_size
)
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import cohort_dir
from scripts.lib.progress import bar, counter
from scripts.lib.schema import ProteinMeta
from scripts.lib.synapse_stream import list_folder, stream_ukbppp_protein
from scripts.lib.cis_extract import RAW_CIS_WINDOW_KB, run_extraction

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


def build_protein_list(
    manifest: list[tuple[str, str, str]],
    build: str = BUILD,
) -> tuple[list[ProteinMeta], dict[str, str]]:
    """Returns (proteins, entity_id_map {seqid: entity_id})."""
    entity_map: dict[str, str] = {}
    proteins = []

    cohort_path = cohort_dir(COHORT)
    tss_cache_path = cohort_path / "_tss_hg19.tsv"
    tss_cache = _load_tss_cache(tss_cache_path)

    new_rows = []
    unresolved_rows = []
    for entity_id, protein_name, gene in bar(manifest, desc="UKB-PPP TSS lookup"):
        if gene not in tss_cache:
            r = resolve_tss(gene, build)
            if r.resolved:
                tss_cache[gene] = (r.chrom, r.tss)
                new_rows.append({
                    "gene": gene,
                    "chrom": r.chrom,
                    "tss": r.tss,
                    "resolved_symbol": r.resolved_symbol,
                    "tier": r.tier,
                    "source": r.source,
                })
            else:
                if r.transient:
                    log.warning(f"Transient TSS lookup failure for {gene!r}; will retry next run")
                else:
                    log.debug(f"TSS not found for UKB-PPP gene {gene}")
                    unresolved_rows.append({
                        "gene": gene,
                        "build": r.build,
                        "attempts": "|".join(r.attempts),
                    })
                continue

        chrom, tss = tss_cache[gene]
        protein = ProteinMeta(
            seqid=protein_name, gene=gene, uniprot=protein_name.split("_")[-1] if "_" in protein_name else "",
            chrom=str(chrom), tss=tss, build=build, source_cohort=COHORT,
        )
        proteins.append(protein)
        entity_map[protein_name] = entity_id

    if new_rows:
        _save_tss_cache(tss_cache_path, tss_cache, new_rows)
    _append_unresolved(cohort_path, unresolved_rows)

    return proteins, entity_map


def normalize_ukbppp_rows(rows: list[dict], n_default: int | None = UKB_N) -> pd.DataFrame | None:
    """Normalize streamed UKB-PPP rows into the shared extraction schema."""
    if not rows:
        return None
    df = pd.DataFrame(rows)
    if df.empty:
        return None
    df = df.rename(
        columns={
            "ALLELE1": "EA",
            "ALLELE0": "OA",
            "A1FREQ": "EAF",
            "BETA": "beta",
            "SE": "se",
        }
    )
    id_parts = df["ID"].str.split(":", expand=True)
    df["chrom"] = id_parts[0].str.replace(r"^chr", "", regex=True)
    df["pos"] = pd.to_numeric(id_parts[1], errors="coerce").astype("Int64")
    df["pval"] = 10 ** (-pd.to_numeric(df["LOG10P"], errors="coerce"))
    df["rsid"] = "."
    if "N" in df.columns:
        n = pd.to_numeric(df["N"], errors="coerce")
        if n_default is None and n.isna().any():
            raise ValueError("UKB-PPP row missing N and no configured sample_size")
        df["N"] = n.astype(int) if n_default is None else n.fillna(n_default).astype(int)
    else:
        if n_default is None:
            raise ValueError("UKB-PPP rows missing N column and no configured sample_size")
        df["N"] = n_default
    return df


def build_read_fn(
    entity_map: dict[str, str],
    window_kb: int,
    n_default: int | None = UKB_N,
    stream_fn: Callable[[str, str, int, int, Path], list[dict]] = stream_ukbppp_protein,
) -> Callable[[ProteinMeta], pd.DataFrame | None]:
    """Construct the per-protein reader used by the cohort-agnostic extraction loop."""
    from scripts.lib.cis import cis_window_bounds

    def read_fn(protein: ProteinMeta) -> pd.DataFrame | None:
        start, end = cis_window_bounds(protein.tss, kb=window_kb)
        eid = entity_map.get(protein.seqid)
        if not eid:
            return None
        with tempfile.TemporaryDirectory(prefix=f"ukb_{protein.seqid}_") as tmp:
            rows = stream_fn(eid, protein.chrom, start, end, Path(tmp))
        return normalize_ukbppp_rows(rows, n_default=n_default)

    return read_fn


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract UKB-PPP cis-pQTLs")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int, default=None)
    add_config_arg(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    cis_cfg = get_section(cfg, "cis_extract")
    build = get_cohort_build(cfg, COHORT)
    n_default = get_cohort_sample_size(cfg, COHORT)
    window_kb = RAW_CIS_WINDOW_KB

    with RunManifest("02_cis_pqtl_extract/ukbppp.py") as manifest:
        manifest_list = load_ukbppp_manifest()
        proteins, entity_map = build_protein_list(manifest_list, build=build)

        read_fn = build_read_fn(entity_map=entity_map, window_kb=window_kb, n_default=n_default)

        n = run_extraction(COHORT, proteins[:args.limit] if args.limit else proteins,
                           read_fn, workers=args.workers, cfg=cis_cfg)
        manifest.n_units = n


if __name__ == "__main__":
    main()
