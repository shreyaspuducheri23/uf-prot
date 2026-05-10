#!/usr/bin/env python3
"""
02_cis_pqtl_extract/fenland.py
Extract cis-pQTLs from Fenland via Synapse (syn51824537).

Each aptamer has 2 files (.txt.gz). Downloads sequentially (Synapse throttles).

Usage:
  python scripts/02_cis_pqtl_extract/fenland.py [--limit N]
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
from scripts.lib.progress import bar
from scripts.lib.schema import ProteinMeta
from scripts.lib.synapse_stream import list_folder, stream_fenland_protein
from scripts.lib.cis_extract import run_extraction

log = setup_logger("02_fenland")

COHORT = "Fenland"
BUILD = "hg19"
SYNAPSE_FOLDER = "syn51824537"
FENLAND_N = 10_708


def load_fenland_manifest() -> dict[str, list[tuple[str, str]]]:
    """
    Return {protein_name: [(entity_id, filename), ...]} — typically 2 files per protein.
    """
    entities = list_folder(SYNAPSE_FOLDER)
    protein_files: dict[str, list[tuple[str, str]]] = {}
    for e in entities:
        name = e["name"]  # e.g. 'PROTEIN.txt.gz' or 'PROTEIN_2.txt.gz'
        base = name.split("_")[0].removesuffix(".txt.gz")
        protein_files.setdefault(base, []).append((e["id"], name))
    log.info(f"Fenland: {len(protein_files)} proteins ({len(entities)} files)")
    return protein_files


def build_protein_list(protein_files: dict) -> tuple[list[ProteinMeta], dict]:
    tss_cache_path = cohort_dir(COHORT) / "_tss_hg19.tsv"
    tss_cache: dict[str, tuple[str, int]] = {}
    if tss_cache_path.exists():
        df = pd.read_csv(tss_cache_path, sep="\t", dtype=str)
        for _, row in df.iterrows():
            try:
                tss_cache[row["gene"]] = (row["chrom"], int(row["tss"]))
            except (ValueError, KeyError):
                pass

    proteins = []
    entity_map: dict[str, list] = {}
    new_rows = []

    for gene, files in bar(protein_files.items(), desc="Fenland TSS lookup"):
        if gene not in tss_cache:
            result = tss_from_ensembl(gene, BUILD)
            if result:
                tss_cache[gene] = result
                new_rows.append({"gene": gene, "chrom": result[0], "tss": result[1]})
            else:
                log.debug(f"TSS not found for Fenland gene {gene}")
                continue

        chrom, tss = tss_cache[gene]
        protein = ProteinMeta(
            seqid=gene, gene=gene, uniprot="",
            chrom=str(chrom), tss=tss, build=BUILD, source_cohort=COHORT,
        )
        proteins.append(protein)
        entity_map[gene] = files

    if new_rows:
        existing = pd.read_csv(tss_cache_path, sep="\t") if tss_cache_path.exists() else pd.DataFrame()
        updated = pd.concat([existing, pd.DataFrame(new_rows)], ignore_index=True)
        updated.drop_duplicates("gene").to_csv(tss_cache_path, sep="\t", index=False)

    return proteins, entity_map


def read_fenland_protein(protein: ProteinMeta, entity_map: dict,
                          cis_start: int, cis_end: int) -> pd.DataFrame | None:
    files = entity_map.get(protein.seqid, [])
    all_rows = []
    for entity_id, _fname in files:
        with tempfile.TemporaryDirectory(prefix=f"fenland_{protein.seqid}_") as tmp:
            rows = stream_fenland_protein(entity_id, protein.chrom, cis_start, cis_end, Path(tmp))
        all_rows.extend(rows)

    if not all_rows:
        return None

    df = pd.DataFrame(all_rows)
    # Fenland columns probed at runtime — expected: CHR, POS, SNPID/rsid, EA, OA, EAF, BETA, SE, P, N
    rename = {}
    col_map = {"CHR": "chrom", "POS": "pos", "EA": "EA", "OA": "OA",
               "EAF": "EAF", "BETA": "beta", "SE": "se", "P": "pval", "N": "N"}
    for src, dst in col_map.items():
        if src in df.columns:
            rename[src] = dst
    # Try alternative names
    for alt, dst in [("SNPID", "rsid"), ("rsid", "rsid"), ("SNP", "rsid")]:
        if alt in df.columns and "rsid" not in rename.values():
            rename[alt] = "rsid"
    df = df.rename(columns=rename)
    if "rsid" not in df.columns:
        df["rsid"] = "."

    df["chrom"] = df["chrom"].astype(str).str.lstrip("chr")
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    df["se"] = pd.to_numeric(df["se"], errors="coerce")
    df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
    df["EAF"] = pd.to_numeric(df["EAF"], errors="coerce")
    if "N" in df.columns:
        df["N"] = pd.to_numeric(df["N"], errors="coerce").fillna(FENLAND_N).astype(int)
    else:
        df["N"] = FENLAND_N

    return df


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Fenland cis-pQTLs")
    parser.add_argument("--limit", type=int, default=None)
    add_config_arg(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    cis_cfg = get_section(cfg, "cis_extract")
    window_kb = cis_cfg["window_kb"]

    with RunManifest("02_cis_pqtl_extract/fenland.py") as manifest:
        protein_files = load_fenland_manifest()
        proteins, entity_map = build_protein_list(protein_files)
        log.info(f"Fenland: {len(proteins)} proteins with TSS")

        from scripts.lib.cis import cis_window_bounds

        def read_fn(protein: ProteinMeta) -> pd.DataFrame | None:
            start, end = cis_window_bounds(protein.tss, kb=window_kb)
            return read_fenland_protein(protein, entity_map, start, end)

        n = run_extraction(COHORT, proteins, read_fn, limit=args.limit, cfg=cis_cfg)
        manifest.n_units = n


if __name__ == "__main__":
    main()
