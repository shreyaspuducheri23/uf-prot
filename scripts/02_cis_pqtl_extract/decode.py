#!/usr/bin/env python3
"""
02_cis_pqtl_extract/decode.py
Extract cis-pQTLs from deCODE per-aptamer .txt.gz files via signed HTTP URLs.

Downloads sequentially (deCODE throttles parallel connections).
Caches EAF from assocvariants.annotated.txt.gz on first run.

Usage:
  python scripts/02_cis_pqtl_extract/decode.py [--limit N]
"""
import argparse
import gzip
import io
import logging

import pandas as pd

from scripts.lib.cis import tss_from_ensembl
from scripts.lib.config import add_config_arg, load_config, get_section
from scripts.lib.decode_stream import iter_decode_rows, parse_bulk_urls
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import DECODE_ANNOTATED, DECODE_URLS, cohort_dir
from scripts.lib.progress import bar
from scripts.lib.schema import ProteinMeta
from scripts.lib.cis_extract import run_extraction

log = setup_logger("02_decode")

COHORT = "deCODE"
BUILD = "hg38"


def load_eaf_dict() -> dict[str, float]:
    """Load EAF for all deCODE variants from the annotation file."""
    eaf_cache = cohort_dir(COHORT) / "_eaf_cache.pkl"
    if eaf_cache.exists():
        import pickle
        with open(eaf_cache, "rb") as fh:
            return pickle.load(fh)

    log.info(f"Loading EAF from {DECODE_ANNOTATED}...")
    eaf: dict[str, float] = {}
    with gzip.open(DECODE_ANNOTATED, "rt") as fh:
        header = fh.readline().strip().split("\t")
        name_i = header.index("Name")
        eaf_i = header.index("effectAlleleFreq")
        for line in bar(fh, desc="EAF load", total=None):
            parts = line.strip().split("\t")
            if len(parts) > eaf_i:
                try:
                    eaf[parts[name_i]] = float(parts[eaf_i])
                except ValueError:
                    pass

    import pickle
    eaf_cache.parent.mkdir(parents=True, exist_ok=True)
    with open(eaf_cache, "wb") as fh:
        pickle.dump(eaf, fh)
    log.info(f"Loaded EAF for {len(eaf):,} variants → cached at {eaf_cache}")
    return eaf


def build_protein_list(urls: list[tuple[str, str]]) -> list[ProteinMeta]:
    """
    Convert (protein_name, url) list to ProteinMeta objects.
    protein_name format: '<id>_<sub>_<gene>_<protein>' (e.g. '10000_28_CRYBB2_CRBB2').
    TSS fetched from Ensembl hg38 REST (cached via @lru_cache).
    """
    tss_cache_path = cohort_dir(COHORT) / "_tss_hg38.tsv"
    # Load existing cache
    tss_cache: dict[str, tuple[str, int]] = {}
    if tss_cache_path.exists():
        df = pd.read_csv(tss_cache_path, sep="\t", dtype=str)
        for _, row in df.iterrows():
            try:
                tss_cache[row["gene"]] = (row["chrom"], int(row["tss"]))
            except (ValueError, KeyError):
                pass

    proteins = []
    new_cache_rows = []

    for protein_name, url in bar(urls, desc="Build deCODE protein list"):
        parts = protein_name.split("_")
        if len(parts) < 3:
            continue
        # Gene is the 3rd component (index 2)
        seqid_base = "_".join(parts[:2])
        gene = parts[2]
        seqid = protein_name

        if gene not in tss_cache:
            result = tss_from_ensembl(gene, BUILD)
            if result:
                tss_cache[gene] = result
                new_cache_rows.append({"gene": gene, "chrom": result[0], "tss": result[1]})
            else:
                log.debug(f"TSS not found for {gene}")
                continue

        chrom, tss = tss_cache[gene]
        proteins.append(ProteinMeta(
            seqid=seqid, gene=gene, uniprot="",
            chrom=str(chrom), tss=tss, build=BUILD, source_cohort=COHORT,
        ))

    # Persist any new TSS lookups
    if new_cache_rows:
        existing = pd.read_csv(tss_cache_path, sep="\t") if tss_cache_path.exists() else pd.DataFrame()
        updated = pd.concat([existing, pd.DataFrame(new_cache_rows)], ignore_index=True)
        updated.drop_duplicates("gene").to_csv(tss_cache_path, sep="\t", index=False)

    return proteins


# Map protein_name → URL for use in read_fn
_url_map: dict[str, str] = {}


def read_decode_protein(protein: ProteinMeta, eaf_dict: dict[str, float]) -> pd.DataFrame | None:
    url = _url_map.get(protein.seqid)
    if not url:
        return None

    rows = list(iter_decode_rows(url))
    if not rows:
        return None

    df = pd.DataFrame(rows)
    # deCODE columns: Chrom, Pos(hg38), Name, rsids, effectAllele, otherAllele, Beta, Pval, SE, N, ImpMAF
    df = df.rename(columns={
        "Chrom": "chrom",
        "Pos(hg38)": "pos",
        "rsids": "rsid",
        "effectAllele": "EA",
        "otherAllele": "OA",
        "Beta": "beta",
        "Pval": "pval",
        "SE": "se",
        "N": "N",
        "Name": "variant_name",
    })
    df["chrom"] = df["chrom"].astype(str).str.lstrip("chr")
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    df["se"] = pd.to_numeric(df["se"], errors="coerce")
    df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
    if "N" in df.columns:
        df["N"] = pd.to_numeric(df["N"], errors="coerce").fillna(35_000).astype(int)
    else:
        df["N"] = 35_000

    n_in = len(df)
    df["EAF"] = df["variant_name"].map(eaf_dict)
    df = df.dropna(subset=["EAF", "pos", "pval", "beta", "se"])
    n_missing_eaf = n_in - len(df)
    if n_missing_eaf:
        pct = 100 * n_missing_eaf / n_in if n_in else 0
        log.info(
            f"{protein.seqid}: EAF lookup — {n_missing_eaf}/{n_in} ({pct:.1f}%) "
            f"variants dropped (not in EAF cache)"
        )

    return df


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract deCODE cis-pQTLs")
    parser.add_argument("--limit", type=int, default=None)
    add_config_arg(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    cis_cfg = get_section(cfg, "cis_extract")

    with RunManifest("02_cis_pqtl_extract/decode.py") as manifest:
        global _url_map
        urls = parse_bulk_urls(DECODE_URLS)
        _url_map = {name: url for name, url in urls}

        eaf_dict = load_eaf_dict()
        proteins = build_protein_list(urls)

        log.info(f"deCODE: {len(proteins)} proteins")

        def read_fn(protein: ProteinMeta) -> pd.DataFrame | None:
            return read_decode_protein(protein, eaf_dict)

        n = run_extraction(COHORT, proteins, read_fn, limit=args.limit, cfg=cis_cfg)
        manifest.n_units = n


if __name__ == "__main__":
    main()
