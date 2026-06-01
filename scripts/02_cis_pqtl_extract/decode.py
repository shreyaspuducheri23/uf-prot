#!/usr/bin/env python3
"""
02_cis_pqtl_extract/decode.py
Extract cis-pQTLs from deCODE per-aptamer .txt.gz files via S3 streaming.

Files are streamed from s3-ext.decode.is:10443 with early abort once the
cis window on the target chromosome is passed — avoiding ~99% of bytes.

EAF comes directly from the ImpMAF column in each per-protein file.
The 2023 deCODE S3 files consistently encode effectAllele as the minor allele,
so ImpMAF == EAF exactly (confirmed empirically: ImpMAF never exceeds 0.5).

Supports threaded extraction workers; tune with --workers.

Usage:
  python scripts/02_cis_pqtl_extract/decode.py [--normalization raw|smp] [--limit N] [--workers N]
"""
import argparse
import json
import re
import logging
import time
from collections.abc import Iterable

import pandas as pd

from scripts.lib.cis import _append_unresolved, _load_tss_cache, _save_tss_cache, resolve_tss
from scripts.lib.somascan import load_seqid_map
from scripts.lib.config import (
    add_config_arg, load_config, get_section, get_cohort_build, get_cohort_sample_size
)
from scripts.lib.decode_stream import (
    _get_s3_client, stream_s3_cis_rows,
    DECODE_S3_ENDPOINT, DECODE_S3_BUCKET, DECODE_S3_ACCESS_KEY, DECODE_S3_SECRET_KEY,
)
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import cohort_dir
from scripts.lib.progress import bar
from scripts.lib.schema import ProteinMeta
from scripts.lib.cis_extract import RAW_CIS_WINDOW_KB, run_extraction

log = setup_logger("02_decode")

BUILD = "hg38"
_DEFAULT_N = 35_559
_REQUIRED_COLS = {
    "Chrom", "Pos", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval", "SE", "N",
    "ImpMAF",
}
_FAST_PARSER_COLS = frozenset(_REQUIRED_COLS)
_prefilter_window_bp = RAW_CIS_WINDOW_KB * 1_000

_S3_ENDPOINT   = DECODE_S3_ENDPOINT
_S3_BUCKET     = DECODE_S3_BUCKET
_S3_ACCESS_KEY = DECODE_S3_ACCESS_KEY
_S3_SECRET_KEY = DECODE_S3_SECRET_KEY

_NORM_CONFIG = {
    "raw": {
        "prefix":  "final_somascan_raw",
        "pattern": r"Proteomics_PC0_(.+)_\d{8}\.txt\.gz",
        "cohort":  "deCODE",
    },
    "smp": {
        "prefix":  "final_somascan_smp",
        "pattern": r"Proteomics_SMP_PC0_(.+)_\d{8}\.txt\.gz",
        "cohort":  "deCODE_smp",
    },
}

# Set in main() based on --normalization
COHORT = "deCODE"

# Populated at startup by _build_s3_key_index; maps core_name → full S3 key
_s3_key_map: dict[str, str] = {}



def _build_s3_key_index(prefix: str, pattern: str) -> dict[str, str]:
    """Build a {core_name: full_s3_key} index by listing the S3 prefix.

    Result is cached to disk so subsequent runs skip the listing.
    """
    cache_path = cohort_dir(COHORT) / "_s3_key_index.json"
    if cache_path.exists():
        with open(cache_path) as fh:
            return json.load(fh)

    log.info(f"Building S3 key index for s3://{_S3_BUCKET}/{prefix}/ ...")
    s3 = _get_s3_client(_S3_ENDPOINT, _S3_ACCESS_KEY, _S3_SECRET_KEY)
    regex = re.compile(pattern)
    index: dict[str, str] = {}

    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=_S3_BUCKET, Prefix=f"{prefix}/"):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            fname = key.split("/")[-1]
            m = regex.match(fname)
            if m:
                core = m.group(1)
                index[core] = key

    log.info(f"S3 key index: {len(index):,} proteins")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as fh:
        json.dump(index, fh)
    return index


def build_protein_list(protein_names: Iterable[str], build: str = BUILD) -> list[ProteinMeta]:
    """
    Convert deCODE protein names to ProteinMeta objects.
    protein_name format: '<id>_<sub>_<gene>_<protein>' (e.g. '10000_28_CRYBB2_CRBB2').
    TSS fetched from Ensembl for the requested build (cached via @lru_cache).
    """
    cohort_path = cohort_dir(COHORT)
    tss_cache_path = cohort_path / "_tss_hg38.tsv"
    # All deCODE cohorts use identical genes; share the raw-cohort TSS cache when
    # the cohort-specific file doesn't exist yet (avoids ~20 min Ensembl re-lookup).
    _tss_seed = cohort_dir("deCODE") / "_tss_hg38.tsv"
    if not tss_cache_path.exists() and _tss_seed.exists() and _tss_seed != tss_cache_path:
        import shutil
        tss_cache_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(_tss_seed, tss_cache_path)
        log.info(f"Seeded TSS cache from deCODE → {tss_cache_path}")
    tss_cache = _load_tss_cache(tss_cache_path)

    seqid_map = load_seqid_map()
    proteins = []
    new_cache_rows = []
    unresolved_rows = []

    for protein_name in bar(protein_names, desc="Build deCODE protein list"):
        parts = protein_name.split("_")
        if len(parts) < 3:
            continue
        seqid_key = f"{parts[0]}_{parts[1]}"
        gene = seqid_map.get(seqid_key, parts[2])
        seqid = protein_name

        if gene not in tss_cache:
            r = resolve_tss(gene, build)
            if r.resolved:
                tss_cache[gene] = (r.chrom, r.tss)
                new_cache_rows.append({
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
                    log.debug(f"TSS not found for {gene}")
                    unresolved_rows.append({
                        "gene": gene,
                        "build": r.build,
                        "attempts": "|".join(r.attempts),
                    })
                continue

        chrom, tss = tss_cache[gene]
        proteins.append(ProteinMeta(
            seqid=seqid, gene=gene, uniprot="",
            chrom=str(chrom), tss=tss, build=build, source_cohort=COHORT,
        ))

    if new_cache_rows:
        _save_tss_cache(tss_cache_path, tss_cache, new_cache_rows)
    _append_unresolved(cohort_path, unresolved_rows)

    return proteins


_S3_MAX_RETRIES = 5
_S3_RETRY_BASE  = 10.0  # seconds; doubles each attempt (10, 20, 40, 80, 160)


def _load_decode_raw_df(protein: ProteinMeta) -> pd.DataFrame | None:
    """Stream cis-window rows for one protein from S3, with retries on transient errors."""
    key = _s3_key_map.get(protein.seqid)
    if not key:
        log.warning(f"{protein.seqid}: no S3 key found in index — skipping")
        return None

    s3 = _get_s3_client(_S3_ENDPOINT, _S3_ACCESS_KEY, _S3_SECRET_KEY)
    last_exc: Exception | None = None
    for attempt in range(1, _S3_MAX_RETRIES + 1):
        try:
            rows = list(stream_s3_cis_rows(
                s3, _S3_BUCKET, key,
                protein.chrom, protein.tss,
                _prefilter_window_bp, _FAST_PARSER_COLS,
            ))
            return pd.DataFrame(rows) if rows else None
        except Exception as exc:
            last_exc = exc
            if attempt < _S3_MAX_RETRIES:
                wait = _S3_RETRY_BASE * (2 ** (attempt - 1))
                log.warning(
                    f"{protein.seqid}: S3 attempt {attempt}/{_S3_MAX_RETRIES} failed "
                    f"({type(exc).__name__}), retrying in {wait:.0f}s"
                )
                time.sleep(wait)
    raise last_exc


def read_decode_protein(protein: ProteinMeta, n_default: int | None = _DEFAULT_N) -> pd.DataFrame | None:
    df = _load_decode_raw_df(protein)
    if df is None:
        return None

    df = df.rename(columns={
        "Chrom": "chrom",
        "Pos": "pos",
        "rsids": "rsid",
        "effectAllele": "EA",
        "otherAllele": "OA",
        "Beta": "beta",
        "Pval": "pval",
        "SE": "se",
        "N": "N",
        "Name": "variant_name",
        "ImpMAF": "EAF",
    })
    df["chrom"] = df["chrom"].astype(str).str.lstrip("chr")
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    df["se"] = pd.to_numeric(df["se"], errors="coerce")
    df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
    df["EAF"] = pd.to_numeric(df["EAF"], errors="coerce")

    if "N" in df.columns:
        n = pd.to_numeric(df["N"], errors="coerce")
        if n_default is None and n.isna().any():
            raise ValueError(f"{protein.seqid}: missing N and no configured sample_size")
        df["N"] = n.astype(int) if n_default is None else n.fillna(n_default).astype(int)
    else:
        if n_default is None:
            raise ValueError(f"{protein.seqid}: missing N column and no configured sample_size")
        df["N"] = n_default

    df = df.dropna(subset=["EAF", "pos", "pval", "beta", "se"])
    if df.empty:
        return None
    df["pos"] = df["pos"].astype(int)
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract deCODE cis-pQTLs")
    parser.add_argument("--normalization", choices=["raw", "smp"], default="raw",
                        help="Which deCODE normalization to use (default: raw)")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--workers", type=int, default=4,
                        help="Thread workers for per-protein downloads (default: 4)")
    add_config_arg(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    cis_cfg = get_section(cfg, "cis_extract")
    workers = max(1, args.workers)

    norm_cfg = _NORM_CONFIG[args.normalization]

    with RunManifest("02_cis_pqtl_extract/decode.py") as manifest:
        global COHORT, _s3_key_map, _prefilter_window_bp
        COHORT = norm_cfg["cohort"]
        build = get_cohort_build(cfg, COHORT) if COHORT in cfg["cohorts"] else BUILD
        n_default = get_cohort_sample_size(cfg, COHORT) if COHORT in cfg["cohorts"] else _DEFAULT_N
        _prefilter_window_bp = RAW_CIS_WINDOW_KB * 1_000

        _s3_key_map = _build_s3_key_index(norm_cfg["prefix"], norm_cfg["pattern"])

        proteins = build_protein_list(_s3_key_map.keys(), build=build)

        log.info(f"{COHORT}: {len(proteins)} proteins")

        def read_fn(protein: ProteinMeta) -> pd.DataFrame | None:
            return read_decode_protein(protein, n_default=n_default)

        n = run_extraction(
            COHORT,
            proteins,
            read_fn,
            workers=workers,
            limit=args.limit,
            cfg=cis_cfg,
        )
        manifest.n_units = n


if __name__ == "__main__":
    main()
