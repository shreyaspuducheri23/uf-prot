#!/usr/bin/env python3
"""
02_cis_pqtl_extract/protonexus_unpack.py
Phase 1: Unpack ProteoNexus tar archives → per-gene raw cis-window TSVs.

Iterates the 26 alphabetical tars from /Volumes/Extreme SSD/ProteoNexus/
sequentially (no SSD contention), extracts each gene's
summ_female2.assoc.txt.gz, filters to the ±1 Mb cis window in-memory,
and writes a plain TSV to processed_data/UKB_female/cis_raw_1000kb/{GENE}.tsv.

Per-gene checkpointing makes it fully resumable.

Usage:
  python scripts/02_cis_pqtl_extract/protonexus_unpack.py [--limit N] [--config PATH]
"""
import argparse
import logging
import tarfile
from pathlib import Path

import pandas as pd

from scripts.lib.checkpoint import Checkpoint
from scripts.lib.cis import cis_window_bounds, tss_from_ensembl
from scripts.lib.config import add_config_arg, load_config, get_section, get_cohort_build
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.paths import UKB_FEMALE_DIR, UKB_FEMALE_CIS_RAW, cohort_dir
from scripts.lib.progress import bar
from scripts.lib.cis_extract import RAW_CIS_WINDOW_KB

log = setup_logger("02e_prep_protonexus_unpack")

COHORT = "UKB_female"
BUILD = "hg19"
# Member path pattern inside each tar: <gene>/output/summ_female2.assoc.txt.gz
_MEMBER_SUFFIX = "/output/summ_female2.assoc.txt.gz"


def _load_tss_cache(cache_path: Path) -> dict[str, tuple[str, int]]:
    """Load TSS cache from disk; returns {gene_upper: (chrom, tss)}."""
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
        log.warning(f"TSS cache read error ({cache_path}): {exc}")
        return {}


def _save_tss_cache(cache_path: Path, cache: dict[str, tuple[str, int]]) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    rows = [{"gene": g, "chrom": c, "tss": t} for g, (c, t) in sorted(cache.items())]
    pd.DataFrame(rows).to_csv(cache_path, sep="\t", index=False)


def run_unpack(limit: int | None, window_kb: int, build: str = BUILD) -> int:
    """
    Main unpack loop.  Returns count of genes successfully written.
    """
    UKB_FEMALE_CIS_RAW.mkdir(parents=True, exist_ok=True)
    cohort_dir(COHORT).mkdir(parents=True, exist_ok=True)

    cp = Checkpoint(cohort_dir(COHORT) / f"_state_02_unpack_{window_kb}kb.json")
    tss_cache_path = cohort_dir(COHORT) / "_tss_hg19.tsv"
    tss_cache = _load_tss_cache(tss_cache_path)

    tars = sorted(UKB_FEMALE_DIR.glob("ProteoNexus_pQTL_protein_*.tar"))
    if not tars:
        log.warning(f"No ProteoNexus tar files found in {UKB_FEMALE_DIR}")
        return 0

    log.info(f"Found {len(tars)} ProteoNexus tar archives in {UKB_FEMALE_DIR}")

    n_ok = 0
    n_skip_tss = 0
    n_done_already = 0
    tss_cache_dirty = False
    genes_processed: list[str] = []

    for tar_path in bar(tars, desc="ProteoNexus tars"):
        if limit is not None and n_ok >= limit:
            break

        log.debug(f"Opening {tar_path.name}")
        try:
            with tarfile.open(tar_path, "r") as tf:
                members = tf.getmembers()
                for member in members:
                    if limit is not None and n_ok >= limit:
                        break

                    if not member.name.endswith(_MEMBER_SUFFIX):
                        continue

                    # gene directory is the first path component, always lowercase in ProteoNexus
                    gene_lower = member.name.split("/")[0]
                    gene = gene_lower.upper()

                    if cp.is_done(gene):
                        n_done_already += 1
                        continue

                    # TSS lookup (cached)
                    if gene not in tss_cache:
                        result = tss_from_ensembl(gene, build)
                        if result is None:
                            log.warning(f"TSS not found for {gene} — skipping")
                            n_skip_tss += 1
                            cp.mark_done(gene)  # mark so we don't retry on every run
                            continue
                        tss_cache[gene] = result
                        tss_cache_dirty = True

                    chrom, tss = tss_cache[gene]
                    start, end = cis_window_bounds(tss, kb=window_kb)

                    # Extract and filter
                    try:
                        fobj = tf.extractfile(member)
                        if fobj is None:
                            log.warning(f"{gene}: extractfile returned None — skipping")
                            continue

                        df = pd.read_csv(fobj, sep="\t", compression="gzip",
                                         dtype={"chr": str})
                    except Exception as exc:
                        log.warning(f"{gene}: failed to read gz — {exc}")
                        cp.mark_failed(gene, str(exc))
                        continue

                    if df.empty:
                        log.debug(f"{gene}: empty dataframe")
                        cp.mark_done(gene)
                        continue

                    # Normalise chrom column for comparison
                    chrom_col = df["chr"].astype(str).str.replace(r"^chr", "", regex=True)
                    mask = (
                        (chrom_col == str(chrom))
                        & (df["ps"] >= start)
                        & (df["ps"] <= end)
                    )
                    cis_df = df[mask].copy()

                    if cis_df.empty:
                        log.debug(f"{gene}: 0 rows in cis window ({chrom}:{start}-{end})")
                        cp.mark_done(gene)
                        continue

                    out_path = UKB_FEMALE_CIS_RAW / f"{gene}.tsv"
                    cis_df.to_csv(out_path, sep="\t", index=False)
                    cp.mark_done(gene)
                    n_ok += 1
                    genes_processed.append(gene)
                    log.debug(f"{gene}: {len(cis_df)} cis rows → {out_path.name}")

        except tarfile.TarError as exc:
            log.error(f"Failed to open {tar_path.name}: {exc}")
            continue

        # Flush TSS cache after each tar so progress is saved incrementally
        if tss_cache_dirty:
            _save_tss_cache(tss_cache_path, tss_cache)
            tss_cache_dirty = False

    # Final cache flush
    if tss_cache_dirty:
        _save_tss_cache(tss_cache_path, tss_cache)

    log.info(
        f"Unpack complete: {n_ok} genes written | "
        f"{n_done_already} already done | "
        f"{n_skip_tss} skipped (no TSS)"
    )
    return n_ok


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Unpack ProteoNexus tars → per-gene cis TSVs (Phase 1)"
    )
    parser.add_argument("--limit", type=int, default=None,
                        help="Stop after N genes (for testing)")
    add_config_arg(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    cis_cfg = get_section(cfg, "cis_extract")
    build = get_cohort_build(cfg, COHORT)
    window_kb = RAW_CIS_WINDOW_KB

    with RunManifest("02_cis_pqtl_extract/protonexus_unpack.py") as manifest:
        n = run_unpack(limit=args.limit, window_kb=window_kb, build=build)
        manifest.n_units = n


if __name__ == "__main__":
    main()
