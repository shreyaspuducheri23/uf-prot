"""
Cohort-agnostic cis-pQTL extraction pipeline.
Each cohort script provides iter_proteins() and read_protein_sumstats() callbacks;
this module provides the shared filter pipeline, output schema normalization, and loop.
"""
import logging
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Iterator

import pandas as pd

from scripts.lib import filters as F
from scripts.lib.checkpoint import Checkpoint, output_exists
from scripts.lib.cis import cis_window_bounds
from scripts.lib.fstat import add_fstat
from scripts.lib.paths import cis_sumstats_dir, cohort_dir
from scripts.lib.progress import bar
from scripts.lib.schema import ProteinMeta
from scripts.lib.sumstats_io import write_norm

log = logging.getLogger(__name__)
_CP_FLUSH_EVERY = 50

# Normalized output columns for cis_sumstats TSVs
OUTPUT_COLS = [
    "seqid", "gene", "uniprot",
    "chrom", "pos", "rsid", "EA", "OA", "EAF",
    "beta", "se", "pval", "N", "build",
]


def run_extraction(
    cohort: str,
    proteins: list[ProteinMeta],
    read_fn: Callable[[ProteinMeta], pd.DataFrame | None],
    workers: int = 1,
    limit: int | None = None,
    cfg: dict | None = None,
) -> int:
    """
    Main extraction loop shared by all four cohort scripts.

    read_fn(protein) -> DataFrame with raw columns for that protein, or None to skip.
    workers > 1 enables concurrent ThreadPoolExecutor (good for IO-bound downloads).
    Returns count of proteins successfully processed.
    """
    out_dir = cis_sumstats_dir(cohort)
    out_dir.mkdir(parents=True, exist_ok=True)

    state_path = cohort_dir(cohort) / "_state_02.json"
    cp = Checkpoint(state_path)

    todo = cp.remaining(proteins, key=lambda p: p.seqid)
    if limit:
        todo = todo[:limit]

    log.info(f"{cohort}: {len(proteins)} proteins total, {len(todo)} remaining "
             f"({cp.n_done} already done)")

    if workers > 1:
        n_ok = _run_parallel(cohort, todo, read_fn, out_dir, cp, workers, cfg)
    else:
        n_ok = _run_sequential(cohort, todo, read_fn, out_dir, cp, cfg)

    # Keep protein_index.tsv in sync for both sequential and parallel paths.
    _write_index(cohort, proteins)
    return n_ok


def _run_sequential(cohort: str, todo: list[ProteinMeta],
                    read_fn: Callable, out_dir: Path, cp: Checkpoint,
                    cfg: dict | None = None) -> int:
    n_ok = 0
    n_empty = 0
    n_fail = 0
    n_since_flush = 0

    for protein in bar(todo, desc=f"{cohort} extract"):
        out_path = out_dir / f"{protein.seqid}.tsv"
        if output_exists(out_path, required_cols=OUTPUT_COLS, min_rows=1):
            cp.mark_done(protein.seqid, save=False)
            n_since_flush += 1
            n_ok += 1
            if n_since_flush >= _CP_FLUSH_EVERY:
                cp.flush()
                n_since_flush = 0
            continue

        try:
            raw = read_fn(protein)
        except Exception as exc:
            log.warning(f"{protein.seqid}: read failed — {exc}")
            cp.mark_failed(protein.seqid, str(exc), save=False)
            n_fail += 1
            n_since_flush += 1
            if n_since_flush >= _CP_FLUSH_EVERY:
                cp.flush()
                n_since_flush = 0
            continue

        if raw is None or raw.empty:
            n_empty += 1
            cp.mark_done(protein.seqid, save=False)
            n_since_flush += 1
            if n_since_flush >= _CP_FLUSH_EVERY:
                cp.flush()
                n_since_flush = 0
            continue

        filtered = _apply_filters(raw, protein, cfg)
        if filtered.empty:
            log.debug(f"{protein.seqid}: 0 variants after filters")
            n_empty += 1
            cp.mark_done(protein.seqid, save=False)
            n_since_flush += 1
            if n_since_flush >= _CP_FLUSH_EVERY:
                cp.flush()
                n_since_flush = 0
            continue

        normalized = _normalize(filtered, protein)
        write_norm(normalized, out_path)
        cp.mark_done(protein.seqid, save=False)
        n_since_flush += 1
        n_ok += 1
        if n_since_flush >= _CP_FLUSH_EVERY:
            cp.flush()
            n_since_flush = 0

    cp.flush()
    log.info(f"{cohort}: done. {n_ok} with cis-pQTLs | {n_empty} empty | {n_fail} read failures")
    return n_ok


def _run_parallel(cohort: str, todo: list[ProteinMeta],
                  read_fn: Callable, out_dir: Path, cp: Checkpoint,
                  workers: int, cfg: dict | None = None) -> int:
    """Concurrent extraction using threads (safe for IO-bound downloads)."""
    lock = threading.Lock()
    n_ok = 0
    n_empty = 0
    n_fail = 0
    n_since_flush = 0

    def process_one(protein: ProteinMeta) -> tuple[str, bool]:
        nonlocal n_ok, n_empty, n_fail, n_since_flush
        out_path = out_dir / f"{protein.seqid}.tsv"
        if output_exists(out_path, required_cols=OUTPUT_COLS, min_rows=1):
            with lock:
                cp.mark_done(protein.seqid, save=False)
                n_ok += 1
                n_since_flush += 1
                if n_since_flush >= _CP_FLUSH_EVERY:
                    cp.flush()
                    n_since_flush = 0
            return protein.seqid, True

        try:
            raw = read_fn(protein)
        except Exception as exc:
            log.warning(f"{protein.seqid}: read failed — {exc}")
            with lock:
                cp.mark_failed(protein.seqid, str(exc), save=False)
                n_fail += 1
                n_since_flush += 1
                if n_since_flush >= _CP_FLUSH_EVERY:
                    cp.flush()
                    n_since_flush = 0
            return protein.seqid, False

        if raw is None or raw.empty:
            with lock:
                n_empty += 1
                cp.mark_done(protein.seqid, save=False)
                n_since_flush += 1
                if n_since_flush >= _CP_FLUSH_EVERY:
                    cp.flush()
                    n_since_flush = 0
            return protein.seqid, True

        filtered = _apply_filters(raw, protein, cfg)
        if filtered.empty:
            log.debug(f"{protein.seqid}: 0 variants after filters")
            with lock:
                n_empty += 1
                cp.mark_done(protein.seqid, save=False)
                n_since_flush += 1
                if n_since_flush >= _CP_FLUSH_EVERY:
                    cp.flush()
                    n_since_flush = 0
            return protein.seqid, True

        normalized = _normalize(filtered, protein)
        with lock:
            write_norm(normalized, out_path)
            cp.mark_done(protein.seqid, save=False)
            n_ok += 1
            n_since_flush += 1
            if n_since_flush >= _CP_FLUSH_EVERY:
                cp.flush()
                n_since_flush = 0
        return protein.seqid, True

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(process_one, p): p for p in todo}
        for _ in bar(as_completed(futures), total=len(futures), desc=f"{cohort} extract"):
            pass  # progress driven by completion; exceptions raised below

    # Re-raise any worker exceptions
    for future in futures:
        future.result()

    cp.flush()
    log.info(f"{cohort}: done (parallel). {n_ok} with cis-pQTLs | {n_empty} empty | {n_fail} read failures")
    return n_ok


def _apply_filters(df: pd.DataFrame, protein: ProteinMeta,
                   cfg: dict | None = None) -> pd.DataFrame:
    """Apply the 6-stage cis-pQTL filter pipeline.

    cfg: cis_extract section from pipeline.json; falls back to hardcoded defaults.
    """
    window_kb = cfg["window_kb"] if cfg else 500
    pval_gw   = cfg["pval_gw"]   if cfg else 5e-8
    maf_min   = cfg["maf_min"]   if cfg else 0.01
    pal_maf   = cfg["palindrome_maf_max"] if cfg else 0.42

    def _step(label: str, result: pd.DataFrame, n_before: int) -> pd.DataFrame:
        n_after = len(result)
        if n_before != n_after:
            log.debug(f"{protein.seqid} [{label}]: {n_before} → {n_after} variants")
        return result

    n = len(df)
    df = _step("cis_window",
               F.cis_window(df, protein.tss, protein.chrom, protein.build,
                            kb=window_kb, chrom_col="chrom", pos_col="pos"), n)
    n = len(df)
    df = _step("gw_significant",
               F.gw_significant(df, p=pval_gw, pval_col="pval"), n)
    n = len(df)
    df = _step("maf_above",
               F.maf_above(df, threshold=maf_min, eaf_col="EAF"), n)
    n = len(df)
    df = _step("exclude_mhc",
               F.exclude_mhc(df, protein.build, chrom_col="chrom", pos_col="pos"), n)
    n = len(df)
    df = _step("drop_ambig_palindromes",
               F.drop_ambig_palindromes(df, maf_threshold=pal_maf,
                                        ea_col="EA", oa_col="OA", eaf_col="EAF"), n)
    return df


def _normalize(df: pd.DataFrame, protein: ProteinMeta) -> pd.DataFrame:
    df = df.copy()
    df["seqid"] = protein.seqid
    df["gene"] = protein.gene
    df["uniprot"] = protein.uniprot
    df["build"] = protein.build
    # Ensure canonical column presence
    for col in OUTPUT_COLS:
        if col not in df.columns:
            df[col] = None
    return df[OUTPUT_COLS]


def _write_index(cohort: str, proteins: list[ProteinMeta]) -> None:
    """Write/update a per-cohort index of all proteins (for downstream joins)."""
    index_path = cohort_dir(cohort) / "protein_index.tsv"
    rows = [
        {"seqid": p.seqid, "gene": p.gene, "uniprot": p.uniprot,
         "chrom": p.chrom, "tss": p.tss, "build": p.build}
        for p in proteins
    ]
    pd.DataFrame(rows).to_csv(index_path, sep="\t", index=False)
