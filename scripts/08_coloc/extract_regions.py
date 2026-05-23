#!/usr/bin/env python3
"""
08_coloc/extract_regions.py
Extract ±1 Mb cis regions for proteins passing MR + sensitivity.
Writes per-protein exposure and outcome TSVs for colocalization.

Usage:
  python scripts/08_coloc/extract_regions.py [--cohort ARIC_EA] [--limit N]
"""
import argparse

import pandas as pd

from scripts.lib.checkpoint import Checkpoint, output_exists
from scripts.lib.cis import cis_window_bounds
from scripts.lib.decode_stream import (
    _get_s3_client, stream_s3_cis_rows,
    DECODE_S3_ENDPOINT, DECODE_S3_BUCKET, DECODE_S3_ACCESS_KEY, DECODE_S3_SECRET_KEY,
)
from scripts.lib.liftover import lift_position, lift_table
from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.outcome import OutcomeLookup
from scripts.lib.paths import (
    ARIC_EA_DIR, COHORTS, COLOC_REGIONS_DIR,
    cohort_dir, cis_sumstats_hg38_dir,
)
from scripts.lib.progress import bar
from scripts.lib.sumstats_io import read_norm, write_norm

log = setup_logger("08_extract_regions")
EXPOSURE_REGION_REQUIRED_COLS = ["chrom", "pos"]
OUTCOME_REGION_REQUIRED_COLS = ["chromosome", "base_pair_location"]


def load_candidates(cohort: str) -> list[str]:
    """Return seqids that passed MR FDR and sensitivity."""
    mr_path   = cohort_dir(cohort) / "mr_results.tsv"
    sens_path = cohort_dir(cohort) / "sensitivity.tsv"

    if not mr_path.exists():
        log.warning(f"{cohort}: no mr_results.tsv")
        return []

    mr = pd.read_csv(mr_path, sep="\t")
    candidates = mr[mr["fdr_pass"] == True]["seqid"].tolist()

    if sens_path.exists():
        sens = pd.read_csv(sens_path, sep="\t")
        pass_sens = set(sens[sens["passes_sensitivity"] == True]["seqid"])
        candidates = [s for s in candidates if s in pass_sens or
                      s not in sens["seqid"].values]  # single-SNP proteins have no sensitivity

    log.info(f"{cohort}: {len(candidates)} candidates for colocalization")
    return candidates


def extract_aric_region(seqid: str, chrom: str, start: int, end: int) -> pd.DataFrame | None:
    """Extract a 1Mb region from ARIC EA local .glm.linear file."""
    import glob
    pattern = str(ARIC_EA_DIR / f"{seqid}.PHENO1.glm.linear")
    matches = glob.glob(pattern)
    if not matches:
        return None

    df = pd.read_csv(matches[0], sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    df = df[df["TEST"] == "ADD"]
    # plink2 .glm.linear: A1 is the tested (effect) allele; it can be REF or ALT.
    # OA must be whichever of REF/ALT is NOT A1 — not simply REF.
    import numpy as np
    df["OA"] = np.where(df["A1"] == df["REF"], df["ALT"], df["REF"])
    df = df.rename(columns={
        "CHROM": "chrom", "POS": "pos", "ID": "rsid",
        "A1": "EA", "A1_FREQ": "EAF",
        "BETA": "beta", "SE": "se", "P": "pval", "OBS_CT": "N",
    })
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = df["pos"].astype(int)
    df = df[(df["chrom"] == chrom) & (df["pos"] >= start) & (df["pos"] <= end)]
    # ARIC .glm.linear positions are already hg38 (the GWAS was run on GRCh38 despite
    # the protein_index build tag of hg19) — no liftover needed.
    return df if not df.empty else None


_DECODE_CIS_USECOLS = frozenset({
    "Chrom", "Pos", "rsids", "effectAllele", "otherAllele",
    "Beta", "SE", "Pval", "N", "ImpMAF",
})


def _load_decode_s3_key_index(cohort_d) -> dict[str, str]:
    """Load the cached S3 key index written by step 2 (decode.py)."""
    import json
    cache_path = cohort_d / "_s3_key_index.json"
    if not cache_path.exists():
        log.warning(f"deCODE S3 key index not found at {cache_path}")
        return {}
    with open(cache_path) as fh:
        return json.load(fh)


def extract_decode_region_s3(
    seqid: str,
    s3_client,
    s3_key_map: dict[str, str],
    chrom: str,
    tss: int,
    window_bp: int = 1_000_000,
) -> pd.DataFrame | None:
    """Extract a cis region for a deCODE protein via S3 streaming (v2 dataset)."""
    key = s3_key_map.get(seqid)
    if not key:
        log.warning(f"deCODE {seqid}: no S3 key in index — skipping")
        return None

    rows = list(stream_s3_cis_rows(
        s3_client, DECODE_S3_BUCKET, key,
        target_chrom=chrom, tss=tss, window_bp=window_bp,
        usecols=_DECODE_CIS_USECOLS,
    ))
    if not rows:
        return None

    df = pd.DataFrame(rows)
    df = df.rename(columns={
        "Chrom": "chrom", "Pos": "pos", "rsids": "rsid",
        "effectAllele": "EA", "otherAllele": "OA", "ImpMAF": "EAF",
        "Beta": "beta", "SE": "se", "Pval": "pval", "N": "N",
    })
    df["chrom"] = df["chrom"].astype(str).str.lstrip("chr")
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    return df if not df.empty else None


def _to_hg38_region(chrom: str, start: int, end: int, build: str) -> tuple[str, int, int] | None:
    build_norm = str(build).lower()
    if build_norm in {"hg38", "grch38"}:
        return chrom, start, end
    if build_norm not in {"hg19", "grch37"}:
        return None

    left = lift_position(chrom, start)
    right = lift_position(chrom, end)
    if left is None or right is None:
        return None
    if left[0] != right[0]:
        return None

    return left[0], min(left[1], right[1]), max(left[1], right[1])


def extract_cohort_regions(
    cohort: str,
    limit: int | None = None,
    retry_failed: bool = False,
) -> int:
    candidates = load_candidates(cohort)
    if not candidates:
        return 0
    if limit:
        candidates = candidates[:limit]

    # Load protein index for TSS info
    index_path = cohort_dir(cohort) / "protein_index.tsv"
    if not index_path.exists():
        log.warning(f"{cohort}: protein_index.tsv not found")
        return 0
    index = pd.read_csv(index_path, sep="\t", dtype=str)
    tss_map = {row["seqid"]: (str(row["chrom"]), int(row["tss"]), row["build"])
               for _, row in index.iterrows()}

    cp = Checkpoint(cohort_dir(cohort) / "_state_08_regions.json")
    todo = cp.remaining(candidates, include_failed=retry_failed)

    # Pre-load S3 resources for deCODE
    s3_key_map: dict[str, str] = {}
    s3_client = None
    if cohort == "deCODE":
        s3_key_map = _load_decode_s3_key_index(cohort_dir(cohort))
        s3_client = _get_s3_client(DECODE_S3_ENDPOINT, DECODE_S3_ACCESS_KEY, DECODE_S3_SECRET_KEY)

    n_ok = 0
    with OutcomeLookup() as outcome:
        for seqid in bar(todo, desc=f"{cohort} coloc regions"):
            try:
                if seqid not in tss_map:
                    cp.mark_failed(seqid, "seqid_not_in_protein_index")
                    continue

                chrom, tss, build = tss_map[seqid]
                start, end = cis_window_bounds(tss, kb=1000)

                out_dir = COLOC_REGIONS_DIR / cohort / seqid
                exp_path = out_dir / "exposure.tsv"
                out_path = out_dir / "outcome.tsv"

                if (
                    output_exists(exp_path, required_cols=EXPOSURE_REGION_REQUIRED_COLS, min_rows=1)
                    and output_exists(out_path, required_cols=OUTCOME_REGION_REQUIRED_COLS, min_rows=1)
                ):
                    cp.mark_done(seqid)
                    n_ok += 1
                    continue

                # Extract exposure cis region
                if cohort == "ARIC_EA":
                    exp_df = extract_aric_region(seqid, chrom, start, end)
                elif cohort == "deCODE":
                    exp_df = extract_decode_region_s3(seqid, s3_client, s3_key_map, chrom, tss, window_bp=1_000_000)
                else:
                    log.warning(f"Region re-extraction for {cohort} not yet implemented (use cached sumstats)")
                    # Fall back to pre-filtered cis_sumstats_hg38 (positions already in hg38)
                    cis_path = cis_sumstats_hg38_dir(cohort) / f"{seqid}.tsv"
                    exp_df = read_norm(cis_path) if cis_path.exists() else None

                if exp_df is None or exp_df.empty:
                    cp.mark_failed(seqid, "no_exposure_region_variants")
                    continue

                query_bounds = _to_hg38_region(chrom, start, end, build)
                if query_bounds is None:
                    cp.mark_failed(seqid, f"liftover_failed_for_{build}_region")
                    continue
                out_chrom, out_start, out_end = query_bounds

                out_df = outcome.fetch_region(out_chrom, out_start, out_end)
                if out_df.empty:
                    cp.mark_failed(seqid, "no_outcome_variants_in_hg38_region")
                    continue

                out_dir.mkdir(parents=True, exist_ok=True)
                write_norm(exp_df, exp_path)
                out_df.to_csv(out_path, sep="\t", index=False)
                cp.mark_done(seqid)
                n_ok += 1
            except Exception as exc:
                cp.mark_failed(seqid, f"exception:{exc.__class__.__name__}")
                log.warning(f"{cohort} {seqid}: region extraction failed — {exc}")

    log.info(
        f"{cohort}: {n_ok}/{len(candidates)} regions extracted "
        f"({cp.n_failed} failed in checkpoint)"
    )
    return n_ok


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract ±1Mb coloc regions")
    parser.add_argument("--cohort", choices=COHORTS + ["all"], default="all")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument(
        "--retry-failed",
        action="store_true",
        help="Include previously failed seqids from checkpoint.",
    )
    args = parser.parse_args()

    cohorts = COHORTS if args.cohort == "all" else [args.cohort]

    with RunManifest("08_coloc/extract_regions.py", args=str(args)) as manifest:
        total = sum(
            extract_cohort_regions(c, limit=args.limit, retry_failed=args.retry_failed)
            for c in cohorts
        )
        manifest.n_units = total


if __name__ == "__main__":
    main()
