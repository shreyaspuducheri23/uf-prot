"""cis-window and TSS lookup helpers."""
import logging
from functools import lru_cache
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

log = logging.getLogger(__name__)

_ENSEMBL_HG19_REST = "https://grch37.rest.ensembl.org"
_ENSEMBL_HG38_REST = "https://rest.ensembl.org"
_HEADERS = {"Content-Type": "application/json"}


@lru_cache(maxsize=10_000)
def tss_from_ensembl(gene_symbol: str, build: str) -> Optional[tuple[str, int]]:
    """
    Fetch TSS for a gene symbol from Ensembl REST API.
    Returns (chrom, tss_1based) or None on failure.
    build: 'hg19' or 'hg38'
    """
    base = _ENSEMBL_HG19_REST if build == "hg19" else _ENSEMBL_HG38_REST
    url = f"{base}/lookup/symbol/homo_sapiens/{gene_symbol}?expand=0"
    try:
        resp = requests.get(url, headers=_HEADERS, timeout=15)
        resp.raise_for_status()
        data = resp.json()
        chrom = str(data["seq_region_name"])
        strand = data["strand"]
        start = int(data["start"])
        end = int(data["end"])
        # Normalise strand to +1 / -1; API returns int 1/-1 or string "+"/"-"
        if strand in (1, "+", "1"):
            tss = start
        elif strand in (-1, "-", "-1"):
            tss = end
        else:
            raise ValueError(
                f"Ensembl returned unexpected strand {strand!r} for {gene_symbol!r}; "
                f"expected 1 or -1"
            )
        return chrom, tss
    except Exception as exc:
        log.debug(f"Ensembl TSS lookup failed for {gene_symbol!r} ({build}): {exc}")
        return None


def load_aric_tss(seqid_path: Path) -> dict[str, tuple[str, int, str, str]]:
    """
    Load ARIC seqid.txt: {seqid: (chrom, tss_hg19, uniprot, gene)}.
    """
    df = pd.read_csv(seqid_path, sep="\t", dtype=str)
    # Columns: seqid_in_sample, uniprot_id, entrezgenesymbol, chromosome_name, transcription_start_site
    result: dict[str, tuple[str, int, str, str]] = {}
    n_skipped = 0
    for _, row in df.iterrows():
        try:
            tss = int(row["transcription_start_site"])
        except (ValueError, TypeError):
            log.warning(
                f"load_aric_tss: skipping {row.get('seqid_in_sample', '?')!r} — "
                f"invalid TSS {row.get('transcription_start_site')!r}"
            )
            n_skipped += 1
            continue
        result[row["seqid_in_sample"]] = (
            str(row["chromosome_name"]),
            tss,
            row["uniprot_id"],
            row["entrezgenesymbol"],
        )
    if n_skipped:
        log.warning(f"load_aric_tss: {n_skipped} proteins skipped due to invalid TSS values")
    return result


def cis_window_bounds(tss: int, kb: int) -> tuple[int, int]:
    """Return (start, end) of a ±kb window around TSS (1-based, clamped at 0)."""
    flank = kb * 1_000
    return max(1, tss - flank), tss + flank
