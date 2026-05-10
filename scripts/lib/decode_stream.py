"""deCODE HTTP streaming helpers.

deCODE files are accessed via signed URLs from data/raw/deCODE/bulk_urls.txt.
Server throttles at ~4 MB/s and drops parallel connections — download sequentially.
"""
import gzip
import io
import logging
import time
from pathlib import Path
from typing import Iterator

import httpx

log = logging.getLogger(__name__)

_TIMEOUT = httpx.Timeout(connect=30.0, read=300.0, write=10.0, pool=10.0)
_MAX_RETRIES = 3
_RETRY_SLEEP = 5.0


def download_bytes(url: str, retries: int = _MAX_RETRIES) -> bytes:
    """Download a URL to memory; retry on transient errors."""
    for attempt in range(1, retries + 1):
        try:
            with httpx.Client(timeout=_TIMEOUT, follow_redirects=True) as client:
                resp = client.get(url)
                resp.raise_for_status()
                return resp.content
        except (httpx.HTTPError, httpx.TimeoutException) as exc:
            if attempt == retries:
                raise
            log.warning(f"Download failed (attempt {attempt}/{retries}): {exc}. Retrying in {_RETRY_SLEEP}s...")
            time.sleep(_RETRY_SLEEP)
    raise RuntimeError("unreachable")


def iter_decode_rows(url: str, sep: str = "\t") -> Iterator[dict]:
    """Download a deCODE .txt.gz file and yield rows as dicts."""
    raw = download_bytes(url)
    with gzip.open(io.BytesIO(raw)) as fh:
        header = fh.readline().decode().strip().split(sep)
        for line in fh:
            parts = line.decode().strip().split(sep)
            if len(parts) == len(header):
                yield dict(zip(header, parts))


def parse_bulk_urls(urls_file: Path) -> list[tuple[str, str]]:
    """
    Parse deCODE bulk_urls.txt.
    Returns list of (protein_name, url) for per-aptamer .txt.gz files.
    protein_name is the filename without .txt.gz (e.g. '10000_28_CRYBB2_CRBB2').
    Skips annotation/md5sum files.
    """
    proteins = []
    with open(urls_file) as fh:
        for line in fh:
            url = line.strip()
            if not url:
                continue
            fname = url.split("file=")[-1] if "file=" in url else url.split("/")[-1]
            if not fname.endswith(".txt.gz"):
                continue
            # Only base deCODE protein files start with a digit
            if not fname[0].isdigit():
                continue
            protein_name = fname.removesuffix(".txt.gz")
            proteins.append((protein_name, url))
    return proteins
