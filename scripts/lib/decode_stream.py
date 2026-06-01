"""deCODE streaming helpers — legacy HTTP download and S3 direct streaming.

deCODE files are accessed directly from s3-ext.decode.is:10443 in the pipeline.
The HTTP helpers remain for legacy signed-URL fixtures and ad hoc checks.
Callers decide concurrency; this module focuses on robust retried downloads.
"""
# ---------------------------------------------------------------------------
# S3 connection defaults — shared by decode.py and extract_regions.py
# ---------------------------------------------------------------------------
DECODE_S3_ENDPOINT   = "https://s3-ext.decode.is:10443"
DECODE_S3_BUCKET     = "largescaleplasma-2023"
DECODE_S3_ACCESS_KEY = "SE0AV795UKCQ338YKWP4"
DECODE_S3_SECRET_KEY = "/mkkvYtFJkO+NAhxcm3OhNKAdvwQivhbdQRLeJ/c"
DECODE_S3_PREFIX_RAW = "final_somascan_raw"
import gzip
import io
import logging
import time
import atexit
from pathlib import Path
from typing import Iterator

import httpx

log = logging.getLogger(__name__)

_TIMEOUT = httpx.Timeout(connect=30.0, read=300.0, write=10.0, pool=10.0)
_LIMITS = httpx.Limits(max_keepalive_connections=8, max_connections=16)
_MAX_RETRIES = 3
_RETRY_SLEEP = 5.0
_CLIENT: httpx.Client | None = None


def _get_client() -> httpx.Client:
    global _CLIENT
    if _CLIENT is None:
        _CLIENT = httpx.Client(
            timeout=_TIMEOUT,
            follow_redirects=True,
            limits=_LIMITS,
            http2=True,
        )
    return _CLIENT


def _close_client() -> None:
    global _CLIENT
    if _CLIENT is not None:
        _CLIENT.close()
        _CLIENT = None


atexit.register(_close_client)


def download_bytes(url: str, retries: int = _MAX_RETRIES) -> bytes:
    """Download a URL to memory; retry on transient errors."""
    client = _get_client()
    for attempt in range(1, retries + 1):
        try:
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


_S3_CLIENT = None


def _get_s3_client(endpoint_url: str, access_key: str, secret_key: str, region: str = "us-east-1"):
    """Return a cached module-level boto3 S3 client (thread-safe for read-only calls)."""
    global _S3_CLIENT
    if _S3_CLIENT is None:
        import boto3
        _S3_CLIENT = boto3.client(
            "s3",
            endpoint_url=endpoint_url,
            aws_access_key_id=access_key,
            aws_secret_access_key=secret_key,
            region_name=region,
        )
    return _S3_CLIENT


def stream_s3_cis_rows(
    s3_client,
    bucket: str,
    key: str,
    target_chrom: str,
    tss: int,
    window_bp: int,
    usecols: frozenset[str],
) -> Iterator[dict]:
    """Stream a deCODE S3 .txt.gz, yielding only rows in the cis window.

    Files are sorted chr1→chr22 by position, so we can abort once we've
    passed the cis window on the target chromosome — saving 50-95% of bytes.
    """
    body = s3_client.get_object(Bucket=bucket, Key=key)["Body"]
    with gzip.open(io.BufferedReader(body), "rt") as fh:
        header = fh.readline().strip().split("\t")
        col_idx = {c: i for i, c in enumerate(header) if c in usecols}
        chrom_i = header.index("Chrom")
        pos_i = header.index("Pos")

        seen_target = False
        low = tss - window_bp
        high = tss + window_bp
        min_cols = max(col_idx.values()) + 1

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < min_cols:
                continue
            row_chrom = parts[chrom_i].lstrip("chr")
            if row_chrom == target_chrom:
                seen_target = True
                try:
                    pos = int(parts[pos_i])
                except ValueError:
                    continue
                if pos > high:
                    break
                if pos < low:
                    continue
                yield {c: parts[i] for c, i in col_idx.items() if i < len(parts)}
            elif seen_target:
                break


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
