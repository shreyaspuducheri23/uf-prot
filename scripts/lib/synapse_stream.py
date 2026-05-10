"""Synapse streaming helpers for UKB-PPP and Fenland.

Downloads one file at a time, extracts the relevant in-memory, saves matched
rows, then deletes. Peak disk ≈ file size per worker.
Credentials read from ~/.synapseConfig.
"""
import gzip
import io
import logging
import os
import tarfile
from pathlib import Path
from typing import Callable, Iterator

import synapseclient

log = logging.getLogger(__name__)

_syn: synapseclient.Synapse | None = None


def _get_syn() -> synapseclient.Synapse:
    global _syn
    if _syn is None:
        _syn = synapseclient.Synapse()
        _syn.login(silent=True)
    return _syn


def list_folder(folder_id: str) -> list[synapseclient.Entity]:
    """Return all file entities in a Synapse folder (non-recursive)."""
    syn = _get_syn()
    return list(syn.getChildren(folder_id, includeTypes=["file"]))


def download_entity(entity_id: str, dest_dir: Path) -> Path:
    """Download a Synapse entity to dest_dir; return path to downloaded file."""
    syn = _get_syn()
    dest_dir.mkdir(parents=True, exist_ok=True)
    entity = syn.get(entity_id, downloadLocation=str(dest_dir))
    return Path(entity.path)


def iter_gz_rows(gz_bytes: bytes, sep: str = "\t") -> Iterator[dict]:
    """Yield rows from a gzip-compressed TSV/text file given as raw bytes."""
    with gzip.open(io.BytesIO(gz_bytes)) as fh:
        header = fh.readline().decode().strip().split(sep)
        for line in fh:
            parts = line.decode().strip().split(sep)
            if len(parts) == len(header):
                yield dict(zip(header, parts))


def iter_tar_gz_rows(tar_path: Path, member_filter: Callable[[str], bool] | None = None,
                     sep: str = "\t") -> Iterator[tuple[str, dict]]:
    """
    Iterate over rows in all .gz members of a tar archive.
    Yields (member_name, row_dict).
    member_filter: if provided, only process members where filter(name) is True.
    """
    with tarfile.open(tar_path, "r") as tf:
        for member in tf.getmembers():
            if not member.name.endswith(".gz"):
                continue
            if member_filter and not member_filter(member.name):
                continue
            fobj = tf.extractfile(member)
            if fobj is None:
                continue
            raw = fobj.read()
            for row in iter_gz_rows(raw, sep=sep):
                yield member.name, row


def stream_ukbppp_protein(entity_id: str,
                           cis_chrom: str, cis_start: int, cis_end: int,
                           tmp_dir: Path) -> list[dict]:
    """
    Download one UKB-PPP tar (≈550 MB), keep rows within the cis window,
    delete tar. Returns matched rows as list of dicts.

    UKB-PPP ID format (from MR_IA notes): position is 2nd colon-component.
    Columns: ID, ALLELE0, ALLELE1, A1FREQ, BETA, SE, LOG10P, N, ...
    """
    tar_path = download_entity(entity_id, tmp_dir)
    try:
        matched: list[dict] = []
        for _member, row in iter_tar_gz_rows(tar_path):
            try:
                parts = row.get("ID", "").split(":")
                if len(parts) < 2:
                    continue
                chrom = parts[0].lstrip("chr") if parts[0].startswith("chr") else parts[0]
                pos = int(parts[1])
            except (ValueError, IndexError):
                continue
            if chrom == cis_chrom and cis_start <= pos <= cis_end:
                matched.append(row)
    finally:
        tar_path.unlink(missing_ok=True)

    return matched


def stream_fenland_protein(entity_id: str,
                            cis_chrom: str, cis_start: int, cis_end: int,
                            tmp_dir: Path) -> list[dict]:
    """
    Download one Fenland .txt.gz file from Synapse, keep rows in the cis window,
    delete file.
    """
    gz_path = download_entity(entity_id, tmp_dir)
    try:
        matched: list[dict] = []
        with gzip.open(gz_path, "rt") as fh:
            header = fh.readline().strip().split("\t")
            for line in fh:
                parts = line.strip().split("\t")
                row = dict(zip(header, parts))
                try:
                    chrom = str(row.get("CHR", "")).lstrip("chr")
                    pos = int(row.get("POS", 0))
                except ValueError:
                    continue
                if chrom == cis_chrom and cis_start <= pos <= cis_end:
                    matched.append(row)
    finally:
        gz_path.unlink(missing_ok=True)

    return matched
