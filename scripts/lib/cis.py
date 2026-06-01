"""cis-window and TSS lookup helpers."""
import dataclasses
import logging
import time
from functools import lru_cache
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests

log = logging.getLogger(__name__)

_ENSEMBL_HG19_REST = "https://grch37.rest.ensembl.org"
_ENSEMBL_HG38_REST = "https://rest.ensembl.org"
_HEADERS = {"Content-Type": "application/json"}
_HGNC_HEADERS = {"Accept": "application/json"}
_OVERRIDE_PATH = Path(__file__).with_name("tss_overrides.tsv")
_TSS_CACHE_COLUMNS = ["gene", "chrom", "tss", "resolved_symbol", "tier", "source"]
_UNRESOLVED_COLUMNS = ["gene", "build", "attempts"]


class _TransientError(Exception):
    """Raised when a REST lookup exhausted retries on transient errors."""


@dataclasses.dataclass(frozen=True)
class TssResolution:
    resolved: bool
    requested_symbol: str
    build: str
    chrom: str | None = None
    tss: int | None = None
    resolved_symbol: str | None = None
    tier: int = 0
    source: str = ""
    attempts: tuple[str, ...] = ()
    transient: bool = False


def _ensembl_lookup(gene_symbol: str, build: str) -> tuple[str, int] | None:
    """Fetch (chrom, TSS) from Ensembl REST, retrying transient failures."""
    base = _ENSEMBL_HG19_REST if build == "hg19" else _ENSEMBL_HG38_REST
    url = f"{base}/lookup/symbol/homo_sapiens/{gene_symbol}?expand=0"
    last_exc: Exception | None = None
    for attempt in range(3):
        try:
            resp = requests.get(url, headers=_HEADERS, timeout=15)
            if resp.status_code in (400, 404):
                return None
            resp.raise_for_status()
            data = resp.json()
            chrom = str(data["seq_region_name"])
            strand = data["strand"]
            start = int(data["start"])
            end = int(data["end"])
            if strand in (1, "+", "1"):
                return chrom, start
            if strand in (-1, "-", "-1"):
                return chrom, end
            raise ValueError(f"unexpected strand {strand!r} for {gene_symbol!r}")
        except ValueError:
            raise
        except requests.exceptions.HTTPError as exc:
            status = exc.response.status_code if exc.response is not None else None
            if status in (400, 404):
                return None
            last_exc = exc
        except Exception as exc:
            last_exc = exc
        if attempt < 2:
            time.sleep(2 ** attempt)
    raise _TransientError(
        f"Ensembl lookup for {gene_symbol!r} ({build}) failed after 3 attempts: {last_exc}"
    )


def _hgnc_approved_for_field(field: str, gene_symbol: str) -> list[str]:
    """Return approved HGNC symbols matching a prev_symbol or alias_symbol search."""
    url = f"https://rest.genenames.org/search/{field}/{quote(gene_symbol, safe='')}"
    last_exc: Exception | None = None
    for attempt in range(3):
        try:
            resp = requests.get(url, headers=_HGNC_HEADERS, timeout=10)
            if resp.status_code in (400, 404):
                return []
            resp.raise_for_status()
            data = resp.json()
            docs = data.get("response", {}).get("docs", [])
            return [str(doc["symbol"]) for doc in docs if doc.get("symbol")]
        except Exception as exc:
            last_exc = exc
        if attempt < 2:
            time.sleep(2 ** attempt)
    log.debug(f"HGNC lookup failed for {field}={gene_symbol!r}: {last_exc}")
    return []


def _hgnc_doc_for_symbol(gene_symbol: str) -> dict:
    """Return the HGNC document for an approved symbol, or {} on a miss."""
    url = f"https://rest.genenames.org/fetch/symbol/{quote(gene_symbol, safe='')}"
    last_exc: Exception | None = None
    for attempt in range(3):
        try:
            resp = requests.get(url, headers=_HGNC_HEADERS, timeout=10)
            if resp.status_code in (400, 404):
                return {}
            resp.raise_for_status()
            data = resp.json()
            docs = data.get("response", {}).get("docs", [])
            if not docs:
                return {}
            return docs[0]
        except Exception as exc:
            last_exc = exc
        if attempt < 2:
            time.sleep(2 ** attempt)
    log.debug(f"HGNC fetch failed for symbol={gene_symbol!r}: {last_exc}")
    return {}


@lru_cache(maxsize=1)
def _load_overrides() -> dict[tuple[str, str], tuple[str, int, str]]:
    """Load curated TSS overrides keyed by (GENE_SYMBOL_UPPER, build)."""
    if not _OVERRIDE_PATH.exists():
        log.warning(f"TSS override file not found: {_OVERRIDE_PATH}")
        return {}

    overrides: dict[tuple[str, str], tuple[str, int, str]] = {}
    try:
        df = pd.read_csv(_OVERRIDE_PATH, sep="\t", dtype=str)
        for _, row in df.iterrows():
            try:
                gene = str(row["gene_symbol"]).upper()
                build = str(row["build"]).lower()
                chrom = str(row["chrom"])
                tss = int(row["tss"])
                source = str(row.get("source", "tss_overrides.tsv"))
            except (KeyError, TypeError, ValueError):
                continue
            overrides[(gene, build)] = (chrom, tss, source)
    except Exception as exc:
        log.warning(f"TSS override read error ({_OVERRIDE_PATH}): {exc}")
        return {}
    return overrides


@lru_cache(maxsize=10_000)
def resolve_tss(gene_symbol: str, build: str) -> TssResolution:
    """
    Resolve a gene symbol to a TSS through Ensembl, HGNC previous/alias
    symbols, and curated overrides. Always returns a TssResolution.
    """
    requested = str(gene_symbol)
    build = str(build).lower()
    attempts: list[str] = []
    tried: set[str] = set()
    any_transient = False

    def try_ensembl(candidate: str) -> tuple[str, int] | None:
        nonlocal any_transient
        attempts.append(candidate)
        tried.add(candidate)
        try:
            return _ensembl_lookup(candidate, build)
        except _TransientError as exc:
            any_transient = True
            log.debug(str(exc))
            return None

    def resolved(
        candidate: str,
        chrom: str,
        tss: int,
        tier: int,
        source: str,
    ) -> TssResolution:
        log.info(f"Resolved {requested!r} via tier {tier}: {source}")
        return TssResolution(
            resolved=True,
            requested_symbol=requested,
            build=build,
            chrom=chrom,
            tss=tss,
            resolved_symbol=candidate,
            tier=tier,
            source=source,
            attempts=tuple(attempts),
        )

    hit = try_ensembl(requested)
    if hit:
        return resolved(requested, hit[0], hit[1], 1, "Ensembl")

    hgnc_doc = _hgnc_doc_for_symbol(requested)
    for field in ("prev_symbol", "alias_symbol"):
        values = hgnc_doc.get(field, [])
        if isinstance(values, str):
            values = [values]
        for prior_symbol in values:
            prior_symbol = str(prior_symbol)
            if prior_symbol in tried:
                continue
            hit = try_ensembl(prior_symbol)
            if hit:
                return resolved(
                    prior_symbol,
                    hit[0],
                    hit[1],
                    2,
                    f"HGNC current symbol {requested} {field} -> {prior_symbol}",
                )

    for approved in _hgnc_approved_for_field("prev_symbol", requested):
        hit = try_ensembl(approved)
        if hit:
            return resolved(approved, hit[0], hit[1], 2, f"HGNC prev_symbol -> {approved}")

    for approved in _hgnc_approved_for_field("alias_symbol", requested):
        if approved not in tried:
            hit = try_ensembl(approved)
            if hit:
                return resolved(approved, hit[0], hit[1], 3, f"HGNC alias_symbol -> {approved}")
        for prev in _hgnc_approved_for_field("prev_symbol", approved):
            if prev in tried:
                continue
            hit = try_ensembl(prev)
            if hit:
                return resolved(
                    prev,
                    hit[0],
                    hit[1],
                    3,
                    f"HGNC alias_symbol -> {approved}; prev_symbol -> {prev}",
                )

    overrides = _load_overrides()
    override = overrides.get((requested.upper(), build)) or overrides.get((requested.upper(), "any"))
    if override:
        chrom, tss, source = override
        return resolved(requested, chrom, tss, 4, f"tss_overrides.tsv: {source}")

    log.warning(f"Unresolved gene {requested!r} ({build}); tried {attempts}")
    return TssResolution(
        resolved=False,
        requested_symbol=requested,
        build=build,
        attempts=tuple(attempts),
        transient=any_transient,
    )


def _load_tss_cache(cache_path: Path, uppercase: bool = False) -> dict[str, tuple[str, int]]:
    """Load {gene: (chrom, tss)} from a legacy or provenance-aware TSS cache."""
    if not cache_path.exists():
        return {}
    try:
        df = pd.read_csv(cache_path, sep="\t", dtype=str, keep_default_na=False)
    except Exception as exc:
        log.warning(f"TSS cache read error ({cache_path}): {exc}")
        return {}

    result: dict[str, tuple[str, int]] = {}
    for _, row in df.iterrows():
        try:
            gene = str(row["gene"])
            if uppercase:
                gene = gene.upper()
            result[gene] = (str(row["chrom"]), int(row["tss"]))
        except (ValueError, KeyError, TypeError):
            continue
    return result


def _save_tss_cache(
    cache_path: Path,
    cache: dict[str, tuple[str, int]],
    rows: list[dict] | None = None,
) -> None:
    """Write a TSS cache with stable columns while preserving existing provenance."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    rows = rows or []
    existing = (
        pd.read_csv(cache_path, sep="\t", dtype=str, keep_default_na=False)
        if cache_path.exists()
        else pd.DataFrame(columns=_TSS_CACHE_COLUMNS)
    )

    existing_meta: dict[str, dict[str, str]] = {}
    if not existing.empty and "gene" in existing.columns:
        for _, row in existing.iterrows():
            existing_meta[str(row["gene"])] = {
                "resolved_symbol": row.get("resolved_symbol", ""),
                "tier": row.get("tier", ""),
                "source": row.get("source", ""),
            }

    cache_rows = []
    for gene, (chrom, tss) in sorted(cache.items()):
        cache_rows.append({
            "gene": gene,
            "chrom": chrom,
            "tss": tss,
            **existing_meta.get(gene, {}),
        })

    frames = [pd.DataFrame(cache_rows)]
    if rows:
        frames.append(pd.DataFrame(rows))
    updated = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    for col in _TSS_CACHE_COLUMNS:
        if col not in updated.columns:
            updated[col] = pd.NA
    updated = updated[_TSS_CACHE_COLUMNS].drop_duplicates("gene", keep="last")
    updated["tier"] = pd.to_numeric(updated["tier"], errors="coerce").astype(pd.Int64Dtype())
    updated.to_csv(cache_path, sep="\t", index=False)


def _append_unresolved(cohort_path: Path, rows: list[dict]) -> None:
    """Append unresolved TSS attempts to a cohort sidecar."""
    if not rows:
        return
    cohort_path.mkdir(parents=True, exist_ok=True)
    path = cohort_path / "_tss_unresolved.tsv"
    df = pd.DataFrame(rows)
    for col in _UNRESOLVED_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    df = df[_UNRESOLVED_COLUMNS]
    df.to_csv(path, sep="\t", index=False, mode="a", header=not path.exists())


def load_aric_tss(seqid_path: Path) -> dict[str, tuple[str, int, str, str]]:
    """
    Load ARIC seqid.txt: {seqid: (chrom, tss, uniprot, gene)}.
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
