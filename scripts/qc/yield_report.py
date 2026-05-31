"""
Stage-aware yield/QC reporter.

The report audits filesystem artifacts and checkpoints without rerunning any
pipeline step. It tracks unit/protein yield plus row and unique-locus yield for
variant-bearing stages, then records per-unit and dropped-locus detail tables.

Usage:
    uv run python scripts/qc/yield_report.py [--cohort COHORT|all] [--strict]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable

import pandas as pd

from scripts.lib.paths import COHORTS, PROCESSED

WARN_FAIL_FRAC = 0.05
WARN_LIFTOVER_LOCUS_DROP_FRAC = 0.0

DETAIL_DROPPED_LOCUS_KINDS = {"clump", "liftover", "harmonise"}
WARN_DONE_WITHOUT_OUTPUT_STAGES = {"mr", "sensitivity", "sharepro", "coloc_abf"}
LIFTOVER_STAGES = {"instruments_hg38", "filtered_cis_pqtls_hg38", "raw_cis_sumstats_hg38"}


@dataclass(frozen=True)
class CheckpointStats:
    done: set[str] = field(default_factory=set)
    failed: set[str] = field(default_factory=set)
    failed_reasons: dict[str, str] = field(default_factory=dict)


@dataclass
class FileStats:
    seqid: str
    path: Path
    n_rows: int
    loci: set[str] = field(default_factory=set)
    fingerprints: dict[str, str] = field(default_factory=dict)
    mr_keep_rows: int | None = None


@dataclass
class ReportBundle:
    rows: list[dict]
    warnings: list[str]
    unit_rows: list[dict]
    dropped_locus_rows: list[dict]


def _strip_known_suffix(path: Path) -> str:
    name = path.name
    for suffix in (".tsv.gz", ".tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _read_index(cdir: Path) -> pd.DataFrame:
    path = cdir / "protein_index.tsv"
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return pd.DataFrame()


def _existing_dir(cdir: Path, names: Iterable[str]) -> Path | None:
    candidates = [cdir / name for name in names]
    nonempty = [
        path
        for path in candidates
        if path.exists() and (any(path.glob("*.tsv")) or any(path.glob("*.tsv.gz")))
    ]
    if nonempty:
        return nonempty[0]
    existing = [path for path in candidates if path.exists()]
    return existing[0] if existing else None


def _tsv_paths(d: Path | None, *, include_gz: bool = True) -> list[Path]:
    if d is None or not d.exists():
        return []
    paths = list(d.glob("*.tsv"))
    if include_gz:
        paths.extend(d.glob("*.tsv.gz"))
    return sorted(paths)


def _empty_marker_units(d: Path | None) -> set[str]:
    if d is None or not d.exists():
        return set()
    units: set[str] = set()
    for path in d.glob("*.tsv.empty"):
        name = path.name.removesuffix(".tsv.empty")
        units.add(name)
    return units


def _valid_id(value) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip() not in {"", ".", "nan", "None", "<NA>"}


def _first_present(columns: Iterable[str], candidates: Iterable[str]) -> str | None:
    present = set(columns)
    for candidate in candidates:
        if candidate in present:
            return candidate
    return None


def _read_file_stats(path: Path, *, prefer_hg38: bool = False) -> FileStats:
    seqid = _strip_known_suffix(path)
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return FileStats(seqid=seqid, path=path, n_rows=0)

    n_rows = len(df)
    if df.empty:
        return FileStats(seqid=seqid, path=path, n_rows=0)

    if prefer_hg38:
        chrom_candidates = (
            "chrom_hg38", "outcome_chrom_hg38", "chromosome",
            "chr.exposure", "chr.outcome", "chrom",
        )
        pos_candidates = (
            "pos_hg38", "outcome_pos_hg38", "base_pair_location",
            "pos.exposure", "pos.outcome", "pos",
        )
    else:
        chrom_candidates = ("chrom", "chrom_hg38", "chromosome", "chr.exposure", "chr.outcome")
        pos_candidates = ("pos", "pos_hg38", "base_pair_location", "pos.exposure", "pos.outcome")

    chrom_col = _first_present(df.columns, chrom_candidates)
    pos_col = _first_present(df.columns, pos_candidates)
    rsid_col = _first_present(df.columns, ("rsid", "SNP", "rsid.exposure", "outcome_rsid", "variant_id"))
    ea_col = _first_present(df.columns, ("EA", "effect_allele", "EA_exp", "effect_allele.exposure"))
    oa_col = _first_present(df.columns, ("OA", "other_allele", "OA_exp", "other_allele.exposure"))

    loci: set[str] = set()
    fingerprints: dict[str, str] = {}
    for i, row in df.iterrows():
        chrom = None
        pos = None
        locus = None
        if chrom_col and pos_col and _valid_id(row.get(chrom_col)) and _valid_id(row.get(pos_col)):
            chrom = str(row[chrom_col]).strip().removeprefix("chr")
            pos_raw = str(row[pos_col]).strip()
            try:
                pos = str(int(float(pos_raw)))
            except ValueError:
                pos = pos_raw
            locus = f"{chrom}:{pos}"
            loci.add(locus)

        if rsid_col and _valid_id(row.get(rsid_col)):
            fingerprint = str(row[rsid_col]).strip()
        elif locus is not None:
            ea = str(row.get(ea_col, "")).strip().upper() if ea_col else ""
            oa = str(row.get(oa_col, "")).strip().upper() if oa_col else ""
            fingerprint = f"{locus}:{ea}:{oa}" if ea or oa else locus
        else:
            fingerprint = f"row:{i}"

        fingerprints[fingerprint] = locus or fingerprint

    mr_keep_rows = None
    if "mr_keep" in df.columns:
        mr_keep_rows = int(df["mr_keep"].astype(str).str.lower().isin({"true", "1", "t"}).sum())

    return FileStats(
        seqid=seqid,
        path=path,
        n_rows=n_rows,
        loci=loci,
        fingerprints=fingerprints,
        mr_keep_rows=mr_keep_rows,
    )


def _file_stats_map(d: Path | None, *, prefer_hg38: bool = False) -> dict[str, FileStats]:
    return {
        _strip_known_suffix(path): _read_file_stats(path, prefer_hg38=prefer_hg38)
        for path in _tsv_paths(d)
    }


def _read_json_checkpoint(cp_path: Path) -> CheckpointStats:
    if not cp_path.exists():
        return CheckpointStats()
    try:
        data = json.loads(cp_path.read_text())
    except (json.JSONDecodeError, OSError):
        return CheckpointStats()

    done = set(map(str, data.get("done", [])))
    failed: set[str] = set()
    reasons: dict[str, str] = {}
    status = data.get("status", {})
    if isinstance(status, dict):
        for key, payload in status.items():
            if isinstance(payload, dict) and payload.get("state") == "failed":
                failed.add(str(key))
                reasons[str(key)] = str(payload.get("reason", ""))
    return CheckpointStats(done=done, failed=failed, failed_reasons=reasons)


def _read_r_checkpoint(cp_path: Path) -> CheckpointStats:
    if not cp_path.exists():
        return CheckpointStats()

    code = f"""
    x <- tryCatch(readRDS({str(cp_path)!r}), error = function(e) list(done = character(0)))
    done <- if (!is.null(x$done)) as.character(x$done) else character(0)
    failed <- character(0)
    reasons <- character(0)
    if (!is.null(x$failed)) {{
      if (is.list(x$failed)) {{
        failed <- names(x$failed)
        reasons <- vapply(x$failed, function(v) {{
          if (is.list(v) && !is.null(v$reason)) as.character(v$reason) else ""
        }}, character(1))
      }} else {{
        failed <- as.character(x$failed)
        reasons <- rep("", length(failed))
      }}
    }}
    cat("__DONE__\\n")
    if (length(done)) cat(done, sep = "\\n")
    cat("\\n__FAILED__\\n")
    if (length(failed)) cat(paste(failed, reasons, sep = "\\t"), sep = "\\n")
    """
    try:
        res = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
    except OSError:
        return CheckpointStats()
    if res.returncode != 0:
        return CheckpointStats()

    section = None
    done: set[str] = set()
    failed: set[str] = set()
    reasons: dict[str, str] = {}
    for raw in res.stdout.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line == "__DONE__":
            section = "done"
            continue
        if line == "__FAILED__":
            section = "failed"
            continue
        if section == "done":
            done.add(line)
        elif section == "failed":
            parts = line.split("\t", 1)
            failed.add(parts[0])
            reasons[parts[0]] = parts[1] if len(parts) > 1 else ""
    return CheckpointStats(done=done, failed=failed, failed_reasons=reasons)


def _read_checkpoint(cdir: Path, names: Iterable[str]) -> CheckpointStats:
    merged = CheckpointStats()
    done: set[str] = set()
    failed: set[str] = set()
    reasons: dict[str, str] = {}
    for name in names:
        path = cdir / name
        if name.endswith(".rds"):
            cp = _read_r_checkpoint(path)
        else:
            cp = _read_json_checkpoint(path)
        done.update(cp.done)
        failed.update(cp.failed)
        reasons.update(cp.failed_reasons)
    return CheckpointStats(done=done, failed=failed, failed_reasons=reasons)


def _pct(n: int | float | None, d: int | float | None) -> float:
    if not d:
        return 0.0
    return round(float(n or 0) / float(d) * 100, 1)


def _read_table(path: Path, *, cohort: str | None = None) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return pd.DataFrame()
    if cohort and "cohort" in df.columns:
        df = df[df["cohort"] == cohort].copy()
    return df


def _candidate_coloc_seqids(cdir: Path) -> set[str]:
    mr = _read_table(cdir / "mr_results.tsv")
    if mr.empty or "seqid" not in mr.columns or "fdr_pass" not in mr.columns:
        return set()
    candidates = set(mr.loc[mr["fdr_pass"].astype(str).str.lower() == "true", "seqid"].astype(str))
    sens = _read_table(cdir / "sensitivity.tsv")
    if not sens.empty and {"seqid", "passes_sensitivity"}.issubset(sens.columns):
        pass_sens = set(
            sens.loc[sens["passes_sensitivity"].astype(str).str.lower() == "true", "seqid"].astype(str)
        )
        sens_all = set(sens["seqid"].astype(str))
        candidates = {seqid for seqid in candidates if seqid in pass_sens or seqid not in sens_all}
    return candidates


def _make_aggregate_row(
    *,
    ts: str,
    cohort: str,
    stage: str,
    input_units: set[str],
    output_units: set[str],
    input_stats: dict[str, FileStats] | None = None,
    output_stats: dict[str, FileStats] | None = None,
    cp: CheckpointStats | None = None,
    empty_units: set[str] | None = None,
    not_applicable_units: set[str] | None = None,
    kind: str = "unit",
) -> dict:
    cp = cp or CheckpointStats()
    input_stats = input_stats or {}
    output_stats = output_stats or {}
    empty_units = empty_units or set()
    not_applicable_units = not_applicable_units or set()

    rows_input = sum(input_stats[seqid].n_rows for seqid in input_units if seqid in input_stats)
    rows_output = sum(output_stats[seqid].n_rows for seqid in output_units if seqid in output_stats)
    loci_input_set = set().union(*(input_stats[seqid].loci for seqid in input_units if seqid in input_stats))
    loci_output_set = set().union(*(output_stats[seqid].loci for seqid in output_units if seqid in output_stats))
    loci_input = len(loci_input_set)
    loci_output = len(loci_output_set)

    done_without_output_units = cp.done & input_units - output_units - empty_units - not_applicable_units
    not_done_units = input_units - output_units - empty_units - not_applicable_units - cp.failed
    if cp.done:
        not_done_units = input_units - cp.done - cp.failed

    units_input = len(input_units)
    units_output = len(output_units)
    units_not_done = max(0, len(not_done_units))
    units_failed_cp = len(cp.failed & input_units) if input_units else len(cp.failed)
    units_done_cp = len(cp.done & input_units) if input_units else len(cp.done)
    units_done_without_output = max(0, len(done_without_output_units))
    rows_dropped = max(0, rows_input - rows_output)
    loci_dropped = max(0, loci_input - loci_output)

    return {
        "ts": ts,
        "cohort": cohort,
        "stage": stage,
        "kind": kind,
        "units_input": units_input,
        "units_output": units_output,
        "pct_unit_yield": _pct(units_output, units_input),
        "rows_input": rows_input,
        "rows_output": rows_output,
        "pct_row_yield": _pct(rows_output, rows_input),
        "loci_input": loci_input,
        "loci_output": loci_output,
        "pct_locus_yield": _pct(loci_output, loci_input),
        "rows_dropped": rows_dropped,
        "loci_dropped": loci_dropped,
        "units_not_done": units_not_done,
        "units_failed_cp": units_failed_cp,
        "units_done_without_output": units_done_without_output,
        "units_not_applicable": len(not_applicable_units),
        "n_done_cp": units_done_cp,
        # Backward-compatible aliases for older callers/tests.
        "n_input": units_input,
        "n_output": units_output,
        "pct_yield": _pct(units_output, units_input),
        "n_failed": max(0, units_input - units_done_cp) if cp.done or cp.failed else units_not_done,
        "pct_failed": _pct(max(0, units_input - units_done_cp) if cp.done or cp.failed else units_not_done, units_input),
    }


def _make_unit_rows(
    *,
    ts: str,
    cohort: str,
    stage: str,
    input_units: set[str],
    output_units: set[str],
    input_stats: dict[str, FileStats] | None = None,
    output_stats: dict[str, FileStats] | None = None,
    cp: CheckpointStats | None = None,
    empty_units: set[str] | None = None,
    not_applicable_units: set[str] | None = None,
) -> list[dict]:
    input_stats = input_stats or {}
    output_stats = output_stats or {}
    cp = cp or CheckpointStats()
    empty_units = empty_units or set()
    not_applicable_units = not_applicable_units or set()
    rows: list[dict] = []
    all_units = sorted(input_units | output_units | empty_units | not_applicable_units | cp.failed)
    for seqid in all_units:
        in_stat = input_stats.get(seqid)
        out_stat = output_stats.get(seqid)
        if seqid in output_units:
            status = "output"
        elif seqid in empty_units:
            status = "empty"
        elif seqid in not_applicable_units:
            status = "not_applicable"
        elif seqid in cp.failed:
            status = "failed"
        elif seqid in cp.done:
            status = "done_without_output"
        else:
            status = "not_done"
        rows.append({
            "ts": ts,
            "cohort": cohort,
            "stage": stage,
            "seqid": seqid,
            "status": status,
            "input_rows": in_stat.n_rows if in_stat else 0,
            "output_rows": out_stat.n_rows if out_stat else 0,
            "input_loci": len(in_stat.loci) if in_stat else 0,
            "output_loci": len(out_stat.loci) if out_stat else 0,
            "input_path": str(in_stat.path) if in_stat else "",
            "output_path": str(out_stat.path) if out_stat else "",
            "failure_reason": cp.failed_reasons.get(seqid, ""),
        })
    return rows


def _make_dropped_locus_rows(
    *,
    ts: str,
    cohort: str,
    stage: str,
    input_units: set[str],
    output_units: set[str],
    input_stats: dict[str, FileStats],
    output_stats: dict[str, FileStats],
) -> list[dict]:
    rows: list[dict] = []
    for seqid in sorted(input_units & set(input_stats)):
        in_stat = input_stats[seqid]
        out_fps = output_stats.get(seqid, FileStats(seqid, Path(""), 0)).fingerprints
        for fingerprint, locus in sorted(in_stat.fingerprints.items()):
            if seqid not in output_units or fingerprint not in out_fps:
                rows.append({
                    "ts": ts,
                    "cohort": cohort,
                    "stage": stage,
                    "seqid": seqid,
                    "variant_id": fingerprint,
                    "locus": locus,
                    "input_path": str(in_stat.path),
                    "output_path": str(output_stats[seqid].path) if seqid in output_stats else "",
                })
    return rows


def _warn_for_row(
    row: dict,
    *,
    warn_fail_frac: float,
    warn_liftover_locus_drop_frac: float,
) -> list[str]:
    warnings: list[str] = []
    stage = row["stage"]
    units_input = int(row["units_input"])
    if units_input and row["n_failed"] / units_input > warn_fail_frac:
        warnings.append(
            f"[WARN] {row['cohort']}/{stage}: {row['n_failed']}/{units_input} units not done "
            f"({row['pct_failed']:.1f}%) exceeds {warn_fail_frac:.0%} threshold"
        )
    if stage in LIFTOVER_STAGES and row["loci_input"]:
        frac = row["loci_dropped"] / row["loci_input"]
        if row["loci_dropped"] > 0 and frac > warn_liftover_locus_drop_frac:
            warnings.append(
                f"[WARN] {row['cohort']}/{stage}: {row['loci_dropped']}/{row['loci_input']} loci "
                f"dropped ({frac * 100:.2f}%) during liftover"
            )
    if stage in WARN_DONE_WITHOUT_OUTPUT_STAGES and row["units_done_without_output"] > 0:
        warnings.append(
            f"[WARN] {row['cohort']}/{stage}: {row['units_done_without_output']} checkpoint-done units "
            f"have no result row/output"
        )
    return warnings


def _append_stage(
    bundle: ReportBundle,
    *,
    ts: str,
    cohort: str,
    stage: str,
    input_units: set[str],
    output_units: set[str],
    input_stats: dict[str, FileStats] | None = None,
    output_stats: dict[str, FileStats] | None = None,
    cp: CheckpointStats | None = None,
    empty_units: set[str] | None = None,
    not_applicable_units: set[str] | None = None,
    kind: str = "unit",
    warn_fail_frac: float = WARN_FAIL_FRAC,
    warn_liftover_locus_drop_frac: float = WARN_LIFTOVER_LOCUS_DROP_FRAC,
) -> None:
    if not input_units and not output_units and not (cp and (cp.done or cp.failed)):
        return

    row = _make_aggregate_row(
        ts=ts,
        cohort=cohort,
        stage=stage,
        input_units=input_units,
        output_units=output_units,
        input_stats=input_stats,
        output_stats=output_stats,
        cp=cp,
        empty_units=empty_units,
        not_applicable_units=not_applicable_units,
        kind=kind,
    )
    bundle.rows.append(row)
    bundle.warnings.extend(
        _warn_for_row(
            row,
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )
    )
    bundle.unit_rows.extend(
        _make_unit_rows(
            ts=ts,
            cohort=cohort,
            stage=stage,
            input_units=input_units,
            output_units=output_units,
            input_stats=input_stats,
            output_stats=output_stats,
            cp=cp,
            empty_units=empty_units,
            not_applicable_units=not_applicable_units,
        )
    )
    if kind in DETAIL_DROPPED_LOCUS_KINDS and input_stats and output_stats:
        bundle.dropped_locus_rows.extend(
            _make_dropped_locus_rows(
                ts=ts,
                cohort=cohort,
                stage=stage,
                input_units=input_units,
                output_units=output_units,
                input_stats=input_stats,
                output_stats=output_stats,
            )
        )


def report_cohort_details(
    cohort: str,
    processed_dir: Path | None = None,
    *,
    warn_fail_frac: float = WARN_FAIL_FRAC,
    warn_liftover_locus_drop_frac: float = WARN_LIFTOVER_LOCUS_DROP_FRAC,
) -> ReportBundle:
    if processed_dir is None:
        processed_dir = PROCESSED

    cdir = processed_dir / cohort
    index = _read_index(cdir)
    if index.empty or "seqid" not in index.columns:
        return ReportBundle([], [], [], [])

    ts = datetime.now(UTC).isoformat(timespec="seconds")
    bundle = ReportBundle([], [], [], [])
    index_units = set(index["seqid"].dropna().astype(str))

    raw_dir = _existing_dir(cdir, ("raw_cis_sumstats", "cis_raw"))
    filtered_dir = _existing_dir(cdir, ("filtered_cis_pqtls", "cis_sumstats"))
    filtered_hg38_dir = _existing_dir(cdir, ("filtered_cis_pqtls_hg38", "cis_sumstats_hg38"))
    instruments_dir = _existing_dir(cdir, ("instruments",))
    instruments_hg38_dir = _existing_dir(cdir, ("instruments_hg38",))
    raw_hg38_dir = _existing_dir(cdir, ("raw_cis_sumstats_hg38",))
    harmonised_dir = _existing_dir(cdir, ("harmonised",))

    raw_stats = _file_stats_map(raw_dir)
    filtered_stats = _file_stats_map(filtered_dir)
    instruments_stats = _file_stats_map(instruments_dir)
    instruments_hg38_stats = _file_stats_map(instruments_hg38_dir, prefer_hg38=True)
    filtered_hg38_stats = _file_stats_map(filtered_hg38_dir, prefer_hg38=True)
    raw_hg38_stats = _file_stats_map(raw_hg38_dir, prefer_hg38=True)
    harmonised_stats = _file_stats_map(harmonised_dir, prefer_hg38=True)

    cp02 = _read_checkpoint(cdir, ("_state_02.json",))
    cp03 = _read_checkpoint(cdir, ("_state_03.json",))
    cp04 = _read_checkpoint(cdir, ("_state_04.json",))
    cp04_filtered = _read_checkpoint(cdir, ("_state_04_filtered_cis.json", "_state_04_cis.json"))
    cp04_raw = _read_checkpoint(cdir, ("_state_04_raw_cis.json",))
    cp05 = _read_checkpoint(cdir, ("_state_05.json",))

    if raw_stats:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="raw_cis_sumstats",
            input_units=index_units,
            output_units=set(raw_stats),
            output_stats=raw_stats,
            cp=cp02,
            kind="raw",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    filtered_input_units = set(raw_stats) if raw_stats else index_units
    filtered_input_stats = raw_stats if raw_stats else {}
    if filtered_stats or raw_stats or cp02.done or cp02.failed:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="filtered_cis_pqtls",
            input_units=filtered_input_units,
            output_units=set(filtered_stats),
            input_stats=filtered_input_stats,
            output_stats=filtered_stats,
            cp=cp02,
            empty_units=_empty_marker_units(filtered_dir),
            kind="filtered",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    if instruments_stats or cp03.done or cp03.failed:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="instruments",
            input_units=set(filtered_stats),
            output_units=set(instruments_stats),
            input_stats=filtered_stats,
            output_stats=instruments_stats,
            cp=cp03,
            kind="clump",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    if instruments_hg38_stats or cp04.done or cp04.failed:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="instruments_hg38",
            input_units=set(instruments_stats),
            output_units=set(instruments_hg38_stats),
            input_stats=instruments_stats,
            output_stats=instruments_hg38_stats,
            cp=cp04,
            kind="liftover",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    if filtered_hg38_stats or cp04_filtered.done or cp04_filtered.failed:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="filtered_cis_pqtls_hg38",
            input_units=set(filtered_stats),
            output_units=set(filtered_hg38_stats),
            input_stats=filtered_stats,
            output_stats=filtered_hg38_stats,
            cp=cp04_filtered,
            kind="liftover",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    if raw_hg38_stats or cp04_raw.done or cp04_raw.failed:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="raw_cis_sumstats_hg38",
            input_units=set(raw_stats),
            output_units=set(raw_hg38_stats),
            input_stats=raw_stats,
            output_stats=raw_hg38_stats,
            cp=cp04_raw,
            kind="liftover",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    if harmonised_stats or cp05.done or cp05.failed:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="harmonised",
            input_units=set(instruments_hg38_stats),
            output_units=set(harmonised_stats),
            input_stats=instruments_hg38_stats,
            output_stats=harmonised_stats,
            cp=cp05,
            kind="harmonise",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    cp06 = _read_checkpoint(cdir, ("_state_06.rds",))
    mr_path = cdir / "mr_results.tsv"
    mr = _read_table(mr_path)
    mr_output_units = set(mr["seqid"].astype(str)) if "seqid" in mr.columns else set()
    mr_input_stats = {}
    mr_input_units: set[str] = set()
    for seqid, stat in harmonised_stats.items():
        keep = stat.mr_keep_rows if stat.mr_keep_rows is not None else stat.n_rows
        if keep and keep > 0:
            mr_input_units.add(seqid)
            mr_input_stats[seqid] = FileStats(
                seqid=seqid,
                path=stat.path,
                n_rows=keep,
                loci=stat.loci,
                fingerprints=stat.fingerprints,
            )
    mr_output_stats = {}
    if not mr.empty and "seqid" in mr.columns:
        for _, row in mr.iterrows():
            seqid = str(row["seqid"])
            n_snps = int(float(row.get("n_snps", 1))) if pd.notna(row.get("n_snps", 1)) else 1
            mr_output_stats[seqid] = FileStats(seqid=seqid, path=mr_path, n_rows=n_snps)
    if mr_output_units or cp06.done or cp06.failed:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="mr",
            input_units=mr_input_units,
            output_units=mr_output_units,
            input_stats=mr_input_stats,
            output_stats=mr_output_stats,
            cp=cp06,
            kind="analysis",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    cp07 = _read_checkpoint(cdir, ("_state_07.rds",))
    sensitivity_path = cdir / "sensitivity.tsv"
    sens = _read_table(sensitivity_path)
    sensitivity_output_units = set(sens["seqid"].astype(str)) if "seqid" in sens.columns else set()
    sensitivity_input_units: set[str] = set()
    sensitivity_not_applicable: set[str] = set()
    if not mr.empty and {"seqid", "n_snps"}.issubset(mr.columns):
        n_snps = pd.to_numeric(mr["n_snps"], errors="coerce").fillna(0)
        sensitivity_input_units = set(mr.loc[n_snps >= 2, "seqid"].astype(str))
        sensitivity_not_applicable = set(mr.loc[n_snps < 2, "seqid"].astype(str))
    if sensitivity_output_units or cp07.done or cp07.failed:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="sensitivity",
            input_units=sensitivity_input_units,
            output_units=sensitivity_output_units,
            cp=cp07,
            not_applicable_units=sensitivity_not_applicable,
            kind="analysis",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    cp08_regions = _read_checkpoint(cdir, ("_state_08_regions.json",))
    coloc_candidates = _candidate_coloc_seqids(cdir)
    region_base = processed_dir / "coloc" / "regions" / cohort
    region_units: set[str] = set()
    region_output_stats: dict[str, FileStats] = {}
    if region_base.exists():
        for d in sorted(p for p in region_base.iterdir() if p.is_dir()):
            exp = d / "exposure.tsv"
            out = d / "outcome.tsv"
            if exp.exists() and out.exists():
                region_units.add(d.name)
                exp_stat = _read_file_stats(exp, prefer_hg38=True)
                out_stat = _read_file_stats(out, prefer_hg38=True)
                region_output_stats[d.name] = FileStats(
                    seqid=d.name,
                    path=d,
                    n_rows=exp_stat.n_rows + out_stat.n_rows,
                    loci=exp_stat.loci | out_stat.loci,
                    fingerprints={**exp_stat.fingerprints, **out_stat.fingerprints},
                )
    if region_units or cp08_regions.done or cp08_regions.failed:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="coloc_regions",
            input_units=coloc_candidates,
            output_units=region_units,
            output_stats=region_output_stats,
            cp=cp08_regions,
            kind="analysis",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    cp08_sharepro = _read_checkpoint(cdir, ("_state_08_sharepro.json",))
    sharepro = _read_table(processed_dir / "coloc" / "sharepro_results.tsv", cohort=cohort)
    sharepro_units = set(sharepro["seqid"].astype(str)) if "seqid" in sharepro.columns else set()
    if sharepro_units or cp08_sharepro.done or cp08_sharepro.failed:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="sharepro",
            input_units=region_units,
            output_units=sharepro_units,
            cp=cp08_sharepro,
            kind="analysis",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    cp08_abf = _read_checkpoint(cdir, ("_state_08_coloc_abf.rds",))
    coloc_abf = _read_table(processed_dir / "coloc" / f"coloc_abf_{cohort}.tsv")
    coloc_abf_units = set(coloc_abf["seqid"].astype(str)) if "seqid" in coloc_abf.columns else set()
    if coloc_abf_units or cp08_abf.done or cp08_abf.failed:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="coloc_abf",
            input_units=region_units,
            output_units=coloc_abf_units,
            cp=cp08_abf,
            kind="analysis",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    final = _read_table(processed_dir / "final_results.tsv", cohort=cohort)
    final_units = set(final["seqid"].astype(str)) if "seqid" in final.columns else set()
    if final_units:
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="final_results",
            input_units=mr_output_units,
            output_units=final_units,
            kind="analysis",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    gene_summary = _read_table(processed_dir / "gene_summary.tsv")
    if not gene_summary.empty and "gene" in gene_summary.columns:
        final_genes = set(final["gene"].dropna().astype(str)) if "gene" in final.columns else set()
        summary_genes = set(gene_summary["gene"].dropna().astype(str))
        _append_stage(
            bundle,
            ts=ts,
            cohort=cohort,
            stage="gene_summary",
            input_units=final_genes,
            output_units=final_genes & summary_genes,
            kind="analysis",
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )

    return bundle


def report_cohort(cohort: str, processed_dir: Path | None = None) -> tuple[list[dict], list[str]]:
    """
    Compute aggregate yield rows for one cohort.
    Returns (rows, warnings). Each row is a dict with stage metrics.
    """
    bundle = report_cohort_details(cohort, processed_dir=processed_dir)
    return bundle.rows, bundle.warnings


def _print_table(cohort: str, rows: list[dict]) -> None:
    header = f"\n{cohort} yield report"
    col_fmt = "{:<26} {:>7} {:>7} {:>7} {:>9} {:>9} {:>7} {:>7} {:>7} {:>7}"
    print(header)
    print(col_fmt.format("stage", "u_in", "u_out", "u_yld", "row_in", "row_out", "loc_in", "loc_out", "notdone", "done_no"))
    print("-" * 102)
    for r in rows:
        print(col_fmt.format(
            r["stage"],
            r["units_input"],
            r["units_output"],
            f"{r['pct_unit_yield']:.1f}%",
            r["rows_input"],
            r["rows_output"],
            r["loci_input"],
            r["loci_output"],
            r["n_failed"],
            r["units_done_without_output"],
        ))


def _append_dicts(path: Path, rows: list[dict], cols: list[str]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists()
    with path.open("a") as fh:
        if write_header:
            fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(col, "")) for col in cols) + "\n")


def _append_outputs(bundle: ReportBundle, processed_dir: Path) -> None:
    aggregate_cols = [
        "ts", "cohort", "stage", "kind",
        "units_input", "units_output", "pct_unit_yield",
        "rows_input", "rows_output", "pct_row_yield",
        "loci_input", "loci_output", "pct_locus_yield",
        "rows_dropped", "loci_dropped",
        "units_not_done", "units_failed_cp", "units_done_without_output", "units_not_applicable",
        "n_done_cp", "n_input", "n_output", "pct_yield", "n_failed", "pct_failed",
    ]
    unit_cols = [
        "ts", "cohort", "stage", "seqid", "status",
        "input_rows", "output_rows", "input_loci", "output_loci",
        "input_path", "output_path", "failure_reason",
    ]
    dropped_cols = [
        "ts", "cohort", "stage", "seqid", "variant_id", "locus", "input_path", "output_path",
    ]
    _append_dicts(processed_dir / "_yield_report.tsv", bundle.rows, aggregate_cols)
    _append_dicts(processed_dir / "_yield_report_by_unit.tsv", bundle.unit_rows, unit_cols)
    _append_dicts(processed_dir / "_yield_report_dropped_loci.tsv", bundle.dropped_locus_rows, dropped_cols)


def run_report(
    cohorts: list[str],
    strict: bool = False,
    processed_dir: Path | None = None,
    warn_fail_frac: float = WARN_FAIL_FRAC,
    warn_liftover_locus_drop_frac: float = WARN_LIFTOVER_LOCUS_DROP_FRAC,
) -> bool:
    """
    Run yield report for the given cohorts.
    Returns True if any warnings were triggered (useful for --strict logic).
    """
    if processed_dir is None:
        processed_dir = PROCESSED

    any_warn = False
    for cohort in cohorts:
        bundle = report_cohort_details(
            cohort,
            processed_dir=processed_dir,
            warn_fail_frac=warn_fail_frac,
            warn_liftover_locus_drop_frac=warn_liftover_locus_drop_frac,
        )
        if not bundle.rows:
            print(f"\n{cohort}: no protein_index.tsv found — skipping", flush=True)
            continue
        _print_table(cohort, bundle.rows)
        _append_outputs(bundle, processed_dir)
        for warning in bundle.warnings:
            print(warning, file=sys.stderr)
            any_warn = True

    return any_warn


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Stage-aware yield/QC report")
    parser.add_argument("--cohort", default="all",
                        help="Cohort name or 'all' (default: all)")
    parser.add_argument("--strict", action="store_true",
                        help="Exit 1 if any warning threshold exceeded")
    parser.add_argument("--warn-protein-fail-frac", type=float, default=WARN_FAIL_FRAC,
                        help="Warn if unit/protein not-done fraction exceeds this value (default: 0.05)")
    parser.add_argument("--warn-liftover-locus-drop-frac", type=float,
                        default=WARN_LIFTOVER_LOCUS_DROP_FRAC,
                        help="Warn if liftover locus drop fraction exceeds this value (default: 0.0)")
    args = parser.parse_args(argv)

    cohorts = COHORTS if args.cohort == "all" else [args.cohort]
    any_warn = run_report(
        cohorts,
        strict=args.strict,
        warn_fail_frac=args.warn_protein_fail_frac,
        warn_liftover_locus_drop_frac=args.warn_liftover_locus_drop_frac,
    )

    if args.strict and any_warn:
        sys.exit(1)


if __name__ == "__main__":
    main()
