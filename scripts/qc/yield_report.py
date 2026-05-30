"""
Yield-funnel reporter: reads filesystem state and prints a per-stage table.

Usage:
    uv run python scripts/qc/yield_report.py [--cohort COHORT|all] [--strict]
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

from scripts.lib.paths import COHORTS, PROCESSED

# (stage_name, output_subfolder, checkpoint_filename)
STAGE_DEFS: list[tuple[str, str, str]] = [
    ("filtered_cis_pqtls", "filtered_cis_pqtls", "_state_02.json"),
    ("instruments",        "instruments",        "_state_03.json"),
    ("instruments_hg38",   "instruments_hg38",   "_state_04.json"),
    ("harmonised",         "harmonised",         "_state_05.json"),
]

WARN_FAIL_FRAC = 0.05


def _count_tsv(d: Path) -> int:
    if not d.exists():
        return 0
    return sum(1 for f in d.iterdir() if f.suffix == ".tsv")


def _read_checkpoint(cp_path: Path) -> tuple[int, int]:
    """Return (n_done, n_failed) from a checkpoint JSON."""
    if not cp_path.exists():
        return 0, 0
    try:
        data = json.loads(cp_path.read_text())
        n_done = len(data.get("done", []))
        status = data.get("status", {})
        n_failed = sum(
            1 for v in status.values()
            if isinstance(v, dict) and v.get("state") == "failed"
        )
        return n_done, n_failed
    except (json.JSONDecodeError, OSError):
        return 0, 0


def report_cohort(cohort: str, processed_dir: Path | None = None) -> tuple[list[dict], list[str]]:
    """
    Compute yield rows for one cohort.
    Returns (rows, warnings). Each row is a dict with stage metrics.
    """
    if processed_dir is None:
        processed_dir = PROCESSED

    cdir = processed_dir / cohort
    index_path = cdir / "protein_index.tsv"

    if not index_path.exists():
        return [], []

    try:
        lines = index_path.read_text().splitlines()
        n_index = max(0, len([l for l in lines if l.strip()]) - 1)
    except OSError:
        n_index = 0

    rows: list[dict] = []
    warnings: list[str] = []
    n_input = n_index
    ts = datetime.now(UTC).isoformat(timespec="seconds")

    for stage_name, subfolder, state_file in STAGE_DEFS:
        cp_path = cdir / state_file
        if not cp_path.exists():
            break

        out_dir = cdir / subfolder
        n_output = _count_tsv(out_dir)
        n_done, _ = _read_checkpoint(cp_path)
        # n_failed = anything not marked done: covers both explicit mark_failed entries
        # and proteins silently abandoned by pre-fix runs (neither done nor failed in cp)
        n_failed = max(0, n_input - n_done)

        pct_yield = (n_output / n_input * 100) if n_input > 0 else 0.0
        pct_failed = (n_failed / n_input * 100) if n_input > 0 else 0.0

        row = {
            "ts": ts,
            "cohort": cohort,
            "stage": stage_name,
            "n_input": n_input,
            "n_output": n_output,
            "pct_yield": round(pct_yield, 1),
            "n_failed": n_failed,
            "pct_failed": round(pct_failed, 1),
            "n_done_cp": n_done,
        }
        rows.append(row)

        if n_input > 0 and n_failed / n_input > WARN_FAIL_FRAC:
            warnings.append(
                f"[WARN] {cohort}/{stage_name}: {n_failed}/{n_input} proteins not done "
                f"({pct_failed:.1f}%) exceeds 5% threshold"
            )

        n_input = n_output

    return rows, warnings


def _print_table(cohort: str, rows: list[dict]) -> None:
    header = f"\n{cohort} yield report"
    col_fmt = "{:<20} {:>8} {:>9} {:>10} {:>9} {:>10}"
    print(header)
    print(col_fmt.format("stage", "n_input", "n_output", "pct_yield", "n_failed", "pct_failed"))
    print("-" * 72)
    for r in rows:
        print(col_fmt.format(
            r["stage"],
            r["n_input"],
            r["n_output"],
            f"{r['pct_yield']:.1f}%",
            r["n_failed"],
            f"{r['pct_failed']:.1f}%",
        ))


def _append_tsv(rows: list[dict], processed_dir: Path) -> None:
    tsv_path = processed_dir / "_yield_report.tsv"
    cols = ["ts", "cohort", "stage", "n_input", "n_output", "pct_yield", "n_failed", "pct_failed"]
    write_header = not tsv_path.exists()
    with tsv_path.open("a") as fh:
        if write_header:
            fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


def run_report(
    cohorts: list[str],
    strict: bool = False,
    processed_dir: Path | None = None,
) -> bool:
    """
    Run yield report for the given cohorts.
    Returns True if any warnings were triggered (useful for --strict logic).
    """
    if processed_dir is None:
        processed_dir = PROCESSED

    any_warn = False
    for cohort in cohorts:
        rows, warnings = report_cohort(cohort, processed_dir=processed_dir)
        if not rows:
            print(f"\n{cohort}: no protein_index.tsv found — skipping", flush=True)
            continue
        _print_table(cohort, rows)
        _append_tsv(rows, processed_dir)
        for w in warnings:
            print(w, file=sys.stderr)
            any_warn = True

    return any_warn


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Yield-funnel report for pipeline QC")
    parser.add_argument("--cohort", default="all",
                        help="Cohort name or 'all' (default: all)")
    parser.add_argument("--strict", action="store_true",
                        help="Exit 1 if any warning threshold exceeded")
    args = parser.parse_args(argv)

    if args.cohort == "all":
        cohorts = COHORTS
    else:
        cohorts = [args.cohort]

    any_warn = run_report(cohorts, strict=args.strict)

    if args.strict and any_warn:
        sys.exit(1)


if __name__ == "__main__":
    main()
