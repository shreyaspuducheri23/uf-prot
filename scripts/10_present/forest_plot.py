#!/usr/bin/env python3
"""Generate a forest plot for UKB_female discoveries and cohort estimates.

Usage:
  uv run python scripts/10_present/forest_plot.py [--out figures/]
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from scripts.lib.paths import GENE_SUMMARY, ROOT

DISPLAY_COHORTS = ["UKB_PPP", "ARIC_EA", "deCODE", "Fenland"]
REPLICATION_COHORTS = {"ARIC_EA", "deCODE", "Fenland"}


def _require_columns(df: pd.DataFrame, cols: list[str], source: Path) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{source} is missing required columns: {', '.join(missing)}")


def _is_valid_ci(or_val: object, lo: object, hi: object) -> bool:
    vals = pd.to_numeric(pd.Series([or_val, lo, hi]), errors="coerce")
    return bool(vals.notna().all() and (vals > 0).all())


def _fmt_p(value: object) -> str:
    p = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(p):
        return "NA"
    return f"{p:.1e}" if p < 0.001 else f"{p:.3f}"


def _fmt_or(or_val: float, lo: float, hi: float) -> str:
    return f"{or_val:.2f} ({lo:.2f}-{hi:.2f})"


def _build_plot_rows(df: pd.DataFrame) -> list[dict]:
    rows: list[dict] = []
    df = df.sort_values("primary_pval", na_position="last").reset_index(drop=True)

    for _, gene_row in df.iterrows():
        gene = str(gene_row["gene"])
        group_start = len(rows)

        if _is_valid_ci(
            gene_row["primary_OR"],
            gene_row["primary_OR_lo95"],
            gene_row["primary_OR_hi95"],
        ):
            rows.append(
                {
                    "gene": gene,
                    "cohort": "UKB_female",
                    "or": float(gene_row["primary_OR"]),
                    "lo": float(gene_row["primary_OR_lo95"]),
                    "hi": float(gene_row["primary_OR_hi95"]),
                    "p": gene_row["primary_pval"],
                    "kind": "primary",
                    "show_gene": True,
                }
            )

        for cohort in DISPLAY_COHORTS:
            or_col = f"OR_{cohort}"
            lo_col = f"OR_lo95_{cohort}"
            hi_col = f"OR_hi95_{cohort}"
            p_col = f"pval_{cohort}"
            if _is_valid_ci(gene_row[or_col], gene_row[lo_col], gene_row[hi_col]):
                kind = "replication" if cohort in REPLICATION_COHORTS else "reference"
                rows.append(
                    {
                        "gene": gene,
                        "cohort": cohort,
                        "or": float(gene_row[or_col]),
                        "lo": float(gene_row[lo_col]),
                        "hi": float(gene_row[hi_col]),
                        "p": gene_row[p_col],
                        "kind": kind,
                        "show_gene": len(rows) == group_start,
                    }
                )

        if _is_valid_ci(
            gene_row["replication_meta_OR"],
            gene_row["replication_meta_OR_lo95"],
            gene_row["replication_meta_OR_hi95"],
        ):
            rows.append(
                {
                    "gene": gene,
                    "cohort": "Replication meta",
                    "or": float(gene_row["replication_meta_OR"]),
                    "lo": float(gene_row["replication_meta_OR_lo95"]),
                    "hi": float(gene_row["replication_meta_OR_hi95"]),
                    "p": gene_row["replication_meta_pval"],
                    "kind": "meta",
                    "show_gene": len(rows) == group_start,
                }
            )

        for row in rows[group_start + 1 :]:
            row["show_gene"] = False
        if rows[group_start:]:
            rows[group_start]["group_start"] = True
            rows[-1]["group_end"] = True

    return rows


def make_forest_plot(out_dir: Path) -> tuple[Path, Path]:
    if not GENE_SUMMARY.exists():
        raise FileNotFoundError(f"Missing {GENE_SUMMARY}; run scripts/09_assemble/cross_cohort.py first")

    df = pd.read_csv(GENE_SUMMARY, sep="\t")
    required = [
        "gene",
        "primary_OR",
        "primary_OR_lo95",
        "primary_OR_hi95",
        "primary_pval",
        "primary_fdr_q",
        "replication_meta_OR",
        "replication_meta_OR_lo95",
        "replication_meta_OR_hi95",
        "replication_meta_pval",
    ]
    for cohort in DISPLAY_COHORTS:
        required.extend([f"OR_{cohort}", f"OR_lo95_{cohort}", f"OR_hi95_{cohort}", f"pval_{cohort}"])
    _require_columns(df, required, GENE_SUMMARY)

    df["primary_fdr_q"] = pd.to_numeric(df["primary_fdr_q"], errors="coerce")
    df = df[df["primary_fdr_q"] <= 0.05].copy()
    if df.empty:
        raise ValueError(f"{GENE_SUMMARY} contains no UKB_female FDR-passing genes")

    rows = _build_plot_rows(df)
    if not rows:
        raise ValueError(f"{GENE_SUMMARY} contains no plottable OR confidence intervals")

    n_rows = len(rows)
    height = max(5.0, 0.32 * n_rows + 1.6)
    fig = plt.figure(figsize=(11.5, height))
    grid = fig.add_gridspec(1, 4, width_ratios=[1.2, 1.7, 4.2, 2.3], wspace=0.03)
    ax_gene = fig.add_subplot(grid[0, 0])
    ax_cohort = fig.add_subplot(grid[0, 1], sharey=ax_gene)
    ax_plot = fig.add_subplot(grid[0, 2], sharey=ax_gene)
    ax_stats = fig.add_subplot(grid[0, 3], sharey=ax_gene)

    y_positions = list(range(n_rows, 0, -1))
    lo_vals = [r["lo"] for r in rows]
    hi_vals = [r["hi"] for r in rows]
    x_min = min(lo_vals) / 1.25
    x_max = max(hi_vals) * 1.25
    x_min = max(x_min, 0.01)

    for i, (row, y) in enumerate(zip(rows, y_positions)):
        if row.get("show_gene"):
            ax_gene.text(0.98, y, row["gene"], ha="right", va="center", fontsize=9, fontweight="bold")
        ax_cohort.text(0.02, y, row["cohort"], ha="left", va="center", fontsize=9)
        ax_stats.text(
            0.02,
            y,
            f"{_fmt_or(row['or'], row['lo'], row['hi'])}    {_fmt_p(row['p'])}",
            ha="left",
            va="center",
            fontsize=9,
        )

        if row["kind"] == "primary":
            marker = "D"
            color = "#c43c2f"
            face = color
            size = 7
            lw = 1.8
        elif row["kind"] == "meta":
            marker = "D"
            color = "#222222"
            face = "white"
            size = 7
            lw = 2.0
        elif row["kind"] == "reference":
            marker = "s"
            color = "#9a9a9a"
            face = "white"
            size = 5
            lw = 1.1
        else:
            marker = "o"
            color = "#6f6f6f"
            face = color
            size = 5
            lw = 1.2

        ax_plot.errorbar(
            row["or"],
            y,
            xerr=[[row["or"] - row["lo"]], [row["hi"] - row["or"]]],
            fmt=marker,
            markersize=size,
            color=color,
            markerfacecolor=face,
            markeredgecolor=color,
            elinewidth=lw,
            capsize=0,
            zorder=3,
        )

        if row.get("group_end") and i < n_rows - 1:
            sep_y = y - 0.5
            for ax in (ax_gene, ax_cohort, ax_plot, ax_stats):
                ax.axhline(sep_y, color="#dddddd", linewidth=0.7, zorder=1)

    ax_plot.axvline(1.0, color="#333333", linestyle="--", linewidth=1.0, zorder=2)
    ax_plot.set_xscale("log")
    ax_plot.set_xlim(x_min, x_max)
    ax_plot.set_xlabel("Odds ratio (log scale)")
    ax_plot.set_yticks([])
    ax_plot.grid(axis="x", color="#eeeeee", linewidth=0.8)

    for ax in (ax_gene, ax_cohort, ax_stats):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, n_rows + 1)
        ax.axis("off")

    ax_gene.text(0.98, n_rows + 0.6, "Gene", ha="right", va="bottom", fontsize=10, fontweight="bold")
    ax_cohort.text(0.02, n_rows + 0.6, "Cohort", ha="left", va="bottom", fontsize=10, fontweight="bold")
    ax_stats.text(0.02, n_rows + 0.6, "OR (95% CI)    p", ha="left", va="bottom", fontsize=10, fontweight="bold")
    ax_plot.set_ylim(0, n_rows + 1)
    ax_plot.set_title("UKB_female discoveries with cohort estimates", fontsize=12, pad=14)

    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "forest_plot.pdf"
    png_path = out_dir / "forest_plot.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "figures", help="Output directory")
    args = parser.parse_args()
    pdf_path, png_path = make_forest_plot(args.out)
    print(f"Wrote {pdf_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
