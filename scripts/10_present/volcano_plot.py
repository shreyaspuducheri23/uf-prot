#!/usr/bin/env python3
"""Generate a volcano plot for the full UKB_female MR landscape.

Usage:
  uv run python scripts/10_present/volcano_plot.py [--out figures/]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

from scripts.lib.paths import FINAL_RESULTS, ROOT


def _require_columns(df: pd.DataFrame, cols: list[str], source: Path) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{source} is missing required columns: {', '.join(missing)}")


def _prepare_ukb_female(df: pd.DataFrame) -> pd.DataFrame:
    _require_columns(df, ["cohort", "seqid", "gene", "beta", "pval", "fdr_q"], FINAL_RESULTS)
    ukb = df[df["cohort"] == "UKB_female"].copy()
    if ukb.empty:
        raise ValueError(f"{FINAL_RESULTS} contains no UKB_female rows")

    for col in ["sharepro_PP_H4", "coloc_abf_PP_H4"]:
        if col not in ukb.columns:
            ukb[col] = np.nan

    ukb["pval"] = pd.to_numeric(ukb["pval"], errors="coerce")
    ukb["beta"] = pd.to_numeric(ukb["beta"], errors="coerce")
    ukb["fdr_q"] = pd.to_numeric(ukb["fdr_q"], errors="coerce")
    ukb["sharepro_PP_H4"] = pd.to_numeric(ukb["sharepro_PP_H4"], errors="coerce")
    ukb["coloc_abf_PP_H4"] = pd.to_numeric(ukb["coloc_abf_PP_H4"], errors="coerce")
    ukb = ukb.dropna(subset=["beta", "pval", "fdr_q"])
    ukb = ukb[ukb["pval"] > 0]
    if ukb.empty:
        raise ValueError("UKB_female rows have no plottable beta/p-value/FDR values")

    ukb["_gene_key"] = ukb["gene"].fillna(ukb["seqid"]).astype(str)
    ukb = ukb.sort_values("pval", na_position="last").drop_duplicates("_gene_key", keep="first")
    ukb["neg_log10_p"] = -np.log10(ukb["pval"].clip(lower=np.nextafter(0, 1)))
    ukb["fdr_sig"] = ukb["fdr_q"] <= 0.05
    ukb["coloc_positive"] = (ukb["sharepro_PP_H4"] >= 0.8) | (ukb["coloc_abf_PP_H4"] >= 0.8)
    return ukb


def make_volcano_plot(out_dir: Path) -> tuple[Path, Path]:
    if not FINAL_RESULTS.exists():
        raise FileNotFoundError(f"Missing {FINAL_RESULTS}; run scripts/09_assemble/assemble.py first")

    df = pd.read_csv(FINAL_RESULTS, sep="\t")
    ukb = _prepare_ukb_female(df)

    fig, ax = plt.subplots(figsize=(9, 7))

    groups = [
        ("Not significant", ~(ukb["fdr_sig"]), "#d7d7d7", 20, 0.75),
        ("FDR < 0.05, beta < 0", ukb["fdr_sig"] & (ukb["beta"] < 0), "#DC3220", 42, 0.95),
        ("FDR < 0.05, beta > 0", ukb["fdr_sig"] & (ukb["beta"] > 0), "#005AB5", 42, 0.95),
    ]
    for label, mask, color, size, alpha in groups:
        subset = ukb[mask]
        if subset.empty:
            continue
        ax.scatter(
            subset["beta"],
            subset["neg_log10_p"],
            s=size,
            c=color,
            alpha=alpha,
            edgecolors="white",
            linewidths=0.4,
            label=label,
        )

    ax.set_xlabel("UKB_female log(OR)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("UKB_female MR results", fontsize=13, pad=12)
    ax.grid(color="#eeeeee", linewidth=0.8)
    ax.legend(frameon=False, loc="best")

    texts = []
    for _, row in ukb[ukb["fdr_sig"]].iterrows():
        texts.append(
            ax.text(
                row["beta"],
                row["neg_log10_p"],
                str(row["_gene_key"]),
                fontsize=8,
                ha="center",
                va="bottom",
            )
        )
    if texts:
        adjust_text(
            texts,
            ax=ax,
            expand_points=(1.2, 1.4),
            expand_text=(1.1, 1.2),
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "volcano_plot.pdf"
    png_path = out_dir / "volcano_plot.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "figures", help="Output directory")
    args = parser.parse_args()
    pdf_path, png_path = make_volcano_plot(args.out)
    print(f"Wrote {pdf_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
