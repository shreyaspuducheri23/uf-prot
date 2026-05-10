"""Canonical column names for normalized summary-statistic tables."""
import logging
from dataclasses import dataclass
from typing import Optional
import pandas as pd

log = logging.getLogger(__name__)

# Normalized column names used in all intermediate TSVs
NORM_COLS = [
    "seqid",        # protein identifier (SeqId string)
    "gene",         # gene symbol
    "uniprot",      # UniProt ID
    "chrom",        # chromosome (str, no "chr" prefix)
    "pos",          # position in native build (int)
    "rsid",         # rs identifier or "." if absent
    "EA",           # effect allele
    "OA",           # other allele
    "EAF",          # effect allele frequency
    "beta",         # effect estimate
    "se",           # standard error
    "pval",         # p-value
    "N",            # sample size
    "build",        # "hg19" or "hg38"
]

# Lifted instrument columns (step 04 output)
INSTRUMENT_LIFTED_COLS = NORM_COLS + [
    "chrom_hg38",   # after liftover (same as chrom for hg38 cohorts)
    "pos_hg38",     # position in GRCh38
    "F_stat",       # (beta/se)^2
]

# Harmonised output columns (step 05)
HARMONISED_COLS = [
    "seqid", "gene", "uniprot",
    "chrom_hg19", "pos_hg19", "chrom_hg38", "pos_hg38", "rsid",
    "EA_exp", "OA_exp", "EAF_exp", "beta_exp", "se_exp", "pval_exp", "N_exp", "F_stat",
    "EA_out", "OA_out", "EAF_out", "beta_out", "se_out", "pval_out", "N_out",
    "proxy_rsid",   # non-empty if a proxy SNP was used
    "proxy_r2",
]

CHROM_ORDER = [str(c) for c in range(1, 23)] + ["X"]


@dataclass
class ProteinMeta:
    seqid: str
    gene: str
    uniprot: str
    chrom: str        # chromosome of protein-coding gene
    tss: int          # TSS in native build
    build: str        # "hg19" or "hg38"
    source_cohort: str


_VALID_CHROMS = {str(c) for c in range(1, 23)} | {"X"}


def validate_norm_df(df: pd.DataFrame, *, where: str = "") -> None:
    """
    Assert invariants that must hold for every normalized sumstats DataFrame.
    Raises ValueError with a descriptive message if any invariant is violated.
    'where' is included in error messages to identify the call site.
    """
    loc = f" [{where}]" if where else ""

    missing = [c for c in NORM_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"validate_norm_df{loc}: missing columns: {missing}")

    # chrom: string type, values in valid set (no "chr" prefix)
    chrom_vals = df["chrom"].dropna().unique()
    bad_chrom = [c for c in chrom_vals if str(c) not in _VALID_CHROMS]
    if bad_chrom:
        raise ValueError(
            f"validate_norm_df{loc}: invalid chrom values (expected '1'–'22', 'X', no 'chr' prefix): {bad_chrom[:5]}"
        )

    # pos: must be int64 (not nullable Int64)
    pos_dtype = df["pos"].dtype
    if hasattr(pos_dtype, "numpy_dtype"):
        # pandas ExtensionDtype (e.g. Int64) — not acceptable
        raise ValueError(
            f"validate_norm_df{loc}: 'pos' column has nullable dtype {pos_dtype!r}; use int64 instead"
        )
    if pd.api.types.pandas_dtype(pos_dtype) != pd.api.types.pandas_dtype("int64"):
        raise ValueError(
            f"validate_norm_df{loc}: 'pos' dtype is {pos_dtype!r}, expected int64"
        )

    # EA and OA: uppercase non-empty strings
    for col in ("EA", "OA"):
        bad = df[col].dropna()
        if (bad == "").any():
            raise ValueError(f"validate_norm_df{loc}: '{col}' contains empty strings")
        if (bad != bad.str.upper()).any():
            raise ValueError(f"validate_norm_df{loc}: '{col}' contains lowercase alleles")

    # Numeric columns: must be float64, no infinities
    for col in ("pval", "beta", "se", "EAF"):
        if col in df.columns:
            s = df[col].dropna()
            if not pd.api.types.is_float_dtype(s):
                log.warning(f"validate_norm_df{loc}: '{col}' is not float dtype ({s.dtype})")
            n_inf = s.isin([float("inf"), float("-inf")]).sum()
            if n_inf:
                log.warning(f"validate_norm_df{loc}: '{col}' has {n_inf} infinite values")


def read_norm_tsv(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"chrom": str, "rsid": str})


def write_norm_tsv(df: pd.DataFrame, path) -> None:
    from pathlib import Path
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
