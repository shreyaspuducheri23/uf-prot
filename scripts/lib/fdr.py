"""Benjamini-Hochberg FDR correction."""
import numpy as np
import pandas as pd


def bh_fdr(pvals: pd.Series, alpha: float = 0.05) -> pd.Series:
    """
    Benjamini-Hochberg FDR. Returns adjusted p-values (q-values).
    NaN inputs are preserved as NaN.
    """
    n = len(pvals)
    if n == 0:
        return pvals.copy()

    nan_mask = pvals.isna()
    valid_idx = pvals[~nan_mask].index
    valid_p = pvals[~nan_mask].values.astype(float)

    order = np.argsort(valid_p)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(valid_p) + 1)

    adjusted = np.minimum(1.0, valid_p * n / ranks)
    # Enforce monotonicity (from largest to smallest rank)
    for i in range(len(adjusted) - 2, -1, -1):
        adjusted[order[i]] = min(adjusted[order[i]], adjusted[order[i + 1]])

    result = pvals.copy().astype(float)
    result[valid_idx] = adjusted
    return result


def add_fdr(df: pd.DataFrame, pval_col: str = "pval",
            out_col: str = "fdr_q", alpha: float = 0.05) -> pd.DataFrame:
    df = df.copy()
    df[out_col] = bh_fdr(df[pval_col], alpha=alpha)
    df["fdr_pass"] = df[out_col] <= alpha
    return df
