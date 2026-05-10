"""F-statistic computation for instrument strength."""
import pandas as pd
import numpy as np


def compute_fstat(df: pd.DataFrame,
                  beta_col: str = "beta", se_col: str = "se") -> pd.Series:
    """Return F = (beta/se)^2 for each row."""
    beta = df[beta_col].astype(float)
    se = df[se_col].astype(float)
    return (beta / se) ** 2


def add_fstat(df: pd.DataFrame,
              beta_col: str = "beta", se_col: str = "se",
              out_col: str = "F_stat") -> pd.DataFrame:
    df = df.copy()
    df[out_col] = compute_fstat(df, beta_col, se_col)
    return df


WEAK_INSTRUMENT_THRESHOLD = 10.0
