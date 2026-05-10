"""Read/write normalized summary-statistic TSVs."""
from pathlib import Path
import pandas as pd


DTYPE_MAP = {
    "chrom": str,
    "rsid": str,
    "EA": str,
    "OA": str,
    "seqid": str,
    "gene": str,
    "uniprot": str,
    "build": str,
}


def read_norm(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=DTYPE_MAP)


def write_norm(df: pd.DataFrame, path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, sep="\t", index=False)


def read_instruments(path: str | Path) -> pd.DataFrame:
    extra = {
        "chrom_hg38": str,
        "chrom_hg19": str,
        "proxy_rsid": str,
    }
    return pd.read_csv(path, sep="\t", dtype={**DTYPE_MAP, **extra})
