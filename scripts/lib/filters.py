"""Standard variant filters applied during cis-pQTL extraction."""
import pandas as pd

# MHC region boundaries (inclusive)
MHC_HG19 = ("6", 25_000_000, 34_000_000)
MHC_HG38 = ("6", 28_500_000, 33_500_000)


def _mhc_bounds(build: str) -> tuple[str, int, int]:
    if build == "hg19":
        return MHC_HG19
    if build == "hg38":
        return MHC_HG38
    raise ValueError(f"Unknown build: {build!r}. Expected 'hg19' or 'hg38'.")


def maf_above(df: pd.DataFrame, threshold: float = 0.01,
              eaf_col: str = "EAF") -> pd.DataFrame:
    """Keep rows where MAF (min(EAF, 1-EAF)) >= threshold."""
    eaf = df[eaf_col].astype(float)
    maf = eaf.where(eaf <= 0.5, 1 - eaf)
    return df[maf >= threshold].copy()


def gw_significant(df: pd.DataFrame, p: float = 5e-8,
                   pval_col: str = "pval") -> pd.DataFrame:
    return df[df[pval_col].astype(float) < p].copy()


def exclude_mhc(df: pd.DataFrame, build: str,
                chrom_col: str = "chrom", pos_col: str = "pos") -> pd.DataFrame:
    chrom, start, end = _mhc_bounds(build)
    in_mhc = (df[chrom_col].astype(str) == chrom) & \
              (df[pos_col].astype(int) >= start) & \
              (df[pos_col].astype(int) <= end)
    return df[~in_mhc].copy()


def cis_window(df: pd.DataFrame, tss: int, gene_chrom: str,
               build: str, kb: int = 500,
               chrom_col: str = "chrom", pos_col: str = "pos") -> pd.DataFrame:
    """Restrict to variants within ±kb of TSS on the same chromosome."""
    flank = kb * 1_000
    same_chrom = df[chrom_col].astype(str) == str(gene_chrom)
    near_tss = (df[pos_col].astype(int) - tss).abs() <= flank
    return df[same_chrom & near_tss].copy()


def drop_ambig_palindromes(df: pd.DataFrame, maf_threshold: float = 0.42,
                           ea_col: str = "EA", oa_col: str = "OA",
                           eaf_col: str = "EAF") -> pd.DataFrame:
    """
    Drop A/T or C/G variants with MAF > maf_threshold (ambiguous strand).
    Variants with MAF <= maf_threshold are kept (strand can be inferred from freq).
    """
    ea = df[ea_col].astype(str).str.upper()
    oa = df[oa_col].astype(str).str.upper()
    is_palindrome = ((ea == "A") & (oa == "T")) | ((ea == "T") & (oa == "A")) | \
                    ((ea == "C") & (oa == "G")) | ((ea == "G") & (oa == "C"))
    eaf = df[eaf_col].astype(float)
    maf = eaf.where(eaf <= 0.5, 1 - eaf)
    ambig = is_palindrome & (maf > maf_threshold)
    return df[~ambig].copy()
