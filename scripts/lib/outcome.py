"""Kim et al. 2025 fibroid GWAS outcome lookup via tabix."""
import logging
from functools import lru_cache
from pathlib import Path
from typing import Optional

import pandas as pd
import pysam

from scripts.lib.paths import KIM_GWAS, OUTCOME_DIR

log = logging.getLogger(__name__)

# Kim GWAS column names (harmonised GWAS-SSF format)
KIM_COLS = [
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error", "effect_allele_frequency", "p_value",
    "rsid", "rs_id", "hm_coordinate_conversion", "hm_code", "variant_id",
]
KIM_N = 434_152
KIM_BUILD = "hg38"


class OutcomeLookup:
    """
    Tabix-backed lookup for the Kim fibroid GWAS.
    Positions must be in GRCh38 (1-based).
    """

    def __init__(self, gwas_path: Path = KIM_GWAS):
        self._tbx = pysam.TabixFile(str(gwas_path))
        self._header = KIM_COLS

    def fetch_region(self, chrom: str, start: int, end: int) -> pd.DataFrame:
        """Fetch all variants in a GRCh38 region (1-based, inclusive)."""
        chrom_str = str(chrom).lstrip("chr")
        rows = []
        try:
            for rec in self._tbx.fetch(chrom_str, start - 1, end):
                rows.append(rec.split("\t"))
        except ValueError:
            # Contig not in index — return empty
            pass
        if not rows:
            return pd.DataFrame(columns=self._header)
        df = pd.DataFrame(rows, columns=self._header)
        df["base_pair_location"] = df["base_pair_location"].astype(int)
        df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
        df["standard_error"] = pd.to_numeric(df["standard_error"], errors="coerce")
        df["effect_allele_frequency"] = pd.to_numeric(df["effect_allele_frequency"], errors="coerce")
        df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
        df["N"] = KIM_N
        return df

    def fetch_snps(self, positions: list[tuple[str, int]]) -> pd.DataFrame:
        """
        Fetch specific (chrom, pos_hg38) pairs.
        Returns rows matching those positions (one row per SNP or none if absent).
        """
        frames = []
        for chrom, pos in positions:
            df = self.fetch_region(chrom, pos, pos)
            frames.append(df)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=self._header)

    def fetch_by_rsid(self, rsids: list[str]) -> pd.DataFrame:
        """
        Fetch variants by rsID by scanning the whole file.
        Slow — only use for small lists where position is unavailable.
        """
        target = set(rsids)
        rows = []
        for rec in self._tbx.fetch():
            parts = rec.split("\t")
            if parts[8] in target:  # rsid column
                rows.append(parts)
        if not rows:
            return pd.DataFrame(columns=self._header)
        return pd.DataFrame(rows, columns=self._header)

    def close(self) -> None:
        self._tbx.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def normalize_outcome_row(row: pd.Series) -> dict:
    """Convert a Kim GWAS row to the harmonised output schema."""
    return {
        "chrom_hg38": str(row["chromosome"]),
        "pos_hg38": int(row["base_pair_location"]),
        "rsid": row["rsid"],
        "EA_out": row["effect_allele"].upper(),
        "OA_out": row["other_allele"].upper(),
        "EAF_out": float(row["effect_allele_frequency"]) if pd.notna(row["effect_allele_frequency"]) else None,
        "beta_out": float(row["beta"]) if pd.notna(row["beta"]) else None,
        "se_out": float(row["standard_error"]) if pd.notna(row["standard_error"]) else None,
        "pval_out": float(row["p_value"]) if pd.notna(row["p_value"]) else None,
        "N_out": KIM_N,
    }
