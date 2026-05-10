"""hg19 ↔ hg38 coordinate liftover using pyliftover."""
import logging
from pathlib import Path

import pandas as pd
from pyliftover import LiftOver

from scripts.lib.paths import CHAIN_HG19_TO_HG38

log = logging.getLogger(__name__)

_cache: dict[str, LiftOver] = {}


def _get_lo(chain_path: Path) -> LiftOver:
    key = str(chain_path)
    if key not in _cache:
        _cache[key] = LiftOver(str(chain_path))
    return _cache[key]


def lift_position(chrom: str, pos: int, chain_path: Path = CHAIN_HG19_TO_HG38
                  ) -> tuple[str, int] | None:
    """
    Lift a single 1-based position. Returns (chrom, pos_hg38) or None on failure.
    Chromosome must be without 'chr' prefix (e.g. '6').
    """
    lo = _get_lo(chain_path)
    result = lo.convert_coordinate(f"chr{chrom}", pos - 1)  # pyliftover is 0-based
    if not result:
        return None
    lifted_chrom = result[0][0].lstrip("chr")
    lifted_pos = result[0][1] + 1  # back to 1-based
    return lifted_chrom, lifted_pos


def lift_table(df: pd.DataFrame,
               chrom_col: str = "chrom", pos_col: str = "pos",
               out_chrom_col: str = "chrom_hg38", out_pos_col: str = "pos_hg38",
               chain_path: Path = CHAIN_HG19_TO_HG38,
               ) -> pd.DataFrame:
    """
    Add hg38 columns to a DataFrame. Rows that fail to lift or change chromosome are dropped.
    Returns the filtered DataFrame with new columns added.
    """
    df = df.copy()
    hg38_chroms: list[str | None] = []
    hg38_pos: list[int | None] = []
    n_no_result = 0
    n_chrom_change = 0

    for _, row in df.iterrows():
        result = lift_position(str(row[chrom_col]), int(row[pos_col]), chain_path)
        if result is None:
            n_no_result += 1
            hg38_chroms.append(None)
            hg38_pos.append(None)
        elif result[0] != str(row[chrom_col]):
            n_chrom_change += 1
            hg38_chroms.append(None)
            hg38_pos.append(None)
        else:
            hg38_chroms.append(result[0])
            hg38_pos.append(result[1])

    df[out_chrom_col] = hg38_chroms
    df[out_pos_col] = hg38_pos

    n_total = len(df)
    if n_no_result and n_total:
        pct = 100 * n_no_result / n_total
        log.warning(f"Liftover: {n_no_result}/{n_total} ({pct:.1f}%) positions failed to lift — dropped")
    if n_chrom_change and n_total:
        pct = 100 * n_chrom_change / n_total
        log.warning(f"Liftover: {n_chrom_change}/{n_total} ({pct:.1f}%) positions changed chromosome — dropped")

    return df.dropna(subset=[out_pos_col]).copy()
