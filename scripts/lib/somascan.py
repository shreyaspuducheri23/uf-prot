"""SomaScan platform utilities shared across cohorts."""
from functools import lru_cache
from pathlib import Path
import pandas as pd

_MENU = Path(__file__).parents[2] / "data/raw/somascan/analyte_menu.tsv"


@lru_cache(maxsize=1)
def load_seqid_map() -> dict[str, str]:
    """Return {seqid_key: gene_symbol} from the official SomaScan.db-derived menu.

    seqid_key format: '14157_21' (underscore-separated, no SeqId_ prefix).
    """
    df = pd.read_csv(_MENU, sep="\t", dtype=str)
    return dict(zip(df["seqid_key"], df["gene_symbol"]))


def parse_fenland_seqid(filename: str) -> str | None:
    """Extract the SeqId key from a Fenland filename without touching the display-name.

    SomaLogic filenames end with _{4-5digit SeqId}_{2digit version}[.txt.gz].
    Splits from the right so multi-word display names (14_3_3E, ADAMTS_4) are irrelevant.

    Returns '14157_21' style key, or None if the filename doesn't match the pattern.
    """
    stem = filename.removesuffix(".txt.gz")
    parts = stem.rsplit("_", maxsplit=2)
    if len(parts) == 3:
        _, major, minor = parts
        if major.isdigit() and len(major) in (4, 5) and minor.isdigit() and len(minor) == 2:
            return f"{major}_{minor}"
    return None
