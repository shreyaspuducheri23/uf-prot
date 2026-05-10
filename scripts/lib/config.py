"""Pipeline configuration: load config/pipeline.json and provide typed access."""
import argparse
import json
from functools import lru_cache
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = _ROOT / "config" / "pipeline.json"

_REQUIRED: dict[str, set] = {
    "cis_extract": {"window_kb", "pval_gw", "maf_min", "palindrome_maf_max"},
    "clump":       {"window_kb", "r2", "p1"},
    "fstat":       {"weak_threshold"},
    "harmonise":   {"maf_proxy_max", "proxy_r2_min"},
    "outcome":     {"kim_N"},
    "mhc":         {"hg19", "hg38"},
    "cohorts":     set(),
}


def _validate(cfg: dict) -> None:
    for section, required_keys in _REQUIRED.items():
        if section not in cfg:
            raise ValueError(f"Config missing required section: {section!r}")
        for key in required_keys:
            if key not in cfg[section]:
                raise ValueError(f"Config [{section}] missing required key: {key!r}")
    ce = cfg["cis_extract"]
    if not (0 < ce["pval_gw"] < 1):
        raise ValueError(f"cis_extract.pval_gw must be in (0, 1), got {ce['pval_gw']}")
    if ce["window_kb"] <= 0:
        raise ValueError(f"cis_extract.window_kb must be > 0, got {ce['window_kb']}")
    if not (0 < ce["maf_min"] < 1):
        raise ValueError(f"cis_extract.maf_min must be in (0, 1), got {ce['maf_min']}")


@lru_cache(maxsize=8)
def load_config(path: str | None = None) -> dict:
    """Load and validate pipeline.json. Results are cached per resolved path."""
    resolved = Path(path).resolve() if path else DEFAULT_CONFIG_PATH
    if not resolved.exists():
        example = _ROOT / "config" / "pipeline.example.json"
        raise FileNotFoundError(
            f"Pipeline config not found: {resolved}\n"
            f"Copy {example} → {DEFAULT_CONFIG_PATH} to get started."
        )
    with open(resolved) as fh:
        cfg = json.load(fh)
    _validate(cfg)
    return cfg


def get_section(cfg: dict, name: str) -> dict:
    if name not in cfg:
        raise KeyError(
            f"Config section {name!r} not found. Available: {sorted(k for k in cfg if not k.startswith('_'))}"
        )
    return cfg[name]


def add_config_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--config", default=None, metavar="PATH",
        help=f"Path to pipeline.json (default: {DEFAULT_CONFIG_PATH})",
    )
