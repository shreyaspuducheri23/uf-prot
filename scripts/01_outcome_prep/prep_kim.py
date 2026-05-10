#!/usr/bin/env python3
"""
01_outcome_prep/prep_kim.py
Parse Kim 2025 fibroid GWAS metadata and smoke-test the tabix outcome lookup.

Usage:
  python scripts/01_outcome_prep/prep_kim.py
"""
import json
import sys

from scripts.lib.logging import setup_logger, RunManifest
from scripts.lib.outcome import OutcomeLookup, KIM_N, KIM_BUILD
from scripts.lib.paths import KIM_META, OUTCOME_DIR, ensure_dirs

log = setup_logger("01_outcome_prep")


def parse_meta(meta_path) -> dict:
    import yaml  # optional dep; use basic parsing if absent
    try:
        import yaml
        with open(meta_path) as fh:
            return yaml.safe_load(fh)
    except ImportError:
        meta = {"N": KIM_N, "build": KIM_BUILD, "source": str(meta_path)}
        log.warning("PyYAML not installed — using hardcoded Kim meta values")
        return meta


def smoke_test(lookup: OutcomeLookup) -> None:
    # Test with a known fibroid-associated locus (chr22 q area near WT1 — example)
    # Use chr1 common variant as a connectivity check
    df = lookup.fetch_region("1", 1_000_000, 1_005_000)
    log.info(f"Smoke test fetch chr1:1M-1.005M → {len(df)} rows")
    if df.empty:
        log.warning("Smoke test returned 0 rows — check tabix index or file path")
    else:
        log.info(f"First row: {df.iloc[0].to_dict()}")


def main() -> None:
    ensure_dirs(OUTCOME_DIR)

    with RunManifest("01_outcome_prep/prep_kim.py") as manifest:
        meta = parse_meta(KIM_META)
        meta_out = {
            "N": KIM_N,
            "build": KIM_BUILD,
            "gwas_id": "GCST90461958",
            "trait": "uterine_fibroids",
            "ancestry": "EUR",
            "raw_meta": meta,
        }
        out_path = OUTCOME_DIR / "kim_meta.json"
        out_path.write_text(json.dumps(meta_out, indent=2, default=str))
        log.info(f"Wrote {out_path}")

        with OutcomeLookup() as lookup:
            smoke_test(lookup)

        manifest.n_units = 1
        log.info("Done.")


if __name__ == "__main__":
    main()
