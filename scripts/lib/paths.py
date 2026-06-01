"""Canonical paths — single source of truth for all scripts."""
import os
from pathlib import Path


def _env_path(name: str, default: Path) -> Path:
    raw = os.getenv(name)
    return Path(raw).expanduser().resolve() if raw else default


ROOT = _env_path("LEIO_ROOT", Path(__file__).resolve().parents[2])

# Raw data
RAW = _env_path("LEIO_RAW_DIR", ROOT / "data" / "raw")
ARIC_EA_DIR = RAW / "ARIC" / "EA"
ARIC_AA_DIR = RAW / "ARIC" / "AA"
ARIC_SEQID = RAW / "ARIC" / "seqid.txt"
DECODE_DIR = RAW / "deCODE"
DECODE_ANNOTATED = DECODE_DIR / "assocvariants.annotated.txt.gz"
FENLAND_DIR = RAW / "fenland"
UKBPPP_DIR = RAW / "ukb_ppp"
UKB_FEMALE_DIR = _env_path("LEIO_UKB_FEMALE_DIR",
                            Path("/Volumes/Extreme SSD/ProteoNexus"))
KIM_GWAS = RAW / "kim_fibroid_gwas" / "GCST90461958.h.tsv.gz"
KIM_META = RAW / "kim_fibroid_gwas" / "GCST90461958.h.tsv.gz-meta.yaml"

# Reference data
REF = _env_path("LEIO_REF_DIR", ROOT / "data" / "ref")
CHAIN_HG19_TO_HG38 = REF / "hg19ToHg38.over.chain.gz"
LD_REF_DIR = _env_path("LEIO_LD_REF_DIR", ROOT / "data" / "ld_ref" / "ld_files")
LD_REF_PREFIX = LD_REF_DIR / "data_maf0.01_rs"  # .bed/.bim/.fam
PLINK2 = "plink2"

# Tools
TOOLS = _env_path("LEIO_TOOLS_DIR", ROOT / "tools")
SHAREPRO_DIR = TOOLS / "SharePro_coloc"
SHAREPRO_SCRIPT = SHAREPRO_DIR / "src" / "SharePro" / "sharepro_coloc.py"

# Processed data per cohort
PROCESSED = _env_path("LEIO_PROCESSED_DIR", ROOT / "processed_data")
UKB_FEMALE_CIS_RAW = _env_path("LEIO_UKB_FEMALE_CIS_RAW",
                                PROCESSED / "UKB_female" / "cis_raw_1000kb")
OUTCOME_DIR = PROCESSED / "outcome"

COHORTS = ["ARIC_EA", "deCODE", "UKB_PPP", "Fenland", "UKB_female"]

def cohort_dir(cohort: str) -> Path:
    return PROCESSED / cohort

def raw_cis_sumstats_dir(cohort: str) -> Path:
    return cohort_dir(cohort) / "raw_cis_sumstats"

def raw_cis_sumstats_hg38_dir(cohort: str) -> Path:
    return cohort_dir(cohort) / "raw_cis_sumstats_hg38"

def filtered_cis_pqtls_dir(cohort: str) -> Path:
    return cohort_dir(cohort) / "filtered_cis_pqtls"

def filtered_cis_pqtls_hg38_dir(cohort: str) -> Path:
    return cohort_dir(cohort) / "filtered_cis_pqtls_hg38"

def cis_sumstats_dir(cohort: str) -> Path:
    """Compatibility alias for the filtered MR-ready cis-pQTL product."""
    return filtered_cis_pqtls_dir(cohort)

def cis_sumstats_hg38_dir(cohort: str) -> Path:
    """Compatibility alias for the lifted filtered MR-ready cis-pQTL product."""
    return filtered_cis_pqtls_hg38_dir(cohort)

def instruments_dir(cohort: str) -> Path:
    return cohort_dir(cohort) / "instruments"

def instruments_hg38_dir(cohort: str) -> Path:
    return cohort_dir(cohort) / "instruments_hg38"

def harmonised_dir(cohort: str) -> Path:
    return cohort_dir(cohort) / "harmonised"

def mr_results_path(cohort: str) -> Path:
    return cohort_dir(cohort) / "mr_results.tsv"

def sensitivity_path(cohort: str) -> Path:
    return cohort_dir(cohort) / "sensitivity.tsv"

COLOC_DIR = PROCESSED / "coloc"
COLOC_REGIONS_DIR = COLOC_DIR / "regions"
MR_ALL_COHORTS = PROCESSED / "mr_all_cohorts.tsv"
FINAL_RESULTS = PROCESSED / "final_results.tsv"
GENE_SUMMARY = PROCESSED / "gene_summary.tsv"

LOGS = _env_path("LEIO_LOGS_DIR", ROOT / "logs")

def ensure_dirs(*dirs: Path) -> None:
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
