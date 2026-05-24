"""Shared test fixtures."""
import json
import pytest

PIPELINE_CFG = {
    "_meta": {"version": "1.0"},
    "cis_extract": {
        "window_kb": 500,
        "pval_gw": 5e-8,
        "maf_min": 0.01,
        "palindrome_maf_max": 0.42,
    },
    "clump": {"window_kb": 1000, "r2": 0.001, "p1": 5e-8},
    "fstat": {"weak_threshold": 10.0},
    "harmonise": {"maf_proxy_max": 0.42, "proxy_r2_min": 0.8},
    "outcome": {"kim_N": 434152},
    "mhc": {"hg19": [25000000, 34000000], "hg38": [28500000, 33500000]},
    "cohorts": {
        "ARIC_EA":    {"N": None, "build": "hg19"},
        "Fenland":    {"N": 10708, "build": "hg19"},
        "deCODE":     {"N_default": 35000, "build": "hg38"},
        "UKB_PPP":    {"N": None, "build": "hg19"},
        "UKB_female": {"N": None, "build": "hg19"},
    },
}


@pytest.fixture
def pipeline_cfg():
    import copy
    return copy.deepcopy(PIPELINE_CFG)


@pytest.fixture
def tmp_config_file(tmp_path):
    path = tmp_path / "pipeline.json"
    path.write_text(json.dumps(PIPELINE_CFG))
    return path
