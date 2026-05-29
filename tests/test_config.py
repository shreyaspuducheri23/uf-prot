"""Tests for scripts.lib.config."""
import copy
import json
import pytest

from scripts.lib.config import (
    load_config, get_section, get_cohort_config, get_cohort_build, _validate
)
from tests.conftest import PIPELINE_CFG


@pytest.fixture(autouse=True)
def clear_lru_cache():
    load_config.cache_clear()
    yield
    load_config.cache_clear()


class TestLoadConfig:
    def test_happy_path(self, tmp_config_file):
        cfg = load_config(str(tmp_config_file))
        assert cfg["cis_extract"]["window_kb"] == 500
        assert cfg["outcome"]["kim_N"] == 434152

    def test_missing_file_error_mentions_example(self, tmp_path):
        missing = str(tmp_path / "nonexistent.json")
        with pytest.raises(FileNotFoundError, match="pipeline.example.json"):
            load_config(missing)

    def test_missing_section_raises(self, tmp_path):
        path = tmp_path / "bad.json"
        path.write_text(json.dumps({"_meta": {"version": "1.0"}}))
        with pytest.raises(ValueError, match="missing required section"):
            load_config(str(path))

    def test_missing_key_raises(self, tmp_path, pipeline_cfg):
        del pipeline_cfg["cis_extract"]["window_kb"]
        path = tmp_path / "bad.json"
        path.write_text(json.dumps(pipeline_cfg))
        with pytest.raises(ValueError, match="missing required key"):
            load_config(str(path))

    def test_lru_cache_same_path_returns_same_object(self, tmp_config_file):
        cfg1 = load_config(str(tmp_config_file))
        cfg2 = load_config(str(tmp_config_file))
        assert cfg1 is cfg2

    def test_lru_cache_different_paths_are_different_entries(self, tmp_path, pipeline_cfg):
        p1 = tmp_path / "a.json"
        p2 = tmp_path / "b.json"
        content = json.dumps(pipeline_cfg)
        p1.write_text(content)
        p2.write_text(content)
        cfg1 = load_config(str(p1))
        cfg2 = load_config(str(p2))
        assert cfg1 is not cfg2


class TestValidate:
    def test_invalid_pval_gw_above_one(self, pipeline_cfg):
        pipeline_cfg["cis_extract"]["pval_gw"] = 1.5
        with pytest.raises(ValueError, match="pval_gw"):
            _validate(pipeline_cfg)

    def test_invalid_pval_gw_zero(self, pipeline_cfg):
        pipeline_cfg["cis_extract"]["pval_gw"] = 0.0
        with pytest.raises(ValueError, match="pval_gw"):
            _validate(pipeline_cfg)

    def test_invalid_window_kb_zero(self, pipeline_cfg):
        pipeline_cfg["cis_extract"]["window_kb"] = 0
        with pytest.raises(ValueError, match="window_kb"):
            _validate(pipeline_cfg)

    def test_invalid_window_kb_negative(self, pipeline_cfg):
        pipeline_cfg["cis_extract"]["window_kb"] = -100
        with pytest.raises(ValueError, match="window_kb"):
            _validate(pipeline_cfg)

    def test_invalid_maf_min_zero(self, pipeline_cfg):
        pipeline_cfg["cis_extract"]["maf_min"] = 0.0
        with pytest.raises(ValueError, match="maf_min"):
            _validate(pipeline_cfg)

    def test_valid_config_passes(self, pipeline_cfg):
        _validate(pipeline_cfg)  # should not raise

    def test_missing_configured_cohort_raises(self, pipeline_cfg):
        del pipeline_cfg["cohorts"]["UKB_PPP"]
        with pytest.raises(ValueError, match="missing required cohort"):
            _validate(pipeline_cfg)

    def test_invalid_cohort_build_raises(self, pipeline_cfg):
        pipeline_cfg["cohorts"]["UKB_PPP"]["build"] = "b36"
        with pytest.raises(ValueError, match="UKB_PPP"):
            _validate(pipeline_cfg)


class TestGetSection:
    def test_returns_section(self, pipeline_cfg):
        section = get_section(pipeline_cfg, "cis_extract")
        assert section["window_kb"] == 500

    def test_missing_section_raises_key_error(self, pipeline_cfg):
        with pytest.raises(KeyError, match="nonexistent"):
            get_section(pipeline_cfg, "nonexistent")


class TestGetCohortConfig:
    def test_returns_cohort_config(self, pipeline_cfg):
        cfg = get_cohort_config(pipeline_cfg, "UKB_PPP")
        assert cfg["build"] == "hg19"

    def test_returns_cohort_build(self, pipeline_cfg):
        assert get_cohort_build(pipeline_cfg, "ARIC_EA") == "hg38"

    def test_missing_cohort_raises_key_error(self, pipeline_cfg):
        with pytest.raises(KeyError, match="missing cohort"):
            get_cohort_config(pipeline_cfg, "NOPE")
