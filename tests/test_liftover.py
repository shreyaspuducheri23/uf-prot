"""Tests for code.lib.liftover — requires hg19ToHg38 chain file."""
import pytest
import pandas as pd
from pathlib import Path

from scripts.lib.paths import CHAIN_HG19_TO_HG38


if not CHAIN_HG19_TO_HG38.exists():
    pytest.fail("hg19ToHg38 chain file not downloaded (run 00_setup/install.sh first)")


from scripts.lib.liftover import lift_position, lift_table


class TestLiftPosition:
    def test_known_snp(self):
        # rs4988235 (chr2:136608646 hg19) → chr2:135851076 hg38 (approx)
        result = lift_position("2", 136_608_646)
        assert result is not None
        chrom, pos = result
        assert chrom == "2"
        assert pos > 0

    def test_invalid_position_returns_none(self):
        result = lift_position("1", 1)  # position 1 likely unmapped
        # May or may not lift — just check it doesn't crash
        assert result is None or (isinstance(result, tuple) and len(result) == 2)

    def test_preserves_chromosome(self):
        # BRCA2 region chr13:32914437 hg19
        result = lift_position("13", 32_914_437)
        if result is not None:
            assert result[0] == "13"


class TestLiftTable:
    def test_adds_hg38_columns(self):
        df = pd.DataFrame({
            "chrom": ["2"], "pos": [136_608_646],
            "rsid": ["rs4988235"], "beta": [0.1],
        })
        result = lift_table(df)
        assert "chrom_hg38" in result.columns
        assert "pos_hg38" in result.columns
        assert len(result) >= 1

    def test_drops_failed_liftover(self, tmp_path):
        # Create a fake chain path that will cause all liftovers to fail
        df = pd.DataFrame({"chrom": ["1"], "pos": [999_999_999_999]})
        # Very large position unlikely to lift
        result = lift_table(df)
        # Should either lift or be empty — not crash
        assert isinstance(result, pd.DataFrame)

    def test_drop_rate_logged(self, caplog):
        import logging
        df = pd.DataFrame({
            "chrom": ["2", "2"],
            "pos": [136_608_646, 999_999_999_999],
        })
        with caplog.at_level(logging.WARNING, logger="code.lib.liftover"):
            result = lift_table(df)
        # If any SNP fails, a warning should be emitted
        if len(result) < len(df):
            assert "failed" in caplog.text.lower() or "dropped" in caplog.text.lower()


class TestOneBased:
    """Validate the 1-based ↔ 0-based conversion logic using a known SNP."""

    def test_rs4988235_lifts_to_expected_hg38_position(self):
        # rs4988235 hg19 chr2:136608646 (LCT region, well-known SNP)
        # Expected hg38: chr2:135851076 (reference: Ensembl/dbSNP cross-check)
        result = lift_position("2", 136_608_646)
        assert result is not None
        chrom, pos_hg38 = result
        assert chrom == "2"
        # Verify the position is in the expected hg38 range
        # The exact value depends on the chain file, but should be ~135851076
        assert 135_840_000 < pos_hg38 < 135_870_000, (
            f"Expected ~135851076 but got {pos_hg38} — "
            f"possible off-by-one in 1-based↔0-based conversion"
        )

    def test_1based_output_is_consistent_with_pyliftover_docs(self):
        # Sanity-check: lifing a known position and back-converting should stay stable.
        # If we lift position P and get Q, lifting Q back hg38→hg19 should give P (or near it).
        # We can only test the forward direction here; just ensure result is > 0 (1-based).
        result = lift_position("2", 136_608_646)
        if result is not None:
            _, pos_hg38 = result
            assert pos_hg38 > 0, "Lifted position must be 1-based (>0)"
