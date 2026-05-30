"""
Tests for UKB_female (ProteoNexus) extraction pipeline.

Tests:
  1. normalize_protonexus_rows — column mapping (EA=allele1, OA=allele0, etc.)
  2. normalize_protonexus_rows — N is per-SNP from n_obs (not a fixed constant)
  3. build_read_fn — returns None for a gene not in cis_raw
  4. build_read_fn — reads a plain TSV and returns normalized DataFrame
  5. protonexus_unpack — filters rows to the requested cis window from a mock tar
  6. Integration — real SSD data (skipped if /Volumes/Extreme SSD not mounted)
"""
import gzip
import importlib
import io
import tarfile
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest

# Use importlib because directory names starting with digits aren't valid Python identifiers
_ukb_mod = importlib.import_module("scripts.02_cis_pqtl_extract.ukb_female")
_unpack_mod = importlib.import_module("scripts.02_cis_pqtl_extract.protonexus_unpack")

normalize_protonexus_rows = _ukb_mod.normalize_protonexus_rows
build_read_fn = _ukb_mod.build_read_fn
run_unpack = _unpack_mod.run_unpack


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_gemma_rows(
    chrom="1",
    positions=(1_000_000, 2_000_000, 3_000_000),
    n_obs_values=(17988, 17950, 17960),
) -> list[dict]:
    rows = []
    for i, (pos, n_obs) in enumerate(zip(positions, n_obs_values)):
        rows.append({
            "chr":       chrom,
            "rs":        f"rs{1000 + i}",
            "ps":        pos,
            "n_mis":     5,
            "n_obs":     n_obs,
            "allele1":   "A",
            "allele0":   "G",
            "af":        0.3,
            "beta":      0.05,
            "se":        0.01,
            "p_wald":    1e-10,
            "pip_susie": 0.9,
            "fwer":      0.001,
        })
    return rows


def _make_cis_raw_tsv(rows: list[dict]) -> str:
    """Render rows as TSV text (header + data)."""
    return pd.DataFrame(rows).to_csv(sep="\t", index=False)


def _make_tar_with_gene(tmp_path: Path, gene: str, rows: list[dict]) -> Path:
    """Create a minimal ProteoNexus tar containing <gene>/output/summ_female2.assoc.txt.gz."""
    gene_lower = gene.lower()
    tar_path = tmp_path / "ProteoNexus_pQTL_protein_test.tar"

    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write(pd.DataFrame(rows).to_csv(sep="\t", index=False).encode())
    gz_bytes = buf.getvalue()

    with tarfile.open(tar_path, "w") as tf:
        info = tarfile.TarInfo(name=f"{gene_lower}/output/summ_female2.assoc.txt.gz")
        info.size = len(gz_bytes)
        tf.addfile(info, io.BytesIO(gz_bytes))

    return tar_path


# ---------------------------------------------------------------------------
# Test 1: column mapping
# ---------------------------------------------------------------------------

class TestNormalizeRows:
    def test_column_mapping(self):
        rows = _make_gemma_rows(chrom="chr1")  # with "chr" prefix

        result = normalize_protonexus_rows(rows)

        assert result is not None
        assert not result.empty

        # allele1 → EA, allele0 → OA
        assert list(result["EA"].unique()) == ["A"]
        assert list(result["OA"].unique()) == ["G"]

        # p_wald → pval
        assert "pval" in result.columns
        assert (result["pval"] == 1e-10).all()

        # rs → rsid
        assert result["rsid"].iloc[0] == "rs1000"

        # chr with "chr" prefix stripped → chrom
        assert (result["chrom"] == "1").all()
        assert "chr" not in result.columns

        # ps → pos (integer)
        assert result["pos"].dtype == "int64"

        # n_obs → N
        assert "N" in result.columns
        assert result["N"].iloc[0] == 17988

        # Unused GEMMA cols dropped
        for col in ("chr", "ps", "n_mis", "pip_susie", "fwer", "n_obs"):
            assert col not in result.columns

    def test_chrom_without_prefix(self):
        rows = _make_gemma_rows(chrom="2")  # no "chr" prefix
        result = normalize_protonexus_rows(rows)
        assert result is not None
        assert (result["chrom"] == "2").all()

    # ---------------------------------------------------------------------------
    # Test 2: N is per-SNP from n_obs
    # ---------------------------------------------------------------------------

    def test_n_is_per_snp(self):
        rows = _make_gemma_rows(n_obs_values=(17988, 17500, 17000))
        result = normalize_protonexus_rows(rows)

        assert result is not None
        ns = list(result["N"])
        # N must vary per row — comes from n_obs, not a fixed cohort constant
        assert ns == [17988, 17500, 17000], f"Expected per-SNP N, got {ns}"

    def test_empty_rows_returns_none(self):
        assert normalize_protonexus_rows([]) is None


# ---------------------------------------------------------------------------
# Tests 3–4: build_read_fn
# ---------------------------------------------------------------------------

class TestReadFn:
    def test_returns_none_for_missing_file(self, tmp_path):
        from scripts.lib.schema import ProteinMeta

        protein = ProteinMeta(
            seqid="NOTHERE", gene="NOTHERE", uniprot="",
            chrom="1", tss=1_000_000, build="hg19", source_cohort="UKB_female",
        )

        with mock.patch.object(_ukb_mod, "UKB_FEMALE_CIS_RAW", tmp_path):
            read_fn = build_read_fn()
            result = read_fn(protein)

        assert result is None

    def test_reads_plain_tsv(self, tmp_path):
        from scripts.lib.schema import ProteinMeta

        gene = "GMPR2"
        rows = _make_gemma_rows(chrom="1")
        (tmp_path / f"{gene}.tsv").write_text(_make_cis_raw_tsv(rows))

        protein = ProteinMeta(
            seqid=gene, gene=gene, uniprot="",
            chrom="1", tss=2_000_000, build="hg19", source_cohort="UKB_female",
        )

        with mock.patch.object(_ukb_mod, "UKB_FEMALE_CIS_RAW", tmp_path):
            read_fn = build_read_fn()
            result = read_fn(protein)

        assert result is not None
        assert not result.empty
        for col in ("chrom", "pos", "rsid", "EA", "OA", "EAF", "beta", "se", "pval", "N"):
            assert col in result.columns, f"Missing column: {col}"
        assert result["EA"].iloc[0] == "A"
        assert result["OA"].iloc[0] == "G"
        assert result["N"].iloc[0] == 17988


# ---------------------------------------------------------------------------
# Test 5: protonexus_unpack cis-window filtering
# ---------------------------------------------------------------------------

class TestProtonexusUnpack:
    def test_filters_to_cis_window(self, tmp_path):
        """Only rows within the requested TSS window should be written to cis_raw."""
        gene = "GMPR2"
        tss = 2_000_000
        window_kb = 500
        flank = window_kb * 1_000

        rows = [
            # inside window
            {"chr": "1", "rs": "rs1", "ps": tss,
             "n_mis": 0, "n_obs": 17988, "allele1": "A", "allele0": "G",
             "af": 0.3, "beta": 0.05, "se": 0.01, "p_wald": 1e-10,
             "pip_susie": 0.9, "fwer": 0.001},
            # outside — too far left
            {"chr": "1", "rs": "rs2", "ps": tss - flank - 1,
             "n_mis": 0, "n_obs": 17988, "allele1": "A", "allele0": "G",
             "af": 0.3, "beta": 0.05, "se": 0.01, "p_wald": 1e-10,
             "pip_susie": 0.9, "fwer": 0.001},
            # outside — too far right
            {"chr": "1", "rs": "rs3", "ps": tss + flank + 1,
             "n_mis": 0, "n_obs": 17988, "allele1": "A", "allele0": "G",
             "af": 0.3, "beta": 0.05, "se": 0.01, "p_wald": 1e-10,
             "pip_susie": 0.9, "fwer": 0.001},
        ]

        _make_tar_with_gene(tmp_path, gene, rows)

        cis_raw_dir = tmp_path / "cis_raw"
        cis_raw_dir.mkdir()
        cohort_base = tmp_path / "cohort"
        cohort_base.mkdir()

        # Pre-populate TSS cache so Ensembl is not contacted
        tss_cache_path = cohort_base / "_tss_hg19.tsv"
        pd.DataFrame([{"gene": gene, "chrom": "1", "tss": tss}]).to_csv(
            tss_cache_path, sep="\t", index=False
        )

        with (
            mock.patch.object(_unpack_mod, "UKB_FEMALE_DIR", tmp_path),
            mock.patch.object(_unpack_mod, "UKB_FEMALE_CIS_RAW", cis_raw_dir),
            mock.patch.object(_unpack_mod, "cohort_dir", return_value=cohort_base),
        ):
            n = run_unpack(limit=None, window_kb=window_kb)

        assert n == 1, f"Expected 1 gene written, got {n}"
        out_path = cis_raw_dir / f"{gene}.tsv"
        assert out_path.exists(), f"Expected {out_path} to exist"

        result = pd.read_csv(out_path, sep="\t")
        assert len(result) == 1, f"Expected 1 row (in-window only), got {len(result)}"
        assert result["rs"].iloc[0] == "rs1"

    def test_window_scoped_checkpoint_ignores_old_unpack_state(self, tmp_path):
        """Old ±500 kb unpack checkpoints must not skip the new ±1 Mb cache."""
        gene = "GMPR2"
        tss = 2_000_000
        window_kb = 1000
        rows = [
            {"chr": "1", "rs": "rs_new_window", "ps": tss + 750_000,
             "n_mis": 0, "n_obs": 17988, "allele1": "A", "allele0": "G",
             "af": 0.3, "beta": 0.05, "se": 0.01, "p_wald": 1e-10,
             "pip_susie": 0.9, "fwer": 0.001},
        ]
        _make_tar_with_gene(tmp_path, gene, rows)

        cis_raw_dir = tmp_path / "cis_raw_1000kb"
        cis_raw_dir.mkdir()
        cohort_base = tmp_path / "cohort"
        cohort_base.mkdir()
        pd.DataFrame([{"gene": gene, "chrom": "1", "tss": tss}]).to_csv(
            cohort_base / "_tss_hg19.tsv", sep="\t", index=False
        )
        # Simulate the pre-migration unpack state. It should not suppress a 1000 kb run.
        (cohort_base / "_state_02_unpack.json").write_text(
            '{"done": ["GMPR2"], "status": {"GMPR2": {"state": "success"}}}'
        )

        with (
            mock.patch.object(_unpack_mod, "UKB_FEMALE_DIR", tmp_path),
            mock.patch.object(_unpack_mod, "UKB_FEMALE_CIS_RAW", cis_raw_dir),
            mock.patch.object(_unpack_mod, "cohort_dir", return_value=cohort_base),
        ):
            n = run_unpack(limit=None, window_kb=window_kb)

        assert n == 1
        result = pd.read_csv(cis_raw_dir / f"{gene}.tsv", sep="\t")
        assert result["rs"].tolist() == ["rs_new_window"]
        assert (cohort_base / "_state_02_unpack_1000kb.json").exists()


# ---------------------------------------------------------------------------
# Test 6: Integration — real SSD data (skipped if not mounted)
# ---------------------------------------------------------------------------

_SSD_MOUNTED = Path("/Volumes/Extreme SSD/ProteoNexus").exists()
_SKIP_REASON = "ProteoNexus SSD not mounted at /Volumes/Extreme SSD/ProteoNexus"


@pytest.mark.skipif(not _SSD_MOUNTED, reason=_SKIP_REASON)
def test_integration_real_data():
    """
    End-to-end: unpack 3 genes from the real SSD, then extract cis-sumstats.
    Verifies output schema and file existence.
    """
    import subprocess
    import sys

    result_unpack = subprocess.run(
        [sys.executable, "-m",
         "scripts.02_cis_pqtl_extract.protonexus_unpack",
         "--limit", "3"],
        capture_output=True, text=True,
    )
    assert result_unpack.returncode == 0, (
        f"protonexus_unpack.py failed:\n{result_unpack.stderr}"
    )

    from scripts.lib.paths import UKB_FEMALE_CIS_RAW
    tsv_files = list(UKB_FEMALE_CIS_RAW.glob("*.tsv"))
    assert len(tsv_files) >= 1, "Expected at least 1 cis_raw TSV after unpack"

    df = pd.read_csv(tsv_files[0], sep="\t")
    for col in ("chr", "rs", "ps", "allele1", "allele0", "af", "beta", "se", "p_wald", "n_obs"):
        assert col in df.columns, f"Expected GEMMA column {col!r} in cis_raw TSV"

    result_extract = subprocess.run(
        [sys.executable, "-m", "scripts.02_cis_pqtl_extract.ukb_female",
         "--limit", "3", "--workers", "1"],
        capture_output=True, text=True,
    )
    assert result_extract.returncode == 0, (
        f"ukb_female.py failed:\n{result_extract.stderr}"
    )

    from scripts.lib.paths import filtered_cis_pqtls_dir, raw_cis_sumstats_dir
    out_files = list(filtered_cis_pqtls_dir("UKB_female").glob("*.tsv"))
    raw_files = list(raw_cis_sumstats_dir("UKB_female").glob("*.tsv.gz"))
    assert len(out_files) >= 1, "Expected at least 1 filtered_cis_pqtls TSV after extraction"
    assert len(raw_files) >= 1, "Expected at least 1 raw_cis_sumstats TSV.GZ after extraction"

    out_df = pd.read_csv(out_files[0], sep="\t")
    for col in ("seqid", "gene", "chrom", "pos", "rsid", "EA", "OA",
                "EAF", "beta", "se", "pval", "N", "build"):
        assert col in out_df.columns, f"Missing output column: {col}"
