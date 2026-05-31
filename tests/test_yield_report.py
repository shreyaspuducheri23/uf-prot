"""Tests for scripts/qc/yield_report.py."""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from scripts.qc.yield_report import CheckpointStats, report_cohort, run_report, main


# ── Fixtures ──────────────────────────────────────────────────────────────────

def _make_index(cdir: Path, n: int) -> None:
    rows = [{"seqid": f"SeqId_{i}", "gene": f"G{i}", "uniprot": f"U{i}",
             "chrom": "1", "tss": 1_000_000, "build": "hg19"}
            for i in range(n)]
    pd.DataFrame(rows).to_csv(cdir / "protein_index.tsv", sep="\t", index=False)


def _make_checkpoint(cdir: Path, state_file: str, done: list[str], failed: dict[str, str]) -> None:
    status = {}
    for k in done:
        status[k] = {"state": "success", "reason": "", "updated_at": ""}
    for k, reason in failed.items():
        status[k] = {"state": "failed", "reason": reason, "updated_at": "2026-01-01T00:00:00+00:00"}
    data = {"done": sorted(done), "status": status}
    (cdir / state_file).write_text(json.dumps(data, indent=2))


def _make_tsv_files(out_dir: Path, n: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for i in range(n):
        (out_dir / f"SeqId_{i}.tsv").write_text("seqid\tgene\n")


def _write_variant_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _variant_rows(seqid: str, positions: list[int]) -> list[dict]:
    return [
        {
            "seqid": seqid,
            "gene": "GENE",
            "uniprot": "U",
            "chrom": "1",
            "pos": pos,
            "rsid": f"rs{pos}",
            "EA": "A",
            "OA": "G",
            "EAF": 0.25,
            "beta": 0.1,
            "se": 0.01,
            "pval": 1e-9,
            "N": 1000,
            "build": "hg19",
        }
        for pos in positions
    ]


# ── report_cohort ─────────────────────────────────────────────────────────────

class TestReportCohort:
    def test_no_index_returns_empty(self, tmp_path):
        rows, warns = report_cohort("deCODE", processed_dir=tmp_path)
        assert rows == []
        assert warns == []

    def test_stage_not_run_stops_iteration(self, tmp_path):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 10)
        # No checkpoint → no rows
        rows, warns = report_cohort("deCODE", processed_dir=tmp_path)
        assert rows == []

    def test_single_stage_counts(self, tmp_path):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 10)
        _make_checkpoint(cdir, "_state_02.json",
                         done=[f"SeqId_{i}" for i in range(8)],
                         failed={"SeqId_8": "timeout", "SeqId_9": "timeout"})
        _make_tsv_files(cdir / "filtered_cis_pqtls", 6)

        rows, warns = report_cohort("deCODE", processed_dir=tmp_path)
        assert len(rows) == 1
        r = rows[0]
        assert r["n_input"] == 10
        assert r["n_output"] == 6
        assert r["n_failed"] == 2
        assert r["pct_yield"] == pytest.approx(60.0, abs=0.1)
        assert r["pct_failed"] == pytest.approx(20.0, abs=0.1)

    def test_multi_stage_n_input_chains(self, tmp_path):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 100)
        _make_checkpoint(cdir, "_state_02.json",
                         done=[f"SeqId_{i}" for i in range(80)], failed={})
        _make_tsv_files(cdir / "filtered_cis_pqtls", 70)
        _make_checkpoint(cdir, "_state_03.json",
                         done=[f"SeqId_{i}" for i in range(60)], failed={})
        _make_tsv_files(cdir / "instruments", 55)

        rows, _ = report_cohort("deCODE", processed_dir=tmp_path)
        assert len(rows) == 2
        assert rows[0]["n_input"] == 100
        assert rows[1]["n_input"] == 70  # previous stage n_output

    def test_warn_on_high_failure_rate(self, tmp_path):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 100)
        # 10 not done (n_input - n_done = 100 - 90 = 10, 10% > 5%)
        _make_checkpoint(cdir, "_state_02.json",
                         done=[f"SeqId_{i}" for i in range(90)], failed={})
        _make_tsv_files(cdir / "filtered_cis_pqtls", 80)

        _, warns = report_cohort("deCODE", processed_dir=tmp_path)
        assert any("not done" in w for w in warns)

    def test_warn_on_silent_abandonment(self, tmp_path):
        """Proteins absent from checkpoint entirely (pre-fix silent failures) are counted."""
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 100)
        # Only 50 marked done, none marked failed — 50 were silently abandoned
        _make_checkpoint(cdir, "_state_02.json",
                         done=[f"SeqId_{i}" for i in range(50)], failed={})
        _make_tsv_files(cdir / "filtered_cis_pqtls", 50)

        rows, warns = report_cohort("deCODE", processed_dir=tmp_path)
        assert rows[0]["n_failed"] == 50  # n_input(100) - n_done(50)
        assert any("not done" in w for w in warns)

    def test_no_warn_when_thresholds_met(self, tmp_path):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 100)
        # 98 done → n_failed = 2 (2% < 5%)
        _make_checkpoint(cdir, "_state_02.json",
                         done=[f"SeqId_{i}" for i in range(98)], failed={})
        _make_tsv_files(cdir / "filtered_cis_pqtls", 96)

        _, warns = report_cohort("deCODE", processed_dir=tmp_path)
        assert warns == []

    def test_tsv_appended(self, tmp_path):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 10)
        _make_checkpoint(cdir, "_state_02.json",
                         done=[f"SeqId_{i}" for i in range(9)], failed={})
        _make_tsv_files(cdir / "filtered_cis_pqtls", 8)

        run_report(["deCODE"], processed_dir=tmp_path)

        tsv_path = tmp_path / "_yield_report.tsv"
        assert tsv_path.exists()
        df = pd.read_csv(tsv_path, sep="\t")
        assert len(df) == 1
        assert df.loc[0, "cohort"] == "deCODE"
        assert df.loc[0, "stage"] == "filtered_cis_pqtls"
        assert df.loc[0, "n_input"] == 10
        assert df.loc[0, "n_output"] == 8


# ── --strict exit code ────────────────────────────────────────────────────────

class TestStrictMode:
    def test_strict_exits_1_on_warn(self, tmp_path, monkeypatch):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 100)
        # 20% failure rate — triggers warning
        _make_checkpoint(cdir, "_state_02.json",
                         done=[f"SeqId_{i}" for i in range(80)],
                         failed={f"SeqId_{i}": "err" for i in range(80, 100)})
        _make_tsv_files(cdir / "filtered_cis_pqtls", 70)

        monkeypatch.setattr("scripts.qc.yield_report.PROCESSED", tmp_path)

        with pytest.raises(SystemExit) as exc_info:
            main(["--cohort", "deCODE", "--strict"])

        assert exc_info.value.code == 1

    def test_strict_exits_0_on_no_warn(self, tmp_path, capsys):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 100)
        _make_checkpoint(cdir, "_state_02.json",
                         done=[f"SeqId_{i}" for i in range(99)],
                         failed={f"SeqId_99": "err"})
        _make_tsv_files(cdir / "filtered_cis_pqtls", 98)

        any_warn = run_report(["deCODE"], strict=True, processed_dir=tmp_path)
        assert not any_warn

    def test_run_report_returns_true_on_warn(self, tmp_path):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 100)
        _make_checkpoint(cdir, "_state_02.json",
                         done=[f"SeqId_{i}" for i in range(50)], failed={})
        _make_tsv_files(cdir / "filtered_cis_pqtls", 45)

        any_warn = run_report(["deCODE"], processed_dir=tmp_path)
        assert any_warn

    def test_run_report_returns_false_on_clean(self, tmp_path):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 100)
        _make_checkpoint(cdir, "_state_02.json",
                         done=[f"SeqId_{i}" for i in range(98)],
                         failed={f"SeqId_{i}": "err" for i in range(98, 100)})
        _make_tsv_files(cdir / "filtered_cis_pqtls", 96)

        any_warn = run_report(["deCODE"], processed_dir=tmp_path)
        assert not any_warn


# ── stage-aware row/locus accounting ─────────────────────────────────────────

class TestComprehensiveStages:
    def test_liftover_locus_loss_visible_when_protein_survives(self, tmp_path):
        cdir = tmp_path / "UKB_female"
        cdir.mkdir()
        _make_index(cdir, 1)
        _write_variant_tsv(cdir / "instruments" / "SeqId_0.tsv", _variant_rows("SeqId_0", [100, 200, 300]))
        lifted = _variant_rows("SeqId_0", [100, 200])
        for row in lifted:
            row["chrom_hg38"] = row["chrom"]
            row["pos_hg38"] = row["pos"] + 1000
        _write_variant_tsv(cdir / "instruments_hg38" / "SeqId_0.tsv", lifted)
        _make_checkpoint(cdir, "_state_04.json", done=["SeqId_0"], failed={})

        rows, warns = report_cohort("UKB_female", processed_dir=tmp_path)
        liftover = next(r for r in rows if r["stage"] == "instruments_hg38")
        assert liftover["units_input"] == 1
        assert liftover["units_output"] == 1
        assert liftover["pct_unit_yield"] == 100.0
        assert liftover["rows_input"] == 3
        assert liftover["rows_output"] == 2
        assert liftover["loci_input"] == 3
        assert liftover["loci_output"] == 2
        assert any("dropped" in w and "liftover" in w for w in warns)

    def test_row_and_unique_locus_yield_can_differ(self, tmp_path):
        cdir = tmp_path / "UKB_female"
        cdir.mkdir()
        _make_index(cdir, 1)
        rows = _variant_rows("SeqId_0", [100, 100])
        rows[0]["rsid"] = "rsA"
        rows[1]["rsid"] = "rsB"
        _write_variant_tsv(cdir / "instruments" / "SeqId_0.tsv", rows)
        lifted = [dict(rows[0], chrom_hg38="1", pos_hg38=1100)]
        _write_variant_tsv(cdir / "instruments_hg38" / "SeqId_0.tsv", lifted)
        _make_checkpoint(cdir, "_state_04.json", done=["SeqId_0"], failed={})

        rows, _ = report_cohort("UKB_female", processed_dir=tmp_path)
        liftover = next(r for r in rows if r["stage"] == "instruments_hg38")
        assert liftover["rows_input"] == 2
        assert liftover["rows_output"] == 1
        assert liftover["loci_input"] == 1
        assert liftover["loci_output"] == 1
        assert liftover["pct_row_yield"] == 50.0
        assert liftover["pct_locus_yield"] == 100.0

    def test_dropped_locus_detail_tsv_contains_missing_variant(self, tmp_path):
        cdir = tmp_path / "UKB_female"
        cdir.mkdir()
        _make_index(cdir, 1)
        _write_variant_tsv(cdir / "instruments" / "SeqId_0.tsv", _variant_rows("SeqId_0", [100, 200]))
        lifted = _variant_rows("SeqId_0", [100])
        lifted[0]["chrom_hg38"] = "1"
        lifted[0]["pos_hg38"] = 1100
        _write_variant_tsv(cdir / "instruments_hg38" / "SeqId_0.tsv", lifted)
        _make_checkpoint(cdir, "_state_04.json", done=["SeqId_0"], failed={})

        run_report(["UKB_female"], processed_dir=tmp_path)
        dropped = pd.read_csv(tmp_path / "_yield_report_dropped_loci.tsv", sep="\t")
        assert "rs200" in set(dropped["variant_id"])
        assert "1:200" in set(dropped["locus"])

    def test_mr_checkpoint_done_without_result_is_flagged(self, tmp_path, monkeypatch):
        cdir = tmp_path / "ARIC_EA"
        cdir.mkdir()
        _make_index(cdir, 2)
        for seqid in ("SeqId_0", "SeqId_1"):
            rows = _variant_rows(seqid, [100])
            rows[0]["mr_keep"] = True
            _write_variant_tsv(cdir / "harmonised" / f"{seqid}.tsv", rows)
        pd.DataFrame([{"cohort": "ARIC_EA", "seqid": "SeqId_0", "n_snps": 1}]).to_csv(
            cdir / "mr_results.tsv", sep="\t", index=False
        )

        monkeypatch.setattr(
            "scripts.qc.yield_report._read_r_checkpoint",
            lambda path: CheckpointStats(done={"SeqId_0", "SeqId_1"}),
        )
        rows, warns = report_cohort("ARIC_EA", processed_dir=tmp_path)
        mr = next(r for r in rows if r["stage"] == "mr")
        assert mr["units_done_without_output"] == 1
        assert any("checkpoint-done" in w for w in warns)

    def test_sensitivity_single_snp_proteins_are_not_applicable(self, tmp_path):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 2)
        pd.DataFrame([
            {"cohort": "deCODE", "seqid": "SeqId_0", "n_snps": 1},
            {"cohort": "deCODE", "seqid": "SeqId_1", "n_snps": 2},
        ]).to_csv(cdir / "mr_results.tsv", sep="\t", index=False)
        pd.DataFrame([{"seqid": "SeqId_1", "passes_sensitivity": True}]).to_csv(
            cdir / "sensitivity.tsv", sep="\t", index=False
        )

        rows, warns = report_cohort("deCODE", processed_dir=tmp_path)
        sensitivity = next(r for r in rows if r["stage"] == "sensitivity")
        assert sensitivity["units_input"] == 1
        assert sensitivity["units_output"] == 1
        assert sensitivity["units_not_applicable"] == 1
        assert warns == []

    def test_coloc_sharepro_checkpoint_failure_reported(self, tmp_path):
        cdir = tmp_path / "Fenland"
        cdir.mkdir()
        _make_index(cdir, 1)
        pd.DataFrame([{
            "cohort": "Fenland",
            "seqid": "SeqId_0",
            "n_snps": 2,
            "fdr_pass": True,
        }]).to_csv(cdir / "mr_results.tsv", sep="\t", index=False)
        pd.DataFrame([{"seqid": "SeqId_0", "passes_sensitivity": True}]).to_csv(
            cdir / "sensitivity.tsv", sep="\t", index=False
        )
        region = tmp_path / "coloc" / "regions" / "Fenland" / "SeqId_0"
        _write_variant_tsv(region / "exposure.tsv", _variant_rows("SeqId_0", [100]))
        pd.DataFrame([{
            "chromosome": "1",
            "base_pair_location": 100,
            "rsid": "rs100",
            "effect_allele": "A",
            "other_allele": "G",
        }]).to_csv(region / "outcome.tsv", sep="\t", index=False)
        _make_checkpoint(cdir, "_state_08_sharepro.json", done=[], failed={"SeqId_0": "insufficient_common_snps"})

        rows, _ = report_cohort("Fenland", processed_dir=tmp_path)
        sharepro = next(r for r in rows if r["stage"] == "sharepro")
        assert sharepro["units_input"] == 1
        assert sharepro["units_output"] == 0
        assert sharepro["units_failed_cp"] == 1

    def test_legacy_cis_sumstats_fallback(self, tmp_path):
        cdir = tmp_path / "deCODE"
        cdir.mkdir()
        _make_index(cdir, 2)
        _make_checkpoint(cdir, "_state_02.json", done=["SeqId_0"], failed={})
        _write_variant_tsv(cdir / "cis_sumstats" / "SeqId_0.tsv", _variant_rows("SeqId_0", [100]))

        rows, _ = report_cohort("deCODE", processed_dir=tmp_path)
        filtered = next(r for r in rows if r["stage"] == "filtered_cis_pqtls")
        assert filtered["units_input"] == 2
        assert filtered["units_output"] == 1
        assert filtered["rows_output"] == 1
