"""Tests for code.lib.checkpoint"""
import json
import tempfile
from pathlib import Path

import pytest

from scripts.lib.checkpoint import Checkpoint, output_exists


class TestCheckpoint:
    def test_empty_on_first_use(self, tmp_path):
        cp = Checkpoint(tmp_path / "state.json")
        assert cp.n_done == 0

    def test_mark_done_persists(self, tmp_path):
        state = tmp_path / "state.json"
        cp = Checkpoint(state)
        cp.mark_done("protein_A")

        cp2 = Checkpoint(state)
        assert cp2.is_done("protein_A")
        assert cp2.n_done == 1

    def test_remaining_excludes_done(self, tmp_path):
        cp = Checkpoint(tmp_path / "state.json")
        cp.mark_done("A")
        result = cp.remaining(["A", "B", "C"])
        assert result == ["B", "C"]

    def test_remaining_with_key_fn(self, tmp_path):
        cp = Checkpoint(tmp_path / "state.json")
        cp.mark_done("A")
        items = [("A", 1), ("B", 2)]
        result = cp.remaining(items, key=lambda x: x[0])
        assert result == [("B", 2)]

    def test_idempotent_mark(self, tmp_path):
        cp = Checkpoint(tmp_path / "state.json")
        cp.mark_done("X")
        cp.mark_done("X")
        assert cp.n_done == 1

    def test_corrupted_state_resets(self, tmp_path):
        state = tmp_path / "state.json"
        state.write_text("NOT VALID JSON")
        cp = Checkpoint(state)
        assert cp.n_done == 0

    def test_saves_sorted(self, tmp_path):
        state = tmp_path / "state.json"
        cp = Checkpoint(state)
        cp.mark_done("B")
        cp.mark_done("A")
        data = json.loads(state.read_text())
        assert data["done"] == ["A", "B"]

    def test_mark_failed_persists_and_excluded_when_requested(self, tmp_path):
        state = tmp_path / "state.json"
        cp = Checkpoint(state)
        cp.mark_failed("protein_Z", "boom")

        cp2 = Checkpoint(state)
        assert cp2.is_failed("protein_Z")
        assert cp2.failure_reason("protein_Z") == "boom"
        assert cp2.remaining(["protein_Z", "protein_A"], include_failed=False) == ["protein_A"]

    def test_mark_done_clears_failed_state(self, tmp_path):
        cp = Checkpoint(tmp_path / "state.json")
        cp.mark_failed("protein_X", "temporary")
        cp.mark_done("protein_X")
        assert cp.is_done("protein_X")
        assert not cp.is_failed("protein_X")


class TestOutputExists:
    def test_existing_nonempty_file(self, tmp_path):
        f = tmp_path / "out.tsv"
        f.write_text("data")
        assert output_exists(f)

    def test_nonexistent_file(self, tmp_path):
        assert not output_exists(tmp_path / "missing.tsv")

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.tsv"
        f.touch()
        assert not output_exists(f)
