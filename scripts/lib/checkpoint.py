"""Per-unit resumable checkpointing."""
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable, TypeVar

T = TypeVar("T")


class Checkpoint:
    """
    Tracks completed units via a JSON ledger at `state_path`.

    Usage:
        cp = Checkpoint(processed_data/ARIC_EA/_state_02.json)
        for protein in cp.remaining(all_proteins, key=lambda p: p.seqid):
            ... process protein ...
            cp.mark_done(protein.seqid)
    """

    def __init__(self, state_path: Path):
        self._path = state_path
        self._done: set[str] = set()
        self._failed: dict[str, dict] = {}
        self._dirty = False
        self._load()

    def _load(self) -> None:
        if self._path.exists():
            try:
                data = json.loads(self._path.read_text())
                self._done = set(data.get("done", []))
                status = data.get("status", {})
                if isinstance(status, dict):
                    self._failed = {
                        key: payload
                        for key, payload in status.items()
                        if isinstance(payload, dict) and payload.get("state") == "failed"
                    }
                else:
                    self._failed = {}
            except (json.JSONDecodeError, OSError):
                self._done = set()
                self._failed = {}

    def _save(self) -> None:
        self._path.parent.mkdir(parents=True, exist_ok=True)
        status = {
            key: {
                "state": "failed",
                "reason": payload.get("reason", ""),
                "updated_at": payload.get("updated_at", ""),
            }
            for key, payload in sorted(self._failed.items())
        }
        for key in sorted(self._done):
            status[key] = {
                "state": "success",
                "reason": "",
                "updated_at": "",
            }
        self._path.write_text(
            json.dumps({"done": sorted(self._done), "status": status}, indent=2)
        )
        self._dirty = False

    def is_done(self, key: str) -> bool:
        return key in self._done

    def is_failed(self, key: str) -> bool:
        return key in self._failed

    def failure_reason(self, key: str) -> str | None:
        if key not in self._failed:
            return None
        return str(self._failed[key].get("reason", ""))

    def mark_done(self, key: str, save: bool = True) -> None:
        self._done.add(key)
        self._failed.pop(key, None)
        self._dirty = True
        if save:
            self._save()

    def mark_failed(self, key: str, reason: str, save: bool = True) -> None:
        self._done.discard(key)
        self._failed[key] = {
            "reason": reason,
            "updated_at": _utc_now_iso(),
        }
        self._dirty = True
        if save:
            self._save()

    def flush(self) -> None:
        if self._dirty:
            self._save()

    def remaining(self, items: Iterable[T], key=lambda x: x, include_failed: bool = True) -> Iterable[T]:
        """Yield items whose key is not already in the done set."""
        if include_failed:
            return [item for item in items if key(item) not in self._done]
        return [
            item for item in items
            if key(item) not in self._done and key(item) not in self._failed
        ]

    @property
    def n_done(self) -> int:
        return len(self._done)

    @property
    def n_failed(self) -> int:
        return len(self._failed)


def output_exists(
    path: Path,
    required_cols: Iterable[str] | None = None,
    min_rows: int = 0,
) -> bool:
    """
    Check whether a TSV-like output file exists and is minimally valid.

    - If ``required_cols`` is provided, the header must contain all required columns.
    - If ``min_rows`` > 0, the file must contain at least that many non-empty data rows.
    """
    if not path.exists() or path.stat().st_size <= 0:
        return False

    if required_cols is None and min_rows <= 0:
        return True

    required = list(required_cols or [])
    try:
        with path.open("r", encoding="utf-8") as fh:
            header_line = fh.readline()
            if not header_line:
                return False
            header = header_line.rstrip("\r\n").split("\t")
            if required and any(col not in header for col in required):
                return False
            if min_rows <= 0:
                return True

            n_rows = 0
            for line in fh:
                if line.strip():
                    n_rows += 1
                    if n_rows >= min_rows:
                        return True
            return False
    except OSError:
        return False


def _utc_now_iso() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")
