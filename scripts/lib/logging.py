"""Timestamped logging setup: stdout + log file + manifest."""
import logging
import sys
import csv
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from scripts.lib.paths import LOGS


def setup_logger(step_name: str, level: int = logging.INFO) -> logging.Logger:
    LOGS.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    log_path = LOGS / f"{step_name}_{ts}.log"

    fmt = logging.Formatter(
        fmt="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    logger = logging.getLogger(step_name)
    logger.setLevel(level)
    logger.handlers.clear()

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(fmt)
    logger.addHandler(file_handler)

    logger.info(f"Log file: {log_path}")
    return logger


def _git_sha() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        return result.stdout.strip() if result.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


class RunManifest:
    """Append a single-line record to logs/_manifest.tsv on context exit."""

    FIELDS = ["script", "args", "git_sha", "start_utc", "end_utc", "n_units", "status"]
    PATH = LOGS / "_manifest.tsv"

    def __init__(self, script: str, args: str = ""):
        self.script = script
        self.args = args
        self.git_sha = _git_sha()
        self._start = datetime.now(timezone.utc)
        self.n_units: int = 0
        self.status: str = "error"

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.status = "ok" if exc_type is None else f"error:{exc_type.__name__}"
        end = datetime.now(timezone.utc)
        row = {
            "script": self.script,
            "args": self.args,
            "git_sha": self.git_sha,
            "start_utc": self._start.isoformat(),
            "end_utc": end.isoformat(),
            "n_units": self.n_units,
            "status": self.status,
        }
        LOGS.mkdir(parents=True, exist_ok=True)
        write_header = not self.PATH.exists()
        with open(self.PATH, "a", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=self.FIELDS, delimiter="\t")
            if write_header:
                writer.writeheader()
            writer.writerow(row)
        return False  # do not suppress exceptions
