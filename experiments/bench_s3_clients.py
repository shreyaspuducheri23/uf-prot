#!/usr/bin/env python3
"""
experiments/bench_s3_clients.py

Tests two hypotheses about deCODE S3 download throughput:
  1. Does per-thread boto3 client vs a shared client matter?
  2. What worker count maximizes aggregate throughput from Iceland S3?

Strategy: stream a fixed number of raw compressed bytes from N protein files
in parallel, then abort. This measures network throughput directly without
waiting for full GWAS files (~13 min each) to complete.

Run from the repo root:
    uv run python experiments/bench_s3_clients.py [--probe-mb 8] [--n 16]
"""
import argparse
import json
import random
import statistics
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import boto3

# ── S3 connection constants ────────────────────────────────────────────────────
S3_ENDPOINT   = "https://s3-ext.decode.is:10443"
S3_BUCKET     = "largescaleplasma-2023"
S3_ACCESS_KEY = "SE0AV795UKCQ338YKWP4"
S3_SECRET_KEY = "/mkkvYtFJkO+NAhxcm3OhNKAdvwQivhbdQRLeJ/c"

KEY_INDEX_PATH = Path("processed_data/deCODE/_s3_key_index.json")
READ_CHUNK = 128 * 1024  # 128 KB per body.read() call

# ── Client factories ───────────────────────────────────────────────────────────
_shared_client = None
_shared_lock = threading.Lock()


def _make_client() -> boto3.client:
    return boto3.client(
        "s3",
        endpoint_url=S3_ENDPOINT,
        aws_access_key_id=S3_ACCESS_KEY,
        aws_secret_access_key=S3_SECRET_KEY,
        region_name="us-east-1",
    )


def get_shared_client() -> boto3.client:
    global _shared_client
    if _shared_client is None:
        with _shared_lock:
            if _shared_client is None:
                _shared_client = _make_client()
    return _shared_client


_thread_local = threading.local()


def get_thread_client() -> boto3.client:
    if not hasattr(_thread_local, "client"):
        _thread_local.client = _make_client()
    return _thread_local.client


CLIENT_MODES = {
    "shared":     get_shared_client,
    "per_thread": get_thread_client,
}

# ── Core probe function ────────────────────────────────────────────────────────

def probe_key(key: str, client_fn, probe_bytes: int) -> dict:
    """
    Open an S3 object, read probe_bytes of raw (compressed) bytes, then close.
    Returns timing and throughput info.
    """
    t0 = time.monotonic()
    client = client_fn()

    body = client.get_object(Bucket=S3_BUCKET, Key=key)["Body"]
    total = 0
    try:
        while total < probe_bytes:
            want = min(READ_CHUNK, probe_bytes - total)
            chunk = body.read(want)
            if not chunk:
                break
            total += len(chunk)
    finally:
        body.close()

    elapsed = time.monotonic() - t0
    return {
        "key":     key,
        "bytes":   total,
        "elapsed": elapsed,
        "mbps":    total / elapsed / 1e6 if elapsed > 0 else 0.0,
    }


# ── Trial runner ───────────────────────────────────────────────────────────────

def run_trial(keys: list[str], workers: int, client_fn, probe_bytes: int) -> dict:
    """
    Submit probe_key for every key using ThreadPoolExecutor(workers).
    Returns aggregate stats.
    """
    t0 = time.monotonic()
    per_protein: list[dict] = []
    errors: list[str] = []

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(probe_key, k, client_fn, probe_bytes): k for k in keys}
        for fut in as_completed(futures):
            try:
                per_protein.append(fut.result())
            except Exception as exc:
                errors.append(f"{futures[fut]}: {exc}")

    wall = time.monotonic() - t0
    total_bytes = sum(r["bytes"] for r in per_protein)
    mbps_values = [r["mbps"] for r in per_protein]
    per_file_elapsed = [r["elapsed"] for r in per_protein]

    return {
        "n":                  len(per_protein),
        "workers":            workers,
        "wall_s":             wall,
        "total_mb":           total_bytes / 1e6,
        "aggregate_mbps":     total_bytes / wall / 1e6 if wall > 0 else 0.0,
        "median_per_file_s":  statistics.median(per_file_elapsed) if per_file_elapsed else 0.0,
        "median_per_file_mbps": statistics.median(mbps_values) if mbps_values else 0.0,
        "errors":             errors,
    }


# ── Calibration ───────────────────────────────────────────────────────────────

def calibrate_probe_bytes(key: str, target_seconds: float = 20.0) -> int:
    """
    Download 1 MB from a single key to estimate bandwidth, then return
    probe_bytes that would take ~target_seconds at that bandwidth.
    """
    print(f"Calibrating with {key.split('/')[-1]} ...", flush=True)
    r = probe_key(key, get_shared_client, 1 * 1024 * 1024)
    bw = r["mbps"]
    print(f"  1 MB in {r['elapsed']:.1f}s → {bw:.3f} MB/s per thread")
    target_bytes = int(bw * target_seconds * 1e6)
    target_bytes = max(2 * 1024 * 1024, min(target_bytes, 30 * 1024 * 1024))
    print(f"  Targeting {target_seconds:.0f}s/probe → probe_bytes = {target_bytes/1e6:.1f} MB")
    return target_bytes


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Benchmark deCODE S3 client strategies")
    parser.add_argument("--probe-mb",  type=float, default=None,
                        help="MB to read per probe (default: auto-calibrate for ~20s/thread)")
    parser.add_argument("--n",         type=int,   default=16,
                        help="Number of proteins per trial (default: 16)")
    parser.add_argument("--workers",   type=int,   nargs="+", default=[1, 4, 8, 12, 16],
                        help="Worker counts to test (default: 1 4 8 12 16)")
    parser.add_argument("--modes",     nargs="+", default=["shared", "per_thread"],
                        choices=list(CLIENT_MODES.keys()),
                        help="Client modes to test (default: shared per_thread)")
    parser.add_argument("--seed",      type=int,   default=42)
    args = parser.parse_args()

    # Load key index
    if not KEY_INDEX_PATH.exists():
        sys.exit(f"Key index not found: {KEY_INDEX_PATH}")
    with open(KEY_INDEX_PATH) as fh:
        key_index = json.load(fh)

    all_keys = list(key_index.values())
    rng = random.Random(args.seed)
    # Oversample so we can draw fresh keys per trial; keeps all trials independent
    probe_keys = rng.sample(all_keys, min(args.n, len(all_keys)))

    print(f"\n{'='*60}")
    print(f"deCODE S3 client benchmark  (n={len(probe_keys)} proteins)")
    print(f"Endpoint: {S3_ENDPOINT}")
    print(f"{'='*60}\n")

    # Probe size
    if args.probe_mb is not None:
        probe_bytes = int(args.probe_mb * 1e6)
        print(f"Probe bytes: {probe_bytes/1e6:.1f} MB (manual)")
    else:
        probe_bytes = calibrate_probe_bytes(probe_keys[0])
    print()

    # Run grid
    results = []
    for workers in args.workers:
        for mode in args.modes:
            client_fn = CLIENT_MODES[mode]
            tag = f"workers={workers:>2d}  mode={mode}"
            print(f"Running {tag} ...", flush=True)
            r = run_trial(probe_keys, workers, client_fn, probe_bytes)
            r["mode"] = mode
            results.append(r)
            err_str = f"  ⚠ {len(r['errors'])} errors" if r["errors"] else ""
            print(
                f"  wall={r['wall_s']:5.1f}s  "
                f"aggregate={r['aggregate_mbps']:.3f} MB/s  "
                f"median_per_file={r['median_per_file_s']:.1f}s "
                f"({r['median_per_file_mbps']:.3f} MB/s/thread)"
                + err_str
            )
            if r["errors"]:
                for e in r["errors"][:3]:
                    print(f"    {e}")

    # ── Summary table ────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    header = f"{'workers':>8}  {'mode':>12}  {'wall_s':>8}  {'agg_MB/s':>10}  {'med_file_s':>11}  {'med_MB/s/thr':>13}"
    print(header)
    print("-" * len(header))
    for r in results:
        print(
            f"{r['workers']:>8}  {r['mode']:>12}  {r['wall_s']:>8.1f}  "
            f"{r['aggregate_mbps']:>10.3f}  {r['median_per_file_s']:>11.1f}  "
            f"{r['median_per_file_mbps']:>13.3f}"
        )

    # ── Interpretation ───────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("INTERPRETATION")
    print(f"{'='*60}")

    # Find best
    best = max(results, key=lambda r: r["aggregate_mbps"])
    w1_shared = next((r for r in results if r["workers"] == 1 and r["mode"] == "shared"), None)
    w4_shared = next((r for r in results if r["workers"] == 4 and r["mode"] == "shared"), None)

    if w1_shared:
        per_protein_bw = w1_shared["aggregate_mbps"]
        print(f"\nSingle-thread bandwidth:    {per_protein_bw:.3f} MB/s")
        if w4_shared:
            ideal_4x = per_protein_bw * 4
            actual_4x = w4_shared["aggregate_mbps"]
            efficiency = actual_4x / ideal_4x * 100
            print(f"w=4 shared (actual/ideal):  {actual_4x:.3f} / {ideal_4x:.3f} MB/s  ({efficiency:.0f}% efficient)")

    print(f"\nBest config: workers={best['workers']} mode={best['mode']} → {best['aggregate_mbps']:.3f} MB/s")

    # Shared vs per_thread at each worker count
    print("\nShared vs per_thread delta (positive = per_thread faster):")
    for w in args.workers:
        s = next((r for r in results if r["workers"] == w and r["mode"] == "shared"), None)
        p = next((r for r in results if r["workers"] == w and r["mode"] == "per_thread"), None)
        if s and p:
            delta_pct = (p["aggregate_mbps"] - s["aggregate_mbps"]) / s["aggregate_mbps"] * 100
            print(f"  w={w:>2d}:  {delta_pct:+.1f}%  ({s['aggregate_mbps']:.3f} vs {p['aggregate_mbps']:.3f} MB/s)")

    # Estimate for current run
    if w4_shared:
        # Current run has 5147 proteins remaining; assume average file is entirely streamed
        # (worst case: gene on chr22 → stream whole file)
        # Use median_per_file_s scaled to full-file size as estimate
        median_file_s = w4_shared["median_per_file_s"]
        probe_fraction = probe_bytes / (median_file_s * w4_shared["median_per_file_mbps"] * 1e6 + 1e-9)
        # safer: just use observed rate
        proteins_remaining = 5147
        est_hours = (proteins_remaining * median_file_s) / best["workers"] / 3600
        print(f"\nProjected remaining time at best config ({best['workers']} workers):")
        print(f"  {proteins_remaining} proteins × {median_file_s:.0f}s (median_per_file) / {best['workers']} workers")
        print(f"  ≈ {est_hours:.1f} hours  (assumes current bandwidth holds)")


if __name__ == "__main__":
    main()
