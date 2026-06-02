#!/usr/bin/env bash
# Parallel variant of run_pipeline.sh — runs cohorts concurrently at the bash level.
#
# Key differences:
#   - Cohort extraction (step 2) runs all 5 cohorts simultaneously; each cohort
#     hits a different data source so bandwidth pools don't compete.
#   - Steps 3, 4, 5, 6, 7, 8a, 8b, 8c split "--cohort all" into 5 per-cohort
#     background jobs.
#   - --workers (default 1) governs network-bound download threads per cohort.
#     For local-file steps, --local-workers (default: nproc) is used instead.
#   - Per-cohort checkpoint files (step_3_ARIC_EA.done etc.) let a failed cohort
#     resume without re-running the others.
#
# Backward-compat note: old global step_3.done / step_4.done / step_5.done /
# step_6.done / step_7.done / step_8a.done / step_8b.done / step_8c.done files are NOT
# recognized here.  To adopt this runner on an existing partial run: delete old
# global .done files and rely on inner checkpoints for resumption, OR manually
# create per-cohort variants (e.g. step_3_ARIC_EA.done) for completed cohorts.
#
# Usage:
#   bash run_pipeline_parallel.sh [--workers N] [--local-workers N] [--force]
#                                  [--skip-step N[,N,...]] [--strict]
#
# Options:
#   --workers N          Download threads per network-bound cohort (default: 1)
#   --local-workers N    Worker threads for local-file steps (default: nproc)
#   --force              Re-run all steps even if already completed
#   --skip-step N,...    Comma-separated step IDs to skip (e.g. --skip-step 2b,2c)
#   --strict             Fail pipeline if yield-report detects suspicious yield drops

set -euo pipefail

# ── Config ──────────────────────────────────────────────────────────────────
WORKERS=1
LOCAL_WORKERS=$(nproc 2>/dev/null || echo 8)
FORCE=0
SKIP_STEPS=()
STRICT_FLAG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workers)       WORKERS="$2";       shift 2 ;;
    --local-workers) LOCAL_WORKERS="$2"; shift 2 ;;
    --force)         FORCE=1;            shift ;;
    --strict)        STRICT_FLAG="--strict"; shift ;;
    --skip-step)
      IFS=',' read -ra _skip <<< "$2"
      SKIP_STEPS+=("${_skip[@]}")
      shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

TIMESTAMP="$(date +%Y%m%dT%H%M%S)"
LOG_DIR="$REPO_ROOT/logs"
CKPT_DIR="$LOG_DIR/.pipeline_steps"
PIPELINE_LOG="$LOG_DIR/pipeline_parallel_${TIMESTAMP}.log"

mkdir -p "$LOG_DIR" "$CKPT_DIR"

# Temp dir used by the parallel helpers (bash 3.2-compatible substitute for
# associative arrays, which require bash 4+).
_PAR_DIR=$(mktemp -d)
trap 'rm -rf "$_PAR_DIR"' EXIT

# ── Helpers ──────────────────────────────────────────────────────────────────
log() { local msg="[$(date +%Y-%m-%dT%H:%M:%S)] $*"; echo "$msg" | tee -a "$PIPELINE_LOG"; }

is_skipped() {
  local step="$1"
  local base="$step"
  case "$step" in
    3_*|4_*|5_*|6_*|7_*|8a_*|8b_*|8c_*) base="${step%%_*}" ;;
  esac
  for s in "${SKIP_STEPS[@]+"${SKIP_STEPS[@]}"}"; do
    [[ "$s" == "$step" || "$s" == "$base" ]] && return 0
  done
  return 1
}

run_step() {
  local step="$1"; shift
  local desc="$1"; shift
  local cmd=("$@")

  if is_skipped "$step"; then
    log "SKIP  step $step ($desc) — explicitly skipped"
    return 0
  fi

  local ckpt="$CKPT_DIR/step_${step}.done"
  if [[ $FORCE -eq 0 && -f "$ckpt" ]]; then
    log "SKIP  step $step ($desc) — already completed (delete $ckpt or use --force to re-run)"
    return 0
  fi

  log "START step $step ($desc)"
  local t0=$SECONDS
  if "${cmd[@]}" 2>&1 | tee -a "$PIPELINE_LOG"; then
    touch "$ckpt"
    log "DONE  step $step ($desc) — $((SECONDS - t0))s"
  else
    log "FAIL  step $step ($desc) — exit $?"
    exit 1
  fi
}

run_step_always() {
  local step="$1"; shift
  local desc="$1"; shift
  local cmd=("$@")

  if is_skipped "$step"; then
    log "SKIP  step $step ($desc) — explicitly skipped"
    return 0
  fi

  local ckpt="$CKPT_DIR/step_${step}.done"
  log "START step $step ($desc)"
  local t0=$SECONDS
  if "${cmd[@]}" 2>&1 | tee -a "$PIPELINE_LOG"; then
    touch "$ckpt"
    log "DONE  step $step ($desc) — $((SECONDS - t0))s"
  else
    log "FAIL  step $step ($desc) — exit $?"
    exit 1
  fi
}

# ── Parallel step helpers ────────────────────────────────────────────────────
# Uses $_PAR_DIR (a temp dir) rather than associative arrays so the script runs
# under macOS system bash 3.2 (which predates bash 4 assoc-array support).

# Queue a step as a background job.  Skipping and checkpoint logic mirrors
# run_step; each job writes output to its own log file so lines don't interleave.
queue_step() {
  local step="$1"; shift
  local desc="$1"; shift
  local cmd=("$@")

  if is_skipped "$step"; then
    log "SKIP  step $step ($desc) — explicitly skipped"
    return 0
  fi

  local ckpt="$CKPT_DIR/step_${step}.done"
  if [[ $FORCE -eq 0 && -f "$ckpt" ]]; then
    log "SKIP  step $step ($desc) — already completed"
    return 0
  fi

  local slog="$LOG_DIR/step_${step}_${TIMESTAMP}.log"
  (
    local t0=$SECONDS
    if "${cmd[@]}" >"$slog" 2>&1; then
      touch "$ckpt"
      log "DONE  step $step ($desc) — $((SECONDS - t0))s"
    else
      log "FAIL  step $step ($desc) — see $slog"
      exit 1
    fi
  ) &
  local pid="$!"
  log "QUEUE step $step ($desc) pid=$pid → $slog"
  echo "$pid" > "$_PAR_DIR/pid_${step}"
  echo "$slog" > "$_PAR_DIR/log_${step}"
  echo "$desc" > "$_PAR_DIR/desc_${step}"
}

# Wait for every queued background step; abort if any failed.
flush_parallel() {
  local any_fail=0 total=0 n_done=0 n_fail=0
  local pid_file step pid slog desc
  local t_group=$SECONDS

  # Count and announce what's in flight.
  local running=""
  for pid_file in "$_PAR_DIR"/pid_*; do
    [[ -f "$pid_file" ]] || continue
    step="${pid_file##*/pid_}"
    pid=$(cat "$pid_file")
    desc=$(cat "$_PAR_DIR/desc_${step}")
    running="$running $step(pid=$pid,$desc)"
    total=$((total + 1))
  done
  if [[ $total -eq 0 ]]; then return 0; fi
  log "WAIT  $total parallel jobs running:$running"

  for pid_file in "$_PAR_DIR"/pid_*; do
    [[ -f "$pid_file" ]] || continue
    step="${pid_file##*/pid_}"
    pid=$(cat "$pid_file")
    slog=$(cat "$_PAR_DIR/log_${step}")
    if wait "$pid"; then
      n_done=$((n_done + 1))
      log "      ($n_done/$total) step $step finished"
    else
      n_fail=$((n_fail + 1))
      log "FAIL  step $step — see $slog"
      any_fail=1
    fi
  done
  rm -f "$_PAR_DIR"/pid_* "$_PAR_DIR"/log_* "$_PAR_DIR"/desc_*

  local elapsed=$((SECONDS - t_group))
  if [[ $any_fail -eq 0 ]]; then
    log "GROUP all $total steps done — ${elapsed}s"
  else
    log "GROUP $n_fail/$total steps FAILED — aborting after ${elapsed}s"
    return 1
  fi
}

# ── Pipeline ─────────────────────────────────────────────────────────────────
log "========================================================"
log "Leiomyoma proteomics pipeline (parallel)  |  $TIMESTAMP"
log "workers=${WORKERS} (net)  local_workers=${LOCAL_WORKERS}  force=$FORCE  strict=${STRICT_FLAG}  skip=(${SKIP_STEPS[*]+"${SKIP_STEPS[*]}"})"
log "log → $PIPELINE_LOG"
log "========================================================"

run_step "1" "outcome prep (Kim GWAS)" \
  uv run python scripts/01_outcome_prep/prep_kim.py

# ── Step 2: cis-pQTL extraction ──────────────────────────────────────────────
# Phase A: all five sources are bandwidth-independent — run simultaneously.
#   Network-bound cohorts use --workers 1 (S3/Synapse bandwidth is the bottleneck;
#   extra threads just queue more requests against the same pipe).
#   Local-file cohorts use LOCAL_WORKERS (CPU/disk-bound, benefits from concurrency).

queue_step "2a" "cis-pQTL extract: ARIC (local)" \
  uv run python scripts/02_cis_pqtl_extract/aric.py
  # Local parquet files; argparse doesn't expose --workers despite docstring

queue_step "2b" "cis-pQTL extract: deCODE (S3, workers=${WORKERS})" \
  uv run python scripts/02_cis_pqtl_extract/decode.py --workers "$WORKERS"

queue_step "2c" "cis-pQTL extract: UKB-PPP (Synapse, workers=${WORKERS})" \
  uv run python scripts/02_cis_pqtl_extract/ukbppp.py --workers "$WORKERS"

queue_step "2d" "cis-pQTL extract: Fenland (Synapse)" \
  uv run python scripts/02_cis_pqtl_extract/fenland.py
  # No --workers arg; Synapse-bound regardless

queue_step "2e_prep" "unpack ProteoNexus tars" \
  uv run python scripts/02_cis_pqtl_extract/protonexus_unpack.py

flush_parallel

# Phase B: UKB_female reads the unpacked local files — must follow 2e_prep.
run_step "2e" "cis-pQTL extract: UKB-female (local, workers=${LOCAL_WORKERS})" \
  uv run python scripts/02_cis_pqtl_extract/ukb_female.py --workers "$LOCAL_WORKERS"

for cohort in ARIC_EA deCODE UKB_PPP Fenland UKB_female; do
  uv run python scripts/qc/yield_report.py --cohort "$cohort" $STRICT_FLAG
done

# ── Step 3: LD clumping ───────────────────────────────────────────────────────
for cohort in ARIC_EA deCODE UKB_PPP Fenland UKB_female; do
  queue_step "3_${cohort}" "LD clump: ${cohort}" \
    uv run python scripts/03_clump/clump.py --cohort "$cohort"
done
flush_parallel
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

# ── Step 4: liftover hg19 → GRCh38 ──────────────────────────────────────────
for cohort in ARIC_EA deCODE UKB_PPP Fenland UKB_female; do
  queue_step "4_${cohort}" "liftover: ${cohort}" \
    uv run python scripts/04_liftover/instruments_to_hg38.py --cohort "$cohort"
done
flush_parallel
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

# ── Step 5: harmonise with Kim outcome ───────────────────────────────────────
# ~308 m sequentially → ~62 m in parallel (biggest win after step 2)
for cohort in ARIC_EA deCODE UKB_PPP Fenland UKB_female; do
  queue_step "5_${cohort}" "harmonise: ${cohort}" \
    uv run python scripts/05_harmonise/harmonise.py --cohort "$cohort"
done
flush_parallel
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

# ── Steps 6–7: MR + sensitivity ──────────────────────────────────────────────
for cohort in ARIC_EA deCODE UKB_PPP Fenland UKB_female; do
  queue_step "6_${cohort}" "two-sample MR + BH-FDR: ${cohort}" \
    Rscript scripts/06_mr/run_mr.R --cohort "$cohort"
done
flush_parallel
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

for cohort in ARIC_EA deCODE UKB_PPP Fenland UKB_female; do
  queue_step "7_${cohort}" "sensitivity analyses: ${cohort}" \
    Rscript scripts/07_sensitivity/run_sensitivity.R --cohort "$cohort"
done
flush_parallel
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

# ── Step 8: colocalization ────────────────────────────────────────────────────
for cohort in ARIC_EA deCODE UKB_PPP Fenland UKB_female; do
  queue_step "8a_${cohort}" "coloc extract regions: ${cohort}" \
    uv run python scripts/08_coloc/extract_regions.py --cohort "$cohort"
done
flush_parallel
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

for cohort in ARIC_EA deCODE UKB_PPP Fenland UKB_female; do
  queue_step "8b_${cohort}" "coloc SharePro: ${cohort}" \
    uv run python scripts/08_coloc/sharepro.py --cohort "$cohort"
done
flush_parallel
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

for cohort in ARIC_EA deCODE UKB_PPP Fenland UKB_female; do
  queue_step "8c_${cohort}" "coloc: coloc.abf sensitivity: ${cohort}" \
    Rscript scripts/08_coloc/coloc_abf.R --cohort "$cohort"
done
flush_parallel
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

# ── Steps 9–9b: assemble results ─────────────────────────────────────────────
run_step_always "9" "assemble final results table" \
  uv run python scripts/09_assemble/assemble.py
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

run_step_always "9b" "cross-cohort gene-level summary" \
  uv run python scripts/09_assemble/cross_cohort.py
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

log "========================================================"
log "Pipeline complete. Results: processed_data/final_results.tsv"
log "========================================================"
