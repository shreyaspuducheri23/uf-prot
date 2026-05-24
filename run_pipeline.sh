#!/usr/bin/env bash
# Full end-to-end pipeline runner for leiomyoma proteomics MR analysis.
# Each step is resumable (inner scripts checkpoint per-protein); this wrapper
# also checkpoints at the step level — re-run after a failure to pick up where
# you left off. Pass --force to ignore step-level checkpoints and re-run all.
#
# Usage:
#   bash run_pipeline.sh [--workers N] [--force] [--skip-step N[,N,...]] [--strict]
#
# Options:
#   --workers N        Parallel workers for UKB-PPP Synapse streaming (default: 4)
#   --force            Re-run all steps even if already completed
#   --skip-step N,...  Comma-separated step numbers to skip (e.g. --skip-step 2b,2c)
#   --strict           Fail pipeline if yield-report detects suspicious yield drops

set -euo pipefail

# ── Config ──────────────────────────────────────────────────────────────────
WORKERS=4
FORCE=0
SKIP_STEPS=()
STRICT_FLAG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workers)  WORKERS="$2"; shift 2 ;;
    --force)    FORCE=1; shift ;;
    --strict)   STRICT_FLAG="--strict"; shift ;;
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
PIPELINE_LOG="$LOG_DIR/pipeline_${TIMESTAMP}.log"

mkdir -p "$LOG_DIR" "$CKPT_DIR"

# ── Helpers ──────────────────────────────────────────────────────────────────
log() { local msg="[$(date +%Y-%m-%dT%H:%M:%S)] $*"; echo "$msg" | tee -a "$PIPELINE_LOG"; }

is_skipped() {
  local step="$1"
  for s in "${SKIP_STEPS[@]+"${SKIP_STEPS[@]}"}"; do [[ "$s" == "$step" ]] && return 0; done
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

# ── Pipeline ─────────────────────────────────────────────────────────────────
log "========================================================"
log "Leiomyoma proteomics pipeline  |  $TIMESTAMP"
log "workers=$WORKERS  force=$FORCE  strict=${STRICT_FLAG}  skip=(${SKIP_STEPS[*]+"${SKIP_STEPS[*]}"})"
log "log → $PIPELINE_LOG"
log "========================================================"

run_step "1"  "outcome prep (Kim GWAS)" \
  uv run python scripts/01_outcome_prep/prep_kim.py

run_step "2a" "cis-pQTL extract: ARIC" \
  uv run python scripts/02_cis_pqtl_extract/aric.py
uv run python scripts/qc/yield_report.py --cohort ARIC_EA $STRICT_FLAG

run_step "2b" "cis-pQTL extract: deCODE (${WORKERS} workers)" \
  uv run python scripts/02_cis_pqtl_extract/decode.py --workers "$WORKERS"
uv run python scripts/qc/yield_report.py --cohort deCODE $STRICT_FLAG

run_step "2c" "cis-pQTL extract: UKB-PPP (Synapse, ${WORKERS} workers)" \
  uv run python scripts/02_cis_pqtl_extract/ukbppp.py --workers "$WORKERS"
uv run python scripts/qc/yield_report.py --cohort UKB_PPP $STRICT_FLAG

run_step "2d" "cis-pQTL extract: Fenland (Synapse)" \
  uv run python scripts/02_cis_pqtl_extract/fenland.py
uv run python scripts/qc/yield_report.py --cohort Fenland $STRICT_FLAG

run_step "2e_prep" "unpack ProteoNexus tars → cis TSVs" \
  uv run python scripts/02_cis_pqtl_extract/protonexus_unpack.py

run_step "2e" "cis-pQTL extract: UKB-female (ProteoNexus)" \
  uv run python scripts/02_cis_pqtl_extract/ukb_female.py --workers "$WORKERS"
uv run python scripts/qc/yield_report.py --cohort UKB_female $STRICT_FLAG

run_step "3"  "LD clumping (all cohorts)" \
  uv run python scripts/03_clump/clump.py --cohort all
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

run_step "4"  "liftover hg19 → GRCh38 (all cohorts)" \
  uv run python scripts/04_liftover/instruments_to_hg38.py --cohort all
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

run_step "5"  "harmonise with Kim outcome (all cohorts)" \
  uv run python scripts/05_harmonise/harmonise.py --cohort all
uv run python scripts/qc/yield_report.py --cohort all $STRICT_FLAG

run_step "6"  "two-sample MR + BH-FDR" \
  Rscript scripts/06_mr/run_mr.R

run_step "7"  "sensitivity analyses" \
  Rscript scripts/07_sensitivity/run_sensitivity.R

run_step "8a" "coloc: extract ±1 Mb regions (all cohorts)" \
  uv run python scripts/08_coloc/extract_regions.py --cohort all

run_step "8b" "coloc: SharePro (all cohorts)" \
  uv run python scripts/08_coloc/sharepro.py --cohort all

run_step "8c" "coloc: coloc.abf sensitivity (R)" \
  Rscript scripts/08_coloc/coloc_abf.R

run_step "9"  "assemble final results table" \
  uv run python scripts/09_assemble/assemble.py

run_step "9b" "cross-cohort gene-level summary" \
  uv run python scripts/09_assemble/cross_cohort.py

log "========================================================"
log "Pipeline complete. Results: processed_data/final_results.tsv"
log "========================================================"
