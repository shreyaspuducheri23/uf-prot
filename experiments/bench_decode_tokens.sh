set -euo pipefail
cd /Users/spuduch/Research/leiomyoma_proteomics

LIMIT=2
BENCH_ROOT="bench_decode_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BENCH_ROOT"

SRC_DECODE_DIR="$(realpath data/raw/deCODE)"
URLS1="$(realpath data/raw/deCODE/bulk_urls.txt)"   # token 1
URLS2="$(realpath data/bulk_urls_2.txt)"            # token 2
EAF_CACHE="$(realpath processed_data/deCODE/_eaf_cache.pkl)"
TSS_CACHE="$(realpath processed_data/deCODE/_tss_hg38.tsv)"

setup_raw () {
  raw_dir="$1"; urls="$2"
  mkdir -p "$raw_dir/deCODE"
  ln -s "$SRC_DECODE_DIR/assocvariants.annotated.txt.gz" "$raw_dir/deCODE/assocvariants.annotated.txt.gz"
  ln -s "$urls" "$raw_dir/deCODE/bulk_urls.txt"
}

setup_proc_cache () {
  proc_dir="$1"
  mkdir -p "$proc_dir/deCODE"
  cp "$EAF_CACHE" "$proc_dir/deCODE/_eaf_cache.pkl"
  cp "$TSS_CACHE" "$proc_dir/deCODE/_tss_hg38.tsv"
}

run_case () {
  label="$1"; workers="$2"; raw_dir="$3"
  proc_dir="$BENCH_ROOT/processed_$label"
  logs_dir="$BENCH_ROOT/logs_$label"
  out_log="$BENCH_ROOT/${label}.log"

  setup_proc_cache "$proc_dir"
  mkdir -p "$logs_dir"

  echo "=== $label (workers=$workers, limit=$LIMIT) ==="
  (
    export LEIO_RAW_DIR="$raw_dir"
    export LEIO_PROCESSED_DIR="$proc_dir"
    export LEIO_LOGS_DIR="$logs_dir"
    /usr/bin/time -p uv run python scripts/02_cis_pqtl_extract/decode.py --workers "$workers" --limit "$LIMIT"
  ) >"$out_log" 2>&1 || true

  echo "real: $(rg '^real ' "$out_log" -m1 || echo NA)"
  echo "done: $(rg 'deCODE: done' "$out_log" -m1 || echo NA)"
  echo "retry_warn: $(rg -c 'Download failed \\(attempt' "$out_log" || true)"
  echo "rate_limit_hits: $(rg -c '429|Too Many Requests|Retry-After' "$out_log" || true)"
  echo
}

RAW1="$BENCH_ROOT/raw_token1"
RAW2="$BENCH_ROOT/raw_token2"
setup_raw "$RAW1" "$URLS1"
setup_raw "$RAW2" "$URLS2"

run_case token1_w4 4 "$RAW1"
run_case token2_w4 4 "$RAW2"
run_case token1_w12 12 "$RAW1"
run_case token2_w12 12 "$RAW2"

echo "logs in: $BENCH_ROOT"
