#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."  # repo root

echo "=== Installing Python dependencies ==="
uv pip install -e .

echo "=== Installing R packages ==="
Rscript scripts/00_setup/install_packages.R

echo "=== Cloning SharePro ==="
bash scripts/00_setup/install_sharepro.sh

echo "=== Fetching hg19→hg38 liftover chain ==="
mkdir -p data/ref
CHAIN=data/ref/hg19ToHg38.over.chain.gz
if [ ! -f "$CHAIN" ]; then
  curl -L -o "$CHAIN" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
  echo "  Downloaded: $CHAIN"
else
  echo "  Already exists: $CHAIN"
fi

echo "=== Checking plink2 executable ==="
if ! command -v plink2 >/dev/null 2>&1; then
  echo "ERROR: plink2 not found on PATH."
  echo "  Install PLINK2 and ensure the 'plink2' command is available."
  exit 1
fi
echo "  plink2: OK ($(plink2 --version 2>&1 | head -1))"

echo "=== Checking Synapse credentials ==="
if [ ! -f "$HOME/.synapseConfig" ]; then
  echo "WARNING: ~/.synapseConfig not found."
  echo "  Synapse streaming (UKB-PPP, Fenland) will fail."
  echo "  Create it with: python -c \"import synapseclient; synapseclient.Synapse().login()\""
else
  echo "  ~/.synapseConfig: OK"
fi

echo "=== Setup complete ==="
