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

echo "=== Checking plink1 binary ==="
PLINK1=data/ld_ref/plink1.90
if [ ! -x "$PLINK1" ]; then
  echo "ERROR: $PLINK1 not found or not executable."
  echo "  The plink1.90 binary should be available via the ld_ref symlink."
  exit 1
fi
echo "  plink1: OK ($("$PLINK1" --version 2>&1 | head -1))"

echo "=== Checking Synapse credentials ==="
if [ ! -f "$HOME/.synapseConfig" ]; then
  echo "WARNING: ~/.synapseConfig not found."
  echo "  Synapse streaming (UKB-PPP, Fenland) will fail."
  echo "  Create it with: python -c \"import synapseclient; synapseclient.Synapse().login()\""
else
  echo "  ~/.synapseConfig: OK"
fi

echo "=== Setup complete ==="
