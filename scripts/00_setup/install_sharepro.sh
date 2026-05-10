#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."  # repo root

SHAREPRO_DIR=tools/SharePro_coloc
SHAREPRO_COMMIT="main"  # pin to a specific commit SHA after first clone

mkdir -p tools

if [ -d "$SHAREPRO_DIR/.git" ]; then
  echo "  SharePro already cloned at $SHAREPRO_DIR"
else
  git clone https://github.com/zhwm/SharePro_coloc.git "$SHAREPRO_DIR"
  echo "  Cloned SharePro to $SHAREPRO_DIR"
fi

# Install SharePro Python dependencies into current venv
if [ -f "$SHAREPRO_DIR/requirements.txt" ]; then
  uv pip install -r "$SHAREPRO_DIR/requirements.txt"
  echo "  Installed SharePro requirements"
fi

echo "  SharePro: OK"
