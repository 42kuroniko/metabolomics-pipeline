#!/usr/bin/env bash
# Author: Rong Wu
# Setup conda env + msconvert for metabolomics pipeline.
# Usage: bash scripts/setup_env.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_DIR="$PROJECT_DIR/envs"

echo "=========================================="
echo "metabolomics-pipeline setup"
echo "=========================================="
echo ""

if ! command -v conda &> /dev/null; then
    if command -v mamba &> /dev/null; then
        CONDA_CMD="mamba"
    else
        echo "[ERROR] conda or mamba not found"
        echo "Install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
else
    CONDA_CMD="conda"
fi
echo "[OK] conda/mamba: $CONDA_CMD"

ENV_FILE="$ENV_DIR/metabolomics.yaml"
if conda env list | grep -q "^metabolomics "; then
    echo "[SKIP] conda env 'metabolomics' already exists"
else
    echo ""
    echo "[INFO] Creating conda env (may take 5-10 min)..."
    $CONDA_CMD env create -f "$ENV_FILE"
fi

echo ""
echo "=========================================="
echo "Setup complete."
echo "=========================================="
echo ""
echo "Activate:  conda activate metabolomics"
echo "msconvert: requires Wine on macOS"
echo "Run:       snakemake -np  # preview"
echo "           snakemake --use-conda -j 4"
