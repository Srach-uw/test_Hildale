#!/usr/bin/env bash
set -euo pipefail
export TARGET_CSV="targets_ld_reference_validation.csv"
export CATALOG_SOURCE="sagear_ld_reference_catalog.csv"
export CATALOG_NAME="sagear_validation_catalog.csv"
export RUN_ID="sagear_validation_reference_lc"
export PROJECT_DIR="$PWD/projects/reference_lc"
export ALDERAAN_REPO="${ALDERAAN_REPO:-$HOME/alderaan_sagear_pinned}"
export CADENCE_MODE="long"
export SEED_OFFSET="0"
export JOBS="${JOBS:-6}"
bash run_batch.sh
