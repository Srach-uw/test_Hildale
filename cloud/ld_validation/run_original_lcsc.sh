#!/usr/bin/env bash
set -euo pipefail
export TARGET_CSV="targets_short_cadence_validation.csv"
export CATALOG_SOURCE="sagear_original_full_system_catalog.csv"
export CATALOG_NAME="sagear_validation_catalog.csv"
export RUN_ID="sagear_validation_original_lcsc"
export PROJECT_DIR="$PWD/projects/original_lcsc"
export ALDERAAN_REPO="${ALDERAAN_REPO:-$HOME/alderaan_sagear_pinned}"
export CADENCE_MODE="both"
export SEED_OFFSET="0"
export JOBS="${JOBS:-6}"
bash run_batch.sh
