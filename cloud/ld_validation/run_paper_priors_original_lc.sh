#!/usr/bin/env bash
set -euo pipefail
export TARGET_CSV="targets_repeatability_validation.csv"
export CATALOG_SOURCE="sagear_original_full_system_catalog.csv"
export CATALOG_NAME="sagear_validation_catalog.csv"
export RUN_ID="sagear_validation_paper_priors_original_lc"
export PROJECT_DIR="$PWD/projects/paper_priors_original_lc"
export ALDERAAN_REPO="${ALDERAAN_PAPER_PRIOR_REPO:-$HOME/alderaan_sagear_paper_priors}"
export CADENCE_MODE="long"
export SEED_OFFSET="0"
export JOBS="${JOBS:-6}"
bash run_batch.sh
