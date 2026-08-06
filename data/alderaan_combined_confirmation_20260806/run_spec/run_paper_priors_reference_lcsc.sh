#!/usr/bin/env bash
set -euo pipefail

# Confirmation arm for the only untested factorial interaction: the reference
# limb-darkening catalog and printed Table 1 transit priors with LC+SC data.
# This is intentionally separate from the completed 82-fit matrix.
export TARGET_CSV="targets_short_cadence_validation.csv"
export CATALOG_SOURCE="sagear_ld_reference_catalog.csv"
export CATALOG_NAME="sagear_validation_catalog.csv"
export RUN_ID="sagear_validation_paper_priors_reference_lcsc"
export PROJECT_DIR="$PWD/projects/paper_priors_reference_lcsc"
export ALDERAAN_REPO="${ALDERAAN_PAPER_PRIOR_REPO:-$HOME/alderaan_sagear_paper_priors}"
export CADENCE_MODE="both"
export SEED_OFFSET="0"
export JOBS="${JOBS:-6}"

bash run_batch.sh
