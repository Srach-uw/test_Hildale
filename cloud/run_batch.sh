#!/usr/bin/env bash
set -euo pipefail

# Conda env activation is per-shell-session and does not persist across a
# fresh SSH login (e.g. after a VM restart) - self-activate here so a
# forgotten manual `conda activate alderaan` doesn't kill the whole batch
# with a cryptic "python: command not found". Found live, 2026-07-03.
if [ "${CONDA_DEFAULT_ENV:-}" != "alderaan" ]; then
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
  conda activate alderaan
fi

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

JOBS="${JOBS:-30}"
RUN_ID="${RUN_ID:-sagear_missing}"
TARGET_CSV="${TARGET_CSV:-targets_missing_launchable.csv}"
PROJECT_DIR="${PROJECT_DIR:-$PWD/alderaan_project}"
ALDERAAN_REPO="${ALDERAAN_REPO:-$HOME/alderaan}"

mkdir -p "$PROJECT_DIR/Catalogs" "$PROJECT_DIR/Data" "$PROJECT_DIR/Results" "$PROJECT_DIR/Figures" "$PROJECT_DIR/status" logs
cp sagear_missing_catalog.csv "$PROJECT_DIR/Catalogs/sagear_missing_catalog.csv"
# bin/detrend_and_estimate_ttvs.py loads Catalogs/holczer_2016_kepler_ttvs.txt
# from PROJECT_DIR unconditionally; ships in the ALDERAAN repo, not our bundle.
cp "$ALDERAAN_REPO/Catalogs/holczer_2016_kepler_ttvs.txt" "$PROJECT_DIR/Catalogs/holczer_2016_kepler_ttvs.txt"
cp "$TARGET_CSV" "$PROJECT_DIR/targets.csv"
chmod +x run_one_target.sh

echo "Running $TARGET_CSV with JOBS=$JOBS RUN_ID=$RUN_ID"
python validate_bundle.py --targets "$TARGET_CSV" --catalog sagear_missing_catalog.csv

tail -n +2 "$TARGET_CSV" | awk -F, '{print $2","$3}' \
  | parallel --colsep ',' -j "$JOBS" --joblog "logs/parallel_joblog_${RUN_ID}.tsv" './run_one_target.sh {1} {2}' \
  > "logs/batch_${RUN_ID}_stdout.log" 2> "logs/batch_${RUN_ID}_stderr.log"

bash summarize_progress.sh
echo "Batch complete. Results are under $PROJECT_DIR/Results/$RUN_ID"
