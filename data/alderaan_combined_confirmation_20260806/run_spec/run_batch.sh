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
ALDERAAN_REPO="${ALDERAAN_REPO:-$HOME/alderaan_sagear_pinned}"
CADENCE_MODE="${CADENCE_MODE:-both}"
SEED_OFFSET="${SEED_OFFSET:-0}"
CATALOG_SOURCE="${CATALOG_SOURCE:-sagear_missing_catalog_FIXED.csv}"
if [ ! -s "$CATALOG_SOURCE" ]; then
  CATALOG_SOURCE="sagear_missing_catalog.csv"
fi
CATALOG_NAME="${CATALOG_NAME:-sagear_missing_catalog.csv}"
export RUN_ID PROJECT_DIR ALDERAAN_REPO CATALOG_NAME CADENCE_MODE SEED_OFFSET
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$BUNDLE_DIR/logs"

mkdir -p "$PROJECT_DIR/Catalogs" "$PROJECT_DIR/Data" "$PROJECT_DIR/Results" "$PROJECT_DIR/Figures" "$PROJECT_DIR/status" "$LOG_DIR"
cp "$CATALOG_SOURCE" "$PROJECT_DIR/Catalogs/$CATALOG_NAME"
# bin/detrend_and_estimate_ttvs.py loads Catalogs/holczer_2016_kepler_ttvs.txt
# from PROJECT_DIR unconditionally; ships in the ALDERAAN repo, not our bundle.
cp "$ALDERAAN_REPO/Catalogs/holczer_2016_kepler_ttvs.txt" "$PROJECT_DIR/Catalogs/holczer_2016_kepler_ttvs.txt"
cp "$TARGET_CSV" "$PROJECT_DIR/targets.csv"
chmod +x "$BUNDLE_DIR/run_one_target.sh"

echo "Running $TARGET_CSV with JOBS=$JOBS RUN_ID=$RUN_ID CATALOG_SOURCE=$CATALOG_SOURCE"
echo "ALDERAAN_REPO=$ALDERAAN_REPO COMMIT=$(git -C "$ALDERAAN_REPO" rev-parse HEAD) CADENCE_MODE=$CADENCE_MODE SEED_OFFSET=$SEED_OFFSET"
python validate_bundle.py --targets "$TARGET_CSV" --catalog "$CATALOG_SOURCE"

# {%} = GNU parallel's job-slot number (1..JOBS). Passed through so each
# concurrent slot gets its own private Theano compiledir - see the
# THEANO_FLAGS export in run_one_target.sh. Without this, all JOBS-many
# processes share one compiledir lock; Theano's lock is a polling wait
# (wchan=do_select), not a native futex, so at high concurrency with mostly
# cold cache it serializes almost everything - observed live 2026-07-04:
# JOBS=30 for 12h+ produced 0/592 results, load average ~2.5 on 32 cores,
# ~90% CPU idle, ~27 of 30 processes sleeping on the lock at any moment.
# Build a fresh pending list so resumptions preserve existing FITS files and
# spend compute only on targets that still lack a valid result.
PENDING_FILE="$PROJECT_DIR/pending_${RUN_ID}.csv"
: > "$PENDING_FILE"
while IFS=, read -r target_index target kepid rest; do
  if [ "$target_index" = "target_index" ]; then
    continue
  fi
  target="${target//$'\r'/}"
  kepid="${kepid//$'\r'/}"
  [ -n "$target" ] && [ -n "$kepid" ] || continue
  result="$PROJECT_DIR/Results/$RUN_ID/$target/$target-results.fits"
  if [ ! -s "$result" ]; then
    printf '%s,%s\n' "$target" "$kepid" >> "$PENDING_FILE"
  fi
done < "$TARGET_CSV"

PENDING_COUNT="$(wc -l < "$PENDING_FILE")"
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
JOBLOG="$LOG_DIR/parallel_joblog_${RUN_ID}_${RUN_STAMP}.tsv"
BATCH_STDOUT="$LOG_DIR/batch_${RUN_ID}_${RUN_STAMP}_stdout.log"
BATCH_STDERR="$LOG_DIR/batch_${RUN_ID}_${RUN_STAMP}_stderr.log"

echo "Pending targets: $PENDING_COUNT"
PARALLEL_RC=0
if [ "$PENDING_COUNT" -gt 0 ]; then
  set +e
  parallel --colsep ',' -j "$JOBS" --joblog "$JOBLOG" "$BUNDLE_DIR/run_one_target.sh {1} {2} {%}" \
    < "$PENDING_FILE" > "$BATCH_STDOUT" 2> "$BATCH_STDERR"
  PARALLEL_RC=$?
  set -e
fi

set +e
bash "$BUNDLE_DIR/summarize_progress.sh"
SUMMARY_RC=$?
set -e

if [ "$PARALLEL_RC" -ne 0 ] || [ "$SUMMARY_RC" -ne 0 ]; then
  echo "Batch incomplete. Inspect $PROJECT_DIR/status_manifest_${RUN_ID}.csv and $LOG_DIR/targets/$RUN_ID/." >&2
  exit 1
fi

echo "Batch complete: all requested targets have result FITS under $PROJECT_DIR/Results/$RUN_ID"
