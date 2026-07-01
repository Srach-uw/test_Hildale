#!/usr/bin/env bash
set -euo pipefail
RUN_ID="${RUN_ID:-sagear_missing}"
PROJECT_DIR="${PROJECT_DIR:-$PWD/alderaan_project}"
TARGET_CSV="${TARGET_CSV:-targets_missing_launchable.csv}"
TOTAL=$(($(wc -l < "$TARGET_CSV") - 1))
COMPLETE=$(find "$PROJECT_DIR/Results/$RUN_ID" -name '*-results.fits' 2>/dev/null | wc -l || true)
FAILED=$(find "$PROJECT_DIR/status" -name '*.status' -print0 2>/dev/null | xargs -0 grep -l -v '^complete$' 2>/dev/null | wc -l || true)
echo "targets_total,$TOTAL" > progress_summary.csv
echo "results_fits,$COMPLETE" >> progress_summary.csv
echo "failed_or_incomplete_status,$FAILED" >> progress_summary.csv
cat progress_summary.csv
