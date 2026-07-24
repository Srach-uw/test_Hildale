#!/usr/bin/env bash
set -euo pipefail

RUN_ID="${RUN_ID:-sagear_missing}"
PROJECT_DIR="${PROJECT_DIR:-$PWD/alderaan_project}"
TARGET_CSV="${TARGET_CSV:-targets_missing_launchable.csv}"
SUMMARY_FILE="$PROJECT_DIR/progress_summary_${RUN_ID}.csv"
MANIFEST_FILE="$PROJECT_DIR/status_manifest_${RUN_ID}.csv"

if [ ! -s "$TARGET_CSV" ]; then
  echo "Target CSV is missing or empty: $TARGET_CSV" >&2
  exit 2
fi

mkdir -p "$PROJECT_DIR"
printf 'target,result_present,status,outcome\n' > "$MANIFEST_FILE"

TOTAL=0
COMPLETE=0
NEVER_STARTED=0
STALE_RUNNING=0
INTERRUPTED=0
FAILED=0
MISSING_LIGHTCURVE=0
MISSING_RESULTS=0
INCONSISTENT=0
OTHER=0

while IFS=, read -r target_index target kepid rest; do
  if [ "$target_index" = "target_index" ]; then
    continue
  fi
  target="${target//$'\r'/}"
  [ -n "$target" ] || continue
  TOTAL=$((TOTAL + 1))
  result="$PROJECT_DIR/Results/$RUN_ID/$target/$target-results.fits"
  status_file="$PROJECT_DIR/status/$target.status"
  status=""
  if [ -f "$status_file" ]; then
    status="$(<"$status_file")"
  fi

  if [ -s "$result" ]; then
    COMPLETE=$((COMPLETE + 1))
    outcome="complete"
    result_present="yes"
  else
    result_present="no"
    case "$status" in
      '') outcome="never_started"; NEVER_STARTED=$((NEVER_STARTED + 1)) ;;
      running) outcome="stale_running"; STALE_RUNNING=$((STALE_RUNNING + 1)) ;;
      interrupted_*) outcome="$status"; INTERRUPTED=$((INTERRUPTED + 1)) ;;
      failed_exit_*|failed_*_exit_*) outcome="$status"; FAILED=$((FAILED + 1)) ;;
      missing_lightcurve) outcome="$status"; MISSING_LIGHTCURVE=$((MISSING_LIGHTCURVE + 1)) ;;
      missing_results) outcome="$status"; MISSING_RESULTS=$((MISSING_RESULTS + 1)) ;;
      complete) outcome="status_complete_but_result_missing"; INCONSISTENT=$((INCONSISTENT + 1)) ;;
      *) outcome="other_$status"; OTHER=$((OTHER + 1)) ;;
    esac
  fi
  printf '%s,%s,%s,%s\n' "$target" "$result_present" "${status:-none}" "$outcome" >> "$MANIFEST_FILE"
done < "$TARGET_CSV"

INCOMPLETE=$((TOTAL - COMPLETE))
{
  echo "targets_total,$TOTAL"
  echo "results_fits,$COMPLETE"
  echo "incomplete_total,$INCOMPLETE"
  echo "never_started,$NEVER_STARTED"
  echo "stale_running,$STALE_RUNNING"
  echo "interrupted,$INTERRUPTED"
  echo "failed_exit,$FAILED"
  echo "missing_lightcurve,$MISSING_LIGHTCURVE"
  echo "missing_results,$MISSING_RESULTS"
  echo "status_result_inconsistent,$INCONSISTENT"
  echo "other_status,$OTHER"
  echo "manifest,$MANIFEST_FILE"
} > "$SUMMARY_FILE"
cat "$SUMMARY_FILE"

if [ "$INCOMPLETE" -ne 0 ]; then
  exit 1
fi
