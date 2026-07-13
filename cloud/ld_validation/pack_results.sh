#!/usr/bin/env bash
set -euo pipefail
RUN_ID="${RUN_ID:-sagear_missing}"
PROJECT_DIR="${PROJECT_DIR:-$PWD/alderaan_project}"
OUT="alderaan_results_${RUN_ID}_$(date +%Y%m%d_%H%M%S).tar.gz"
bash summarize_progress.sh

# Only tar paths that actually exist so a missing optional file can't abort the
# pack after a long run (set -e would otherwise fail the whole tarball).
project_items=()
for item in "Results/$RUN_ID" status; do
  [ -e "$PROJECT_DIR/$item" ] && project_items+=("$item")
done
pwd_items=()
for item in logs progress_summary.csv targets_missing_launchable.csv sagear_missing_catalog.csv cloud_missing_manifest.md; do
  [ -e "$PWD/$item" ] && pwd_items+=("$item")
done

tar -czf "$OUT" \
  ${project_items:+-C "$PROJECT_DIR" "${project_items[@]}"} \
  ${pwd_items:+-C "$PWD" "${pwd_items[@]}"}
echo "$OUT"
