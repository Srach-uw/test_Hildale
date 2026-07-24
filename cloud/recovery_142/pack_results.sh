#!/usr/bin/env bash
set -euo pipefail
RUN_ID="${RUN_ID:-sagear_missing}"
PROJECT_DIR="${PROJECT_DIR:-$PWD/alderaan_project}"
OUT="alderaan_results_${RUN_ID}_$(date +%Y%m%d_%H%M%S).tar.gz"
set +e
bash summarize_progress.sh
SUMMARY_RC=$?
set -e

# Only tar paths that actually exist so a missing optional file can't abort the
# pack after a long run (set -e would otherwise fail the whole tarball).
project_items=()
for item in \
  "Results/$RUN_ID" \
  status \
  "progress_summary_${RUN_ID}.csv" \
  "status_manifest_${RUN_ID}.csv"
do
  [ -e "$PROJECT_DIR/$item" ] && project_items+=("$item")
done
pwd_items=()
for item in \
  logs \
  provenance \
  targets_missing_launchable.csv \
  sagear_missing_catalog.csv \
  published_inventory_missing_population_rows.csv \
  published_inventory_missing_full_system_inventory.csv \
  bundle_manifest_sha256.csv \
  target_shards \
  README.md
do
  [ -e "$PWD/$item" ] && pwd_items+=("$item")
done
[ -e "$PWD/shard_id.txt" ] && pwd_items+=("shard_id.txt")

tar -czf "$OUT" \
  ${project_items:+-C "$PROJECT_DIR" "${project_items[@]}"} \
  ${pwd_items:+-C "$PWD" "${pwd_items[@]}"}
echo "$OUT"
if [ "$SUMMARY_RC" -ne 0 ]; then
  echo "Archive created with incomplete targets; inspect the bundled status manifest and logs." >&2
fi
