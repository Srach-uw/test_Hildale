#!/usr/bin/env bash
set -euo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

cat > "$TMP_DIR/targets.csv" <<'EOF'
target_index,koi_target,kepid
0,K00001,1
1,K00002,2
2,K00003,3
EOF
mkdir -p "$TMP_DIR/project/Results/test/K00001" "$TMP_DIR/project/status"
printf 'fits\n' > "$TMP_DIR/project/Results/test/K00001/K00001-results.fits"
printf 'running\n' > "$TMP_DIR/project/status/K00002.status"

set +e
RUN_ID=test PROJECT_DIR="$TMP_DIR/project" TARGET_CSV="$TMP_DIR/targets.csv" \
  bash "$BUNDLE_DIR/summarize_progress.sh" > "$TMP_DIR/summary.log"
RC=$?
set -e

[ "$RC" -eq 1 ]
grep -qx 'targets_total,3' "$TMP_DIR/project/progress_summary_test.csv"
grep -qx 'results_fits,1' "$TMP_DIR/project/progress_summary_test.csv"
grep -qx 'incomplete_total,2' "$TMP_DIR/project/progress_summary_test.csv"
grep -qx 'stale_running,1' "$TMP_DIR/project/progress_summary_test.csv"
grep -qx 'never_started,1' "$TMP_DIR/project/progress_summary_test.csv"
grep -qx 'K00002,no,running,stale_running' "$TMP_DIR/project/status_manifest_test.csv"
grep -qx 'K00003,no,none,never_started' "$TMP_DIR/project/status_manifest_test.csv"

echo 'runner reliability tests passed'
