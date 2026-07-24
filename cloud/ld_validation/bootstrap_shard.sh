#!/usr/bin/env bash
set -euo pipefail

SHARD_ID="${1:?usage: bash bootstrap_shard.sh SHARD_ID [JOBS]}"
JOBS="${2:-20}"
printf -v SHARD_PADDED "%02d" "$((10#$SHARD_ID))"
LOG="$PWD/bootstrap_shard_${SHARD_PADDED}.log"
FINAL_ARCHIVE="$HOME/alderaan_shard_${SHARD_PADDED}_results.tar.gz"

exec > >(tee -a "$LOG") 2>&1

publish_status() {
  curl -sf -X PUT \
    -H "Metadata-Flavor: Google" \
    --data-binary "$1" \
    "http://metadata.google.internal/computeMetadata/v1/instance/guest-attributes/alderaan/status" \
    || true
}

stop_vm() {
  sync
  echo "Stopping VM; persistent boot disk and partial results are retained."
  sudo shutdown -h now || true
}
trap stop_vm EXIT

echo "$SHARD_PADDED" > shard_id.txt
echo "Bootstrap started at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
publish_status "setup"

bash setup_vm.sh
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate alderaan
publish_status "running"

set +e
bash run_shard.sh "$SHARD_PADDED" "$JOBS"
RUN_RC=$?
set -e

LATEST_ARCHIVE="$(ls -1t alderaan_results_sagear_published_inventory_missing_*.tar.gz 2>/dev/null | head -1 || true)"
if [ -z "$LATEST_ARCHIVE" ] || [ ! -s "$LATEST_ARCHIVE" ]; then
  echo "No packaged shard archive was produced; inspect $LOG on the retained disk." >&2
  exit 3
fi

cp "$LATEST_ARCHIVE" "$FINAL_ARCHIVE"
sha256sum "$FINAL_ARCHIVE" | tee "$FINAL_ARCHIVE.sha256"
publish_status "archived_run_rc_${RUN_RC}"
echo "SHARD_RUN_RC=$RUN_RC"
echo "FINAL_ARCHIVE=$FINAL_ARCHIVE"
echo "Bootstrap finished at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
