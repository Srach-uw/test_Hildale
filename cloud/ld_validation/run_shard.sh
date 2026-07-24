#!/usr/bin/env bash
set -euo pipefail

SHARD_ID="${1:?usage: bash run_shard.sh SHARD_ID [JOBS]}"
JOBS="${2:-28}"
printf -v SHARD_PADDED "%02d" "$((10#$SHARD_ID))"
TARGET_CSV="target_shards/targets_shard_${SHARD_PADDED}.csv"

if [ ! -s "$TARGET_CSV" ]; then
  echo "Missing shard target list: $TARGET_CSV" >&2
  exit 2
fi

export TARGET_CSV JOBS
export RUN_ID="sagear_published_inventory_missing"
export PROJECT_DIR="$PWD/alderaan_project"
export CADENCE_MODE="both"
export SHARD_ID="$SHARD_PADDED"

echo "Starting shard $SHARD_PADDED with $(($(wc -l < "$TARGET_CSV") - 1)) targets and JOBS=$JOBS"
bash run_recovery_two_pass.sh
