#!/usr/bin/env bash
set -euo pipefail

RUN_ID="${RUN_ID:-sagear_published_inventory_missing}"
JOBS="${JOBS:-20}"
CADENCE_MODE="${CADENCE_MODE:-both}"
export RUN_ID JOBS CADENCE_MODE

echo "Canonical pass: RUN_ID=$RUN_ID JOBS=$JOBS CADENCE_MODE=$CADENCE_MODE SEED_OFFSET=0"
SEED_OFFSET=0 bash run_batch.sh || true

echo "Changed-seed retry for targets still lacking a nonempty result FITS"
SEED_OFFSET=1 bash run_batch.sh || true

echo "Packaging every result, status, provenance record, and target log"
bash pack_results.sh
