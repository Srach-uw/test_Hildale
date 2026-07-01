#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_ID:?Set PROJECT_ID first}"
ZONE="${ZONE:-us-central1-a}"
INSTANCE="${INSTANCE:-alderaan-missing-e2-32}"
MACHINE="${MACHINE:-e2-standard-32}"
DISK_SIZE="${DISK_SIZE:-150GB}"
MAX_RUN_DURATION="${MAX_RUN_DURATION:-20h}"
IMAGE_FAMILY="${IMAGE_FAMILY:-ubuntu-2204-lts}"
IMAGE_PROJECT="${IMAGE_PROJECT:-ubuntu-os-cloud}"

gcloud config set project "$PROJECT_ID"
gcloud compute instances create "$INSTANCE"   --zone "$ZONE"   --machine-type "$MACHINE"   --provisioning-model=SPOT   --instance-termination-action=STOP   --max-run-duration "$MAX_RUN_DURATION"   --boot-disk-size "$DISK_SIZE"   --boot-disk-type pd-ssd   --image-family "$IMAGE_FAMILY"   --image-project "$IMAGE_PROJECT"   --metadata=google-logging-enabled=true

echo "Created $INSTANCE in $ZONE"
echo "Runtime cap: $MAX_RUN_DURATION, then VM stops automatically. Delete it after copying results."
echo "Copy bundle with:"
echo "gcloud compute scp --recurse . $INSTANCE:~/sagear_cloud_missing --zone $ZONE"
