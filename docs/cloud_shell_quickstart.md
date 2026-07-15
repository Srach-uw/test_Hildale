# Cloud Shell Quickstart For Missing ALDERAAN Runs

Use this route first. It avoids local Windows `gcloud` install/auth issues.

Before creating the VM:

1. In Google Cloud Console, open Billing and confirm you are using an active Free Trial / credit-covered billing account if you want no out-of-pocket spend.
2. Create a budget alert for the project, for example `$5` with alerts at 50%, 90%, and 100%.
3. Do not continue if the billing account says it is upgraded to paid and you are not comfortable with real charges.

Local file to upload:

`cloud_missing_batch_ready_for_gcp.zip`

In Google Cloud Shell:

```bash
unzip cloud_missing_batch_ready_for_gcp.zip -d sagear_cloud_missing
cd sagear_cloud_missing
python3 validate_bundle.py

export PROJECT_ID=your-gcp-project
export ZONE=us-central1-a
export MAX_RUN_DURATION=20h
bash create_gcp_spot_vm.sh

gcloud compute scp --recurse . alderaan-missing-e2-32:~/sagear_cloud_missing --zone $ZONE
gcloud compute ssh alderaan-missing-e2-32 --zone $ZONE
```

On the VM, run a test shard first:

```bash
cd ~/sagear_cloud_missing
bash setup_vm.sh
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alderaan
python validate_bundle.py
TARGET_CSV=shards/targets_shard_000.csv JOBS=8 bash run_batch.sh
bash pack_results.sh
```

Copy the test tarball back, postprocess locally, and inspect coverage. If it works, return to the VM and run the full queue:

```bash
cd ~/sagear_cloud_missing
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alderaan
JOBS=30 bash run_batch.sh
bash pack_results.sh
```

After results are copied out, delete the VM:

```bash
gcloud compute instances delete alderaan-missing-e2-32 --zone $ZONE
```

Cost-control rule: do not start the 320 existing-but-flagged posterior reruns until the 718 truly missing launch-ready posteriors have been filled and merged.
