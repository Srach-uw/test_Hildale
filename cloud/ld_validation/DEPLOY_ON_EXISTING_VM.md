# Deploy On The Existing GCP VM

The existing `alderaan-factorial` VM is stopped in `us-central1-b`. Its
60 GB persistent disk already contains the working ALDERAAN environment and
the completed factorial validation. The completion batch uses a separate
project directory and does not overwrite those results.

Starting or resizing the VM is billable. Confirm the remaining promotional
credit in the Google Cloud Billing Credits page before running these commands.
The 24-hour maximum runtime is a hard compute stop, unlike a budget alert.

## 1. Local Windows PowerShell

```powershell
gcloud.cmd config set project project-7f7ff467-5d61-4072-8f4

gcloud.cmd compute instances set-machine-type alderaan-factorial `
  --zone=us-central1-b `
  --machine-type=e2-standard-32

gcloud.cmd compute instances set-scheduling alderaan-factorial `
  --zone=us-central1-b `
  --max-run-duration=24h `
  --instance-termination-action=STOP

gcloud.cmd compute instances start alderaan-factorial `
  --zone=us-central1-b

Start-Sleep -Seconds 75

gcloud.cmd compute scp `
  "<local_bundle_path>\Hildale_ALDERAAN_Published_Inventory_Missing_142_20260724_v4.zip" `
  "shreshth_rach1@alderaan-factorial:~/" `
  --zone=us-central1-b

gcloud.cmd compute ssh shreshth_rach1@alderaan-factorial `
  --zone=us-central1-b
```

## 2. Inside The VM

```bash
cd ~
unzip -q -o Hildale_ALDERAAN_Published_Inventory_Missing_142_20260724_v4.zip
cd ~/cloud_published_inventory_missing_batch

source ~/miniforge3/etc/profile.d/conda.sh
conda activate alderaan

python validate_bundle.py \
  --targets targets_missing_launchable.csv \
  --catalog sagear_missing_catalog.csv

RUN_ID=sagear_published_inventory_missing \
CADENCE_MODE=both \
JOBS=20 \
nohup bash run_recovery_two_pass.sh > published_inventory_missing.log 2>&1 &

echo "STARTED_PID=$!"
sleep 10
tail -n 40 published_inventory_missing.log
```

Expected validation:

```text
VALIDATION OK: 142 targets, 181 catalog rows
```

## 3. Reconnect And Check

From Cloud Shell:

```bash
gcloud config set project project-7f7ff467-5d61-4072-8f4
gcloud compute ssh shreshth_rach1@alderaan-factorial \
  --zone=us-central1-b
```

Inside the VM:

```bash
cd ~/cloud_published_inventory_missing_batch

echo "ACTIVE:"
pgrep -af "run_recovery_two_pass|run_batch|run_one_target|parallel|detrend|analyze_autocorrelated|fit_transit" || true

echo
RUN_ID=sagear_published_inventory_missing \
bash summarize_progress.sh || true

echo
tail -n 50 published_inventory_missing.log
```

The wrapper runs the canonical deterministic seed first, retries only unresolved
targets with one changed seed, and then packages all available results, statuses,
logs, and provenance. The runner skips every existing nonempty result FITS, so
restarting after an automatic stop resumes rather than repeating completed
targets. Packaging still completes when scientifically difficult targets remain
unresolved; inspect the bundled status manifest and target logs before deciding
whether any sensitivity rerun is justified.
