# Deploy the reliability repair

This revision is designed for the existing `alderaan-factorial` VM. It does
not replace `projects/`, downloaded light curves, or any valid ALDERAAN FITS.
It updates only the runner scripts and documentation.

## 1. Extend the automatic stop window

In Cloud Shell, before starting new work:

```bash
gcloud config set project <gcp-project-id>
gcloud compute instances set-scheduling alderaan-factorial \
  --zone us-central1-b \
  --max-run-duration=48h \
  --instance-termination-action=STOP
```

## 2. Upload the archive

In local PowerShell:

```powershell
gcloud.cmd compute scp "<local_bundle_path>\Hildale_ALDERAAN_Factorial_Validation_Reliability_Repair_20260712.zip" <vm-linux-user>@alderaan-factorial:~ --zone us-central1-b
```

## 3. Install only the repaired files

SSH to the VM, then run this whole block. It backs up the live scripts before
copying the repaired versions and leaves `projects/` unchanged.

```bash
cd ~
LIVE="$HOME/sagear_ld_validation_batch"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
REPAIR_ROOT="$HOME/sagear_ld_validation_batch_repair_$STAMP"
unzip -q Hildale_ALDERAAN_Factorial_Validation_Reliability_Repair_20260712.zip -d "$REPAIR_ROOT"
REPAIR="$REPAIR_ROOT/sagear_ld_validation_batch"
BACKUP="$LIVE/runner_backup_$STAMP"
mkdir -p "$BACKUP"

for f in run_batch.sh run_one_target.sh summarize_progress.sh run_validation_matrix.sh README_LD_VALIDATION.md test_runner_reliability.sh DEPLOY_RELIABILITY_REPAIR.md
do
  [ -f "$LIVE/$f" ] && cp -a "$LIVE/$f" "$BACKUP/$f"
  cp -a "$REPAIR/$f" "$LIVE/$f"
done

cd "$LIVE"
chmod +x run_batch.sh run_one_target.sh summarize_progress.sh run_validation_matrix.sh test_runner_reliability.sh
bash -n run_batch.sh run_one_target.sh summarize_progress.sh run_validation_matrix.sh test_runner_reliability.sh
bash test_runner_reliability.sh
```

Expected final line: `runner reliability tests passed`.

## 4. Inspect The Existing Partial Arms

```bash
cd ~/sagear_ld_validation_batch

RUN_ID=sagear_validation_original_lc \
PROJECT_DIR="$PWD/projects/original_lc" \
TARGET_CSV=targets_ld_reference_validation.csv \
bash summarize_progress.sh || true

RUN_ID=sagear_validation_reference_lc \
PROJECT_DIR="$PWD/projects/reference_lc" \
TARGET_CSV=targets_ld_reference_validation.csv \
bash summarize_progress.sh || true
```

The old partial state is expected to be incomplete. The new manifests identify
every target as complete, stale-running, or never-started.

## 5. Resume The Full Matrix

```bash
cd ~/sagear_ld_validation_batch
JOBS=6 nohup bash run_validation_matrix.sh > validation_matrix_repaired.log 2>&1 &
echo "Started PID: $!"
```

The runner will skip targets that already have a valid result FITS and run only
unresolved targets. It will continue into later arms even if a target fails,
then exit nonzero at the end if any arm remains incomplete.

## 6. Monitor

```bash
cd ~/sagear_ld_validation_batch
tail -n 80 validation_matrix_repaired.log
find projects -name '*-results.fits' | wc -l
find projects -name 'status_manifest_*.csv' -print
```

For a specific numerical failure, inspect:

```bash
find logs/targets -name '*.stderr.log' -print
```

Do not delete old result FITS or status files. A rerun overwrites a stale
`running` status only for the target being resumed, and preserves all valid
results.
