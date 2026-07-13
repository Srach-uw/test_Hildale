# GCP ALDERAAN Missing-Posterior Bundle

This bundle is for the true missing-posterior queue only.

- Run id: `sagear_missing`
- Runnable target systems: `592`
- Missing planet rows covered by those targets: `718`
- ALDERAAN catalog rows after whole-system expansion: `767`
- Default parallel jobs: `30`
- VM runtime cap in `create_gcp_spot_vm.sh`: `20h` by default

## Local Sanity Check

```bash
python validate_bundle.py
```

## Billing Safety First

Do not create the VM until you have checked billing.

- If your account is an active Google Cloud Free Trial, Compute Engine usage should be charged against the trial credit.
- If your account has been upgraded to paid billing, the VM can charge your card.
- Budgets/alerts are useful warnings, but they are not a hard spending cap.
- The VM script includes `--max-run-duration 20h`, so compute stops automatically after 20 hours if you forget. The disk can still incur a small charge until the VM is deleted.
- Delete the VM after copying out results.

## Recommended Route: Google Cloud Shell

Use this if `gcloud` is not installed or authenticated on your laptop. Cloud Shell already has `gcloud`, so it avoids local Windows auth/PATH problems.

1. Open Google Cloud Console, select your billing-enabled project, then open Cloud Shell.
2. Upload `cloud_missing_batch_ready_for_gcp.zip` to Cloud Shell.
3. In Cloud Shell:

```bash
unzip cloud_missing_batch_ready_for_gcp.zip -d sagear_cloud_missing
cd sagear_cloud_missing
python3 validate_bundle.py
export PROJECT_ID=your-gcp-project
export ZONE=us-central1-a
bash create_gcp_spot_vm.sh
gcloud compute scp --recurse . alderaan-missing-e2-32:~/sagear_cloud_missing --zone $ZONE
gcloud compute ssh alderaan-missing-e2-32 --zone $ZONE
```

## Create VM

```bash
export PROJECT_ID=your-gcp-project
export ZONE=us-central1-a
bash create_gcp_spot_vm.sh
gcloud compute scp --recurse . alderaan-missing-e2-32:~/sagear_cloud_missing --zone $ZONE
gcloud compute ssh alderaan-missing-e2-32 --zone $ZONE
```

## Run On VM

```bash
cd ~/sagear_cloud_missing
bash setup_vm.sh
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alderaan
python validate_bundle.py
JOBS=30 bash run_batch.sh
bash pack_results.sh
```

For a test slice:

```bash
TARGET_CSV=shards/targets_shard_000.csv JOBS=8 bash run_batch.sh
```

## Retrieve

```bash
gcloud compute scp alderaan-missing-e2-32:~/sagear_cloud_missing/alderaan_results_sagear_missing_*.tar.gz . --zone $ZONE
```

Extract locally into `sagear_reproduction/alderaan_project` or a separate results folder, then run:

```powershell
python sagear_reproduction\extract_eccentricity_posteriors.py `
  --sample sagear_reproduction\outputs\canonical_sample_old_astropy_rawcc.csv `
  --run-id sagear_missing `
  --posterior-subdir eccentricity_posteriors_sagear_missing `
  --summary-out sagear_reproduction\outputs\eccentricity_posterior_summary_sagear_missing.csv `
  --coverage-out sagear_reproduction\outputs\eccentricity_posterior_coverage_sagear_missing.csv

python sagear_reproduction\merge_posterior_summaries.py `
  --new sagear_reproduction\outputs\eccentricity_posterior_summary_sagear_missing.csv `
  --out sagear_reproduction\outputs\eccentricity_posterior_summary_merged_sagear_missing.csv `
  --coverage-out sagear_reproduction\outputs\eccentricity_posterior_coverage_merged_sagear_missing.csv
```

## Stop Costs Afterward

After results are packed and copied out, delete the VM:

```bash
gcloud compute instances delete alderaan-missing-e2-32 --zone $ZONE
```

If you used a temporary Cloud Storage bucket, delete the result tarball there too after downloading it.
