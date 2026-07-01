# Pre-Cloud ALDERAAN Checklist

Current status:

| item | status |
|---|---|
| Existing posterior archive inventoried | done |
| Usable existing posterior baseline separated | done |
| Existing suspicious/refit posterior queue separated | done |
| Missing-posterior queue separated | done |
| Missing-posterior launch-ready subset validated | done |
| Missing-depth blocked subset separated | done |
| GCP bundle for missing launch-ready systems generated | done |
| Local extractor patched for cloud run IDs | done |
| Posterior merge script added | done |
| Local post-cloud script added | done |

Cloud bundle:

`C:\Users\shres\Desktop\HILDALE RESEARCH\sagear_reproduction\cloud_missing_batch`

Ready-to-upload zip:

`C:\Users\shres\Documents\Codex\2026-06-29\8\outputs\cloud_missing_batch_ready_for_gcp.zip`

Bundle contents:

| file | role |
|---|---|
| `targets_missing_launchable.csv` | 592 KOI systems to run |
| `sagear_missing_catalog.csv` | ALDERAAN catalog, 767 seedable catalog rows |
| `shards/targets_shard_*.csv` | four shards for test/resume runs |
| `validate_bundle.py` | local/VM validation |
| `setup_vm.sh` | install Miniforge, ALDERAAN, and environment |
| `run_batch.sh` | parallel ALDERAAN runner |
| `run_one_target.sh` | resumable per-target runner |
| `summarize_progress.sh` | completion counter |
| `pack_results.sh` | tarball packer for return to laptop |
| `README_GCP_MISSING.md` | exact VM/run/retrieve commands |

Counts:

| quantity | count |
|---|---:|
| missing launch-ready planet rows | 718 |
| target systems | 592 |
| ALDERAAN catalog rows after system expansion | 767 |
| thick multi missing planets | 71 |
| thick single missing planets | 120 |
| thin multi missing planets | 193 |
| thin single missing planets | 334 |

Recommended run sequence:

0. Use Google Cloud Shell unless you specifically want to install Cloud SDK locally. Local `gcloud` was not found on this Windows session, while Cloud Shell already includes it.

1. Run one shard first:

```bash
TARGET_CSV=shards/targets_shard_000.csv JOBS=8 bash run_batch.sh
```

2. Retrieve that shard, run local postprocessing, and check posterior extraction.

3. If the shard works, run the full queue:

```bash
JOBS=30 bash run_batch.sh
```

4. Pack results:

```bash
bash pack_results.sh
```

5. Local postprocess after retrieving tarball:

```powershell
cd "C:\Users\shres\Desktop\HILDALE RESEARCH\sagear_reproduction"
.\postprocess_missing_cloud_results.ps1 -TarPath "C:\path\to\alderaan_results_sagear_missing_YYYYMMDD_HHMMSS.tar.gz"
```

Do not run the 320 suspicious existing posteriors first. The priority is to fill the 718 true missing launch-ready rows, then reassess the merged posterior coverage and Rayleigh fits.

After the tarball has been copied out, delete the VM:

```bash
gcloud compute instances delete alderaan-missing-e2-32 --zone $ZONE
```
