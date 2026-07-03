# Full 592-Target ALDERAAN Run — Runbook

Written 2026-07-03, after the 5-target smoke test validated the pipeline end to end
(K00179: all three stages clean, results.fits produced, zero tracebacks after 7 fixes).

## What the smoke test taught us (baked into this plan)

| lesson | consequence here |
|---|---|
| 7 bugs found live (wget redirect, PYTHONPATH, Holczer, validate-import, CRLF, analyze args, conda activation) | use ONLY the current fixed bundle zip — no manual VM patching needed |
| Spot VMs preempted 3× in one evening (us-central1-a stockout, -b churn) | use a STANDARD (on-demand) VM: ~$1.07/hr for e2-standard-32, total ~$50–70 — reliability worth it against $300 credit |
| detrend peaks ~3.8 GB/process | JOBS=24 on 128 GB (not 30) for OOM headroom |
| theano JIT-compiles for ~20–30 min on first run per fresh VM | warm-up: run 1 target alone first, then launch the full queue (finished target auto-skips) |
| stage-3 nested sampling: 30 min – 2 h+ per target, single-core, no resume | preempted/killed targets restart stage 3 from scratch — another reason for STANDARD |
| per-target wall time ~1–2.5 h (K00179 60 min, K00277 2 h+) | 592 × ~1.75 h ÷ 24 ≈ 43 h wall ≈ 2 days; cap at 72h for cost safety |

## Phase 0 — prerequisites (Cloud Shell)

1. Quota check (free trial often caps region vCPUs):
   ```bash
   gcloud compute regions describe us-central1 --format="table(quotas.filter(metric='CPUS'))"
   ```
   - limit ≥ 32 → proceed as written.
   - limit = 8 → either request a quota increase (Console → IAM → Quotas) or fall back to
     e2-standard-8 + JOBS=6 (slow: ~1 week wall; quota increase strongly preferred).
2. Billing sanity: confirm free-trial credit remaining covers ~$80 worst case.
3. Upload the CURRENT fixed bundle (`cloud_missing_batch_ready_for_gcp_FIXED.zip`,
   all 7 fixes baked in) to Cloud Shell; `rm -rf` any stale unzipped copy first.

## Phase 1 — create the VM (STANDARD, not Spot)

Chosen config (2026-07-03): a single **n2-standard-64** (64 vCPU, 256 GB) at JOBS=56 —
one VM to babysit, ~overnight (~18 h), ~$52. Cost is trivial vs the $300 credit; the goal
is minimizing operational hassle (one SSH session, no target-list sharding). Region CPU
quota confirmed at 200 with only 8 in use, so 64 fits easily.

```bash
export ZONE=us-central1-b
export INSTANCE=alderaan-full-n2-64

gcloud compute instances create "$INSTANCE" \
  --zone "$ZONE" \
  --machine-type n2-standard-64 \
  --provisioning-model=STANDARD \
  --max-run-duration=72h \
  --instance-termination-action=STOP \
  --boot-disk-size 200GB \
  --boot-disk-type pd-ssd \
  --image-family ubuntu-2204-lts \
  --image-project ubuntu-os-cloud
```

Zone stockout → try us-central1-c, then -f. The 72h cap is a cost backstop, not a target.
Boot disk 200 GB (vs 150) because 56 concurrent targets stage more light curves at once.

## Phase 2 — stage and set up

```bash
cd ~/sagear_cloud_missing   # fresh unzip of the FIXED bundle
gcloud compute scp --recurse . "$INSTANCE":~/sagear_cloud_missing --zone "$ZONE"
gcloud compute ssh "$INSTANCE" --zone "$ZONE"
# on the VM:
cd ~/sagear_cloud_missing
bash setup_vm.sh            # ~15 min; auto-applies the ALDERAAN validate-import patch
```

## Phase 3 — warm-up (one target, populates theano cache)

```bash
source ~/miniforge3/etc/profile.d/conda.sh && conda activate alderaan
python validate_bundle.py
head -n 2 targets_missing_launchable.csv > warmup_1.csv
TARGET_CSV=warmup_1.csv JOBS=1 nohup bash run_batch.sh > warmup.log 2>&1 &
disown
```

Wait for it to complete (~1–1.5 h; check `tail warmup.log` / joblog). One clean
results.fits here = full-stack confirmation on THIS VM. It auto-skips in the full run.

## Phase 4 — launch the full queue

```bash
cd ~/sagear_cloud_missing
JOBS=56 nohup bash run_batch.sh > full_run.log 2>&1 &
disown
echo "PID $!"
```

(run_batch.sh self-activates conda since the bug-7 fix — a fresh SSH session is fine.)
JOBS=56 on 64 vCPU / 256 GB: leaves 8 cores for OS + GNU parallel overhead; steady-state
memory ~90–120 GB (detrend spikes are the short stage), comfortably under 256 GB. The
Phase-3 warm-up is REQUIRED before this — without a populated theano cache, 56 targets
would cold-compile simultaneously and thrash the compile lock.

## Phase 5 — monitor (paste as one block, ~2–3×/day)

```bash
cd ~/sagear_cloud_missing
date
MAIN_PID=$(pgrep -f "bash run_batch.sh" | head -1)
[ -n "$MAIN_PID" ] && ps -p "$MAIN_PID" -o pid,etime,cmd || echo "BATCH NOT RUNNING"
echo "done/failed so far:"; awk 'NR>1{print $7}' logs/parallel_joblog_sagear_missing.tsv | sort | uniq -c
echo "results: $(find alderaan_project/Results -name '*-results.fits' | wc -l) of 592"
echo "active stages: $(pgrep -fc 'detrend_and_estimate|analyze_autocorrelated|fit_transit_shape')"
echo "tracebacks: $(grep -c Traceback logs/batch_sagear_missing_stderr.log 2>/dev/null)"
```

Expectations:
- Some targets WILL fail (celerite2 LinAlgError etc.) — genuine numerical attrition,
  not infrastructure. Collect, don't chase, until the second pass.
- Nonzero exitvals in the joblog are normal; what matters is the results count climbing.

## Phase 6 — second pass (after first pass drains)

Rerun the exact same command. Completed targets skip instantly (results.fits check);
only failures rerun. Stragglers still failing after pass 2 → export the failed list from
the joblog and triage separately — do not hold the big VM open for a handful of targets.

## Phase 7 — pack, retrieve, DELETE

```bash
bash pack_results.sh          # on the VM
exit
gcloud compute scp "$INSTANCE":~/sagear_cloud_missing/alderaan_results_sagear_missing_*.tar.gz . --zone "$ZONE"
gcloud compute instances delete "$INSTANCE" --zone "$ZONE"
gcloud compute instances list  # verify GONE (stopped VM still bills its 150 GB disk)
```

Download the tarball via Cloud Shell ⋮ → Download, then run local postprocessing
(`postprocess_missing_cloud_results.ps1`) and merge. Flag K01316/K06516 results for
review (incomplete-system fits, see README).

## Recovery drills (any failure, any time)

- Cloud Shell disconnect → reconnect, `gcloud config set project ...`, SSH back. Job unaffected.
- SSH refused → wait 30 s, retry; check `instances describe --format="get(status)"`.
- VM stopped (72h cap) → `instances start`, SSH, rerun Phase 4 command. Completed targets skip.
- Batch not running after reconnect → just rerun the Phase 4 command (self-activating, resume-safe).
