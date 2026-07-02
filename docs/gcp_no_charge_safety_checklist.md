# GCP No-Surprise-Cost Checklist

Use this before launching ALDERAAN on Google Cloud.

## Zero-Out-Of-Pocket Requirement

You should only run the VM if one of these is true:

- You have an active Google Cloud Free Trial with enough remaining credit.
- You are using a university/lab billing account or grant where charges are expected.
- You explicitly accept that the VM is billable and might charge a card.

If none of those is true, do not create the VM.

## Rough Cost Estimate (2026-07-02, verify current prices)

Assumptions: `e2-standard-32` Spot in `us-central1` ≈ $0.30–0.45/hr; 150 GB pd-ssd
≈ $0.9/day; ALDERAAN ≈ 20–90 min per target single-core; `JOBS=30`.

- 592 targets ≈ 15–40 VM-hours → **roughly $6–20 compute + $1–2 disk**.
- One test shard (148 targets, `JOBS=8`) ≈ **$5–10**.
- Well within a $300 free trial, but preemptions can stretch wall time; the
  `--max-run-duration 20h` cap bounds each VM session regardless.

Flag-compatibility note (verified against Google Cloud docs, 2026-07-02):
`--max-run-duration` is supported with `--provisioning-model=SPOT` and
`--instance-termination-action=STOP` (duration 30 s – 120 d). The
`--discard-local-ssds-at-termination-timestamp` requirement applies only to
local SSDs; this VM uses a persistent boot disk, so it does not apply.

## Before Creating The VM

1. Open Google Cloud Console.
2. Go to Billing.
3. Confirm whether the billing account is Free Trial / credit-covered or paid.
4. Check remaining credit.
5. Create a project-specific budget alert, for example `$5`, with alerts at 50%, 90%, and 100%.
6. Confirm the VM script has a runtime cap:

```bash
grep max-run-duration create_gcp_spot_vm.sh
```

It should show:

```text
--max-run-duration "$MAX_RUN_DURATION"
```

## Safer Test Run

Run only one shard first:

```bash
TARGET_CSV=shards/targets_shard_000.csv JOBS=8 bash run_batch.sh
```

Do not run the full queue until the shard returns and the postprocessing works locally.

## After The Run

Copy results out, then delete the VM:

```bash
gcloud compute instances delete alderaan-missing-e2-32 --zone $ZONE
```

Then check Billing again.

## Important Caveat

Budget alerts notify you; they do not guarantee a hard cap. The practical safeguards are:

- use active free credits,
- run a shard first,
- set `MAX_RUN_DURATION`,
- delete the VM,
- do not leave disks or buckets sitting around.
