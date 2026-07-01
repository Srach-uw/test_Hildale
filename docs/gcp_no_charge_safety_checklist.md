# GCP No-Surprise-Cost Checklist

Use this before launching ALDERAAN on Google Cloud.

## Zero-Out-Of-Pocket Requirement

You should only run the VM if one of these is true:

- You have an active Google Cloud Free Trial with enough remaining credit.
- You are using a university/lab billing account or grant where charges are expected.
- You explicitly accept that the VM is billable and might charge a card.

If none of those is true, do not create the VM.

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
