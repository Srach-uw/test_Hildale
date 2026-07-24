# Exact published-inventory ALDERAAN completion batch

This bundle covers every published-host system that still lacks a usable
ALDERAAN result in the local archive plus the completed cloud run.

- target systems: 142
- population planets missing posteriors: 174
- full-system ALDERAAN catalog rows: 181
- full-system inventory rows including unseedable KOIs: 181
- cadence mode: long plus available short cadence
- ALDERAAN commit: 7443dff16b7f9092e14a6f0cc1f8948d457c9e0b
- runner: repaired resumable target-specific logging and deterministic seed

Multiplicity and population labels come from the exact published host
inventory. The transit fit includes every non-false-positive, seedable
KOI companion in each selected system before the 1-100 day population cut.

## Target breakdown

| disk | system | prior launch status | targets |
|---|---|---|---:|
| thick | multi | not_in_previous_592_target_launch | 3 |
| thick | multi | previously_launched_no_usable_result | 2 |
| thick | single | not_in_previous_592_target_launch | 9 |
| thick | single | previously_launched_no_usable_result | 9 |
| thin | multi | not_in_previous_592_target_launch | 24 |
| thin | multi | previously_launched_no_usable_result | 2 |
| thin | single | not_in_previous_592_target_launch | 80 |
| thin | single | previously_launched_no_usable_result | 13 |

## Run

After setup and environment activation:

```bash
JOBS=20 nohup bash run_recovery_two_pass.sh > published_inventory_missing.log 2>&1 &
```

For a multi-VM run, assign exactly one shard to each VM:

```bash
nohup bash run_shard.sh 00 28 > shard_00.log 2>&1 &
```

The target shards are deterministic, mutually exclusive, and their
union is validated against the complete 142-target master list.

This runs the canonical seed, retries only unresolved targets with
`SEED_OFFSET=1`, and packages all results and failure evidence.
The run is resumable. Existing nonempty result FITS are skipped.
Failures are recorded with their active stage (`download`, `detrend`,
`noise`, or `transit_fit`) and are never silently replaced by a
different noise or transit model.

A successful retry records its seed offset in the target provenance.
Do not alter GP jitter, priors, or the transit model in the canonical
branch. Any such numerical rescue belongs in a separate sensitivity
run.

For a manual resume or extra snapshot, package the current results
at any time with:

```bash
RUN_ID=sagear_published_inventory_missing bash pack_results.sh
```

On Windows, process the downloaded archive from the canonical
`sagear_reproduction` directory with:

```powershell
python .\postprocess_published_inventory_recovery.py --archive "C:\path\to\alderaan_results.tar.gz"
```

The postprocessor withholds the final hierarchy until the exact
published-inventory coverage gate is complete.
