# Claude Review Handoff: Sagear / Hilldale Replication

Date: 2026-07-01

Purpose: independently audit Codex's replication pipeline work, catch mistakes, and then assemble only the correct pieces into a new clean GitHub repository.

## Instruction For Claude

This zip is intended to be self-contained. The user should not need to paste any extra prompt. Start here, read the audit context, inspect the `new_repo_candidate/` folder, and then tell the user what is correct, what is wrong, and what should go into a new GitHub repo.

Do not trust Codex's conclusions by default. Rerun/inspect the files in this bundle.

Do not modify any previous GitHub repo from earlier attempts. If GitHub work is requested, use a **new clean repository** based on `new_repo_candidate/`.

## Critical Instruction

Do not treat this as a finished or trusted implementation. Treat it as a research-engineering audit package:

- verify counts from source CSVs,
- rerun validation commands,
- inspect queue logic,
- check cloud/cost safety,
- then decide what belongs in GitHub.

The generated working folder is not a git repo:

`C:\Users\shres\Desktop\HILDALE RESEARCH\sagear_reproduction`

The deliverables folder is:

`C:\Users\shres\Documents\Codex\2026-06-29\8\outputs`

Only ALDERAAN's external dependency clone was found as a git repo under:

`C:\Users\shres\Desktop\HILDALE RESEARCH\external\alderaan`

Therefore, do not modify any previous GitHub repo. Create a clean new repo layout from the approved files after review.

## User Wants From This Review

Please give the user:

1. A concise verdict: correct / partly correct / wrong.
2. Concrete mistakes, stale assumptions, or unverified claims.
3. Which files in `new_repo_candidate/` should be kept.
4. Which files should be removed or ignored.
5. Any code edits needed before making the new repo.
6. Whether the ALDERAAN cloud plan is safe and scientifically sensible.
7. The next safest command to run, if any.

## Current Best Operational State Claimed By Codex

| item | count / state |
|---|---:|
| canonical sample planets | 2474 |
| existing extracted eccentricity posterior summaries | 1729 |
| true missing posterior rows | 745 |
| missing launch-ready posterior rows | 718 |
| missing launch-ready target systems | 592 |
| ALDERAAN catalog rows after whole-system expansion | 767 |
| missing-depth / unseeded blocked rows | 27 |
| existing posterior rows flagged for later review/refit | 320 |

Current category counts claimed from best canonical sample:

| category | current | Sagear comparison target |
|---|---:|---:|
| thin singles | 1109 | 1121 |
| thick singles | 272 | 275 |
| thin multis | 860 | 862 |
| thick multis | 233 | 207 |

Interpretation claimed by Codex:

- The pipeline is close before ALDERAAN, but not "perfect."
- Three groups are very close to Sagear; thick multis remain high.
- The exact total is complicated by Sagear manuscript/macro inconsistency: some places imply 2465 planets, while thin+thick macro totals imply 2474.
- Filling missing ALDERAAN posteriors is now the best next scientific step, but first run one shard only.

## Files To Read First

These are in the review bundle root:

| file | why it matters |
|---|---|
| `report_Hildale.md` | full narrative; final dated section supersedes stale intermediate snapshots |
| `posterior_inventory_and_queues_best.md` | concise posterior coverage and queue logic |
| `pre_cloud_alderaan_checklist.md` | operational checklist |
| `cloud_shell_quickstart.md` | current recommended cloud route |
| `gcp_no_charge_safety_checklist.md` | billing/cost controls |
| `README_GCP_MISSING.md` | cloud bundle command guide |
| `cloud_missing_manifest.md` | exact missing-run bundle manifest |
| `alderaan_needed_validation.md` | validation result for broader needed manifest |
| `alderaan_unseeded_historical_seed_audit.md` | why 27 rows remain blocked |

Read scripts under `source_scripts/` for implementation details.

## Key Scripts To Audit

| script | purpose |
|---|---|
| `source_scripts/build_alderaan_needed_manifest.py` | separates existing/flagged/missing ALDERAAN needs |
| `source_scripts/validate_alderaan_needed_manifest.py` | validates needed manifest coverage and joins |
| `source_scripts/audit_unseeded_transit_recovery.py` | checks whether missing-depth rows can be safely seeded from old catalogs |
| `source_scripts/prepare_gcp_missing_alderaan_bundle.py` | builds cloud bundle for launch-ready missing queue |
| `source_scripts/alderaan_batch.py` | ALDERAAN project/catalog/run helper; patched so missing impact no longer blocks cataloging |
| `source_scripts/extract_eccentricity_posteriors.py` | extracts posterior summaries, now supports run IDs/subdirs |
| `source_scripts/merge_posterior_summaries.py` | merges new cloud summaries into existing archive |
| `source_scripts/postprocess_missing_cloud_results.ps1` | local post-cloud extraction/merge wrapper |
| `cloud_scripts/create_gcp_spot_vm.sh` | creates Spot VM with `--max-run-duration 20h` |
| `cloud_scripts/run_batch.sh` | parallel target runner |
| `cloud_scripts/run_one_target.sh` | resumable per-target runner |

## Commands Codex Already Ran / Expected

From the cloud bundle:

```powershell
cd "C:\Users\shres\Desktop\HILDALE RESEARCH\sagear_reproduction\cloud_missing_batch"
& "C:\Users\shres\anaconda3\python.exe" validate_bundle.py
```

Expected:

```text
VALIDATION OK: 592 targets, 767 catalog rows
```

The VM creation script should include:

```bash
--max-run-duration "$MAX_RUN_DURATION"
```

with:

```bash
MAX_RUN_DURATION="${MAX_RUN_DURATION:-20h}"
```

## Specific Things To Try To Falsify

1. Are the 1729 existing posterior summaries really unique planet-level rows?
2. Are the 718 launch-ready missing rows truly missing, not already represented under alternate KOI/KIC naming?
3. Does system-level batching accidentally include or exclude planets incorrectly?
4. Does the ALDERAAN catalog use seed parameters with the right units and required columns?
5. Does dropping the impact-parameter requirement in `alderaan_batch.py` make sense for ALDERAAN, or does it merely move failure later?
6. Are any of the 27 unseeded planets recoverable with a less strict but still defensible historical catalog rule?
7. Is `--max-run-duration 20h` compatible with Spot VM creation for the selected image/zone/current gcloud version?
8. Do cloud scripts resume cleanly after preemption or failed targets?
9. Should generated CSV outputs be committed, ignored, or stored as release artifacts?
10. Does the final recommended order remain: one shard, local postprocess, full missing queue, then consider 320 flagged refits?

## Suggested New GitHub Repo Strategy

Do not dump everything into GitHub. Create a new clean repository with code, docs, and small reproducibility metadata only.

Recommended commit groups after audit:

1. `scripts/` or `sagear_reproduction/` code only:
   - manifest/validation scripts,
   - posterior extraction/merge scripts,
   - cloud-bundle preparation scripts,
   - README explaining how to regenerate outputs.
2. Small metadata outputs only:
   - summary CSVs,
   - markdown reports/checklists,
   - no large posterior FITS,
   - no downloaded Kepler light curves,
   - no generated ALDERAAN result directories.
3. Add/verify `.gitignore` for:
   - `alderaan_project/`,
   - `cloud_missing_batch/outputs/`,
   - `*.fits`,
   - `*.tar.gz`,
   - large posterior/light-curve folders,
   - `__pycache__/`.

Before committing, rerun validation and include the output in the commit message or repository README.

## Desired Claude Output

Please produce:

1. A concise verdict: correct / partly correct / wrong.
2. A list of concrete mistakes or unverified assumptions.
3. A list of files to put into the new GitHub repo.
4. A list of files to exclude.
5. Any code edits needed before porting.
6. The next safest command to run.
