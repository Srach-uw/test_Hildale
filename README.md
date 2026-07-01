# Hilldale Sagear Replication

Clean repository candidate for reproducing and auditing the Hilldale / Sagear Kepler disk-eccentricity analysis.

This repository candidate should be reviewed before publishing. It contains code, cloud-run helpers, and small reproducibility metadata. It intentionally excludes large local data products such as ALDERAAN posterior FITS files, Kepler light curves, downloaded archives, generated run directories, and cloud result tarballs.

## Current Status

- Canonical pre-ALDERAAN sample and posterior inventory have been rebuilt.
- Existing eccentricity posterior summaries cover 1729 / 2474 planets.
- True missing posterior queue contains 745 planets.
- Launch-ready missing ALDERAAN queue contains 718 planets across 592 KOI systems.
- 27 missing rows remain blocked by missing/unreliable transit seeds.
- 320 existing posterior rows are flagged for later review or possible refit.

## Manifest Scope & Reproducibility

Two different manifests appear in `docs/` and `metadata/`; they are consistent (subset vs.
superset), not contradictory:

- **Missing-only cloud bundle** (`cloud/`, `docs/cloud_missing_manifest.md`): **718 planets /
  592 targets / 767 catalog rows** — only truly missing, launch-ready posteriors.
- **Needed superset** (`docs/alderaan_needed_validation.md`): **1065 planets / 889 targets /
  1225 catalog rows** — the 745 missing plus the 320 existing-but-flagged refit candidates
  (745 + 320 = 1065; minus 27 unseeded = 1038 runnable). Run this only after the missing-only
  queue is filled.

**Reproducibility caveat:** the real data CSVs (`targets_missing_launchable.csv`,
`sagear_missing_catalog.csv`, the extracted posterior archive, and the canonical sample) live
outside this repo and are intentionally git-ignored. The `VALIDATION OK: 592 targets, 767
catalog rows` result is regenerated locally by running `cloud/validate_bundle.py` against those
inputs — it cannot be reproduced from the repo alone.

## Recommended Scientific Order

1. Audit sample construction and disk classification counts.
2. Audit posterior inventory and missing/flagged queue logic.
3. Run one missing-posterior ALDERAAN shard only.
4. Postprocess and merge that shard locally.
5. Run the full missing launch-ready queue only after the shard passes.
6. Reassess hierarchical fits.
7. Only then consider rerunning the 320 existing-but-flagged posterior cases.

## Repository Layout

| path | contents |
|---|---|
| `scripts/` | Python and PowerShell analysis/reproduction scripts |
| `cloud/` | GCP/ALDERAAN run helper scripts |
| `docs/` | handoff reports, checklists, and audit notes |
| `metadata/` | small summary CSVs and validation summaries |

## Important

Cloud VM usage is billable unless covered by free trial, grant, or university credits. Read `docs/gcp_no_charge_safety_checklist.md` before launching anything.

