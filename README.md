# Hilldale Sagear Replication

Repository for reproducing and auditing the Hilldale / Sagear Kepler disk-eccentricity
analysis. It contains code, cloud-run helpers, reference materials, and small
reproducibility metadata. It intentionally excludes large local data products such as
ALDERAAN posterior FITS files, Kepler light curves, downloaded archives, generated run
directories, and cloud result tarballs.

Supersedes the earlier attempt at
https://github.com/Srach-uw/Shreshth_Hildale_Project (see `legacy/README.md`).

## Current Status

- Canonical pre-ALDERAAN sample and posterior inventory have been rebuilt.
- Existing eccentricity posterior summaries cover 1729 / 2474 planets.
- True missing posterior queue contains 745 planets.
- Launch-ready missing ALDERAAN queue contains 718 planets across 592 KOI systems.
- 27 missing rows remain blocked by missing/unreliable transit seeds.
- 320 existing posterior rows are flagged for later review or possible refit.

## Independent Verification (2026-07-02)

All headline counts were re-verified against the real local data CSVs (not just
internal consistency of the summaries):

- `cloud/validate_bundle.py` against the real bundle → `VALIDATION OK: 592 targets, 767 catalog rows`.
- Existing posterior archive: 1729 rows, all unique `kepoi_name` (no double-counting).
- Queues partition exactly: 1409 baseline + 320 flagged = 1729; 718 launchable + 27 unseeded = 745 missing; total 2474.
- Zero overlap between the missing queue and the existing archive, by `kepoi_name`
  **and** by same-KIC/same-period aliasing — the 718 are genuinely missing.
- Catalog units match ALDERAAN's ingestion (`depth` in ppm ×1e-6, `duration` in hours /24).
- Catalog rows are strictly period-ascending per system (ALDERAAN hard-errors otherwise).
- **Postprocess chain smoke-tested end-to-end**: a synthetic ALDERAAN results FITS
  (matching the real format: `NPL` header, `TTIMES_nn` HDUs, `SAMPLES` table) was run
  through the real `extract_eccentricity_posteriors.py` (flat `--results-dir` mode) and
  `merge_posterior_summaries.py` — correct category match, 1729+1=1730 merged rows, zero
  duplicates. The local half of the cloud round-trip works before any cloud spend.
- **Cloud flag compatibility verified against Google docs**: `--max-run-duration` works
  with SPOT + STOP (see `docs/gcp_no_charge_safety_checklist.md`, incl. cost estimate).
- **3 of the 27 unseeded planets are recoverable** via DR25 TCE seeds
  (see `docs/unseeded_dr25_tce_recovery.md`) — recommended for the follow-up batch,
  not the current bundle. K01316.02/K06516.02 are confirmed unrecoverable.

### Known caveats found during verification

1. **Two systems will be fit incomplete.** K01316 and K06516 each have a second
   sample planet (K01316.02, K06516.02) that is unseeded (no valid depth), so the
   catalog carries them as `npl=1`. Their unmodeled sibling transits remain in the
   light curve and may bias those two fits. **Flag the K01316 / K06516 results for
   review after the cloud run** (or exclude the 2 targets; the other 590 are clean).
2. **32 catalog rows have KOI impact seeds > 1** (max 8.7, unphysical). This is
   harmless: ALDERAAN clamps them internally
   (`detrend_and_estimate_ttvs.py`: `if b > 1 - sqrt(depth): b = (1 - sqrt(depth))**2`),
   and the seed only initializes the fit. Similarly, `scripts/alderaan_batch.py`
   seeds missing impacts with a neutral 0.5.

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
inputs — it cannot be reproduced from the repo alone. External catalog inputs (Furlan 2017,
APOGEE DR17 crossmatch, etc.) are regenerated with `scripts/prepare_external_inputs.py`.

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
| `scripts/` | Python and PowerShell analysis/reproduction scripts (incl. diagnostics) |
| `cloud/` | GCP/ALDERAAN run helper scripts and target shards |
| `docs/` | handoff reports, checklists, and audit notes |
| `metadata/` | small summary CSVs and validation summaries |
| `reference/` | Sagear et al. paper, reference figures, original photoeccentric script |
| `legacy/` | superseded results from the earlier 1716-planet attempt (provenance only) |

## Important

Cloud VM usage is billable unless covered by free trial, grant, or university credits. Read
`docs/gcp_no_charge_safety_checklist.md` before launching anything.
