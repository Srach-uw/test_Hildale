# Hilldale Sagear Replication

Repository for reproducing and auditing the Hilldale / Sagear Kepler disk-eccentricity
analysis. It contains code, cloud-run helpers, reference materials, and small
reproducibility metadata. It intentionally excludes large local data products such as
ALDERAAN posterior FITS files, Kepler light curves, downloaded archives, generated run
directories, and cloud result tarballs.

The local analysis suite is tested on Python 3.11. Use `requirements-lock.txt` for the
verified environment or `requirements.txt` when adapting to another supported platform.

Supersedes the earlier attempt at
https://github.com/Srach-uw/Shreshth_Hildale_Project (see `legacy/README.md`).

## Current Status

- The final published article is now the reference: Sagear et al. (2026), AJ 172, 42,
  DOI [10.3847/1538-3881/ae71bf](https://doi.org/10.3847/1538-3881/ae71bf).
- The official machine-readable host table is bundled under `reference/data/` and
  supplies authoritative velocities, `P_thick`, and disk labels for all 1,888 hosts.
- The missing-posterior cloud campaign recovered 547 / 592 target systems (92.4%).
  The remaining 45 systems failed upstream ALDERAAN numerical or transit-quality checks.
- The postprocessed archive currently covers 2,395 / 2,474 planets (96.8%), before
  final publication-label/sample reconciliation and validation-arm selection.
- A factorial ALDERAAN validation is testing original versus published transit priors,
  long versus available short cadence, deterministic repeatability, and limb-darkening
  choices before one arm is declared canonical.
- Current eccentricity population values remain diagnostic until that validation and the
  Berger-density ambiguity described below are resolved.

## Most Important Published-Record Finding (2026-07-13)

Sagear's machine-readable Table 1 contains 1,515 thin and **373 thick** hosts. The
article prose says 378 thick hosts, but that is internally impossible: 1,515 + 373 =
1,888, and 275 thick singles + 98 thick-multi hosts = 373.

Our broad reconstructed sample overlaps 1,764 published hosts, but 272 of those disk
labels disagree (84.58% agreement). This is a major upstream explanation for the
incorrect reconstructed Toomre diagram. Primary replication runs must now use the
published labels; the reconstructed APOGEE/GMM classifier is sensitivity-only.

Run the publication contract and regenerate the authoritative host/Toomre products:

```bash
python scripts/published_sagear_audit.py
pytest -q scripts/test_published_sagear_audit.py
```

The parser asserts 1,888 hosts, 1,515 thin, 373 thick, 585 measured-velocity hosts,
unique KIC identifiers, and exact agreement between labels and `P_thick > 0.5`.

## Canonical Safety Rules

- Preserve paired ALDERAAN `T14`, `Rp/R*`, and `b` samples and nested-sampling weights.
- Use the direct MacDougall/Sagear postmodel importance sampler with `e in [0,0.95]`.
- Require an explicit posterior-summary path for hierarchical inference.
- Exclude `qc_primary_exclude` rows by default; low importance ESS is not publishable.
- Never invent missing stellar-density errors in a canonical run. The 13% fallback is
  available only through an explicitly labeled sensitivity flag.
- Use reciprocal geometric transit-probability correction for Sagear comparison.
- Treat deterministic sigma-grid inference as a numerical cross-check. The exact paper
  reports NumPyro, two chains, 1,000 steps, and `Rhat < 1.05`.

## Unresolved Publication Ambiguity

The final paper explicitly says the eccentricity analysis uses stellar densities from
Berger et al. (2018), but that catalog publishes radii rather than the homogeneous
log-density column used by this pipeline. The available homogeneous density table is
Berger et al. (2020). This cannot be resolved safely by inference; it requires an author
clarification or a labeled Berger-2018-derived sensitivity reconstruction.

## Earlier Independent Verification (2026-07-02)

The following section records the pre-cloud inventory audit for provenance. Its queue
counts are historical and are superseded by the current-status section above.

- `cloud/validate_bundle.py` against the real bundle → `VALIDATION OK: 592 targets, 767 catalog rows`.
- Existing posterior archive: 1729 rows, all unique `kepoi_name` (no double-counting).
- Queues partition exactly: 1409 baseline + 320 flagged = 1729; 718 launchable + 27 unseeded = 745 missing; total 2474.
- Zero overlap between the missing queue and the existing archive, by `kepoi_name`
  **and** by same-KIC/same-period aliasing — the 718 are genuinely missing.
- Catalog units match ALDERAAN's ingestion (`depth` in ppm ×1e-6, `duration` in hours /24).
- Catalog rows are strictly period-ascending per system (ALDERAAN hard-errors otherwise).
- **Postprocess chain validated twice**: first with a synthetic ALDERAAN results FITS
  (before any cloud spend), then for real (2026-07-03) with an actual `K00179-results.fits`
  produced by ALDERAAN on GCP — period recovered to 8 ppm accuracy (20.740124 fit vs
  20.740286 catalog), a physically sensible eccentricity posterior (e50=0.226,
  16th/84th-percentile bounds, 4000 resampled draws), merged into the real 1729-row
  archive with zero duplicates and the value preserved through the merge.
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
