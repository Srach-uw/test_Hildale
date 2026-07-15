# Hilldale Sagear Replication

Research code and audit materials for reproducing the Kepler disk-eccentricity
analysis of Sagear et al. (2026), *The Astronomical Journal*, 172, 42.

Reference article: [doi:10.3847/1538-3881/ae71bf](https://doi.org/10.3847/1538-3881/ae71bf)

## Scientific Status

This is an active replication, not a completed reproduction. The sample construction,
published disk labels, multiplicity bookkeeping, ALDERAAN posterior extraction, and
population model have been audited separately. The current population estimates remain
diagnostic because posterior coverage is incomplete and a controlled ALDERAAN validation
is testing cadence, limb-darkening, prior, and repeatability effects. The complete 82-fit
factorial result release is included through Git LFS; its comparison analysis has not yet
been folded into the population result.

| population | current fit N | current mean eccentricity (16th-84th) | Sagear N | Sagear mean eccentricity (16th-84th) |
|---|---:|---:|---:|---:|
| thin singles | 304 | 0.335 (0.319-0.353) | 1121 | 0.022 (0.017-0.029) |
| thick singles | 108 | 0.288 (0.263-0.315) | 275 | 0.066 (0.045-0.096) |
| thin multis | 222 | 0.122 (0.108-0.136) | 862 | 0.030 (0.023-0.031) |
| thick multis | 69 | 0.117 (0.086-0.154) | 207 | 0.033 (0.015-0.065) |

These values do not reproduce Table 2. They should not be interpreted as a physical
measurement of disk-population eccentricity. The discrepancy is the object of the
remaining validation work.

See [docs/replication_status.md](docs/replication_status.md) for the current evidence,
limitations, and acceptance criteria.

## Findings That Changed The Replication

1. The machine-readable published host table is authoritative for the primary disk
   labels. It contains 1,515 thin and 373 thick hosts. The article's statement of 378
   thick hosts is inconsistent with both the total of 1,888 and subgroup arithmetic.
2. Multiplicity must be assigned from the full eligible planet inventory before
   planet-level quality cuts. Recounting only surviving planets would misclassify 45
   planets in the reconstructed sample. The canonical pipeline now preserves the pre-cut
   system label and fails closed when full-system provenance is unavailable.
3. Paired ALDERAAN samples of transit duration, radius ratio, and impact parameter must
   remain paired during eccentricity extraction. Geometric impact draws are retained only
   as a sensitivity analysis.
4. A mixed archive of heterogeneous posterior constructions is not a valid canonical
   input to the hierarchical fit. Canonical inference requires explicit provenance and
   quality-control fields.
5. The Berger et al. (2018) stellar-density construction used in the article is not fully
   specified by the public catalog. Berger et al. (2020) densities are therefore a labeled
   sensitivity, not a silent substitute.

## Reproduce The Audits

Use Python 3.11 and the pinned analysis dependencies:

```bash
python -m venv .venv
python -m pip install -r requirements-lock.txt
python -m pytest -q
python scripts/published_sagear_audit.py
```

ALDERAAN uses a separate pinned environment and repository commit. Cloud execution is
optional and billable. Read `docs/gcp_no_charge_safety_checklist.md` before creating a VM.

## Repository Layout

| path | contents |
|---|---|
| `scripts/` | sample, posterior, population, diagnostic, and validation code |
| `cloud/` | resumable ALDERAAN and GCP execution helpers |
| `metadata/` | compact derived tables and validation summaries |
| `data/alderaan_factorial_validation_20260715/` | 82 completed factorial ALDERAAN result FITS and their provenance, stored with Git LFS |
| `reference/` | published article, machine-readable tables, and reference material |
| `docs/` | scientific status, methods audits, runbooks, and historical worklog |
| `legacy/` | superseded analysis retained only for provenance |

Large or regenerable products are intentionally excluded: Kepler light curves,
intermediate detrending products, checkpoints, posterior archives, virtual environments,
and cloud result bundles. The 82 curated factorial result FITS are the explicit exception.

After cloning, retrieve the validation data with:

```bash
git lfs install
git lfs pull
```

Read `data/alderaan_factorial_validation_20260715/README.md` before using these outputs.

## Release Check

Before sharing or publishing a revision, run:

```bash
python scripts/check_professor_release.py
python -m pytest -q
git diff --check
```

The release check rejects personal home paths, Unicode em dashes, merge markers, and
common credential patterns in tracked text files.

## Scope

This repository supports reproducibility and diagnosis. It does not claim an independent
astrophysical result until the validation gates in `docs/replication_status.md` pass.
