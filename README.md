# Hildale Sagear Replication

Code, compact evidence tables, and validation data for an independent
replication of Sagear et al. (2026), *The Orbital Eccentricities of Planets in
the Kinematic Thin and Thick Galactic Disks*.

Reference article:
[doi:10.3847/1538-3881/ae71bf](https://doi.org/10.3847/1538-3881/ae71bf)

## Status

This repository does not reproduce the population eccentricities reported in
the paper. It does reproduce the published planet inventory, preserves the
paper's host classifications and pre-cut multiplicity labels, and constructs a
uniform public-data posterior set for all 2,465 planets.

The table below uses the literal public-data sensitivity: paired ALDERAAN
transit-shape samples, dynesty `LN_WT` weights, row-paired period samples,
fixed Berger et al. (2020) stellar densities, and the reciprocal selection rule
preserved in the arXiv v1 source comments. The validated generative default is
the forward-normalized selection model; neither convention resolves the
discrepancy.

| Population | Public-data N | QC N | QC mean eccentricity | Paper N | Paper mean eccentricity |
| --- | ---: | ---: | ---: | ---: | ---: |
| Thin singles | 1,109 | 1,104 | 0.232 (0.225-0.239) | 1,121 | 0.022 (0.017-0.029) |
| Thick singles | 269 | 268 | 0.208 (0.194-0.223) | 275 | 0.066 (0.045-0.096) |
| Thin multis | 878 | 873 | 0.094 (0.088-0.103) | 862 | 0.030 (0.023-0.031) |
| Thick multis | 209 | 209 | 0.133 (0.121-0.147) | 207 | 0.033 (0.015-0.065) |

These values are diagnostic. They are not an independent astrophysical
measurement. The remaining disagreement is upstream of the population grid and
cannot be resolved from the public products alone with the checks completed so
far.

See [docs/replication_status.md](docs/replication_status.md) for the scientific
interpretation and remaining information gaps. Compact evidence for the current
reconstruction is in
[`metadata/public_reconstruction_20260727/`](metadata/public_reconstruction_20260727/).

## Main Findings

1. The machine-readable host table contains 1,888 stars: 1,515 labeled thin and
   373 labeled thick. The text's value of 378 thick hosts is inconsistent with
   the total and subgroup arithmetic.
2. Multiplicity must be assigned from the full eligible planet inventory before
   planet-level quality cuts. A superseded 2,474-row reconstruction moved 45
   planets when this rule was corrected; that count is historical evidence,
   not a claim about the final 2,465-row inventory.
3. ALDERAAN duration, radius-ratio, impact-parameter, and period samples must
   remain paired. The extractor now enforces this contract.
4. Dynesty nested points must be weighted with `LN_WT`. Equal weighting is kept
   only as an explicitly invalid diagnostic because it produced a misleading
   numerical coincidence.
5. The completed 82-fit validation matrix shows that cadence, limb darkening,
   printed transit priors, and sampler seed do not explain the population-wide
   discrepancy.
6. A direct audit of 2,124 planets with available weighted ALDERAAN fits finds
   median absolute circular-density disagreements of 0.154-0.196 dex across the
   four populations. The corresponding posterior widths are 0.598-0.675 dex.
   This localizes the mismatch upstream of the population fit without yet
   identifying a single cause.
7. The remaining high-value unknowns are the exact stellar-density product, the
   final rejected-planet list, the released or intermediate `(e, omega)`
   posteriors, and the population-export convention used for Table 3.

## Reproduce the Checks

Use Python 3.11:

```bash
python -m venv .venv
python -m pip install -r requirements-lock.txt
python -m pytest -q
python scripts/check_professor_release.py
git diff --check
```

ALDERAAN uses a separate pinned environment and repository commit. Cloud
execution is optional and billable. Read
[`docs/gcp_no_charge_safety_checklist.md`](docs/gcp_no_charge_safety_checklist.md)
before creating a VM.

## Repository Layout

| Path | Contents |
| --- | --- |
| `scripts/` | Sample construction, posterior extraction, hierarchy, diagnostics, and tests |
| `cloud/` | Reproducible ALDERAAN execution and recovery helpers |
| `metadata/` | Compact count audits, fit summaries, and provenance tables |
| `data/alderaan_factorial_validation_20260715/` | 82 validation FITS and manifests stored with Git LFS |
| `reference/` | Published article, tables, and source references |
| `docs/` | Current status, focused audits, and runbooks |
| `legacy/` | Superseded early analysis retained for provenance |

Large and regenerable products are excluded: Kepler light curves, posterior
archives, cloud result bundles, virtual environments, checkpoints, temporary
directories, and private chronological worklogs.

Retrieve the validation FITS after cloning:

```bash
git lfs install
git lfs pull
```

## Citation

Citation metadata for this repository and the source article are provided in
[`CITATION.cff`](CITATION.cff). GitHub can render this file as APA or BibTeX.

## Scientific Boundary

The repository supports reproducibility and diagnosis. A match obtained by
relabeling populations, removing influential planets after looking at the
answer, or ignoring nested-sampling weights is not accepted as a replication.
The reported discrepancy remains open pending author products or a new,
falsifiable methodological explanation.
