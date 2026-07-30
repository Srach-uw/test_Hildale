# Replication Status

Updated: 2026-07-31

## Bottom Line

The public-data pipeline now has a usable ALDERAAN result for every planet in
the reconstructed 2,465-planet inventory. The sample, multiplicity, posterior,
and population-fit contracts have been checked independently. The resulting
Rayleigh eccentricities remain substantially above Sagear et al. (2026)
Table 3.

This is a failed public-data replication, not evidence that the paper's
astrophysical conclusion is wrong. Several inputs required for an exact
reproduction are not public.

## Published and Reconstructed Populations

| Population | Reconstructed N | Paper N | Difference |
| --- | ---: | ---: | ---: |
| Thin singles | 1,109 | 1,121 | -12 |
| Thick singles | 269 | 275 | -6 |
| Thin multis | 878 | 862 | +16 |
| Thick multis | 209 | 207 | +2 |
| Total | 2,465 | 2,465 | 0 |

The total is exact, but individual planet membership in the paper's final
post-fit rejection set is unavailable. Count-matched and adversarial-removal
tests do not account for the full eccentricity difference.

## Current Public-Data Fit

The canonical sensitivity uses:

- published host labels;
- multiplicity assigned before planet-level cuts;
- paired ALDERAAN `T14`, `Rp/R*`, impact, and period samples;
- dynesty `LN_WT` weights;
- fixed Berger et al. (2020) density;
- direct MacDougall-style post-model importance sampling;
- manuscript transit-selection normalization;
- deterministic posterior QC.

| Population | All N | All mean eccentricity | QC N | QC mean eccentricity | Paper |
| --- | ---: | ---: | ---: | ---: | ---: |
| Thin singles | 1,109 | 0.249 (0.242-0.257) | 1,104 | 0.232 (0.225-0.239) | 0.022 |
| Thick singles | 269 | 0.234 (0.220-0.249) | 268 | 0.208 (0.194-0.223) | 0.066 |
| Thin multis | 878 | 0.183 (0.176-0.190) | 873 | 0.094 (0.088-0.103) | 0.030 |
| Thick multis | 209 | 0.133 (0.121-0.147) | 209 | 0.133 (0.121-0.147) | 0.033 |

The QC sensitivity is reported because the thin-multi result is strongly
affected by a small number of poor or high-leverage posteriors. Neither table
reproduces the paper.

## Checks That Are Closed

- The published inventory total is recovered exactly.
- Published host labels are used directly for the primary analysis.
- Multiplicity is frozen before planet-level fit cuts.
- ALDERAAN transit-shape samples remain paired.
- Row-paired ALDERAAN periods replace fixed catalog periods. The controlled
  change shifts the population values by at most `1.2e-5`.
- Nested rows use `LN_WT`. Equal raw-row weighting is rejected.
- The hierarchy requires homogeneous provenance and rejects a mixed archive.
- The hierarchy reproduces the archived M-dwarf Rayleigh benchmark.
- Fixed, Gaussian, and split-normal density sensitivities do not recover
  Table 3.
- The 82-fit experiment rules out ordinary cadence, limb-darkening,
  printed-prior, and sampler-seed differences as the main explanation.
- Removing the most influential planets cannot reproduce the thin-single value.

## Remaining Information Gaps

1. The exact stellar-density table and uncertainty representation supplied to
   the eccentricity calculation.
2. The identities of the final rejected planets and systems.
3. The per-planet `(e, omega)` posterior samples used in the paper.
4. The exact sampling-importance-resampling convention applied to ALDERAAN
   nested points.
5. The population-fitting implementation and explicit Table 3 array ordering.

These are appropriate questions for the corresponding author. Further
parameter tuning against the published answer would not be a valid
replication.

## Evidence

- `metadata/public_reconstruction_20260727/`: current compact reconstruction.
- `metadata/recovery_preflight_20260724/`: published-host and recovery audits.
- `metadata/factorial_validation_20260715/`: six-arm ALDERAAN validation.
- `docs/full_factorial_validation_assessment.md`: interpretation of the 82-fit
  experiment.
- `docs/author_clarification_request.md`: concise request for unavailable
  reproduction inputs.
