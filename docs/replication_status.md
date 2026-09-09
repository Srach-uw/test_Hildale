# Replication status

Updated: 2026-08-13

## Bottom line

The public reconstruction recovers the paper's total of 2,465 planets and uses
the published host classifications. The strict source-faithful hierarchy does
not reproduce the four Rayleigh means in Table 3. This is a failed numerical
replication from public inputs, not evidence that the paper's astrophysical
interpretation is wrong.

The remaining discrepancy is already present in the comparison between the
ALDERAAN circular-density posteriors and the adopted stellar-density priors.
Changing the population model cannot remove that upstream disagreement.

## Inventory

| Population | Reconstructed | Paper | Difference |
| --- | ---: | ---: | ---: |
| Thin singles | 1,109 | 1,121 | -12 |
| Thick singles | 269 | 275 | -6 |
| Thin multis | 878 | 862 | +16 |
| Thick multis | 209 | 207 | +2 |
| Total | 2,465 | 2,465 | 0 |

The exact total is recovered. The paper does not publish the final accepted
planet ledger, so the remaining category differences cannot be resolved by ID.
Multiplicity is assigned from the full eligible system before planet-level
quality cuts.

## Current source-faithful result

The comparison below uses paired ALDERAAN transit-shape samples, dynesty
`LN_WT` weights, the documented post-fit quality rules, and the paper's
reciprocal transit-selection prescription as a literal replication branch.
The hierarchy-ready counts are smaller than the inventory because only planets
with valid posterior support and complete quality fields enter the fit.

| Population | Fit N | Reconstructed mean | Paper N | Paper mean |
| --- | ---: | ---: | ---: | ---: |
| Thin singles | 887 | 0.214 | 1,121 | 0.022 |
| Thick singles | 218 | 0.153 | 275 | 0.066 |
| Thin multis | 821 | 0.086 | 862 | 0.030 |
| Thick multis | 193 | 0.114 | 207 | 0.033 |

These values remain diagnostic rather than publishable reproduction results.
The public products do not establish that the same planet-level inputs and
final rejection ledger were used in the paper.

## What has been ruled out

- Published host labels replace reconstructed disk labels in the primary run.
- Multiplicity is frozen before planet-level fit cuts.
- ALDERAAN `T14`, `Rp/R*`, impact, and period samples remain paired.
- Nested samples use `LN_WT`; equal raw-row weighting is invalid.
- Fixed catalog periods and paired fitted periods give indistinguishable
  population results in the controlled test.
- The forward-normalized hierarchy passes synthetic recovery tests.
- Fixed, Gaussian, and split-normal density treatments do not recover Table 3.
- The 82-fit factorial experiment does not identify cadence, limb darkening,
  printed priors, or sampler seed as the principal cause.
- The Gilbert real-data control recovers the expected low-eccentricity scale.
- High-leverage removal does not recover the thin-single value robustly.
- Importance-sampling effective sample sizes are adequate; Monte Carlo scatter
  is much smaller than the discrepancy.
- The ALDERAAN radius-ratio prior difference printed in Table 2 has negligible
  effect on the circular-density width in the direct test.

## Remaining public-data boundary

An exact replication now requires at least one unpublished intermediate:

1. the accepted planet ledger, including the final rejected KOIs;
2. the stellar-density values and uncertainty representation supplied to the
   eccentricity calculation;
3. the final per-planet `(e, omega)` posteriors; or
4. the exact population-model implementation and input ordering.

The article is internally inconsistent about whether Berger et al. (2018) or
Berger et al. (2020) supplied the stellar properties. The public 2018 catalog
does not contain the homogeneous density product needed to reconstruct the
stated prior without additional choices. Those choices materially affect the
inference and should not be guessed to match the published answer.

See `docs/final_public_data_boundary.md` for the complete forensic boundary and
`docs/author_clarification_request.md` for the minimal data request.

## Evidence map

- `metadata/public_reconstruction_20260727/`: compact 2,465-planet inventory.
- `metadata/recovery_preflight_20260724/`: host, multiplicity, and recovery audits.
- `metadata/factorial_validation_20260715/`: six-arm ALDERAAN validation.
- `metadata/combined_confirmation_20260806/`: combined cadence and prior check.
- `metadata/final_forensic_20260808/`: controls, QC sensitivity, and weighting audits.
- `metadata/uncertainty_calibration_20260810/`: density and interval calibration.
- `docs/full_factorial_validation_assessment.md`: interpretation of the 82 fits.
- `docs/uncertainty_and_density_findings_20260810.md`: density diagnostics.
