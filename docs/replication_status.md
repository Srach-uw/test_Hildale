# Replication Status

Updated: 2026-08-08

## Bottom Line

The public-data pipeline now has a usable ALDERAAN result for every planet in
the reconstructed 2,465-planet inventory. The sample, multiplicity, posterior,
and population-fit contracts have been checked independently. A real-data
Gilbert control recovers its published small-planet eccentricity scale, while
the Sagear Rayleigh eccentricities remain substantially above Table 3.

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

## Current Weighted Diagnostic

The primary public-data diagnostic uses:

- published host labels;
- multiplicity assigned before planet-level cuts;
- paired ALDERAAN `T14`, `Rp/R*`, impact, and period samples;
- dynesty `LN_WT` weights;
- fixed Berger et al. (2020) density;
- direct MacDougall-style post-model importance sampling;
- the arXiv v1 reciprocal transit-selection rule, retained only as a literal
  replication sensitivity;
- deterministic posterior QC.

| Population | Reconstructed N | Weighted Rayleigh mean | Paper N | Paper mean |
| --- | ---: | ---: | ---: | ---: |
| Thin singles | 1,109 | 0.246 | 1,121 | 0.022 |
| Thick singles | 269 | 0.234 | 275 | 0.066 |
| Thin multis | 878 | 0.171 | 862 | 0.030 |
| Thick multis | 209 | 0.141 | 207 | 0.033 |

These values are the deterministic 1,000-draw SIR branch used in the final
weighting audit. The earlier grid and QC sensitivities remain in the evidence
tables. None reproduces the paper.

Equal raw-row weighting gives 0.0747, 0.0630, 0.0332, and 0.0444 for thin
singles, thick singles, thin multis, and thick multis. It is not a valid
nested-sampling estimator. More importantly, its BIC comparison favors a
half-Gaussian in all four populations, whereas the paper reports a strong
Rayleigh preference. It cannot be used as an answer-matching shortcut.

## Checks That Are Closed

- The published inventory total is recovered exactly.
- Published host labels are used directly for the primary analysis.
- Multiplicity is frozen before planet-level fit cuts.
- ALDERAAN transit-shape samples remain paired.
- Row-paired ALDERAAN periods replace fixed catalog periods. The controlled
  change shifts the population values by at most `1.2e-5`.
- Nested rows use `LN_WT`. Equal raw-row weighting is rejected.
- The hierarchy requires homogeneous provenance and rejects a mixed archive.
- Nonzero-eccentricity round trips test the MacDougall equation's sign,
  angular units, and velocity factor.
- The generative forward-normalized hierarchy passes synthetic recovery; the
  arXiv v1 reciprocal rule is not used as a validated scientific estimator.
- Fixed, Gaussian, and split-normal density sensitivities do not recover
  Table 3.
- The 82-fit experiment and nine-system combined confirmation show that the
  tested cadence, limb-darkening, printed-prior, and sampler-seed differences
  are not the main explanation. The combined arm changes median eccentricity
  by only +0.00042 across 13 matched planets, although several individual
  systems move more than the small repeat-run threshold.
- Removing the most influential planets cannot reproduce the thin-single value.
- A real-ALDERAAN Gilbert control recovers a Beta mean of 0.0487 for all small
  planets and 0.0684 for observed singles. The same quality-cut ladder applied
  to the Sagear populations does not recover Table 3.
- The weighted and equal-row branches fail in different ways: the weighted
  branch misses both values and model ordering, while the equal-row branch
  partly approaches the values but still misses the ordering.
- The paired ALDERAAN circular-density audit contains all 2,465 planets and places
  the disagreement in the transit-shape and density inputs before hierarchy.
- Impact-parameter uncertainty contributes to broad density posteriors. In an
  exploratory narrower-impact subset, 114 thin singles remain after five
  primary-QC failures are removed and give a forward-normalized hierarchical
  mean eccentricity of 0.333 (0.317-0.350). This selected subset is diagnostic
  rather than an unbiased population estimate.

## Remaining Information Gaps

1. The exact stellar-density and planet-radius tables supplied to the
   eccentricity calculation.
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
- `metadata/final_forensic_20260808/`: Gilbert control, residual attribution,
  full-sample QC sensitivity, and nested-weighting diagnostics.
- `metadata/uncertainty_calibration_20260810/`: bootstrap interval calibration
  and the Berger-2018 density-reconstruction branch.
- `docs/full_factorial_validation_assessment.md`: interpretation of the 82-fit
  experiment.
- `docs/uncertainty_and_density_findings_20260810.md`: identifies stellar-density
  uncertainty as the dominant driver of inferred eccentricity (Spearman +0.52
  against per-planet e50, larger than any other variable tested). Records that
  the Berger-2018 reconstruction carries a defective uncertainty model, bimodal
  at 0.315 and 1.712 fractional with 35 percent of planets having sigma larger
  than rho, and that the defect tracks planet multiplicity, which is physically
  impossible. That confound explains why this reconstruction reports singles as
  twice as eccentric as multis. The canonical Berger-2020 path is verified sound
  (unimodal, median 0.101, no planet with sigma above rho). Also measures the
  per-planet information content directly (median KL from prior of 0.019 nats
  for thick singles, 0.141 for thick multis), shows the published delta-BIC model
  ranking of about 20 to 25 is not reproduced (maximum 6.1), closes five
  hypotheses by measurement, and withdraws the Table 3 ordering lead.
- `docs/author_clarification_request.md`: concise request for unavailable
  reproduction inputs.
