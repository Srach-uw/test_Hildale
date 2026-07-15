# Complete Factorial Validation Metadata

This directory contains compact derived products from the complete 82-fit
ALDERAAN factorial validation. The underlying immutable FITS results and their
run-level provenance live in
`data/alderaan_factorial_validation_20260715/` and are tracked with Git LFS.

## Contents

- `factorial_validation_discovery.csv`: strict discovery audit for all expected
  arm-target result FITS.
- `factorial_validation_direct_exclusions.csv`: direct eccentricity-extraction
  exclusion ledger. It contains no excluded matched planets for this matrix.
- `factorial_validation_arm_planets.csv`: per-planet direct-extraction
  summaries for every arm.
- `factorial_validation_paired_planets.csv`: paired planet-level comparison
  table. This is the required input to
  `scripts/factorial_validation_diagnostics.py`.
- `factorial_validation_paired_metrics.csv`: clustered paired summaries for
  transit-shape, eccentricity, and zeta quantities.
- `factorial_validation_repeatability_thresholds.csv`: observed repeat-seed
  variation used only as a descriptive comparison threshold.
- `factorial_validation_arm_summaries.csv` and
  `factorial_validation_report.md`: direct-extraction and paired-arm summary.
- `secondary_diagnostics/`: population-stratified descriptive summaries,
  system leverage, leave-system-out robustness, exploratory correlations, and
  a rendered diagnostic plot.

## Interpretation

The panel is deliberately targeted, not a random draw from the Kepler planet
population. It can measure whether a configuration change moves the same
posterior, but it cannot by itself revise Sagear et al.'s four-population
eccentricity result. Read
`docs/full_factorial_validation_assessment.md` before interpreting these files.

## Path Sanitization

The published copies of `factorial_validation_arm_planets.csv` and
`factorial_validation_discovery.csv` replace machine-local Windows paths with
repository-relative FITS paths. The `posterior_file` values in the arm table
are logical identifiers under `direct_posteriors_not_released/`; the temporary
direct posterior grids were used to create the published summaries but are not
part of this release. No scientific values or QC flags were changed.

## Regeneration

The analysis was run with 150,000 importance proposals per target, 10,000
target-system cluster-bootstrap replicates, and seed `20260715`. The exact
command and inputs are recorded in
`docs/full_factorial_validation_assessment.md`.
