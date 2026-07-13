# Replication Status

Updated: 2026-07-13

## Bottom Line

The repository now reproduces and audits the major data contracts, but it does not yet
reproduce Sagear et al. (2026) Table 2. The remaining discrepancy is too large to describe
as numerical noise. Current evidence points to posterior-construction and sample-provenance
differences rather than a single population-fit bug.

## Published Targets

| population | planets | mean eccentricity | 16th-84th interval |
|---|---:|---:|---:|
| thin singles | 1121 | 0.022 | 0.017-0.029 |
| thick singles | 275 | 0.066 | 0.045-0.096 |
| thin multis | 862 | 0.030 | 0.023-0.031 |
| thick multis | 207 | 0.033 | 0.015-0.065 |

The four planet bins sum to 2,465. The paper also reports disk-total macros that sum to
2,474, so the nine-planet difference is documented as a publication ambiguity rather
than silently forced to match.

## Current Uniform Diagnostic Fit

| population | fit N | mean eccentricity | 16th-84th interval | ratio to published mean |
|---|---:|---:|---:|---:|
| thin singles | 304 | 0.335 | 0.319-0.353 | 15.2 |
| thick singles | 108 | 0.288 | 0.263-0.315 | 4.4 |
| thin multis | 222 | 0.122 | 0.108-0.136 | 4.1 |
| thick multis | 69 | 0.117 | 0.086-0.154 | 3.5 |

This fit uses only the uniformly processed direct subset. Its coverage is incomplete and
non-random, so the table is a diagnostic of the mismatch, not a final measurement.

## Completed Checks

- The published machine-readable table contains 1,888 hosts: 1,515 thin and 373 thick.
- Primary disk labels now come from the published host table. A reconstructed GMM label is
  retained only for sensitivity analysis.
- Multiplicity is assigned before planet-level quality cuts and preserved afterward.
- The corrected 2,474-row reconstructed inventory has zero multiplicity-contract
  mismatches. Recounting after cuts would have moved 45 planets between bins.
- The missing-posterior campaign produced 547 successful target-system fits from 592
  attempted systems. Failures and incomplete systems are represented in QC manifests.
- Eccentricity extraction preserves nested-sampling weights and paired ALDERAAN transit
  parameters. Synthetic geometric impact draws are not canonical.
- Hierarchical inference rejects missing provenance, mixed transit-prior conventions, and
  absent required QC fields.

## Factorial Validation

The controlled validation compares original and reference limb darkening, long cadence
and available short cadence, a repeated nested-sampling run, and the paper's printed
priors. The first 23 matched systems contain 31 matched planets.

For reference minus original limb darkening, the median change in posterior median
eccentricity is +0.00168, with a bootstrap interval of -0.00103 to +0.00523. The median
absolute change is 0.00604. This is much smaller than the full Table 2 discrepancy, but
repeatability, cadence, and prior arms are still required before attributing causes.

## Remaining Gates

1. Finish and analyze the repeatability, cadence, and prior validation arms.
2. Quantify whether between-run stochastic variation exceeds the limb-darkening shift.
3. Freeze one posterior-construction arm and one explicit stellar-density source.
4. Regenerate the deterministic QC manifest and population inputs from that arm.
5. Run leave-10-percent-out and per-planet leverage checks.
6. Compare with Table 2 only after all four populations pass the same provenance rules.

## Known Limitation

The article states that Berger et al. (2018) stellar densities were used, but the public
catalog does not provide the homogeneous density field consumed by the current pipeline.
Berger et al. (2020) densities are therefore labeled as a sensitivity. Resolving the exact
2018 construction requires author clarification or a documented reconstruction.

## Canonical Evidence

- `reference/data/`: published machine-readable tables.
- `metadata/sagear2026_publication_audit_summary.csv`: published host-count contract.
- `metadata/rayleigh_population_fit_transit_selection_manuscript_reciprocal_direct_count_calibrated_q0535.csv`:
  current uniform diagnostic fit.
- `docs/sagear_diagnosis_report.md`: detailed scientific audit.
- `docs/report_Hildale.md`: chronological project worklog; older entries are superseded by
  later dated findings where they conflict.
