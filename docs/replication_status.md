# Replication Status

Updated: 2026-07-15

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

## Complete Factorial Validation

The controlled ALDERAAN matrix is complete: 82 of 82 expected target-system FITS across
six paired arms. It compares original and reference limb darkening, long cadence and
available short cadence, independent nested-sampling seeds, and the paper's printed
transit priors. Direct re-extraction found no matched-matrix QC exclusions.

At the panel level, none of these tested choices produces an eccentricity shift remotely
large enough to explain the current thin-single diagnostic mismatch. Reference minus
original limb darkening has median delta e = +0.00120 across 34 planets in 24 systems;
the corresponding 95% system-bootstrap interval is -0.01034 to +0.00512. The remaining
arms likewise have median effects of order 0.001 in this targeted panel.

Some individual systems are sensitive, especially `K00283` and `K02533`, so target-level
QC remains necessary. The full evidence and decision record are in
`docs/full_factorial_validation_assessment.md` and
`metadata/factorial_validation_20260715/`.

## Remaining Gates

1. Freeze an explicit stellar-density construction. The paper identifies Berger et al.
   (2018), but the exact density calculation and uncertainties are not public.
2. Obtain or reconstruct Sagear's final planet-level inclusion list and visual-fit
   rejection list.
3. Apply target-level visual and posterior-shape QC to the high-leverage systems before
   any full population refit.
4. Run the nine-system combined reference-LD, LC+SC, printed-prior confirmation arm.
5. Only then choose whether a large uniform rerun is justified and compare with Table 2.

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
- `docs/scientific_interim_assessment.md`: current evidence synthesis and robustness results.
- `docs/report_Hildale.md`: chronological project worklog; older entries are superseded by
  later dated findings where they conflict.

## Completed Factorial Result Release

The complete factorial matrix now contains 82 successful target-system result FITS:
24 original-long-cadence, 24 reference-long-cadence, 9 original long-plus-short-cadence,
9 reference long-plus-short-cadence, 8 repeat-seed, and 8 paper-prior runs. The immutable
results, input catalogs, target selections, status manifests, runner specification, and
SHA-256 manifest are published under
`data/alderaan_factorial_validation_20260715/` through Git LFS.

The 82-file release is a reproducibility artifact, not a final population result. Its
matched target and planet analysis is complete and rules out several easy explanations;
it does not remove the unresolved density and final-QC provenance limitations.
