# Complete Factorial Validation Assessment

Updated: 2026-07-15

## Purpose

This controlled ALDERAAN experiment tests whether a small set of documented
implementation differences can plausibly explain the eccentricity mismatch in
the current Sagear replication. It is a paired target-level validation, not a
population measurement. Every comparison keeps the same target and planet
where possible, and bootstrap resampling treats a host system as the
independent unit.

The immutable result FITS and their run provenance are in
`data/alderaan_factorial_validation_20260715/`. Compact analysis tables are in
`metadata/factorial_validation_20260715/`.

## Completed Matrix

| arm | target-system FITS | purpose |
|---|---:|---|
| `original_lc` | 24 | Baseline public ALDERAAN configuration with long cadence |
| `reference_lc` | 24 | Reference limb-darkening centers with long cadence |
| `original_lcsc` | 9 | Original limb darkening with available long and short cadence |
| `reference_lcsc` | 9 | Reference limb darkening with available long and short cadence |
| `original_lc_repeat` | 8 | Independent repeat seed for sampler variability |
| `paper_priors_original_lc` | 8 | Printed Table 1 priors versus public defaults |

All 82 expected FITS were discovered and directly re-extracted. No planet in
the matched matrix was excluded by the direct eccentricity QC.

## What the Matrix Shows

The values below are paired differences in posterior median eccentricity,
comparison minus baseline. They are descriptive effects in a selected panel,
not estimates for the full Kepler population.

| comparison | matched planets / systems | median delta e | 95% system-bootstrap interval | median absolute delta e | interpretation |
|---|---:|---:|---:|---:|---|
| Reference minus original limb darkening | 34 / 24 | +0.00120 | -0.01034 to +0.00512 | 0.01545 | No coherent panel-wide shift |
| Long plus short cadence minus long, original LD | 13 / 9 | +0.00101 | -0.00059 to +0.01280 | 0.00682 | No coherent panel-wide shift |
| Long plus short cadence minus long, reference LD | 13 / 9 | -0.00081 | -0.04771 to +0.00263 | 0.00431 | Median stable; individual systems sensitive |
| Printed priors minus public defaults | 9 / 8 | +0.00014 | -0.00053 to +0.00065 | 0.00065 | No material median effect in this panel |
| Same configuration, different sampler seed | 9 / 8 | +0.00071 | -0.00108 to +0.00276 | 0.00108 | Empirical repeatability reference |

For comparison, the current thin-single diagnostic mean eccentricity is 0.335
versus 0.022 in Sagear Table 2. A panel-wide change at the scale measured here
cannot reconcile that difference.

## Individual-System Sensitivity

The lack of a global median shift does not mean all fits are interchangeable.
The largest changes occur in a small number of multi-planet or high-impact
systems:

- `K00283` has the largest reference-LD cadence response, with a maximum
  absolute planet-level eccentricity shift of 0.194.
- `K02533` is the next largest cadence-sensitive system and has large duration
  changes.
- The reference-LD comparison has an exploratory positive association between
  signed eccentricity change and catalog impact parameter (Spearman rho 0.50,
  nominal p 0.012). This was one of many uncorrected exploratory tests and is
  not a population claim.

These systems should receive visual and posterior-shape QC before they are
allowed to have high leverage in a future population fit. They do not justify
discarding a population or changing the canonical configuration by themselves.

## Decision Record

1. Do not claim that limb darkening, cadence treatment, printed priors, or
   nested-sampling randomness resolves the current Table 2 discrepancy.
2. Do not select a new full-run configuration from this targeted panel alone.
   A scientific choice can follow the published methods, but the panel does
   not establish a population-level correction.
3. Preserve the paired ALDERAAN transit-shape extraction, pre-cut multiplicity
   labels, published host disk labels, and explicit QC manifest as mandatory
   contracts.
4. Prioritize resolving the Berger et al. (2018) density construction and the
   final Sagear planet-level inclusion and visual-rejection list. Those inputs
   can change the entire posterior sample and remain unidentifiable from the
   public article.
5. Run the nine-system combined reference-LD, LC+SC, printed-prior confirmation
   arm before assuming these configuration effects are additive. This is a
   narrow interaction check, not a replacement population study. See
   `docs/combined_configuration_confirmation.md`.
6. Any future larger ALDERAAN rerun must use a frozen density source,
   documented cadence and limb-darkening catalog, a target-level exclusion
   ledger, and the same direct paired-impact extractor used here.

## Reproduction

The complete analysis was generated with the strict full-matrix mode:

```powershell
python scripts/compare_factorial_validation.py `
  --validation-root data/alderaan_factorial_validation_20260715 `
  --metadata-root data/alderaan_factorial_validation_20260715/provenance/target_sets `
  --inventory data/alderaan_factorial_validation_20260715/provenance/input_catalogs/full_system_inventory.csv `
  --sample <live-sagear-reproduction>/outputs/canonical_sample_old_astropy_rawcc.csv `
  --config <live-sagear-reproduction>/config.json `
  --output-dir <output-directory> `
  --n-proposals 150000 --bootstrap-replicates 10000 --seed 20260715
```

Run `scripts/factorial_validation_diagnostics.py` on the resulting
`factorial_validation_paired_planets.csv` to reproduce the target-system
leverage and robustness summaries.
