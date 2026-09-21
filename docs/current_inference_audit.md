# Current inference audit

Updated: 2026-09-20

## Scope

This note records the current public-data diagnosis. It does not replace the
published reconstruction outputs or claim a new astrophysical measurement.

## What the reconstruction establishes

- The reconstructed sample contains the paper's reported total of 2,465
  planets and preserves the published host labels.
- Multiplicity is assigned before planet-level fit cuts.
- The primary eccentricity extractor preserves paired ALDERAAN transit-shape
  samples and uses Dynesty `LN_WT` weights.
- The forward-normalized hierarchy passes its algebraic and synthetic recovery
  checks.
- A real-data Gilbert control recovers the expected low-eccentricity scale.

These checks support the mechanics of the reconstruction. They do not show
that every planet-level input and final fit decision matches Sagear et al.

## What remains unresolved

The four reconstructed Rayleigh means remain above the paper's values. Median
density offsets alone do not identify the source of that difference because
the transit-density posteriors are broad. Individual eccentricity summaries
also cannot determine a population mean by themselves.

The 82-fit factorial experiment tested several useful choices, including
limb darkening, selected priors, and sampler seeds. Its intended cadence test
does not close the cadence question: the archived LC+SC runners did not pass
ALDERAAN's `--use_sc True` argument, and their logs report no short-cadence
data. A corrected cadence test is documented in
[`corrected_cadence_validation_plan.md`](corrected_cadence_validation_plan.md).

The public 2018 stellar table does not expose the full density, mass, and
mass-radius uncertainty representation needed to reconstruct the exact
stellar inputs. The final accepted-fit ledger is also not identified by the
available article tables.

## Recent bounded checks

The native likelihood replay for K02712 reproduces stored likelihood values at
archived samples and retains a high-impact preference in fixed-impact profile
tests. This is a numerical-fidelity result for one target, not a population
result or an adequacy test for every noise model.

A bounded saved-posterior comparison found only three valid same-system
matches between the original archive and the recovery archive. Most shared
target identifiers have different jointly fit planet counts, so they cannot
be treated as repeated fits. The six comparable planet slots show small median
shifts in radius ratio and impact. The comparison does not establish identical
light curves, priors, code, or sampler configuration.

## Next discriminating work

1. Test recovery at the published population scales using simulated photometry
   with observed cadence, uncertainty, and selection characteristics.
2. Compare direct eccentric transit fits with the importance-sampling path
   under the same data and priors.
3. Evaluate overlap between full transit-density and stellar-density
   distributions, rather than comparing medians alone.
4. Run a corrected long-plus-short-cadence validation on a predeclared sample.

The author products that would most directly locate any remaining divergence
are the accepted-fit ledger and the stellar-density rows or posterior samples
supplied to the eccentricity calculation. See
[`author_clarification_request.md`](author_clarification_request.md).
