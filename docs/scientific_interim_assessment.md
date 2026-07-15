# Scientific Interim Assessment

Updated: 2026-07-15

## Executive Finding

The replication is not at a dead end, but the current four-population result is not
publishable. The exact published host classifications and the corrected multiplicity
contract have removed two major upstream ambiguities. The remaining mismatch is centered
on the eccentricity posteriors and the target-level quality decisions that created
Sagear et al.'s final posterior sample.

The current uniformly processed 710-planet subset remains much more eccentric than the
published result. The discrepancy is not explained by limb darkening, transit-selection
weighting, the reconstructed disk classifier, or a small number of thin-single systems.
The complete factorial validation now quantifies short-cadence, prior, and
nested-sampling effects. These tested choices do not create a panel-wide shift large
enough to explain the discrepancy.

## Published Ground Truth

The final AJ article provides a machine-readable host table. It contains 1,888 KIC hosts,
with 1,515 labeled thin disk and 373 labeled thick disk. All labels obey
`P_thick > 0.5`. This table is now the primary classification source.

The article's main sample paragraph contains a typographical inconsistency: one sentence
calls all 378 thick-disk hosts "single hosts." The internally consistent breakdown later
in the article is 275 thick singles plus 98 thick-multi hosts, for 373 unique thick hosts.
The four planet bins are:

| population | published planets | published Rayleigh mean e |
|---|---:|---:|
| thin singles | 1,121 | 0.022 (0.017-0.029) |
| thick singles | 275 | 0.066 (0.045-0.096) |
| thin multis | 862 | 0.030 (0.023-0.031) |
| thick multis | 207 | 0.033 (0.015-0.065) |

## Current Uniform Diagnostic

The direct paired-impact subset uses nested-sampling weights, row-paired ALDERAAN
`T14`, `Rp/Rstar`, and impact samples, the exact MacDougall density equation, and a
deterministic QC manifest. Coverage is incomplete and nonrandom.

| population | planets | hosts | current mean e | published mean e | ratio |
|---|---:|---:|---:|---:|---:|
| thin singles | 304 | 304 | 0.335 (0.319-0.353) | 0.022 | 15.2 |
| thick singles | 108 | 108 | 0.288 (0.263-0.315) | 0.066 | 4.4 |
| thin multis | 222 | 101 | 0.122 (0.108-0.136) | 0.030 | 4.1 |
| thick multis | 69 | 28 | 0.117 (0.086-0.154) | 0.033 | 3.5 |

These values reproduce Sagear's printed reciprocal transit-selection convention. A
generative forward-normalized selection model gives 0.342, 0.316, 0.124, and 0.121 for
the same four groups. Transit selection therefore does not explain the discrepancy.

## Host-Clustered Robustness

Planets in a multi share a star, density prior, limb-darkening parameters, and parts of
the transit fit. A new diagnostic resamples and removes complete target systems rather
than treating every planet as independent. Results below use 5,000 trials under the
printed reciprocal convention.

| population | host-bootstrap 16th-84th | leave-10%-hosts within 5% | 95th percentile shift | largest one-host shift |
|---|---:|---:|---:|---:|
| thin singles | 0.313-0.357 | 97.5% | 4.1% | 2.5%, K04921 |
| thick singles | 0.252-0.320 | 83.6% | 8.3% | 8.4%, K04943 |
| thin multis | 0.099-0.139 | 67.6% | 10.7% | 8.8%, K02857 |
| thick multis | 0.072-0.181 | 31.9% | 34.8% | 34.3%, K01992 |

No group passes the strict requirement that every random leave-10% trial move the mean
by less than 5%. Thin singles are nevertheless comparatively stable, so their high value
is not the product of one target. The current thick-multi inference is strongly
host-dominated and cannot support a population claim.

## Complete Factorial Validation

The complete matrix contains 82 successful target-system FITS: 24 original-long-cadence,
24 reference-long-cadence, 9 original long-plus-short-cadence, 9 reference
long-plus-short-cadence, 8 repeat-seed, and 8 printed-prior runs. All expected FITS were
directly re-extracted and all matched planets passed direct extraction QC.

| comparison | planets / systems | median delta e | 95% system-bootstrap interval | median absolute delta e |
|---|---:|---:|---:|---:|
| reference versus original limb darkening | 34 / 24 | +0.00120 | -0.01034 to +0.00512 | 0.01545 |
| long plus short versus long, original LD | 13 / 9 | +0.00101 | -0.00059 to +0.01280 | 0.00682 |
| long plus short versus long, reference LD | 13 / 9 | -0.00081 | -0.04771 to +0.00263 | 0.00431 |
| printed priors versus public defaults | 9 / 8 | +0.00014 | -0.00053 to +0.00065 | 0.00065 |
| same configuration, different sampler seed | 9 / 8 | +0.00071 | -0.00108 to +0.00276 | 0.00108 |

The tested choices cause meaningful changes for a few systems, but not a coherent panel
shift. `K00283` and `K02533` deserve explicit target-level investigation; they must not
be allowed to drive a decision about the full population. See
`docs/full_factorial_validation_assessment.md` for the interpretation boundary, results,
and next decisions.

## Causes Assessed

| candidate cause | current assessment |
|---|---|
| Wrong reconstructed Toomre classifier | Real historical problem; solved for the primary analysis by the published host table. |
| Recounting multiplicity after cuts | Real historical bug; solved. It would mislabel 45 current planet rows. |
| Synthetic geometric impact samples | Real posterior bug; solved for the uniform direct subset by preserving paired ALDERAAN samples. |
| Limb-darkening prior centers | Measurable for individual systems, but the partial median shift is far too small to explain thin singles. |
| Transit-selection convention | Not the solution. The printed reciprocal rule lowers the current means only modestly and fails a generative injection test. |
| Berger 2018 versus 2020 density | Still not exactly reproducible. Radius-based 2018-like shifts are generally too small and inconsistent in direction to remove the high-e tail. |
| A few thin-single outliers | Not sufficient. Thin singles remain high under host-level resampling and top-leverage removal. |
| Hidden convergence or visual vetting | Still plausible and not publicly identifiable because the final planet and rejected-fit tables are unavailable. |
| ALDERAAN cadence, priors, and run stochasticity | Measured in the complete matrix. They affect some individual fits but do not explain the panel-wide discrepancy. |

## Scientific Versus Literal Replication

For literal replication, `manuscript_reciprocal` remains the comparison mode because it
matches the published analysis description. For a generative population analysis,
`legacy_forward_norm` is the defensible transit-selection model and recovers an injected
intrinsic Rayleigh distribution more accurately. Both results must be reported and must
not be blended into one headline number.

## Remaining Decision Gates

1. Inspect the high-leverage systems named in the factorial leverage tables.
2. Freeze one stellar-density construction and record it per planet.
3. Request Sagear's final planet list, visually rejected targets, density-prior table,
   and NumPyro likelihood if an exact numerical match remains impossible.
4. Re-extract and refit only after the density and inclusion contracts are fixed.

The current evidence rules out several easy explanations. It does not justify giving up.
It narrows the unresolved problem to a smaller and testable set of posterior-generation
and unpublished provenance choices.
