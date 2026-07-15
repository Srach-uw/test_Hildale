# Scientific Interim Assessment

Updated: 2026-07-13

## Executive Finding

The replication is not at a dead end, but the current four-population result is not
publishable. The exact published host classifications and the corrected multiplicity
contract have removed two major upstream ambiguities. The remaining mismatch is centered
on the eccentricity posteriors and the target-level quality decisions that created
Sagear et al.'s final posterior sample.

The current uniformly processed 710-planet subset remains much more eccentric than the
published result. The discrepancy is not explained by limb darkening, transit-selection
weighting, the reconstructed disk classifier, or a small number of thin-single systems.
The running factorial validation is still needed to quantify short-cadence, prior, and
nested-sampling effects.

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

## Factorial Validation So Far

The completed partial snapshot contains 24 original-limb-darkening systems and 23
reference-limb-darkening counterparts, giving 31 paired planets in 23 systems. Reference
minus original results are:

| diagnostic | result |
|---|---:|
| median signed change in e50 | +0.00168 |
| system-bootstrap interval | -0.00103 to +0.00523 |
| median absolute change in e50 | 0.01568 |
| 95th percentile absolute change | 0.05548 |
| thin-single median signed change | +0.00168 |
| thin-single median absolute change | 0.00604 |

The global median is stable to removing the highest-leverage system and to impact and
posterior-width sensitivity cuts. Larger eccentricity changes are concentrated in a few
high-impact or duration-sensitive planets. Exploratory correlations suggest that signed
eccentricity change increases with catalog impact parameter, but this is a targeted
sample with multiple tests and is not a population result.

The final reference-LD system, K01127, contains three thick-multi planets. The partial
thick-multi LD result currently has only one system and must not be interpreted until
K01127 and the remaining validation arms are available.

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
| ALDERAAN cadence, priors, and run stochasticity | Actively being measured by the factorial matrix. |

## Scientific Versus Literal Replication

For literal replication, `manuscript_reciprocal` remains the comparison mode because it
matches the published analysis description. For a generative population analysis,
`legacy_forward_norm` is the defensible transit-selection model and recovers an injected
intrinsic Rayleigh distribution more accurately. Both results must be reported and must
not be blended into one headline number.

## Remaining Decision Gates

1. Complete the factorial repeatability, cadence, and printed-prior arms.
2. Require matched target-level comparisons and quantify arm effects relative to repeat
   seed variability.
3. Inspect the high-leverage systems named in the host and planet leverage tables.
4. Freeze one stellar-density construction and record it per planet.
5. Re-extract the selected validation arm with deterministic QC and rerun both literal
   and generative hierarchical fits.
6. Request Sagear's final planet list, visually rejected targets, density-prior table,
   and NumPyro likelihood if an exact numerical match remains impossible.

The current evidence rules out several easy explanations. It does not justify giving up.
It narrows the unresolved problem to a smaller and testable set of posterior-generation
and unpublished provenance choices.
