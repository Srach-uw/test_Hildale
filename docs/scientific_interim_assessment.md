# Scientific Interim Assessment

Updated: 2026-08-08

## Executive Finding

The public-data replication is technically complete but does not reproduce the
four population results. The exact published host classifications and the
corrected multiplicity contract removed two major upstream ambiguities. A
real-data Gilbert control now verifies that the extraction and hierarchy can
recover an external published eccentricity scale. The remaining Sagear
mismatch is centered on unreleased posterior construction, resampling, and
target-level quality decisions.

The July uniformly processed diagnostic began with 710 planets and retained
703 across the four hierarchy cells after its deterministic exclusions. That
historical subset remains much more eccentric than the published result; it is
not the current 2,465-planet reconstruction. The individual factorial arms
quantify short-cadence, limb-darkening, prior, and nested-sampling effects.
None creates a panel-wide shift large enough to explain the discrepancy on its
own. The closest-to-paper combined arm is also complete: its 13 matched planets
have a median eccentricity shift of +0.00042, so it does not supply the missing
population-wide shift.

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

## Current Full-Sample Diagnostic

The full 2,465-planet branch uses dynesty weights, row-paired ALDERAAN `T14`,
`Rp/Rstar`, impact, and period samples, the exact MacDougall density equation,
and deterministic posterior QC.

| population | reconstructed planets | weighted Rayleigh mean e | published planets | published mean e |
|---|---:|---:|---:|---:|
| thin singles | 1,109 | 0.246 | 1,121 | 0.022 |
| thick singles | 269 | 0.234 | 275 | 0.066 |
| thin multis | 878 | 0.171 | 862 | 0.030 |
| thick multis | 209 | 0.141 | 207 | 0.033 |

The total planet count is exact, but the four cell memberships are not. Equal
raw-row weighting gives values closer to the paper in two cells, but its BIC
comparison favors a half-Gaussian in all four populations. Sagear reports a
strong Rayleigh preference. The equal-row branch is therefore a diagnostic of
an unresolved posterior-export convention, not a valid replication result.

The earlier 703-planet uniform subset remains useful as a historical stress
test. It is not the current headline sample.

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
| Limb-darkening prior centers | Measurable for individual systems, but the completed single-factor median shift is far too small to explain thin singles. |
| Transit-selection convention | Not the solution. The printed reciprocal rule lowers the current means only modestly and fails a generative injection test. |
| Berger 2018 versus 2020 density | Still not exactly reproducible. Radius-based 2018-like shifts are generally too small and inconsistent in direction to remove the high-e tail. |
| A few thin-single outliers | Not sufficient. Thin singles remain high under host-level resampling and top-leverage removal. |
| Hidden convergence or visual vetting | Still plausible and not publicly identifiable because the final planet and rejected-fit tables are unavailable. |
| ALDERAAN cadence, priors, and run stochasticity | Measured in the complete matrix and nine-system combined confirmation. They affect some individual fits but do not explain the panel-wide discrepancy. |
| Generic extraction or hierarchy failure | Disfavored by the real-ALDERAAN Gilbert control, which recovers Beta means of 0.0487 for all small planets and 0.0684 for observed singles. |
| Importing Gilbert's full quality cuts | Not justified for the Sagear sample and does not recover Table 3 when applied as a labeled sensitivity. |

## Circular-Density Diagnostic

The July 27 posterior ledger permits a direct check using paired ALDERAAN
transit samples rather than DR25 point estimates. For each planet, the audit
compares the weighted circular-density posterior inferred from `T14`,
`Rp/Rstar`, and impact parameter with the adopted stellar density.

| population | planets | median signed delta log10 rho | median absolute delta | median posterior width |
|---|---:|---:|---:|---:|
| thin singles | 1,109 | +0.085 | 0.225 | 0.698 |
| thick singles | 269 | +0.080 | 0.217 | 0.712 |
| thin multis | 878 | -0.024 | 0.159 | 0.622 |
| thick multis | 209 | -0.064 | 0.183 | 0.611 |

The single-planet samples show a positive signed shift, while both multi-planet
samples are near zero or slightly negative. All four groups have broad
circular-density posteriors. This confirms that the disagreement is already
present in the transit-shape and stellar-density inputs.

The raw weighted FITS allow a second, population-wide test of impact-parameter
uncertainty. Using `b84 - b16 <= 0.4` as an exploratory constrained threshold:

| population | constrained / broad N | constrained median absolute delta | broad median absolute delta | constrained / broad median e50 |
|---|---:|---:|---:|---:|
| thin singles | 119 / 990 | 0.268 | 0.219 | 0.297 / 0.297 |
| thick singles | 18 / 251 | 0.191 | 0.221 | 0.258 / 0.290 |
| thin multis | 66 / 812 | 0.156 | 0.160 | 0.230 / 0.237 |
| thick multis | 14 / 195 | 0.558 | 0.167 | 0.385 / 0.242 |

Broad impact posteriors increase circular-density uncertainty for most groups.
A hierarchy fit to the constrained subsets provides a diagnostic comparison:

| population | constrained N | forward-normalized mean e | arXiv-v1 reciprocal sensitivity | paper mean e |
|---|---:|---:|---:|---:|
| thin singles | 114 | 0.333 (0.317-0.350) | 0.403 (0.375-0.436) | 0.022 |
| thick singles | 18 | 0.312 (0.276-0.356) | 0.349 (0.297-0.424) | 0.066 |
| thin multis | 61 | 0.214 (0.195-0.236) | 0.236 (0.212-0.265) | 0.030 |
| thick multis | 14 | 0.248 (0.193-0.308) | 0.361 (0.259-0.602) | 0.033 |

This subset is selected on a posterior property and is therefore neither an
unbiased population sample nor a replacement population measurement. High
inferred eccentricity persists among planets with narrower impact posteriors,
so simple posterior broadening alone is not supported as the full explanation.
Five thin-single and five thin-multi constrained rows fail the primary
importance-sampling QC and are excluded from the hierarchy table. The selection
can still interact with transit geometry or fit quality, and the thick-multi
cell has only 14 planets and is especially uncertain.

Target-level fit quality, the adopted density construction, and the unpublished
rejection contract therefore remain live explanations.

The compact result and generating script are
`metadata/public_reconstruction_20260727/photoeccentric_density_audit.csv` and
`scripts/photoeccentric_density_audit.py`. The impact-width stratification is
`metadata/public_reconstruction_20260727/photoeccentric_impact_width_audit.csv`.
The corresponding hierarchy sensitivity is
`metadata/public_reconstruction_20260727/rayleigh_impact_constrained_sensitivity.csv`.

## Scientific Versus Literal Replication

For literal replication, `manuscript_reciprocal` remains a labeled sensitivity because
it follows the hierarchy equation preserved in the arXiv v1 source comments. That
equation is not printed in the final journal PDF and fails the repository's
uninformative-data normalization check. For a generative population analysis,
`legacy_forward_norm` is the defensible transit-selection model and recovers an injected
intrinsic Rayleigh distribution more accurately. Both results must be reported and must
not be blended into one headline number.

## Remaining Reproduction Inputs

1. Sagear's final included and rejected KOI list.
2. The stellar-density and planet-radius tables used in the eccentricity step.
3. The intermediate `(e, omega)` posteriors or the exact nested-point SIR rule.
4. The population-fitting implementation and Table 3 array ordering.

The public-data investigation has exhausted the defensible internal tests
identified so far. Further tuning against the published values would weaken,
rather than improve, the replication. The next scientifically useful step is a
focused request for these products, followed by a preregistered rerun of the
affected stages.
