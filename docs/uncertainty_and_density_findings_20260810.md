# Uncertainty calibration and the Berger-2018 density branch

Updated: 2026-08-10

This document records the two findings that most change the interpretation of
the remaining discrepancy, plus five hypotheses closed by measurement. All
numbers were computed against the full 2,465-planet source-matched inventory.

## Scope, stated before the results

**Everything in sections 1 and 2 applies to the equal-raw-row diagnostic branch,
not to the canonical `LN_WT` branch.** This distinction is essential and is not
softened anywhere below.

The canonical dynesty-weighted branch remains discrepant by a wide margin
(thin singles about 0.246 against a published 0.022). No uncertainty
recalibration closes that: a 1.5x wider interval leaves the canonical value far
above the published upper bound of 0.029. The canonical branch is also the
statistically defensible one, and equal raw-row weighting still reverses the
published model comparison toward a half-Gaussian, which is why the repository
retains it only as a diagnostic (see README Status and
`docs/replication_status.md`).

What follows therefore **narrows where the remaining problem lives**. It does
not overturn the repository's headline result.

## Bottom line

Two results, in order of importance.

**1. Within the equal-row diagnostic branch, the residual thin-single difference
is not statistically significant.** The population fit's reported interval
understates its own uncertainty for thin singles by about 1.5x. Bootstrap
resampling over planets gives a 16-84 interval of 0.0175 to 0.0620, against a
fit-reported 0.0308 to 0.0603. With correctly propagated uncertainty the
comparison is **0.83 sigma**, and **27 percent of bootstrap resamples fall
inside the published 16-84 interval**. The apparent factor-of-two residual in
that branch was an artefact of an under-dispersed likelihood interval.

**2. WITHDRAWN, see section 6.** An earlier version of this document reported
that reconstructing stellar densities from Berger et al. (2018) radii lowers the
sixteen-value Table 3 chi-square from 183.2 to 54.0 and called that a genuine
improvement in density provenance. The chi-square measurement is correct, but
the attribution is not: the Berger-2018 reconstruction as implemented carries a
**defective uncertainty model**, and that defect is the strongest single
predictor of inferred eccentricity in the whole dataset. The improvement cannot
currently be credited to better density values.

**3. The finding that replaces it: stellar-density uncertainty is the dominant
driver of inferred eccentricity.** Across 2,465 planets the rank correlation
between fractional density uncertainty and per-planet median eccentricity is
**Spearman +0.52**, larger than any other variable tested in this project by a
wide margin (impact parameter is +0.13; host temperature, metallicity and
evolutionary state are all null). This is the mechanism behind the
singles-versus-multis pattern and behind the low information content of the
per-planet posteriors.

## 1. Uncertainty calibration

Bootstrap over planets, 60 replicates for thin singles, 30 for the others,
refitting with the unmodified production `fit_rayleigh`.

| Population | n | Point | Fit-reported 16-84 | Bootstrap 16-84 | Width ratio |
| --- | ---: | ---: | --- | --- | ---: |
| thin singles | 1109 | 0.0471 | 0.0308-0.0603 | **0.0175-0.0620** | **1.51x** |
| thick singles | 269 | 0.0430 | 0.0174-0.0780 | 0.0096-0.0566 | 0.77x |
| thick multis | 209 | 0.0098 | 0.0060-0.0250 | 0.0062-0.0200 | 0.73x |
| thin multis | 878 | 0.0243 | 0.0186-0.0310 | 0.0199-0.0274 | 0.60x |

The likelihood-based interval is **too narrow for thin singles** and somewhat
conservative for the other three. Thin singles is also the population where the
apparent discrepancy sat, so this matters directly.

Thin singles against the published value, three ways of stating the same
comparison:

| Statistic | Value |
| --- | ---: |
| z, published uncertainty only | +3.63 |
| z, both uncertainties, fit-reported sigma | +1.42 |
| **z, both uncertainties, bootstrap sigma** | **+0.83** |
| P(bootstrap replicate <= 0.022) | 0.267 |
| P(bootstrap replicate inside published 16-84) | 0.267 |
| Bootstrap range over 60 replicates | 0.0111-0.0774 |

**Recommendation:** report bootstrap intervals alongside the likelihood
intervals for any population-level number quoted in the release, and do not
describe the thin-single result as a significant discrepancy without that
qualification.

## 2. The Berger-2018 density branch

Reconstructing `rho_star` from Berger 2018 radii with a catalog surface gravity,
via `rho = 3g / (4 pi G R)`, rather than adopting Berger 2020 densities
directly.

| Population | B20 branch | **B18 branch** | Published |
| --- | ---: | ---: | ---: |
| thick singles | 0.0556 | 0.0431 | 0.066 |
| thin singles | 0.0721 | **0.0474** | 0.022 |
| thick multis | 0.0181 | 0.0099 | 0.033 |
| thin multis | 0.0226 | 0.0244 | 0.030 |

Sixteen-value asymmetric chi-square against all four published model families:

| Branch | chi2 (16 values) | Values at abs(z) > 2 |
| --- | ---: | ---: |
| Berger 2020 | 183.2 | not recorded |
| **Berger 2018 reconstruction** | **54.0** | 3 of 16 |

Caveat: two of the four full-sample fits (thick singles, thick multis) are
boundary-flagged and should be treated as poorly constrained rather than as
agreement.

## 3. Hypotheses closed by measurement

Each of these is closed by a number, not by argument.

| Hypothesis | Test | Result |
| --- | --- | --- |
| Some global convention (density, weighting, normalisation) still reconciles all four populations | Rescale that fixes thin singles is x0.464; apply to the others | thick singles z = -2.19, thin multis z = -2.67. **No global change can work; the residual is population-specific.** |
| The unpublished rejection list explains the thin-single excess | Remove the most eccentric planets and refit | Removing the top **30 percent** moves the median per-planet e50 from 0.4299 to 0.4171. The paper removes about 1 percent. **Closed.** |
| Table 3 rows are mis-ordered or mis-transcribed | Recompute the 24-permutation audit on the B18 branch | Identity-mapping disadvantage falls from about **9x to 1.4-2.1x**, which is what best-of-24 selection yields from noise. **Withdrawn.** The anomaly weakens as the density input improves, the opposite of a real transcription error. |
| Evolved stars drive the thick-over-thin separation | Split by Berger `evol_b18` and compare per-planet posteriors | Host-level: thick hosts are evolved 27.9% vs thin 15.3%, p = 1.4e-07, a real structural difference. Per-planet posteriors: **no significant difference in any population** (Mann-Whitney p = 0.25 to 1.00). Subset population fits differ but are boundary-flagged and unstable. **Not supported.** |
| Distance/brightness selection makes thick hosts noisier | Compare Kepler magnitude and SNR by disk | thick hosts 0.12 mag fainter (p = 0.016); thick singles 14% lower SNR (p = 0.081). Direction is right, **magnitude far too small** to produce a 3x effect. |

## 4. An internal-consistency test worth keeping

Comparing the measured-RV subsample (Figure 9) to the full sample within a
single convention cancels every multiplicative choice and isolates sample
structure:

| Population | Published Fig9/Table3 | This branch, Fig9/full |
| --- | ---: | ---: |
| thick singles | 1.46 | 1.80 |
| thin singles | 0.97 | 0.34 |
| thick multis | 1.50 | 0.69 |
| thin multis | 0.49 | 1.20 |

Three of four disagree and two reverse sign. Any candidate convention should be
required to pass this test before it is treated as an explanation. The
Figure 9 values are digitised from a plot and carry reading error, but the
discrepancies here are factors of two to three.

## 5. What remains genuinely open

Restating the scope: within the **equal-row diagnostic branch** the four
populations are now consistent or marginally consistent with the paper. Within
the **canonical `LN_WT` branch** they are not, and nothing here changes that.

Two things remain open regardless of branch.

**The size of the thick-over-thin separation is not independently confirmed.**
This reconstruction measures thick/thin = 0.91x in single-planet systems where
the paper reports 3.0x, and the two populations are statistically
indistinguishable at the per-planet level (Kolmogorov-Smirnov p = 0.57 on e50,
p = 0.18 on e16). Given the interval calibration in section 1, that non-detection
is best described as **insufficient sensitivity in the public-data
reconstruction** rather than as a contradiction of the published result.

**The weighting question is still the pivot.** The equal-row branch is the one
that approaches the published values, and it is the statistically weaker choice:
the raw ALDERAAN rows are nested-sampling points with a median effective sample
size of about 40 percent of the row count, so treating them equally
over-represents low-likelihood prior volume. Whether the published analysis
consumed exported rows directly or reweighted them is the single highest-value
question remaining, and it is answerable by the authors in one sentence.

Establishing either point to the precision the paper reports would require the
planet-level posterior products that remain unreleased.

## 6. Correction: the Berger-2018 branch has a defective uncertainty model

This section withdraws the density-provenance claim made earlier in this
document and records what replaced it.

### What is wrong

The Berger-2018 reconstruction (`B18_KG_LOGG_FIXED_EQUAL`) produces fractional
stellar-density uncertainties that are **bimodal** and far too large:

| Quantity | B18-KG reconstruction | Canonical Berger 2020 | Berger published |
| --- | --- | --- | --- |
| distribution | **bimodal: 0.315 and 1.712** | unimodal | single value |
| median fractional sigma | 0.476 | **0.101** | 0.13 |
| planets with sigma > rho | **35.2 percent** | **0.0 percent** | none expected |
| singles vs multis median | **1.696 vs 0.320 (5.3x)** | 0.107 vs 0.090 (1.19x) | no dependence expected |

Two things are wrong independently. The high cluster sits at **13.2x** Berger's
published 13 percent uncertainty, and a density prior with sigma larger than rho
is not a constraint at all. And the split tracks **planet multiplicity**, which
is physically impossible: stellar density is a property of the star and cannot
depend on how many planets transit it.

The canonical Berger-2020 path is sound and was verified in source:
`absolute_density_error(x) = 10**x` (`scripts/alderaan_shape_diagnostics.py:207`)
correctly interprets the CDS `E_rho` and `e_rho` columns as log10 of the linear
uncertainty, recovering a median of 0.101 against Berger's published 0.13.

### Why it matters

The defect is not cosmetic. Fractional density uncertainty is the strongest
predictor of inferred eccentricity in the dataset:

| Split | n | median fractional sigma | median per-planet e50 |
| --- | ---: | ---: | ---: |
| low-uncertainty cluster | 1574 | 0.318 | 0.3156 |
| high-uncertainty cluster | 891 | 1.712 | **0.4393** |

Mann-Whitney p = 6e-149; Spearman(fractional sigma, e50) = **+0.52**.

The inflated cluster is concentrated in singles: **58.0 percent of thick singles
and 53.7 percent of thin singles**, against 13.4 and 12.8 percent of the two
multi populations. **That alone explains why this reconstruction reports singles
as roughly twice as eccentric as multis**, without invoking any dynamics.

### Consequences

- The chi-square improvement from 183.2 to 54.0 is a real measurement but
  **cannot be attributed to better density provenance**. It is confounded with a
  four-to-thirteen-fold inflation of the density uncertainties.
- The Berger-2018 **values** may still be the more source-faithful reading, since
  the manuscript states Berger 2018. That part is not withdrawn. Only the
  uncertainty model is defective.
- **Required follow-up:** rebuild the Berger-2018 branch propagating
  uncertainties to Berger's published scale, confirm the bimodality and the
  multiplicity dependence are gone, and only then re-evaluate the chi-square.
  Until that is done, no result from this branch should be quoted.
- The bootstrap calibration in section 1 was computed on this branch and
  inherits the caveat. The finding that the reported interval is under-dispersed
  is a property of the estimator rather than of the density input, so it is
  expected to survive, but it should be recomputed on the corrected branch.

### Information content, and why nothing else worked

Measured directly as the Kullback-Leibler divergence of each per-planet
eccentricity posterior from its uniform proposal prior:

| Population | median KL (nats) |
| --- | ---: |
| thick singles | **0.019** |
| thin singles | **0.043** |
| thick multis | 0.141 |
| thin multis | 0.132 |

A KL of 0.02 nats means the posterior is within a few percent of the prior: that
planet contributes almost nothing. **Singles carry three to seven times less
information than multis**, and thick singles least of all, which is precisely
the population carrying the paper's headline claim.

This is the unifying explanation for the whole investigation. It is why the
models cannot be discriminated (see section 7), why the fit intervals are
under-dispersed, why thin and thick singles are statistically indistinguishable,
and why no global convention could ever have reconciled the four populations.

## 7. The model-comparison constraint, which has not been met

The paper states a hard, quantitative result at `main.tex:215`: the Rayleigh
distribution is the best fit in **both** disks, with Beta and monotonic Beta
disfavoured at **delta-BIC of about 20** and half-Gaussian at **delta-BIC above
25**.

Recomputed on the Berger-2018 branch, using matched files (identical n and
identical likelihood normalisation, both checked):

| Population | dBIC Beta | dBIC mono-Beta | dBIC half-Gaussian | Rayleigh best |
| --- | ---: | ---: | ---: | :---: |
| thick singles | +5.2 | +1.0 | +0.7 | yes |
| thin singles | +3.6 | -2.4 | -1.9 | no |
| thick multis | +5.2 | -1.8 | -1.5 | no |
| thin multis | +6.1 | -0.7 | -0.7 | no |

Maximum absolute delta-BIC is **6.1**, against a published 20 to 25. **These
posteriors cannot discriminate the four model families at all**, and Rayleigh is
preferred in only one of four populations.

This is a demanding test because delta-BIC measures discriminating power rather
than any eccentricity value, so it is immune to every multiplicative convention
that has been explored. It should be adopted as a standard acceptance criterion:
a candidate reconstruction that cannot reproduce the published model ranking has
not reproduced the paper's posteriors, whatever its Table 3 values look like.

An equivalent comparison for the canonical branch could not be made. No matched
pair of files exists (the available Rayleigh and multi-model outputs come from
different runs with different planet counts, 1109 against 1105, and different
reported mean eccentricities), and comparing them would produce a meaningless
number. Generating a matched canonical pair is a worthwhile follow-up.

## Reproducing

```bash
python hierarchical_table3_order_diagnostic.py --summary outputs/eccentricity_posterior_summary_FULL2465_B18_KG_LOGG_FIXED_EQUAL_20260810.csv --tag B18KG_EQUAL_FULL2465 --selection-mode manuscript_reciprocal --allow-non-dynesty-weights
```

Supporting artifacts are in `metadata/uncertainty_calibration_20260810/`.
