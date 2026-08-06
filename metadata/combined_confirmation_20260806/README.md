# Combined Confirmation Analysis

These files regenerate the seven-arm comparison from the immutable 82-fit
factorial release and the nine-system combined confirmation release.

The paired comparison of interest is
`paper_prior_conditional_reference_lcsc`, with `reference_lcsc` as the baseline.
It contains 13 planets in nine systems and no direct-posterior QC exclusions.

| parameter | median signed change | median absolute change | repeatability 95th percentile |
|---|---:|---:|---:|
| eccentricity | +0.00042 | 0.01135 | 0.02476 |
| zeta | -0.00019 | 0.00982 | 0.01295 |
| impact parameter | -0.00144 | 0.00740 | 0.07934 |
| transit duration (hr) | -0.00454 | 0.01110 | 0.02020 |
| radius ratio | +0.000017 | 0.000037 | 0.05727 |

Five of 13 eccentricity shifts exceed the small repeat-run sample's 95th
percentile, but the median shift is near zero and the changes are not coherent
in sign. This targeted arm does not explain the population-wide eccentricity
discrepancy. The repeatability comparison is exploratory, not a p-value.

The output was generated with 150,000 proposals per planet, 10,000
target-system bootstrap replicates, and seed `20260715`.
