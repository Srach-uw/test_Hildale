# Factorial validation secondary diagnostics

These are paired, targeted validation diagnostics. They are not population estimates or
formal discovery tests. Nominal correlation p-values are uncorrected and are reported only
to identify engineering relationships worth checking in the completed matrix.

## Coverage

- QC-passing paired planets: 78
- Paired target systems: 24
- Available comparisons: 5
- Discovered FITS: 82/82

## Population summaries

| comparison_id         | effect                | population   |   n_planets |   n_systems |   median_delta_e |   median_delta_e_ci_low |   median_delta_e_ci_high |   mean_delta_e |   median_abs_delta_e |   p95_abs_delta_e |   fraction_positive |   median_baseline_e_width |   median_abs_delta_over_e_width | interpretation                    |
|:----------------------|:----------------------|:-------------|------------:|------------:|-----------------:|------------------------:|-------------------------:|---------------:|---------------------:|------------------:|--------------------:|--------------------------:|--------------------------------:|:----------------------------------|
| cadence_original_ld   | cadence               | thick_multi  |           2 |           1 |     -0.0393759   |            -0.0393759   |             -0.0393759   |   -0.0393759   |          0.0414073   |       0.0768456   |            0.5      |                  0.503537 |                      0.0813594  | underpowered_fewer_than_5_systems |
| cadence_original_ld   | cadence               | thick_single |           1 |           1 |      0.0128023   |             0.0128023   |              0.0128023   |    0.0128023   |          0.0128023   |       0.0128023   |            1        |                  0.333949 |                      0.0383361  | underpowered_fewer_than_5_systems |
| cadence_original_ld   | cadence               | thin_multi   |           5 |           2 |     -0.00129027  |            -0.00129027  |              0.0074216   |   -0.00599008  |          0.0141889   |       0.0504865   |            0.4      |                  0.452281 |                      0.0300328  | underpowered_fewer_than_5_systems |
| cadence_original_ld   | cadence               | thin_single  |           5 |           5 |      0.00101218  |            -0.000590662 |              0.0155285   |    0.00352041  |          0.00101218  |       0.012697    |            0.8      |                  0.349244 |                      0.00367068 | descriptive_targeted_sample       |
| cadence_reference_ld  | cadence               | thick_multi  |           2 |           1 |     -0.0437969   |            -0.0437969   |             -0.0437969   |   -0.0437969   |          0.0437969   |       0.0809317   |            0        |                  0.503038 |                      0.0843316  | underpowered_fewer_than_5_systems |
| cadence_reference_ld  | cadence               | thick_single |           1 |           1 |      0.00881325  |             0.00881325  |              0.00881325  |    0.00881325  |          0.00881325  |       0.00881325  |            1        |                  0.331705 |                      0.0265695  | underpowered_fewer_than_5_systems |
| cadence_reference_ld  | cadence               | thin_multi   |           5 |           2 |     -0.0477079   |            -0.126677    |              0.00128953  |   -0.0732052   |          0.0477079   |       0.180712    |            0.2      |                  0.464582 |                      0.0987442  | underpowered_fewer_than_5_systems |
| cadence_reference_ld  | cadence               | thin_single  |           5 |           5 |      0.000727055 |            -0.000806506 |              0.0346468   |    0.00709777  |          0.000806506 |       0.0279073   |            0.6      |                  0.357118 |                      0.00257715 | descriptive_targeted_sample       |
| ld_reference_lc       | limb_darkening        | thick_multi  |           5 |           2 |     -0.0152259   |            -0.0152259   |              0.027596    |    0.000400482 |          0.0167259   |       0.0630596   |            0.2      |                  0.49988  |                      0.0328455  | underpowered_fewer_than_5_systems |
| ld_reference_lc       | limb_darkening        | thick_single |           4 |           4 |      0.0119546   |            -0.00102639  |              0.0385964   |    0.0153698   |          0.0119546   |       0.0362319   |            0.75     |                  0.435949 |                      0.0352821  | underpowered_fewer_than_5_systems |
| ld_reference_lc       | limb_darkening        | thin_multi   |          12 |           5 |     -0.00567277  |            -0.0205862   |              0.00512379  |   -0.00588329  |          0.0231824   |       0.0575806   |            0.5      |                  0.472159 |                      0.050221   | descriptive_targeted_sample       |
| ld_reference_lc       | limb_darkening        | thin_single  |          13 |          13 |      0.00168083  |            -0.00748102  |              0.00522757  |   -0.00150763  |          0.00603627  |       0.0311337   |            0.615385 |                  0.349244 |                      0.0149272  | descriptive_targeted_sample       |
| paper_prior_ambiguity | prior_ambiguity       | thick_single |           1 |           1 |      0.00061023  |             0.00061023  |              0.00061023  |    0.00061023  |          0.00061023  |       0.00061023  |            1        |                  0.352381 |                      0.00173173 | underpowered_fewer_than_5_systems |
| paper_prior_ambiguity | prior_ambiguity       | thin_multi   |           2 |           1 |     -0.00111969  |            -0.00111969  |             -0.00111969  |   -0.00111969  |          0.002237    |       0.00324472  |            0.5      |                  0.466388 |                      0.00472814 | underpowered_fewer_than_5_systems |
| paper_prior_ambiguity | prior_ambiguity       | thin_single  |           6 |           6 |      0.000121146 |            -0.000912912 |              0.0238544   |    0.00768756  |          0.000589363 |       0.0356192   |            0.666667 |                  0.380382 |                      0.0014223  | descriptive_targeted_sample       |
| sampler_repeatability | sampler_repeatability | thick_single |           1 |           1 |      0.000886859 |             0.000886859 |              0.000886859 |    0.000886859 |          0.000886859 |       0.000886859 |            1        |                  0.352381 |                      0.00251676 | underpowered_fewer_than_5_systems |
| sampler_repeatability | sampler_repeatability | thin_multi   |           2 |           1 |      0.00840664  |             0.00840664  |              0.00840664  |    0.00840664  |          0.00840664  |       0.0134926   |            1        |                  0.466388 |                      0.0176747  | underpowered_fewer_than_5_systems |
| sampler_repeatability | sampler_repeatability | thin_single  |           6 |           6 |      4.93427e-05 |            -0.00153311  |              0.0174832   |    0.00533313  |          0.000894457 |       0.0261906   |            0.5      |                  0.380382 |                      0.00252327 | descriptive_targeted_sample       |

Any row with fewer than five independent systems is explicitly underpowered.

## Highest-leverage systems

| comparison_id         | koi_target   |   rank_by_max_abs_delta_e |   n_planets |   system_median_delta_e |   system_max_abs_delta_e |   full_median_delta_e |   leave_system_out_median_delta_e |   leave_system_out_shift |
|:----------------------|:-------------|--------------------------:|------------:|------------------------:|-------------------------:|----------------------:|----------------------------------:|-------------------------:|
| cadence_original_ld   | K02533       |                         1 |           2 |            -0.0393759   |              0.0807832   |           0.00101218  |                       0.00101218  |              0           |
| cadence_original_ld   | K00283       |                         2 |           3 |            -0.00129027  |              0.0576922   |           0.00101218  |                       0.00119171  |              0.000179535 |
| cadence_original_ld   | K01001       |                         3 |           2 |             0.0074216   |              0.0216639   |           0.00101218  |                       0.00101218  |              0           |
| cadence_original_ld   | K00716       |                         4 |           1 |             0.0155285   |              0.0155285   |           0.00101218  |                       0.00064649  |             -0.00036569  |
| cadence_original_ld   | K00064       |                         5 |           1 |             0.0128023   |              0.0128023   |           0.00101218  |                       0.00064649  |             -0.00036569  |
| cadence_reference_ld  | K00283       |                         1 |           3 |            -0.126677    |              0.19422     |          -0.000806506 |                       0.000349662 |              0.00115617  |
| cadence_reference_ld  | K02533       |                         2 |           2 |            -0.0437969   |              0.0850578   |          -0.000806506 |                      -2.77308e-05 |              0.000778775 |
| cadence_reference_ld  | K00716       |                         3 |           1 |             0.0346468   |              0.0346468   |          -0.000806506 |                      -0.00126629  |             -0.000459782 |
| cadence_reference_ld  | K00064       |                         4 |           1 |             0.00881325  |              0.00881325  |          -0.000806506 |                      -0.00126629  |             -0.000459782 |
| cadence_reference_ld  | K01001       |                         5 |           2 |             0.00128953  |              0.00430514  |          -0.000806506 |                      -0.000806506 |              0           |
| ld_reference_lc       | K00283       |                         1 |           3 |             0.00512379  |              0.0802484   |           0.00119777  |                       0.00107602  |             -0.000121755 |
| ld_reference_lc       | K02533       |                         2 |           2 |             0.027596    |              0.0719178   |           0.00119777  |                       0.00119777  |              0           |
| ld_reference_lc       | K01001       |                         3 |           2 |             0.0124323   |              0.0390342   |           0.00119777  |                       0.00119777  |              0           |
| ld_reference_lc       | K02109       |                         4 |           1 |             0.0385964   |              0.0385964   |           0.00119777  |                       0.00107602  |             -0.000121755 |
| ld_reference_lc       | K00108       |                         5 |           2 |             0.00248528  |              0.0378228   |           0.00119777  |                       0.00119777  |              0           |
| paper_prior_ambiguity | K00716       |                         1 |           1 |             0.0470604   |              0.0470604   |           0.000139402 |                       0.000121146 |             -1.82562e-05 |
| paper_prior_ambiguity | K01001       |                         2 |           2 |            -0.00111969  |              0.00335669  |           0.000139402 |                       0.000139402 |              0           |
| paper_prior_ambiguity | K00791       |                         3 |           1 |            -0.00129555  |              0.00129555  |           0.000139402 |                       0.000374816 |              0.000235414 |
| paper_prior_ambiguity | K00428       |                         4 |           1 |             0.000648451 |              0.000648451 |           0.000139402 |                       0.000121146 |             -1.82562e-05 |
| paper_prior_ambiguity | K00890       |                         5 |           1 |             0.00061023  |              0.00061023  |           0.000139402 |                       0.000121146 |             -1.82562e-05 |
| sampler_repeatability | K00716       |                         1 |           1 |             0.0342594   |              0.0342594   |           0.00070691  |                       0.000462989 |             -0.000243921 |
| sampler_repeatability | K01001       |                         2 |           2 |             0.00840664  |              0.0140577   |           0.00070691  |                       0.000219069 |             -0.000487841 |
| sampler_repeatability | K00428       |                         3 |           1 |            -0.00198421  |              0.00198421  |           0.00070691  |                       0.000796884 |              8.99745e-05 |
| sampler_repeatability | K01553       |                         4 |           1 |            -0.001082    |              0.001082    |           0.00070691  |                       0.000796884 |              8.99745e-05 |
| sampler_repeatability | K00890       |                         5 |           1 |             0.000886859 |              0.000886859 |           0.00070691  |                       0.000462989 |             -0.000243921 |

## Robustness

| comparison_id         | analysis                        |   n_planets |   n_systems |   median_delta_e |   shift_from_full_median |     shift_p95 |     shift_max |   n_trials |
|:----------------------|:--------------------------------|------------:|------------:|-----------------:|-------------------------:|--------------:|--------------:|-----------:|
| cadence_original_ld   | all_qc_pass                     |          13 |           9 |      0.00101218  |              0           | nan           | nan           |        nan |
| cadence_original_ld   | exclude_highest_leverage_system |          11 |           8 |      0.00101218  |              0           | nan           | nan           |        nan |
| cadence_original_ld   | impact_below_0p8_both_arms      |           9 |           5 |      0.000280801 |             -0.000731379 | nan           | nan           |        nan |
| cadence_original_ld   | e84_below_0p9_both_arms         |          12 |           8 |      0.000826025 |             -0.000186155 | nan           | nan           |        nan |
| cadence_original_ld   | leave_10pct_systems_out_trials  |          13 |           9 |      0.00101218  |              0.000179535 |   0.00036569  |   0.00036569  |       1000 |
| cadence_reference_ld  | all_qc_pass                     |          13 |           9 |     -0.000806506 |              0           | nan           | nan           |        nan |
| cadence_reference_ld  | exclude_highest_leverage_system |          10 |           8 |      0.000349662 |              0.00115617  | nan           | nan           |        nan |
| cadence_reference_ld  | impact_below_0p8_both_arms      |           8 |           5 |     -0.00213103  |             -0.00132453  | nan           | nan           |        nan |
| cadence_reference_ld  | e84_below_0p9_both_arms         |          12 |           8 |     -0.0008769   |             -7.03946e-05 | nan           | nan           |        nan |
| cadence_reference_ld  | leave_10pct_systems_out_trials  |          13 |           9 |     -0.000806506 |              0.000459782 |   0.00115617  |   0.00115617  |       1000 |
| ld_reference_lc       | all_qc_pass                     |          34 |          24 |      0.00119777  |              0           | nan           | nan           |        nan |
| ld_reference_lc       | exclude_highest_leverage_system |          31 |          23 |      0.00107602  |             -0.000121755 | nan           | nan           |        nan |
| ld_reference_lc       | impact_below_0p8_both_arms      |          28 |          18 |     -0.00101357  |             -0.00221134  | nan           | nan           |        nan |
| ld_reference_lc       | e84_below_0p9_both_arms         |          33 |          23 |      0.00107602  |             -0.000121755 | nan           | nan           |        nan |
| ld_reference_lc       | leave_10pct_systems_out_trials  |          34 |          24 |      0.00119777  |              0.000121755 |   0.00219851  |   0.00219851  |       1000 |
| paper_prior_ambiguity | all_qc_pass                     |           9 |           8 |      0.000139402 |              0           | nan           | nan           |        nan |
| paper_prior_ambiguity | exclude_highest_leverage_system |           8 |           7 |      0.000121146 |             -1.82562e-05 | nan           | nan           |        nan |
| paper_prior_ambiguity | impact_below_0p8_both_arms      |           7 |           6 |      0.000139402 |              0           | nan           | nan           |        nan |
| paper_prior_ambiguity | e84_below_0p9_both_arms         |           9 |           8 |      0.000139402 |              0           | nan           | nan           |        nan |
| paper_prior_ambiguity | leave_10pct_systems_out_trials  |           9 |           8 |      0.000139402 |              1.82562e-05 |   0.000235414 |   0.000235414 |       1000 |
| sampler_repeatability | all_qc_pass                     |           9 |           8 |      0.00070691  |              0           | nan           | nan           |        nan |
| sampler_repeatability | exclude_highest_leverage_system |           8 |           7 |      0.000462989 |             -0.000243921 | nan           | nan           |        nan |
| sampler_repeatability | impact_below_0p8_both_arms      |           7 |           6 |      0.00070691  |              0           | nan           | nan           |        nan |
| sampler_repeatability | e84_below_0p9_both_arms         |           9 |           8 |      0.00070691  |              0           | nan           | nan           |        nan |
| sampler_repeatability | leave_10pct_systems_out_trials  |           9 |           8 |      0.00070691  |              8.99745e-05 |   0.000487841 |   0.000487841 |       1000 |

## Largest exploratory correlations

| comparison_id         | response    | feature                     |   n_planets |   spearman_rho |   nominal_p_value | interpretation                           |
|:----------------------|:------------|:----------------------------|------------:|---------------:|------------------:|:-----------------------------------------|
| paper_prior_ambiguity | abs_delta_e | fractional_delta_rp_over_rs |           9 |      -0.733333 |         0.0245542 | exploratory_multiple_tests_not_corrected |
| paper_prior_ambiguity | abs_delta_e | delta_impact                |           9 |      -0.7      |         0.0357696 | exploratory_multiple_tests_not_corrected |
| cadence_original_ld   | delta_e     | catalog_impact              |           9 |       0.6      |         0.0876228 | exploratory_multiple_tests_not_corrected |
| cadence_reference_ld  | delta_e     | delta_impact                |          13 |      -0.538462 |         0.0576344 | exploratory_multiple_tests_not_corrected |
| ld_reference_lc       | delta_e     | catalog_impact              |          24 |       0.501739 |         0.0124872 | exploratory_multiple_tests_not_corrected |
| cadence_reference_ld  | abs_delta_e | ld_abs_max_delta            |           9 |      -0.5      |         0.170471  | exploratory_multiple_tests_not_corrected |
| paper_prior_ambiguity | delta_e     | fractional_delta_rp_over_rs |           9 |      -0.466667 |         0.205386  | exploratory_multiple_tests_not_corrected |
| cadence_reference_ld  | delta_e     | catalog_impact              |           9 |       0.466667 |         0.205386  | exploratory_multiple_tests_not_corrected |
| cadence_original_ld   | delta_e     | delta_t14_hr                |          13 |       0.450549 |         0.122332  | exploratory_multiple_tests_not_corrected |
| cadence_reference_ld  | delta_e     | e_width_baseline            |          13 |      -0.43956  |         0.13286   | exploratory_multiple_tests_not_corrected |
| cadence_original_ld   | abs_delta_e | e_width_baseline            |          13 |       0.417582 |         0.155675  | exploratory_multiple_tests_not_corrected |
| cadence_reference_ld  | abs_delta_e | e_width_baseline            |          13 |       0.417582 |         0.155675  | exploratory_multiple_tests_not_corrected |
