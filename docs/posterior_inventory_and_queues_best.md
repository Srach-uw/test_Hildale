# Posterior Inventory And Operational Queues

This clarifies that the existing ALDERAAN posterior archive already covers most of the current best sample. The ALDERAAN-needed manifest mixes truly missing posteriors with existing posteriors that are suspicious enough to inspect/refit.

## Queue Summary

| queue                                      |   planets |   targets |   thin_single |   thick_single |   thin_multi |   thick_multi |
|:-------------------------------------------|----------:|----------:|--------------:|---------------:|-------------:|--------------:|
| baseline_usable_existing_posteriors        |      1409 |      1069 |           602 |            116 |          558 |           133 |
| existing_posterior_flagged_review_or_refit |       320 |       302 |           154 |             33 |          106 |            27 |
| missing_extracted_posterior_total          |       745 |       617 |           353 |            123 |          196 |            73 |
| missing_launchable_alderaan_now            |       718 |       592 |           334 |            120 |          193 |            71 |
| missing_unseeded_manual_seed_needed        |        27 |        27 |            19 |              3 |            3 |             2 |

## Category Status Versus Sagear

| disk   | system   |   sagear_planets |   current_sample_planets |   delta_current_minus_sagear |   extracted_posteriors |   coverage_fraction |   baseline_usable_existing |   existing_flagged_review_or_refit |   missing_posterior |   missing_launchable_now |   missing_unseeded |
|:-------|:---------|-----------------:|-------------------------:|-----------------------------:|-----------------------:|--------------------:|---------------------------:|-----------------------------------:|--------------------:|-------------------------:|-------------------:|
| thick  | multi    |              207 |                      233 |                           26 |                    160 |               0.687 |                        133 |                                 27 |                  73 |                       71 |                  2 |
| thick  | single   |              275 |                      272 |                           -3 |                    149 |               0.548 |                        116 |                                 33 |                 123 |                      120 |                  3 |
| thin   | multi    |              862 |                      860 |                           -2 |                    664 |               0.772 |                        558 |                                106 |                 196 |                      193 |                  3 |
| thin   | single   |             1121 |                     1109 |                          -12 |                    756 |               0.682 |                        602 |                                154 |                 353 |                      334 |                 19 |

## Existing Posterior Flag Reasons

| flag                        |   flagged_existing_posterior_planets |
|:----------------------------|-------------------------------------:|
| flag_zeta_lt_0p7            |                                  242 |
| flag_e50_gt_0p5             |                                  239 |
| flag_low_snr_lt_20          |                                  178 |
| flag_ror_shift_25pct        |                                   79 |
| flag_duration_shift_25pct   |                                   65 |
| flag_zeta_gt_1p3            |                                   23 |
| flag_impact_gt_1            |                                    8 |
| flag_results_file_bad_match |                                    0 |

## Interpretation

- The strongest current resource is `existing_eccentricity_posteriors_best.csv`: 1729 planets already have extracted e/zeta summaries and posterior FITS references.
- The immediately usable baseline subset is 1409 planets with extracted posteriors and no current automated refit/stress-test flag.
- The truly missing-posterior set is 745 planets. Of those, 718 are launchable now and 27 are blocked by missing transit-depth seeds.
- The 320 existing-but-flagged posteriors should be reviewed/refit after the missing-posterior queue, unless a specific flagged system is driving the thin-single discrepancy.
