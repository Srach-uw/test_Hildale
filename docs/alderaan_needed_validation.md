# ALDERAAN Needed Manifest Validation

## Pass/Fail Checks

| check                                      |   value | passed   | note                                                                         |
|:-------------------------------------------|--------:|:---------|:-----------------------------------------------------------------------------|
| sample_rows                                |    2474 | True     | Canonical best sample rows.                                                  |
| sample_unique_kepoi                        |    2474 | True     | No duplicate planet IDs expected.                                            |
| all_status_rows                            |    2474 | True     | Status table should mirror the sample.                                       |
| all_status_unique_kepoi                    |    2474 | True     | No duplicate status planet IDs.                                              |
| needed_planets_rows                        |    1065 | True     | Needed rows must equal status needs_alderaan count.                          |
| runnable_needed_planets_rows               |    1038 | True     | Runnable needed rows must have valid launch seeds.                           |
| unseeded_needed_planets_rows               |      27 | True     | Unseeded needed rows should be tracked separately.                           |
| needed_unique_kepoi                        |    1065 | True     | No duplicate needed planet IDs.                                              |
| needed_targets_rows                        |     889 | True     | One row per KOI target.                                                      |
| runnable_planets_missing_target_row        |       0 | True     | Every runnable needed planet target must be in target manifest.              |
| target_rows_without_runnable_needed_planet |       0 | True     | Every target row should be justified by at least one runnable needed planet. |
| expanded_catalog_rows                      |    1228 | True     | Expanded rows should include all sample planets in needed systems.           |
| alderaan_catalog_rows                      |    1225 | True     | ALDERAAN catalog should include every target system.                         |
| runnable_planets_missing_from_catalog      |       0 | True     | Every runnable flagged planet should have valid ALDERAAN catalog parameters. |
| catalog_required_field_failures            |       0 | True     | Required catalog fields must be finite and positive where appropriate.       |
| catalog_npl_mismatches                     |       0 | True     | Catalog npl should match rows per target.                                    |
| system_expansion_extra_sibling_planets     |     160 | True     | Non-needed planets added because ALDERAAN fits whole systems.                |

## Summary

- Needed planet rows: 1065
- Runnable needed planet rows: 1038
- Unseeded needed planet rows: 27
- Unique KOI systems to run/refit: 889
- ALDERAAN catalog rows after system expansion: 1225
- Extra sibling planets included by system expansion: 160
- Failed validation checks: 0

## Needed Planets By Tier

|   priority_tier | recommended_action                       |   planets |
|----------------:|:-----------------------------------------|----------:|
|               1 | rerun_results_file_bad_match_or_bad_zeta |        39 |
|               1 | run_missing_no_results_file              |       706 |
|               2 | refit_hard_shape_pathology               |       109 |
|               3 | refit_high_e_zeta_tail                   |       204 |
|               4 | refit_low_snr_high_e_tail                |         7 |

## Needed Targets By Tier

|   min_priority_tier | target_run_type                    |   targets |
|--------------------:|:-----------------------------------|----------:|
|                   1 | required_missing_or_bad_extraction |       593 |
|                   2 | required_refit_hard_shape          |        99 |
|                   3 | recommended_refit_tail             |       191 |
|                   4 | recommended_refit_tail             |         6 |

## Needed Planets By Population

| disk   | system   |   needed_planets |
|:-------|:---------|-----------------:|
| thick  | multi    |              100 |
| thick  | single   |              156 |
| thin   | multi    |              302 |
| thin   | single   |              507 |
