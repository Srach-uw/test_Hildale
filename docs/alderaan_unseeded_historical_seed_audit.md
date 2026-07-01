# Unseeded Transit Seed Recovery Audit

This checks local historical KOI/TCE catalogs for the 27 needed ALDERAAN planets that lack current KOI depth seeds.

## Summary By Catalog

| catalog                      |   matched_rows |   valid_depth_rows |   valid_depth_and_duration_consistent_rows |
|:-----------------------------|---------------:|-------------------:|-------------------------------------------:|
| alderaan_cumulative_20240816 |             27 |                  0 |                                          0 |
| current_cumulative_2026      |             27 |                  0 |                                          0 |
| dr22_mullally_q1_q16         |             16 |                  1 |                                          0 |
| dr24_coughlin_q1_q17         |             11 |                  3 |                                          0 |
| dr25_thompson_q1_q17         |              4 |                  0 |                                          0 |
| merged_planets_full          |             27 |                  0 |                                          0 |

## Rows With Historical Depth

| kepoi_name   | catalog              | candidate_disposition   | candidate_pdisposition   |   current_duration |   candidate_duration |   duration_ratio_candidate_to_current |   candidate_depth |   candidate_impact |   candidate_snr | autofill_recommended   |
|:-------------|:---------------------|:------------------------|:-------------------------|-------------------:|---------------------:|--------------------------------------:|------------------:|-------------------:|----------------:|:-----------------------|
| K03722.01    | dr22_mullally_q1_q16 | NOT DISPOSITIONED       | NOT DISPOSITIONED        |             3.1500 |              26.6680 |                                8.4660 |        40434.0000 |             1.1780 |         85.0000 | False                  |
| K03722.01    | dr24_coughlin_q1_q17 | FALSE POSITIVE          | FALSE POSITIVE           |             3.1500 |              26.6680 |                                8.4660 |        40434.0000 |             1.1780 |         85.0000 | False                  |
| K06698.01    | dr24_coughlin_q1_q17 | FALSE POSITIVE          | FALSE POSITIVE           |             4.5600 |              55.4840 |                               12.1675 |        48875.0000 |             8.4366 |        160.9000 | False                  |
| K07368.01    | dr24_coughlin_q1_q17 | CONFIRMED               | CANDIDATE                |             2.8000 |              12.8700 |                                4.5964 |         1723.0000 |             0.8417 |          4.8000 | False                  |

Conclusion: historical depths exist for a few blocked planets, but none pass the simple automatic-fill rule requiring depth plus duration agreement within 25 percent of the current catalog row. These should stay manual/blocked for now.
