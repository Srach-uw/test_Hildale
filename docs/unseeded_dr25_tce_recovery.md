# Unseeded-Planet Recovery via DR25 TCE Table (2026-07-02)

Follow-up to `alderaan_unseeded_historical_seed_audit.md`, which found no local
historical catalog could seed any of the 27 depth-blocked planets. This audit
queried the NASA Exoplanet Archive **Q1–Q17 DR25 TCE table** (independent
pipeline transit fits) for all 27 KICs, matching by KIC + period within 1%.
Evidence: `metadata/unseeded_dr25_tce_matches.csv`.

## Result

| planet | verdict | DR25 seed depth (ppm) | note |
|---|---|---:|---|
| K07973.01 | **recoverable** | 1063 | duration ratio 1.0005, SNR 37 — textbook match |
| K07368.01 | **recoverable** | 657.7 | duration ratio 1.13, SNR 32 |
| K06698.01 | **recoverable (with caution)** | ~2660 | two consistent TCEs (2686/2640 ppm) at the same period; a third divergent TCE exists |
| K03722.01 | stays blocked | — | two conflicting depths (20690 vs 1300 ppm), both duration ratios > 1.3 |
| other 23 | stay blocked | — | no period-matched DR25 TCE at all |

**K01316.02 and K06516.02 are NOT recoverable** (no DR25 TCE at their periods),
so the incomplete-system caveat for K01316/K06516 in `README.md` stands:
flag those two targets' results for review after the cloud run.

## Recommendation

Do **not** regenerate the current 592-target bundle for these 3 planets — they are
sole-planet targets, so adding them changes the target list and would force
re-validation of an already-validated bundle for a 0.4% gain. Instead:

1. Run the current 592-target missing queue as planned.
2. When building the follow-up batch (the 320 flagged refits), include these 3
   targets with DR25 TCE seed depths (and re-derive `duration`/`epoch` seeds from
   the same TCE rows for internal consistency — do not mix catalogs within a row).
3. Keep K03722.01 and the remaining 23 blocked/manual.

Rule used (stricter than the earlier historical audit in one way — it requires a
*self-consistent* alternative fit, not a depth grafted onto current-catalog rows):
period match < 1%, duration agreement < 25%, unambiguous or mutually consistent
depth, from an independent pipeline fit (DR25 TCE).
