# Real-Data Audit Findings — 2026-07-04

First audit against REAL ALDERAAN output: a snapshot of 48 completed `*-results.fits` from the
live 592-target GCP run (`sagear_missing` batch). Complements the code-only
`deep_audit_findings_2026-07-04.md` — everything here is measured from actual fitted posteriors.
No script fixes applied; the live run was not touched.

Evidence files (committed alongside this doc, `docs/evidence_real48/`):
`summary_48.csv` (67 extracted posterior summaries), `fits_crosscheck_table.csv`
(per-planet fitted-vs-catalog cross-check), `coverage_48.csv`, `merged_coverage_48.csv`.

---

## 1. Extraction and merge results

| Step | Result |
|---|---|
| Results files found | 48/48 (flattened copies of `Results/sagear_missing/<T>/<T>-results.fits`) |
| Sample planets in the 48 systems | 68 |
| Posteriors extracted | **67/68** (`extract_eccentricity_posteriors.py`, default geometric impact mode) |
| Not extracted | **K01316.02** only — see §3; excluded *upstream*, not an extraction failure |
| Merge vs 1729-row archive | 1729 + 67 = **1796 unique planets, zero collisions** (all 67 are net-new) |
| Coverage after merge | thick multi 197/233 (84.5%), thick single 179/272 (65.8%), thin multi 664/860 (77.2%), thin single 756/1109 (68.2%) |

All 48 systems are thick disk (37 multi + 30 single planets), consistent with the missing-batch
composition. All 67 `.npz` e–omega grids pass integrity checks: finite, non-negative, e_pdf
normalized to 1 within 1e-6, and **zero significant bimodality** (no posterior has a secondary
mode >50% of the primary with a real dip; the raw >1-local-max count of 25 is smoothing-kernel
ripple only).

Note: the prompt's `outputs/existing_eccentricity_posteriors_best.csv` does not exist; the merge
used the script's default base `eccentricity_posterior_summary_old_astropy_rawcc.csv` (1729 rows),
which matches the expected archive count.

## 2. Sanity checks on the 67 extracted posteriors

| Check | Result | Verdict |
|---|---|---|
| Period match (ALDERAAN ttimes-derived vs KOI) | median rel diff 1.1e-4, min 3.8e-7, max 1.09e-3 (K02148.03: 3.6146 d vs 3.6107 d) | OK — see caveat below |
| e ordering (e16 ≤ e50 ≤ e84) | 0 violations | clean |
| Grid-boundary pile-up | 0 planets with e50 < 0.01 or > 0.90 | clean |
| Posterior width e84−e16 | 0.27–0.52 (median 0.47); none suspiciously narrow (<0.02) | expected for single-planet photoeccentric |
| e50 range | 0.160–0.837, median 0.254 | physical |
| zeta_median | 0.32–1.38, all inside [0.3, 2.5] | sane |
| n_zeta | 3704–4000 of 4000 draws; worst is K01549.01 (296 draws lost to tcirc=0 when geometric b-draw exceeds the transit chord — the source of the benign `divide by zero` RuntimeWarning) | benign |
| Fitted DUR14 vs koi_duration | ratio 1.00 ± 0.11, range [0.67, 1.36], none outside [0.5, 2] | clean |
| Fitted depth (ROR²) vs koi_depth | median ratio 0.88 (expected: ROR² understates limb-darkened depth by ~10–15%); **1 outlier: K01549.01** (§3) | clean except grazer |
| Sibling mix-up test (15 multi systems with ≥2 extracted planets) | ALDERAAN planet index is monotonic with period in all 15 systems; 1:1 matching; no duplicated/copied posteriors between siblings | clean |
| Distribution vs archive | new-67 e50 deciles {0.18, 0.20, 0.25, 0.34, 0.58} vs archive-thick {0.17, 0.19, 0.25, 0.37, 0.58} | **nearly identical — strong consistency** |
| Paper direction (thick single > thick multi) | e50 mean: single 0.389 vs multi 0.260 | consistent with paper |

Period caveat (honest, not alarming): the earlier single-target check found ppm-level agreement;
across 67 planets the median is ~112 ppm and 1 planet exceeds 1e-3. The ALDERAAN-side period is
derived from the median of noisy fitted transit-time differences (`period_from_ttimes`), so
1e-4–1e-3 scatter for small/TTV planets is expected. It is used **only for matching**
(tol 0.01); `koi_period` is what enters the circular-duration calculation. No science impact.

Spearman(zeta_median, e50) = −0.78 across the 67 — the expected strong (but not deterministic)
duration-ratio→eccentricity mapping; no sign of a pathological one-to-one lookup.

## 3. Grazing-transit cross-reference (the "human decision" evidence)

Of the prior audit's flagged sets, the 48 targets contain: **1 catalog b>1 planet, 1 NaN-impact
planet, and 11 planets with 0.9 < b ≤ 1.**

**K01316.02 (NaN impact, NaN depth) — never fit.** `alderaan_needed_catalog_best.csv` contains
K01316 once with `npl=1`; the FITS header confirms `NPL=1`. The planet failed the
`koi_depth > 0` catalog filter (depth is NaN), so ALDERAAN modeled K01316 as a single-planet
system. Consequence: (a) the NaN-impact worry is moot for this planet — it self-excludes; (b) the
sibling K01316.01 was fit with the 12.4-d sibling's transits unmodeled in the light curve. Its
posterior looks normal (e50 = 0.254) but it has one of the larger period diffs (6.5e-4). Low
severity; worth a flag column ("fit with incomplete system model") rather than action.

**K02775.01 (catalog b = 1.25) — the one true b>1 planet here.** Extracted e50 = 0.770
(e16–e84 = 0.53–0.87), zeta_median = 0.396. Not bimodal, not boundary-piled, but pushed hard to
high e: exactly the predicted failure mode where a grazing-shortened duration is read as
eccentricity when b is marginalized geometrically.

**K01549.01 (catalog b = 0.983) — objectively pathological fit.** ALDERAAN's fitted median
ROR = 0.493 vs catalog 0.169 (implied depth 24.3% vs measured 1.19%, a 20× discrepancy) with
fitted b = 1.37 — the classic runaway ror–b grazing degeneracy. It is the most extreme
eccentricity of all 67 (e50 = 0.837, the only e84 > 0.90) and lost the most zeta draws (3704).
This one planet is *demonstrably* not a credible eccentricity measurement.

**Systematic b–e correlation across the 67 (catalog b bins):**

| b bin | n | e50 mean | e50 median | zeta median |
|---|---|---|---|---|
| <0.5 | 36 | 0.259 | 0.245 | 1.04 |
| 0.5–0.8 | 14 | 0.294 | 0.258 | 0.98 |
| 0.8–0.9 | 5 | 0.270 | 0.254 | 1.21 |
| **>0.9** | **12** | **0.544** | **0.604** | **0.57** |

The jump is confined to b > 0.9 and tracks low zeta, i.e. duration deficit attributed to
eccentricity. Excluding b > 0.9 lowers thick-single mean e50 from 0.389 to 0.309 (multi: 0.260 →
0.242). Because the paper's headline result is *thick singles are more eccentric*, and 6 of the
12 high-b planets here are thick singles, a thin/thick asymmetry in grazing fraction could
partially mimic the signal. **Evidence-based recommendation for the human decision:** carry a
`b_flag` (catalog b > 0.9, plus fitted-ror/catalog-ror > 2 as a runaway-degeneracy tripwire —
which catches K01549.01) into the population fit and run the hierarchical fit with and without
flagged planets as a sensitivity test. K02775.01 and K01549.01 are the strongest candidates for
outright exclusion; do not silently drop the other ten.

## 4. Bugs / anomalies found

| # | Item | Severity | Proposed fix (not applied) |
|---|---|---|---|
| R1 | K01549.01 runaway grazing fit (ROR 0.49 vs 0.17, fitted b 1.37) propagates to the top of the e distribution | Medium (per-planet; population fix via flag) | Add flag: fitted-ror/catalog-ror > 2 or fitted b > 1 → exclude/flag in population fit |
| R2 | b>0.9 planets systematically read as high-e under geometric impact marginalization (12/67 here, e50 mean 0.54 vs 0.26) | Medium (fidelity; affects headline thin-vs-thick singles contrast) | Sensitivity rerun of hierarchical fit excluding b>0.9; deep-audit C4 decision now has quantitative backing |
| R3 | Multi systems with an upstream-excluded sibling (NaN depth) are fit with incomplete system models (K01316) | Low | Add "incomplete system model" flag column when catalog npl < sample npl for a target |
| R4 | `divide by zero` RuntimeWarning in `extract_eccentricity_posteriors.py:178` (tcirc = 0 when geometric b draw exceeds chord) | Cosmetic (draws are filtered) | Optional: `np.errstate` guard around the division |
| R5 | Prompt/docs reference `outputs/existing_eccentricity_posteriors_best.csv`, which does not exist | Cosmetic | Use the script default `eccentricity_posterior_summary_old_astropy_rawcc.csv` |

No evidence of: unit errors, sibling transit mix-ups, boundary pile-ups, bimodal pathologies,
duplicated posteriors, disk/system-correlated artifacts, or extraction/matching bugs.

## 5. Verdict

**Confidence in the pipeline is increased.** 67/68 planets extracted cleanly with the single
miss fully explained upstream; fitted durations agree with the KOI catalog at the 11% level with
zero gross outliers; depths agree at the level expected from limb darkening; sibling matching is
provably correct in all 15 multi systems; and the new posteriors' e50 distribution is
statistically indistinguishable from the 1729-planet archive built from independent earlier fits
— the strongest end-to-end replication signal available at this stage. The two genuine issues
(K01549.01's runaway grazing fit; the b>0.9 high-e systematic) are properties of grazing
geometry, were predicted by the code-only audit (C4), and are now quantified with real numbers —
they need a flagged sensitivity test in the population fit, not a pipeline change. Nothing found
here blocks continued ingestion of results from the live run.
