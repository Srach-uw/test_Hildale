# Replication match ledger

Updated: 2026-09-20

## Reading this ledger

This ledger separates four states that are often blurred in a reproduction
project: matched, tested but unresolved, invalidated, and unavailable from the
public record. A matching total or a passing software test does not by itself
show that every input to the published analysis has been reproduced.

## Published analysis versus current reconstruction

| Stage | Current evidence | Status | Consequence |
| --- | --- | --- | --- |
| Published kinematic host labels | The machine-readable host table supplies 1,888 labels, 1,515 thin and 373 thick. The reconstruction uses those labels in the primary comparison. | Matched | Disk assignment is not the current primary discrepancy. |
| Planet total | The reconstructed inventory contains 2,465 planets, the reported paper total. | Matched at total level | The final accepted planet IDs are not published, so equal totals do not establish equal membership. |
| Single versus multi label | Multiplicity is frozen from the full eligible system before planet-level fit cuts. | Matched method | A planet-level fit failure cannot turn a multi into a single. |
| ALDERAAN transit-shape samples | The extractor retains paired duration, radius-ratio, impact, and period samples and uses Dynesty `LN_WT` weights. | Matched implementation | Equal weighting of nested rows is excluded from primary results. |
| Eccentricity conversion | The current branch uses the finite-duration importance calculation with paired samples. | Matched implementation | It does not identify whether the historical per-planet stellar-density inputs match. |
| Population hierarchy | The forward-normalized hierarchy passes algebraic and synthetic checks. The paper describes a NumPyro implementation; the public evidence does not expose its exact input ordering and configuration. | Tested, not historically identical | The public result is a reconstruction diagnostic, not an exact software replay. |
| Limb darkening, selected priors, sampler seed | The 82-fit matrix contains valid comparisons for these inputs. Their panel-level effects are too small to explain the reported discrepancy. | Tested within a targeted panel | These choices remain part of provenance, but are not leading explanations on their own. |
| Cadence | The 18 archived LC+SC arms omitted `--use_sc True` and processed no short-cadence data. | Invalidated test | Cadence remains open. A new test must actually pass the short-cadence argument. |
| Public ALDERAAN injection branch | The `Kepler-Validation` branch imports a missing `alderaan.validate` module. | Blocked by missing public code | Do not claim a native ALDERAAN injection recovery test until this module is recovered or a separate diagnostic implementation passes tests. |
| Stellar-density prior | The paper cites both Berger 2018 and Berger 2020 in different places. The public 2018 table does not contain the complete density uncertainty representation. | Not uniquely recoverable | Catalog substitutions are sensitivity tests, not a reconstruction of the original density ledger. |
| Final fit decisions | The final accepted/rejected KOI list and rejection reasons are not identified by the article tables. | Not publicly available | Population membership cannot be proven identical from public products. |

## Current numerical comparison

The current public QC branch is recorded in
[`metadata/public_reconstruction_20260727/population_comparison.csv`](../metadata/public_reconstruction_20260727/population_comparison.csv).

| Population | QC planets | Reconstructed Rayleigh mean, 16th-84th | Paper planets | Paper Rayleigh mean, 16th-84th |
| --- | ---: | ---: | ---: | ---: |
| Thin singles | 1,104 | 0.232, 0.225-0.239 | 1,121 | 0.022, 0.017-0.029 |
| Thick singles | 268 | 0.208, 0.194-0.223 | 275 | 0.066, 0.045-0.096 |
| Thin multis | 873 | 0.094, 0.088-0.103 | 862 | 0.030, 0.023-0.031 |
| Thick multis | 209 | 0.133, 0.121-0.147 | 207 | 0.033, 0.015-0.065 |

The result differs from Table 3 in every cell. It should not be interpreted as
a new astrophysical measurement.

## Work in progress

A predeclared six-system pilot now has verified ALDERAAN timing tables and
101 verified long-cadence public MAST files. It can answer a narrower
question: whether a controlled circular signal is recovered under observed
Kepler sampling. The input inventory and the injection-code limitation are
recorded in [pilot_injection_preflight.md](pilot_injection_preflight.md).

The pilot cannot determine the population result by itself. It can only decide
whether the transit-fitting and noise path deserves further investigation
before attributing the discrepancy to stellar densities, selections, or
unreleased author intermediates.

## What would close the remaining historical gap

The smallest author-provided products with the highest diagnostic value are:

1. the final accepted planet ledger with KOI IDs and rejection reasons; and
2. the stellar-density rows or posterior samples used in the eccentricity
   calculation.

The final per-planet `(e, omega)` samples and the exact population-model input
order would allow an even stronger reproduction, but the two products above
would identify most remaining upstream differences.
