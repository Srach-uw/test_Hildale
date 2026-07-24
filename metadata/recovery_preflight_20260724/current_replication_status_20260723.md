# Current Sagear Replication Status

## Bottom line

The large discrepancy is no longer a generic failure of the whole pipeline.
Two specific upstream differences remain capable of producing it:

1. The exact published-host reconstruction revealed 142 systems, containing
   174 population planets, for which no usable ALDERAAN result is available
   locally.
2. The published values are reproduced remarkably well if raw dynesty rows
   are treated equally and the thin/thick population arrays are inverted after
   the published host table. This is a historical-code diagnostic, not a
   scientifically valid relabeling.

The missing systems must be fitted before the inversion hypothesis can be
accepted or rejected. Real observed posterior shapes, used as extreme analogs
for those missing systems, can move every incomplete population through
Sagear's central value.

## Exact sample reconstruction

Joining Sagear's published 1,888-host table to non-false-positive KOIs with
periods from 1 to 100 days, while preserving the pre-cut `koi_count`, gives:

| population | reconstructed before final visual QC | Sagear | difference |
|---|---:|---:|---:|
| thick singles | 275 | 275 | 0 |
| thin singles | 1121 | 1121 | 0 |
| thick multis | 212 | 207 | 5 |
| thin multis | 883 | 862 | 21 |
| total | 2491 | 2465 | 26 |

The professor's multiplicity warning is a valid pipeline safeguard: labels
must be assigned from the full host inventory before planet-level cuts. The
direct audit shows that Sagear's published pre-visual counts are reproduced by
the catalog `koi_count` definition, while recounting only surviving
non-false-positive rows changes 23 hosts and moves away from the published
contract. The remaining 26 unpublished final removals are count-constrained
to 5 thick-multi and 21 thin-multi planets, subject to catalog-epoch and
visual-QC differences.

## Current posterior coverage

The archive and completed cloud run provide 2,317 unique posterior products on
the exact reconstructed inventory. After primary posterior QC, the
Sagear-count coverage used by the equal-row fit is:

| population | usable posteriors | Sagear planets | missing |
|---|---:|---:|---:|
| thick singles | 256 | 275 | 19 |
| thin singles | 1028 | 1121 | 93 |
| thick multis | 204 | 207 | 3 |
| thin multis | 828 | 862 | 34 |

The 142 missing systems divide into 116 systems that were never included in
the previous 592-target cloud launch and 26 systems that were launched but did
not produce a usable result.

The missing set is demonstrably non-random. Its target-level median Berger
`logg` is 4.15 versus 4.43 for covered systems. Forty-nine of the 142 missing
systems have a nonzero Berger evolutionary flag, compared with 275 of 1,746
covered systems. The missing systems also have larger median stellar radii and
lower median density metadata in every population. This is consistent with
the earlier strict dwarf/evolutionary-state filtering and makes the Berger
catalog and stellar-state cuts a primary unresolved cause, especially for thin
singles.

## Weighting and label diagnostics

### Scientifically correct dynesty weights

| population | ours, posterior median and 16th-84th | Sagear |
|---|---:|---:|
| thick singles | 0.195 [0.180, 0.211] | 0.066 [0.045, 0.096] |
| thin singles | 0.225 [0.218, 0.232] | 0.022 [0.017, 0.029] |
| thick multis | 0.132 [0.120, 0.145] | 0.033 [0.015, 0.065] |
| thin multis | 0.081 [0.075, 0.087] | 0.030 [0.023, 0.031] |

### Equal raw rows with published labels

| population | ours, posterior median and 16th-84th | Sagear |
|---|---:|---:|
| thick singles | 0.034 [0.022, 0.051] | 0.066 [0.045, 0.096] |
| thin singles | 0.064 [0.049, 0.080] | 0.022 [0.017, 0.029] |
| thick multis | 0.030 [0.021, 0.041] | 0.033 [0.015, 0.065] |
| thin multis | 0.021 [0.016, 0.027] | 0.030 [0.023, 0.031] |

### Equal raw rows with downstream disk-array inversion

| paper-labeled population | diagnostic value | Sagear | interval overlap |
|---|---:|---:|---|
| thick singles | 0.064 [0.049, 0.080] | 0.066 [0.045, 0.096] | yes |
| thin singles | 0.034 [0.022, 0.051] | 0.022 [0.017, 0.029] | yes |
| thick multis | 0.021 [0.016, 0.027] | 0.033 [0.015, 0.065] | yes |
| thin multis | 0.030 [0.021, 0.041] | 0.030 [0.023, 0.031] | yes |

The equal-row inverted diagnostic has four of four interval overlaps,
RMSE 0.0086, and approximate chi-square 1.00. Correct dynesty weighting misses
all four intervals whether labels are inverted or not.

This does not justify changing the published disk labels. Table 1 explicitly
defines `P_thick` and disk assignment, and the sample counts also support those
labels. If an inversion occurred, it would have to be downstream in the
posterior-array or population-label mapping.

## Why missing coverage still matters

The hierarchical likelihood is a product over planet posteriors, so a modest
number of strongly informative posteriors can move the population scale much
more than their raw percentage suggests. An empirical stress test replaced
all missing rows in each bin with copies of each real observed posterior from
that same bin. The resulting possible ranges were:

| population | empirical analog-imputed range | Sagear |
|---|---:|---:|
| thick singles | 0.015-0.128 | 0.066 |
| thin singles | 0.003-0.212 | 0.022 |
| thick multis | 0.025-0.035 | 0.033 |
| thin multis | 0.004-0.073 | 0.030 |

This is not an estimate of the missing planets. It proves that observed,
non-synthetic posterior shapes are capable of explaining the residual. The
142 systems therefore must be completed before attributing the match to a
label-ordering bug.

## What has been ruled out as the main explanation

- the professor-identified post-cut multiplicity bug, now fixed;
- ordinary sample-count mismatch;
- cadence, limb darkening, printed transit priors, and sampler randomness in
  the complete 82-fit factorial validation;
- paired versus synthetic impact after the extractor correction;
- tested stellar-density error scales and offsets;
- tested transit-selection conventions;
- a small set of isolated high-e outliers.

The missing-coverage selection-bias table is in
`outputs/missing_coverage_selection_bias.csv`; it should be read before any
population interpretation.

## Next action

A validated, resumable completion bundle is ready:

`HILDALE ALDERAN/Hildale_ALDERAAN_Published_Inventory_Missing_142_20260724_v4.zip`

It contains 142 targets and a 181-row full-system ALDERAAN catalog, including
seedable companions outside the 1-100 day population cut. It uses the pinned
Sagear ALDERAAN commit, deterministic target seeds, target-specific logs,
resume-safe result detection, long plus available short cadence, and the
repaired runner.

In parallel, send the author clarification request for the exact 26 final
rejections, dynesty-weight handling, stellar-density source, and the mapping
from Table 1 disk labels to the four population arrays.

## Recovery implementation status

The authoritative bundle is the v4 archive above. Its SHA-256 digest is
`064A02C4BF1AFBA8FC7521C9C7E517062FC22A339DEB66B8CAEBACFCF79B9B46`.
It validates as 142 targets and 181 full-system catalog rows.

The runner now provides deterministic seeds for detrending, noise analysis,
and dynesty; stage-specific failure statuses; resumable retries with a changed
seed offset; and a complete packaged audit trail. Packaging completes even if
some targets remain unresolved, so their statuses and terminal logs are not
lost. The local postprocessor
extracts weighted and equal-row sensitivity branches but withholds final
hierarchical results unless the exact 2,491-planet inventory has complete
posterior coverage.

The full local test suite passes 106 tests. No GitHub files were changed or
pushed.

## Multi-VM recovery preflight

The live project quotas allow 24 E2 CPUs and 24 C3 CPUs per region but no
preemptible CPUs. The fast candidate therefore uses seven standard
`c3-highcpu-22` VMs in seven distinct regions, with deterministic 20-21 target
shards. Their union is exactly 142 systems and no host is duplicated.

Observed individual ALDERAAN validation fits reached 4.9-5.7 hours, so a
guaranteed one- or two-hour completion is not possible merely by adding VMs.
The realistic estimate is 4-7 hours and approximately $25-$55 of welcome
credit, bounded by eight-hour automatic VM deletion.

The sharding, shard runner, partial archive merge, and duplicate/missing-target
guards are implemented locally. A read-only Claude preflight package is ready
at:

`HILDALE ALDERAN/Claude_Final_Preflight_Sagear_20260724_v3.zip`

SHA-256:
`C245E0D33E46E177D299888790FCDEE59ECB010A1C719FD7ACA7DCAE78D32B26`.

The full local suite now passes 113 tests. GitHub remains untouched.

## Independent launch verdict

Claude returned `GO AFTER ONE FIX`. All scientific and shard invariants it
checked passed. Its valid operational finding was that the draft multi-VM
design used termination action `DELETE`, risking loss of an unretrieved shard
disk at the eight-hour limit.

The final design now uses `STOP` only. Shard disks persist, each shard writes a
predictable archive plus SHA-256 sidecar, and no VM is deleted before local
download, hash verification, seven-archive completeness, and duplicate-safe
merge validation.

Final remote bundle:

`HILDALE ALDERAN/Hildale_ALDERAAN_Published_Inventory_Missing_142_20260724_v5_sharded.zip`

SHA-256:
`E17969A916FB3AA8544213A7F58F68C096892476A90BD2A29D4245AA621E0563`.

The complete local suite now passes 116 tests. GitHub remains untouched.
