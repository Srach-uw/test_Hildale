# Pilot injection preflight

Updated: 2026-09-20

## Purpose

This note records the preparation state for a six-system photometric
injection and recovery pilot. The pilot is intended to test whether the
public transit-fitting path can recover known circular transits under realistic
Kepler sampling. It is not a population rerun and does not update any canonical
eccentricity product.

## Frozen pilot systems

| KOI system | KIC | ALDERAAN result source | Long-cadence files | Short-cadence files listed by MAST |
| --- | ---: | --- | ---: | ---: |
| K00367 | 4,815,520 | original archive | 14 | 13 |
| K01868 | 6,773,862 | original archive | 17 | 0 |
| K01474 | 12,365,184 | recovery archive | 18 | 25 |
| K02204 | 8,494,542 | recovery archive | 17 | 0 |
| K00366 | 3,545,478 | published-inventory recovery | 18 | 16 |
| K01852 | 9,763,348 | published-inventory recovery | 18 | 0 |

The six result FITS passed the timing-readiness validator. Each has a matching
`TARGET` header, complete `TTIMES_nn` coverage, finite required timing columns,
and a recorded SHA-256. The downloaded long-cadence files are public MAST
PDCSAP FITS. All 101 files passed KIC identity, observing-mode, finite-PDCSAP,
and SHA-256 checks. These local audit artifacts are intentionally excluded from
the public repository because they contain downloaded inputs and machine-local
paths.

## Public-source support

Sagear's public `photoeccentric` tutorial constructs injected circular transits
on Kepler sampling, fits a circular transit model, and then compares circular
transit density with stellar density to infer eccentricity. It also explicitly
distinguishes Gaussian simulated noise from real Kepler light curves. This
supports the pilot's purpose, but the tutorial is not evidence that it is the
2026 disk-paper implementation.

## Implementation boundary

The public `ssagear/alderaan` release at commit
`7443dff16b7f9092e14a6f0cc1f8948d457c9e0b` contains a
`Kepler-Validation` branch in `detrend_and_estimate_ttvs.py`. That branch calls
`remove_known_transits` and `inject_synthetic_transits` from
`alderaan.validate`. The module is absent from the release, from the packaged
0.1.0 source examined for this project, and from the available Git history for
that path. The public normal-run repair makes the missing import optional only
because normal Kepler runs never enter the validation branch.

Consequently, the validation branch must not be used as though it were a
verified injection engine. Before fitting injected data, the project needs one
of the following:

1. a source-backed copy of the missing validation module, or
2. a separately implemented injection path with unit tests against known
   noiseless transits and an explicit statement that it is a new diagnostic
   implementation.

Either path must preserve the real-transit masking rule, injection parameters,
exposure integration, random seed, selected cadence files, and every rejected
trial. Only then can the pilot distinguish a transit-fitting or noise-treatment
mechanism from an upstream stellar-density or membership difference.

## Next gate

Do not start a broad cloud run. First implement the new diagnostic path in a
separate workspace, run the noiseless recovery test at `e=0` for declared
impact parameters, and retain the generated input and output manifest. If the
noiseless test fails, stop and repair geometry or exposure integration before
adding observed Kepler noise.
