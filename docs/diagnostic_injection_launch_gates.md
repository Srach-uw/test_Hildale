# Circular-injection launch gates

Date: 2026-09-24

## Decision

Do not launch a cloud recovery fit until every gate below has evidence attached
to the diagnostic run directory. A circular injection can test the available
photometry and fitting path. It cannot, on its own, prove that the paper's
stellar densities or sample definition caused the population mismatch.

## Fixture integrity

1. The prepared-input manifest and rendered-injection manifest must identify
   the same NPZ by SHA-256.
2. The injection exposure time must agree with the PDCSAP exposure declaration.
3. Each source file must match the preparation manifest by filename and
   SHA-256. Every retained source row must map to exactly one FITS row.
4. The target and KIC in the preparation manifest must agree with the fixture
   catalog. The first run is restricted to a one-planet system.
5. The fixture must preserve raw files outside retained rows, mask all known
   real-transit rows, and place the injected flux only in eligible rows.

## What ALDERAAN reads

Before sampling, load the fixture with the pinned `alderaan.io.read_mast_files`
path and compare its retained cadence identifiers, flux, errors, and quarters
with the fixture manifest. Confirm that the original transits are absent after
the normal read and that injected ingress and egress cadences remain.

ALDERAAN also imports the Holczer timing catalog during detrending. For a
phase-shifted circular injection, the fixture must contain a private copy with
the synthetic target's historical Holczer rows removed. Otherwise the matching
code may replace the catalog epoch and period with the real system's timing
solution before it searches for the injected transits. Record the source and
derived timing-file hashes, the target-row count removed, and the count of
unrelated rows retained. This is a fixture-only change; it must never alter the
production timing catalog or a canonical result.

## Recovery design

The first run uses a circular injection with free recovery geometry. Its catalog
values define a detection window and a starting configuration. They must not
provide the injected stellar density to the recovery likelihood. Run each
ALDERAAN stage in a fresh Python process with the same explicit NumPy seed and
pass the same value through `ALDERAAN_SEED`, which the pinned nested-sampling
driver reads when it constructs Dynesty's generator. Record both values with
the integration setting, transit count, cadence coverage, and every terminal
failure.

The initial exact-catalog run is a plumbing control only. A second, predeclared
run must perturb the catalog duration, impact parameter, and limb-darkening
starting values within realistic bounds. If the conclusion depends on knowing
the truth in advance, the control has not tested the inference.

## Interpretation matrix

| Result | Most direct implication | What remains unresolved |
| --- | --- | --- |
| Noiseless geometry control fails | Renderer, geometry conversion, or recovery implementation needs repair | Detrending and stellar inputs |
| Noiseless and post-detrending controls pass, full PDCSAP recovery fails | Detrending, timing, or noise treatment is implicated | Historical stellar inputs |
| All controlled recoveries pass | Tested public photometry path can recover a circular density | Whether it matches the paper's per-target inputs and historical choices |
| Controlled recoveries are unstable | The available configuration cannot support a population-level diagnosis | Both upstream inputs and population inference |

## Other high-value checks

The corrected short-cadence experiment remains open because the earlier LC+SC
arms did not enable short cadence. The historical paper-prior sensitivity needs
execution-time provenance, not only a patch file. Direct eccentric-transit
inference and the importance-sampling branch should also be compared under the
same priors and selection treatment. These are separate from the circular
fixture and must not be conflated with it.
