# Direct chord-eccentricity diagnostic

## Question

The production reconstruction first fits ALDERAAN's transit-shape parameters
and then maps those samples into eccentricity space. This diagnostic asks a
narrower question: does the same public, detrended ALDERAAN chord likelihood
give a materially different result when eccentricity and argument of
periastron are sampled directly and duration is derived from a fixed stellar
density?

It is not a population fit. It cannot validate Sagear's unpublished inputs,
and it must not be merged into the canonical eccentricity archive.

## Fixed comparison

For K00367's seeded circular fixture, the two routes must use the same:

- pinned ALDERAAN installation;
- filtered light curve and quick-TTV table;
- native ALDERAAN limb-darkening penalty;
- exposure integration and finite chord renderer;
- radius-ratio, impact-parameter, and ephemeris priors;
- ALDERAAN's finite duration support from one short cadence through three
  times the catalog duration;
- declared stellar density; and
- deterministic random seed, recorded in the output manifest.

The only intentional change is replacing sampled `T14` with the finite chord
duration derived from `(rho_star, e, omega, r, b)`. The public implementation
uses ALDERAAN's `dynesty_helpers.lnlike`; it does not replace that likelihood
with an exact eccentric-orbit light-curve renderer.

## Preflight

In the pinned Linux environment, first produce a new manifest. This hashes the
catalog, filtered light curve, and quick-TTV table and fails if the fixture is
not a single-planet long-cadence target.

```bash
python scripts/direct_chord_eccentricity_diagnostic.py \
  --project-dir /path/to/project_exact_seeded \
  --run-id diagnostic_k00367_exact_seeded \
  --target K00367 \
  --rho-star-solar 1.0198685421655524 \
  --manifest /path/to/k00367_direct_chord_manifest.json
```

Check that the input hashes match the retained fixture record before executing.

## Bounded execution

The first run is intentionally small. It is a compatibility and posterior
geometry check, not a production-quality posterior. Increase `maxcall` only
after inspecting the first output, its evidence trajectory, and its effective
sample size.

```bash
python scripts/direct_chord_eccentricity_diagnostic.py \
  --project-dir /path/to/project_exact_seeded \
  --run-id diagnostic_k00367_exact_seeded \
  --target K00367 \
  --rho-star-solar 1.0198685421655524 \
  --manifest /path/to/k00367_direct_chord_manifest.json \
  --execute --seed 20260925 --maxcall 20000 \
  --output /path/to/k00367_direct_chord_samples.npz
```

## Interpretation gates

1. The preflight manifest must list exactly the expected fixture files.
2. The output hash, row count, and random seed must be recorded.
3. Compare direct samples to duration-reweighted samples only after applying
   `LN_WT` weights to the latter.
4. Treat any disagreement as a parameterization diagnostic until the test is
   repeated with an independent target such as K01868.
5. Do not use either diagnostic result in a population estimate.
