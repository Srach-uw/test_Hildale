# Diagnostic transit injection protocol

This protocol defines a small, standalone test of photometry-to-density
recovery. It is a new diagnostic implementation. It is not a reconstruction of
the unavailable `alderaan.validate` module, and it does not replace, modify, or
reinterpret any ALDERAAN posterior.

## Purpose

The current replication disagreement could originate in the transit fit, the
conversion from transit shape to eccentricity, or the population model. A
circular injection test isolates the first of those stages. It asks whether a
known circular transit, placed into a declared light-curve baseline, returns a
circular density under a specified fitting workflow. It cannot by itself prove
that the published population analysis was reproduced.

## Input contract

`scripts/diagnostic_transit_injection.py` accepts an NPZ file containing:

| Field | Requirement | Reason |
| --- | --- | --- |
| `time` | Finite, one-dimensional float array | Defines the fixed observing times. |
| `flux` | Finite, one-dimensional float array matching `time` | Supplies the untouched baseline. |
| `mask_index` | Non-empty, strictly increasing integer indices | Defines exactly where the synthetic transit may be applied. |

The tool rejects changed array lengths, unordered or non-finite points,
duplicate mask rows, and indices outside the input. It modifies flux only at
`mask_index`; time and the selected index array are preserved in the output and
checked in the tests. In a real pilot, `mask_index` must contain the
injection-eligible cadences after conservative masks have removed every known
transit. The manifest records the complete input hash, separate time, flux,
and selected-row hashes, plus the rendered output hash. The companion
preparation command, `prepare_diagnostic_injection_input.py`, creates this
input from local PDCSAP and ALDERAAN files only after it checks the complete
system catalog, fitted timing tables, Kepler time convention, quality-mask
version, cadence mode, and exposure time. Its output also retains errors,
cadence identifiers, quarters, source row mapping, and `known_transit_mask`.
The recovery stage must exclude `known_transit_mask`; avoiding injection in
those cadences alone would leave the original planets in the fit.

## Circular geometry

The injection specifies period `P`, stellar density `rho_star` in solar units,
radius ratio `Rp/Rstar`, impact parameter `b`, transit epoch, quadratic limb
darkening, cadence exposure, oversampling factor, and seed. It derives
`a/Rstar` from the circular density relation

```
a/Rstar = [G rho_sun rho_star P^2 / (3 pi)]^(1/3).
```

The reciprocal relation is recorded in the manifest as a round-trip check.
The implementation validates a transiting geometry (`0 <= b <= 1 + Rp/Rstar`),
a radius ratio below one, physically bounded quadratic limb darkening, and an
inclination that is defined for the derived `a/Rstar`.

## Rendering and provenance

Rendering requires `batman-package`. When available, the utility uses
`batman.TransitModel` with the stated `exp_time` and `supersample_factor`, so
each synthetic flux point is exposure integrated. It does not add random noise;
the declared seed is still recorded to prevent ambiguity when this controlled
tool is paired with an external noise-realization procedure.

`--dry-run` validates the specification and input contract without importing
`batman`. It writes a deterministic JSON preflight manifest. A rendered run
writes the injected NPZ and a manifest. The output directory must be absent or
empty, and an existing manifest is never replaced. This makes accidental
overwriting visible rather than silently changing provenance.

`prepare_diagnostic_injection_background.py` provides the declared
noise-realization step without changing the cadence contract. It creates a
separate immutable NPZ and manifest for one of three baselines: observed raw
PDCSAP flux, independent draws from `PDCSAP_FLUX_ERR`, or independent draws
with a robust per-source-file PDCSAP residual scale. Synthetic baselines are
centered on the median of injection-eligible cadences in each source file, so
known-transit cadences cannot set the baseline or noise scale. The fixture
builder still masks real transits before ALDERAAN reads the synthetic files.

Example preflight:

```powershell
python scripts/diagnostic_transit_injection.py `
  --input-npz pilot_input.npz `
  --output-dir pilot_preflight `
  --period-days 10.0 `
  --rho-star-solar 1.0 `
  --radius-ratio 0.02 `
  --impact-parameter 0.2 `
  --t0-days 134.5 `
  --limb-darkening 0.30 0.20 `
  --exposure-minutes 29.4 `
  --supersample-factor 15 `
  --seed 20260920 `
  --dry-run
```

## Bounded experiment

The first scientific use should remain narrow:

1. Use K00367 and K01868 as an original-archive high-leverage/comparator pair.
   K01868 is not assumed to be circular.
2. Run noiseless circular cases at `b = 0.2` and `b = 0.8` for each target.
   Verify independent exposure-integrated flux agreement within 1 ppm and
   fixed-geometry density recovery within 1 percent before adding noise.
3. Then use at most 24 noisy recoveries. The observed-PDCSAP background is
   deterministic and needs one realization. Each synthetic background gets two
   fixed realizations. For two targets and two impact parameters, that is a
   20-fit matrix: two targets times two geometries times five distinct
   backgrounds. The backgrounds are observed PDCSAP, quoted-error Gaussian
   noise, and variance-matched independent Gaussian noise. A time-ordered
   PDCSAP-residual background is optional only after its detrending and
   residual-construction rule is specified in advance. Record failed trials
   and selection losses.
4. Record `log10(rho_circ/rho_injected)`, interval coverage, input hashes,
   masks, timing treatment, and software versions. Do not tune the experiment
   against Sagear's Table 3 values.

The immediate warning threshold is a median signed density bias above 0.10 dex
within an eight-case background arm, or three of eight nominal 95 percent
intervals missing the injected truth. Those are diagnostic flags, not
population-level significance tests. A passing circular pilot establishes only
that the declared pipeline can recover the limited injected cases. It does not
resolve cadence, detrending, eccentricity-prior, stellar-density, or accepted
sample differences.

The full photometry test injects into raw, pre-detrending PDCSAP flux, then
uses a declared ALDERAAN recovery workflow. A post-detrending injection is a
separate control for isolating the detrending response. Neither should be
presented as an ALDERAAN validation unless the recovery holds the transit
geometry sufficiently free to infer circular density, reports coverage, and
keeps the injected stellar density out of the transit-shape likelihood.
## ALDERAAN fixture gate

For a single-planet system, `build_diagnostic_alderaan_fixture.py` converts a
rendered injection into a separate set of MAST-format PDCSAP files and a new
catalog row. It checks the source-file hashes and row mapping, masks the known
real transits with `NaN`, and updates the catalog period, epoch, depth,
duration, impact parameter, and limb-darkening values to the injected circular
specification. The resulting files can exercise the public ALDERAAN detrending
path without modifying a source file or any canonical result. A successful
fixture build is still only an input gate. The recovery fit, output validation,
and density comparison remain separate gates.

The single-target recovery runner is
[`../cloud/diagnostic_injection/run_one_fixture.sh`](../cloud/diagnostic_injection/run_one_fixture.sh).
It rejects reused project directories, copies the fixture rather than a source
file, and never calls a downloader. Its result belongs in a separate
diagnostic directory until the recovered density has passed the planned gates.
