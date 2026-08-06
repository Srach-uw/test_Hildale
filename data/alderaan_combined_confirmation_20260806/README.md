# Combined ALDERAAN Confirmation

This release contains the nine-system ALDERAAN arm that combines the closest
available published-method settings:

- Sagear-reference limb-darkening inputs;
- long cadence plus available short cadence;
- the printed-prior sensitivity patch;
- the pinned ALDERAAN commit
  `7443dff16b7f9092e14a6f0cc1f8948d457c9e0b`.

All nine target jobs completed with exit code 0. The result FITS are stored with
Git LFS.

## Scope

The arm is compared only with the matched `reference_lcsc` arm in the original
82-fit validation release. This estimates the effect of the printed-prior patch
conditional on reference limb darkening and LC+SC. It is not a full factorial
difference-in-differences interaction because the long-cadence printed-prior
arm uses a different eight-system target set.

## Contents

- `results/`: nine nonempty ALDERAAN result FITS.
- `status/` and `status_manifest.csv`: terminal per-target accounting.
- `input_catalogs/`: the reference catalog and full-system inventory.
- `run_spec/`: the exact runner scripts.
- `provenance/`: run metadata, Conda package lock, the complete ALDERAAN patch,
  and copies of all modified source files.
- `FITS_SHA256SUMS.txt`: SHA-256 for the nine published result FITS.

The source archive was independently verified after download:

```text
archive: alderaan_combined_confirmation_20260806_PROVENANCE_COMPLETE_20260806T204634Z.tar.gz
sha256: dd9ea39b9715ed8254eea309ebd2d5aa54a5db078ed9900962f59fc28e027655
bytes: 27749588
```

The full archive, including batch and per-target logs, remains in the local
research record. Verbose VM logs are omitted here because the terminal status,
run specification, source patch, environment, and result products provide the
reproducible public record without machine-specific noise.

Compact regenerated comparisons are in
`metadata/combined_confirmation_20260806/`.
