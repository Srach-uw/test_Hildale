# Repository Contents

This repository contains source code, compact metadata, reference material, and audit
documentation needed to inspect the replication. Large generated products are excluded.

## Included

- `scripts/`: analysis, extraction, diagnostics, tests, and release checks.
- `cloud/`: resumable ALDERAAN execution helpers and cloud validation scripts.
- `metadata/`: compact output tables that support documented counts and comparisons.
- `reference/`: the published article, machine-readable tables, and source references.
- `docs/`: current scientific status, audit records, and operational runbooks.
- `legacy/`: superseded work retained for provenance and excluded from canonical results.

## Excluded

- Raw Kepler light curves.
- ALDERAAN result FITS files and posterior sample archives.
- Virtual environments and downloaded dependency repositories.
- Cloud result archives and generated run directories.
- Credentials, private keys, access tokens, and billing information.

The exclusions are enforced by `.gitignore`. The tracked release surface is additionally
checked by `scripts/check_professor_release.py`.
