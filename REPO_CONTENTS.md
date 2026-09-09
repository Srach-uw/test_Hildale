# Repository contents

This repository contains the code and compact evidence needed to inspect the
Sagear et al. (2026) replication.

## Included

- `scripts/`: analysis code, diagnostics, tests, and release checks.
- `cloud/`: ALDERAAN runner, recovery, packaging, and collection utilities.
- `metadata/`: compact derived tables with provenance.
- `data/alderaan_factorial_validation_20260715/`: the curated 82-fit validation
  release stored through Git LFS.
- `reference/`: the article, machine-readable tables, and cited source material.
- `docs/`: current scientific status, focused method audits, and runbooks.
- `legacy/`: superseded early analysis, clearly separated from canonical work.
- `CITATION.cff`: citation metadata for the repository and source article.
- `LICENSE`: MIT license for original repository software.
- `CONTRIBUTING.md`: scientific and repository contribution rules.
- `docs/reproducibility_map.md`: claim-to-code and evidence index.
- `metadata/final_public_boundary_20260813/`: final compact boundary audits.

## Excluded

- Raw Kepler light curves.
- Full eccentricity-posterior archives.
- Cloud recovery archives and generated run directories.
- Virtual environments, caches, and temporary files.
- Credentials, billing data, and personal filesystem paths.
- Private chronological research logs.

The exclusions are enforced by `.gitignore` and
`scripts/check_professor_release.py`.

Third-party data and article files remain subject to their source licenses.
See `docs/data_availability.md` for details.
