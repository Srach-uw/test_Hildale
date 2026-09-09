# Data availability

This repository separates public reference material, curated validation data,
compact derived evidence, and large local products.

## Included data

The `reference/` directory contains the published article, its machine-readable
tables, and archived source material used to define the replication target.
The `data/` directory contains the completed ALDERAAN factorial and combined
configuration validation results. FITS files in these releases use Git LFS.
The `metadata/` directory contains compact tables generated from larger local
analysis products.

Run the following commands after cloning:

```bash
git lfs install
git lfs pull
```

## Excluded data

Raw Kepler light curves, full ALDERAAN recovery archives, full eccentricity
posterior collections, cloud disks, temporary files, and private research logs
are not versioned. These products are large, regenerable, operationally
sensitive, or unsuitable for a public research release. The cloud scripts and
provenance tables document how the retained validation products were made.

## Source licenses

The MIT license covers original software in this repository. Third-party data,
article files, and archived source packages remain subject to their original
licenses and terms of use. Citation metadata for this repository and the Sagear
et al. article is in `CITATION.cff`.

## Missing reproduction products

The public article data do not include the final planet-level fit ledger, the
stellar-density rows or posterior samples used in the eccentricity calculation,
or the individual `(e, omega)` posteriors used for Table 3. The absence of these
products limits exact numerical reproduction. See
`docs/final_public_data_boundary.md` for the tests supporting this conclusion.
