# Published journal asset inventory

Final article: https://doi.org/10.3847/1538-3881/ae71bf

The publisher article and article-data pages were inspected directly on
2026-07-24. The public machine-readable assets are:

| Asset | Local file | Bytes | SHA-256 | Scientific content |
|---|---|---:|---|---|
| Table 1 | `inputs/sagear_2026_published/sagear2026_table1_kinematic_hosts_mrt.txt` | 270060 | `92280ede0c828413abb7c8314ac6f35b0ccc3e68aabbcf6a94075e84b30e76ed` | 1,888 host KICs, Gaia identifiers, Galactocentric coordinates and velocities, `P_thick`, and final disk assignment |
| Table 2 | `inputs/sagear_2026_published/sagear2026_table2_transit_priors_ascii.txt` | 958 | `3cd83106e11820423b207ba93a66c8c76e6dcf2092313211aef42c0efbdce7df` | ALDERAAN transit parameters and printed priors |
| Table 3 | `inputs/sagear_2026_published/sagear2026_table3_population_results_ascii.txt` | 1329 | `56a079c52c804f5021c51b56b1c5dd8b11d873b9c251e4be1520257fb1ea4fa9` | Published Rayleigh, Beta, monotonic-Beta, and half-Gaussian summaries |

The newly retrieved publisher filenames `ajae71bft1_mrt.txt`,
`ajae71bft2_ascii.txt`, and `ajae71bft3_ascii.txt` are byte-identical to the
three canonical local files above.

No public journal asset contains:

- planet-level ALDERAAN FITS or transit-shape posterior samples;
- planet-level `e, omega` posterior samples;
- the final planet inclusion table;
- the 26 visual/nonconvergence exclusions;
- the stellar-density prior table used in post-model importance sampling;
- the NumPyro population-model source or posterior chains.

Crossref has no related-data relation for the DOI. A DataCite search found no
matching deposit. The author GitHub account contains ALDERAAN and
`photoeccentric` source/catalog repositories but no repository for this
paper's posterior archive. The cited Gilbert et al. 2025 data-availability
statement likewise releases the upstream Kepler light curves and Berger
stellar catalog, not its 1,646 planet-level eccentricity posteriors.

Therefore no public posterior product can replace the missing ALDERAAN fits.
The public Table 1 should be treated as authoritative for host disk labels,
Table 2 as the printed transit-prior contract, and Table 3 as the final
population-level acceptance target.
