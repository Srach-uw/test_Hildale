# Berger Density Provenance Audit

## Finding
Sagear's published methods prose names Berger et al. (2018) for the stellar-density prior, while the public source comments and conclusion point to Berger et al. (2020). Berger 2018 Table 1 supplies radius, evolution, and binary flags, but no stellar-density column, so the public materials do not uniquely identify an exact B18 density input.

## Evidence
- Berger 2020 `table2.dat.gz`: 186,301 rows, with mass, logg, radius, and asymmetric log-density columns.
- Berger 2018 table: 177,911 rows, with radius/evolution/binary fields and no density field.
- Exact published-inventory hosts audited: 1,888; B20 density coverage: 1,888; B18 radius coverage: 1,887.
- The operational extractor uses the B20 published density and its asymmetric uncertainties because that is the only public source with the required density posterior inputs and it is supported by the source comments/conclusion; this remains an explicit operational choice, not proof of Sagear's hidden run configuration.

## Diagnostic comparison
For hosts with both sources, a hybrid density using B20 mass and B18 radius differs from B20's published density by median delta_log10(rho) = 0.0010 dex (16th-84th percentile -0.0302 to 0.0302). This is a sensitivity diagnostic, not evidence that Sagear used the hybrid.

## Decision
Keep B20 as the operational branch, but do not describe it as definitively Sagear-equivalent. Treat B18-radius/B20-mass and the KG-RADII-derived B18 diagnostic as explicit sensitivity branches, and request author clarification about the apparent 2018/2020 provenance mismatch before calling the replication exact.
