# Corrected tropical lineage entry/loading replay — 2026-10-03

## Purpose

This replay revisits an already frozen source-matched lineage-representation bridge after
the 24 September 2026 geography repair. It is a **measurement-repair replay**, not a new
assembly-model search.

All historical analysis choices are retained:
- tropical context;
- source-backed `native_nonendemic` stratum;
- evidence scopes `broad` and `broad_direct`;
- source matching `prevalence_richness`;
- minimum five represented genera;
- four frozen source modes;
- equal-island OLS with island area and climate PC1–PC4;
- spatial-block cluster-robust covariance.

Only `log_distance_to_continent_km` is replaced by the corrected spherical-coastline
distance.

The historical response is the frozen floral functional position
`(-large_bee_like + generalized_accessible) / 2`. Therefore this replay is a
supplementary assembly diagnostic and **not a direct decomposition of the current
seven-response H1**.

## Reproduction gate

Before replacing geography, the archived island-level scores reproduce all 24 frozen
slopes (2 evidence scopes × 4 source modes × 3 outcomes) to numerical precision:

- maximum absolute slope difference: **1.15 × 10^-16**
- maximum absolute SE difference: **1.15 × 10^-16**
- maximum absolute P-value difference: **2.17 × 10^-15**

## Corrected result

### Genus entry

The corrected distance slope for `entry_enrichment` is negative and FDR-supported in
all four source modes in both evidence scopes.

- broad: β range **-0.06861 to -0.05552**, maximum q = **0.02925**
- broad_direct: β range **-0.06713 to -0.06269**, maximum q = **0.02059**

Because the historical functional-position score is oriented so that larger values
represent the focal floral position, increasingly isolated tropical islands show a
source-mode-robust shift in **which source-available genera are represented**.

### Species-weighted enrichment

`species_enrichment` is likewise negative and FDR-supported in all four source modes:

- broad: β range **-0.08541 to -0.07380**, maximum q = **0.02860**
- broad_direct: β range **-0.08541 to -0.08153**, maximum q = **0.02606**

### Within-represented-genus loading

`loading_increment = species_enrichment - entry_enrichment` remains negative in every
source mode but is not FDR-supported in either evidence scope:

- broad: β range **-0.02046 to -0.01657**, maximum q = **0.15221**
- broad_direct: β range **-0.02122 to -0.01617**, maximum q = **0.17060**

Thus the supported species-weighted signal is already present at genus entry; there is no
supported additional shift from weighting species within genera that are already
represented.

## Interpretation

For this **historically frozen floral functional-position bridge**, the tropical
source-backed native-nonendemic assembly pattern is most consistent with **genus entry /
representation rather than additional within-genus species loading**.

This is concordant with the current seven-response H1 status-provenance audit, where the
strict tropical native-nonendemic H1 is supported before genus adjustment but is not
supported after genus residualization.

The two analyses are not identical estimands and must not be collapsed. This replay does
not show which ecological process controls genus entry and does not establish dispersal,
colonization, extinction, pollinator loss, or within-lineage evolution.

## Provenance

- validation workflow run: **37097462474**
- artifact: **11264653941**
- artifact digest:
  `sha256:afc980240819a0ce11fe8269fb7ec71d89e001184481dd1f0970e55de37625fb`
