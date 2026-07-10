# Frozen v1 analysis

This directory preserves the original R-based island floral-syndrome
analysis exactly for provenance. It is not maintained as active code.

`r/01`–`r/09` are the original end-to-end scripts, including the
original 20 km² island threshold, local desktop paths, GBIF queries,
Bombus/PDI construction, trait recoding, and final Bayesian path model.

The frozen v1 input and intermediate data retained beside the scripts are:

- `bombus_climate_niche.rds`
- `island_bombus_unsuitability.rds`
- `island_floral_traits_raw.xlsx`
- `island_metadata.csv`
- `island_species_list.csv`
- `island_species_traits_matrix.rds`
- `species_floral_traits.csv`

The active implementation is the repository root v2 pipeline.
