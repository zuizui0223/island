# v1 artifact migration manifest

This manifest was generated from file references in `legacy/v1/r/*.R`.
Only files already tracked in this repository were moved. Local Desktop,
external-drive, or previously unversioned inputs cannot be recovered by Git.

## Moved tracked artifacts

- `data_global_islands/global_islands_over20km2.csv` → `artifacts/data_global_islands/global_islands_over20km2.csv`

## Referenced by v1 but not tracked in this repository

- `README_NorthernTemperate_IslandFloralSyndrome_BayesianPathModelAnalysis.txt`
- `_loo.rds`
- `analysis_metadata.txt`
- `apis_occ_thinned.csv`
- `bombus_niche_5_95.rds`
- `bombus_occ_thinned.csv`
- `component_model_fixed_effects_all_scenarios.csv`
- `diagnostic_analysis_group_counts.csv`
- `diagnostic_flower_color_recoding.csv`
- `diagnostic_flower_form_recoding.csv`
- `diagnostic_self_compatibility_recoding.csv`
- `direct_indirect_and_total_path_effects.csv`
- `env_island_vif.tif`
- `env_stack_vif.tif`
- `external_projection_negative_log_likelihood.csv`
- `global_islands_over20km2.gpkg`
- `island_PDI_bombus_CONTINUOUS.csv`
- `island_PDI_bombus_CONTINUOUS.rds`
- `island_PDI_bombus_FAST.csv`
- `island_flowering_species_facet_progress.rds`
- `island_level_all_regions_trait_counts.csv`
- `island_level_northern_temperate_model_ready.csv`
- `island_zone_table.csv`
- `islandinfo.csv`
- `list.csv`
- `list_with_traits_island_distance_bom.csv`
- `list_with_traits_island_distance_bom_with_PDI.csv`
- `list_with_traits_island_distance_bom_with_PDI_reclassified.csv`
- `list_with_traits_island_distance_bom_with_PDI_reclassified_isolated.csv`
- `list_with_traits_island_distance_bom_with_PDI_reclassified_isolated_animal_trop_temp.csv`
- `regional_island_level_summary.csv`
- `response_specific_LOO_ELPD_by_scenario.csv`
- `scenario_comparison_total_LOO_ELPD.csv`
- `significant_edges_<representative_scenario>.csv`
- `traits.csv`
- `vif_results.csv`

## Protection rule

The migration never moves `data/v2/`, `src/`, `tests/`, `config/`, `docs/`,
or active GitHub Actions workflows. Those remain part of v2.
