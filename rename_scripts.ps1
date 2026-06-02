# Run this script in the root of your island repository
# Usage: powershell -ExecutionPolicy Bypass -File .\rename_scripts.ps1

# R scripts
git mv islandinfo.R              01_compile_island_metadata.R
git mv distance.R                02_calc_continent_distance.R
git mv islandenv.R               03_extract_island_climate.R
git mv islandgbif.R              04_download_plant_traits_gbif.R
git mv getbombusforislands.R     05_get_bombus_occurrence.R
git mv pdi.R                     06_calc_bombus_unsuitability.R
git mv reclassify_isolated.R     07_filter_oceanic_islands.R
git mv animal_trop_temp.R        08_filter_by_pollination_climate.R
git mv binary.R                  09_recode_floral_traits.R
git mv list.R                    10_summarize_island_traits.R
git mv pollpreference.R          11_assign_pollinator_preference.R
git mv "sem(brms).R"             12_fit_bayesian_path_models.R
git mv "sem(lavaan).R"           13_check_sem_lavaan.R
git mv NorthernTemperate_IslandFloralSyndrome_BayesianPathModelAnalysis.R 14_full_analysis_pipeline.R

# Data files
git mv islandinfo.csv                          island_metadata.csv
git mv list.csv                                island_species_list.csv
git mv traits.csv                              species_floral_traits.csv
git mv "島嶼花形質.xlsx"                        island_floral_traits_raw.xlsx
git mv bombus_niche_5_95.rds                   bombus_climate_niche.rds
git mv island_PDI_bombus_CONTINUOUS.rds        island_bombus_unsuitability.rds
git mv island_flowering_species_facet_progress.rds island_species_traits_matrix.rds

git commit -m "rename R scripts and data files to follow analysis pipeline order"
git push
