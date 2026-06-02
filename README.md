# Plant Island Syndrome — Floral Trait Simplification in Global Island Floras

This repository contains all R scripts for the analysis of floral trait simplification across global island floras, testing ecological conditions and pathways underlying plant island syndrome.

---

## Analysis Pipeline

Scripts should be run in the following order:

| Step | Filename | Description |
|------|----------|-------------|
| 1 | `01_compile_island_metadata.R` | Compile island area, geographic coordinates, and climate zone classification |
| 2 | `02_calc_continent_distance.R` | Calculate minimum distance from each island to the nearest continent |
| 3 | `03_extract_island_climate.R` | Extract climate variables for each island |
| 4 | `04_download_plant_traits_gbif.R` | Download and clean GBIF-based plant occurrence and trait data |
| 5 | `05_get_bombus_occurrence.R` | Download and clean *Bombus* spp. occurrence data from GBIF |
| 6 | `06_calc_bombus_unsuitability.R` | Calculate *Bombus* climatic unsuitability (PDI) for each island |
| 7 | `07_filter_oceanic_islands.R` | Filter and reclassify islands to retain truly isolated oceanic islands |
| 8 | `08_filter_by_pollination_climate.R` | Restrict dataset to animal-pollinated species; classify islands by climate zone (north-temperate, tropical, southern non-tropical) |
| 9 | `09_recode_floral_traits.R` | Recode floral color (inconspicuous/conspicuous) and floral form (open/specialized) as binary traits |
| 10 | `10_summarize_island_traits.R` | Aggregate trait proportions at the island level; z-score standardization |
| 11 | `11_assign_pollinator_preference.R` | Assign pollinator preference categories to species |
| 12 | `12_fit_bayesian_path_models.R` | Fit four Bayesian piecewise path-model scenarios using brms; compare with LOO-ELPD |
| 13 | `13_check_sem_lavaan.R` | Confirmatory SEM using lavaan (cross-check) |
| 14 | `14_full_analysis_pipeline.R` | Full integrated analysis: north-temperate model fitting and projection to tropical and southern non-tropical islands |

---

## Data Files

| File | Description |
|------|-------------|
| `islandinfo.csv` | Island metadata (area, coordinates, climate zone) |
| `list.csv` | Island-level plant species list |
| `traits.csv` | Species-level floral and reproductive traits |
| `島嶼花形質.xlsx` | Raw floral trait data sheet |
| `bombus_niche_5_95.rds` | *Bombus* spp. climate niche (5–95th percentile range per variable) |
| `island_PDI_bombus_CONTINUOUS.rds` | Island-level *Bombus* climatic unsuitability (PDI, continuous) |
| `island_flowering_species_facet_progress.rds` | Intermediate output: island × species trait matrix |
| `data_global_islands/` | Supplementary global island data |

---

## Key Variables

| Variable | Description |
|----------|-------------|
| `dist_continent` | Minimum distance to nearest continent (km); proxy for island isolation |
| `bombus_unsuitability` | *Bombus* climatic unsuitability index (PDI); measures how far island climate deviates from the *Bombus* realized niche |
| `prop_inconspicuous` | Proportion of inconspicuously colored flowers (white, green, or brown) per island |
| `prop_open_form` | Proportion of open or brush-type flowers (actinomorphic, brush-type, or composite) per island |
| `prop_selfcompat` | Proportion of self-compatible species per island |

---

## Path Model Scenarios

Four Bayesian piecewise path-model scenarios were compared using LOO-ELPD:

1. **Geographic Filtering** — Island isolation and area directly predict floral traits
2. **Selfing Syndrome** — Isolation promotes self-compatibility, which in turn shapes floral traits
3. **Full Mediation** — *Bombus* climatic unsuitability drives self-compatibility, which then leads to floral simplification
4. **Partial Mediation** *(best-fit model)* — As above, with *Bombus* climatic unsuitability additionally exerting direct selection on floral color independently of self-compatibility

---

## Dependencies

```r
# Core packages
library(brms)       # Bayesian path models
library(lavaan)     # Confirmatory SEM
library(rgbif)      # GBIF data access
library(sf)         # Spatial operations
library(tidyverse)  # Data wrangling
```

---
