# Plant Island Syndrome — Floral Trait Simplification in Global Island Floras

This repository contains R scripts for analyzing floral trait simplification across global island floras. The current main workflow builds the analysis-ready island flora dataset in steps 01–08, then runs the manuscript-level Bayesian path-model analysis in step 14.

---

## Main Analysis Workflow

Run the scripts in the following order for the final manuscript-style analysis:

| Step | Filename | Description |
|------|----------|-------------|
| 1 | `01_compile_island_metadata.R` | Compile global island polygons and retain islands ≥ 20 km² |
| 2 | `02_calc_continent_distance.R` | Join plant trait/island metadata and calculate minimum distance from each island to the nearest continent |
| 3 | `03_extract_island_climate.R` | Build the VIF-filtered environmental raster stack |
| 4 | `04_download_plant_traits_gbif.R` | Download GBIF-based island plant species lists |
| 5 | `05_get_bombus_occurrence.R` | Download, clean, and thin *Bombus* occurrence data from GBIF |
| 6 | `06_calc_bombus_unsuitability.R` | Calculate island-level *Bombus* climatic mismatch / PDI and join it to the plant table |
| 7 | `07_filter_oceanic_islands.R` | Reclassify floral traits and reproductive traits; remove invalid/non-isolated island records |
| 8 | `08_filter_by_pollination_climate.R` | Retain animal-pollinated species and add tropical / temperate climate-zone labels |
| 9 | `14_full_analysis_pipeline.R` | Full integrated analysis: island-level aggregation, Bayesian path-model fitting, LOO scenario comparison, path-effect estimation, figures, and external projection |

In short: **01–08 create the analysis-ready data; 14 runs the final integrated analysis.**

---

## Archived / Supplementary Scripts

The following scripts are useful for checking older analyses or producing supplementary outputs, but they are **not required in the main 01–08 → 14 workflow**.

| Filename | Status | Notes |
|----------|--------|-------|
| `09_recode_floral_traits.R` | Archived / superseded | Earlier SEM-ready binary trait recoding. The same conceptual recoding is now performed inside `14_full_analysis_pipeline.R`. |
| `10_summarize_island_traits.R` | Archived utility | Converts the GBIF facet RDS into an island × species long table. This output is not used by the final workflow. |
| `11_assign_pollinator_preference.R` | Supplementary | Independent pollinator preference / Dirichlet heatmap analysis. Keep only if used as a supplementary figure. |
| `12_fit_bayesian_path_models.R` | Archived / superseded | Earlier brms SEM-style path model using species-level SEM-ready data. Superseded by the island-level Bayesian model comparison in `14_full_analysis_pipeline.R`. |
| `13_check_sem_lavaan.R` | Optional cross-check | lavaan confirmatory SEM check. Not part of the main Bayesian workflow. |

Recommended future folder organization:

```text
archive/
  09_recode_floral_traits.R
  10_summarize_island_traits.R
  12_fit_bayesian_path_models.R

supplementary/
  11_assign_pollinator_preference.R
  13_check_sem_lavaan.R
```

---

## Data Files

| File | Description |
|------|-------------|
| `islandinfo.csv` | Island metadata |
| `list.csv` | Island-level plant species list |
| `traits.csv` | Species-level floral and reproductive traits |
| `島嶼花形質.xlsx` | Raw floral trait data sheet |
| `bombus_niche_5_95.rds` | *Bombus* spp. climate niche, defined by the 5–95th percentile range per environmental variable |
| `island_PDI_bombus_CONTINUOUS.rds` | Island-level *Bombus* climatic mismatch / PDI |
| `island_flowering_species_facet_progress.rds` | Intermediate GBIF facet output from island plant species download |
| `data_global_islands/` | Global island polygon data and derived island tables |

---

## Key Variables

| Variable | Description |
|----------|-------------|
| `dist_to_continent_km` | Minimum distance to nearest continent, in km; proxy for island isolation |
| `mean_PDI` | Island-level mean *Bombus* climatic mismatch / PDI |
| `prop_plain` / `n_plain` | Proportion / count of plain-colored flowers per island, based on white, green-brown-inconspicuous, or other simple colors |
| `prop_generalized` / `n_generalized` | Proportion / count of generalized floral forms per island |
| `prop_self_compatible` / `n_self_compatible` | Proportion / count of self-compatible species per island |

---

## Path Model Scenarios

`14_full_analysis_pipeline.R` compares four Bayesian piecewise path-model scenarios using summed LOO-ELPD:

1. **M0: Geographic filtering** — island isolation and area directly predict floral traits.
2. **M2: Selfing syndrome** — isolation predicts self-compatible flora, which then predicts floral traits.
3. **M3: Full mediation** — isolation predicts *Bombus* climatic mismatch; mismatch predicts self-compatible flora; self-compatible flora predicts floral traits.
4. **M4: Partial mediation** — same as M3, but mismatch can also directly predict floral traits.

LOO-ELPD evaluates predictive consistency among candidate path scenarios. It does not by itself prove causality; biological interpretation should combine scenario ranking, posterior uncertainty, and path-effect estimates.

---

## Dependencies

```r
# Core packages
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(brms)
library(posterior)
library(loo)
library(rgbif)
library(sf)
library(terra)
```

---
