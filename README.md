# From theory to macroecological evidence: a global test of island syndrome in plant floral traits

This repository contains the R workflow for testing whether floral trait simplification in island plants is a universal feature of plant island syndrome or a condition-dependent outcome of pollinator environments and self-compatibility.

The analysis uses a global island-flora dataset integrating island area, distance to continent, GBIF-based *Bombus* environmental mismatch, floral traits, pollination mode, and self-compatibility. The main analysis focuses on animal-pollinated species and asks whether plain-colored flowers and generalized floral forms are best explained by simple geographic filtering or by pathways involving pollinator mismatch and increased self-compatibility.

The central result is that the floral island syndrome is most clearly supported for **northern-temperate islands**. In this region, more isolated islands show greater *Bombus* environmental mismatch, mismatch is associated with higher self-compatibility, and self-compatible floras show higher proportions of plain-colored flowers and generalized floral forms. Projection to tropical and southern non-tropical islands is weaker, suggesting that the syndrome is condition-dependent rather than globally universal.

---

## Conceptual workflow

```text
Global island floras
        ↓
Island area + distance to continent
        ↓
GBIF-based Bombus environmental mismatch / PDI
        ↓
Floral color, floral form, pollination mode, self-compatibility
        ↓
Animal-pollinated species only
        ↓
Island-level proportions:
  - plain-colored flowers
  - generalized floral forms
  - self-compatible species
        ↓
Northern-temperate Bayesian path-model comparison
        ↓
Projection to tropical and southern non-tropical islands
```

---

## Main analysis workflow

Run the scripts in this order for the final analysis:

| Step | Script | Role in the analysis |
|------|--------|----------------------|
| 1 | `01_compile_global_island_polygons.R` | Compile global island polygons and retain islands ≥ 20 km². |
| 2 | `02_join_traits_and_calculate_continent_distance.R` | Join plant trait/island metadata and calculate distance from each island to the nearest continent. |
| 3 | `03_build_vif_filtered_environment_stack.R` | Build the VIF-filtered environmental raster stack used for *Bombus* niche mismatch. |
| 4 | `04_download_gbif_island_plant_species.R` | Download GBIF-based island plant species lists. |
| 5 | `05_download_and_thin_bombus_occurrences.R` | Download, clean, and thin *Bombus* occurrence data from GBIF. |
| 6 | `06_calculate_bombus_environmental_mismatch_pdi.R` | Calculate island-level *Bombus* environmental mismatch / PDI and join it to the plant table. |
| 7 | `07_reclassify_traits_and_filter_isolated_islands.R` | Standardize floral/reproductive traits and remove invalid or non-isolated island records. |
| 8 | `08_filter_animal_pollinated_and_assign_climate_zones.R` | Restrict the dataset to animal-pollinated species and assign tropical / temperate climate-zone labels. |
| 9 | `14_run_final_bayesian_path_model_analysis.R` | Run the final integrated analysis: island-level aggregation, Bayesian path-model comparison, path-effect estimation, figures, and external projection. |

In short: **01–08 build the analysis-ready dataset; 14 performs the final manuscript analysis.**

The final input expected by `14_run_final_bayesian_path_model_analysis.R` is:

```text
list_with_traits_island_distance_bom_with_PDI_reclassified_isolated_animal_trop_temp.csv
```

---

## Final Bayesian path-model analysis

`14_run_final_bayesian_path_model_analysis.R` compares four Bayesian piecewise path-model scenarios for northern-temperate islands using summed LOO-ELPD:

1. **M0: Geographic filtering**  
   Island isolation and island area directly predict floral trait composition.

2. **M2: Selfing syndrome**  
   Island isolation predicts the proportion of self-compatible species, which then predicts floral trait composition.

3. **M3: Full mediation**  
   Island isolation predicts *Bombus* environmental mismatch; mismatch predicts self-compatible flora; self-compatible flora predicts floral trait composition.

4. **M4: Partial mediation**  
   Same as M3, but *Bombus* mismatch can also directly predict floral traits.

The script also projects the representative northern-temperate pathway to tropical and southern non-tropical islands to test whether the same mechanism generalizes outside the northern-temperate region.

---

## Main response variables

At the island level, the final analysis uses:

| Variable | Meaning |
|----------|---------|
| `n_plain` / `n_color_total` | Number of plain-colored animal-pollinated species out of species with usable color data. |
| `n_generalized` / `n_form_total` | Number of species with generalized floral forms out of species with usable floral-form data. |
| `n_self_compatible` / `n_self_total` | Number of self-compatible species out of species with usable self-compatibility data. |
| `mean_PDI` | Island-level mean *Bombus* environmental mismatch / pollinator deficit index. |
| `dist_to_continent_km` | Distance from island to nearest continent, used as an isolation proxy. |
| `area_km2` | Island area, used as a geographic filtering / founder-effect proxy. |
| `climate_zone` | Tropical or temperate classification used for regional comparison and projection. |

Trait coding used in the final model:

```text
plain-colored flowers:
  white, green-brown-inconspicuous, other

conspicuous-colored flowers:
  blue-purple, red-pink, yellow-orange

generalized floral forms:
  actinomorphic/open, brush/puff, composite head, reduced wind-type

specialized floral forms:
  tubular, zygomorphic

self-compatible:
  SC, likely_SC

self-incompatible:
  SI, likely_SI, obligate_SI
```

---

## Main outputs from `14_run_final_bayesian_path_model_analysis.R`

The final script creates a manuscript-style output folder containing:

```text
01_data/
  island_level_all_regions_trait_counts.csv
  island_level_northern_temperate_model_ready.csv
  regional_island_level_summary.csv

03_model_comparison/
  scenario_comparison_total_LOO_ELPD.csv
  response_specific_LOO_ELPD_by_scenario.csv
  component_model_fixed_effects_all_scenarios.csv

04_path_effects/
  direct_indirect_and_total_path_effects.csv
  significant_edges_<representative_scenario>.csv

05_figures/
  Figure_1_scenario_comparison_total_ELPD.png/pdf
  Figure_2_response_specific_ELPD.png/pdf
  Figure_3_path_effect_estimates.png/pdf
  Figure_4_significant_path_diagram_<representative_scenario>.png/pdf/svg
  Figure_S1_external_projection_NLL.png/pdf

06_external_projection/
  external_projection_negative_log_likelihood.csv

07_metadata/
  analysis_metadata.txt
```

---

## Removed exploratory scripts

The previous repository included older exploratory or confirmatory scripts numbered `09`–`13`. These were removed from the main repository because they are not required for the abstract-aligned analysis.

| Removed script | Reason |
|----------------|--------|
| `09_recode_floral_traits.R` | Superseded by trait recoding inside `14_run_final_bayesian_path_model_analysis.R`. |
| `10_summarize_island_traits.R` | Utility RDS-to-long-table conversion; not used by the final workflow. |
| `11_assign_pollinator_preference.R` | Independent pollinator-preference heatmap; not part of the abstract's main analysis. |
| `12_fit_bayesian_path_models.R` | Older species-level / SEM-ready brms path model; superseded by the island-level Bayesian model comparison in `14_run_final_bayesian_path_model_analysis.R`. |
| `13_check_sem_lavaan.R` | Optional lavaan cross-check; not part of the final Bayesian workflow. |

---

## Interpretation

The analysis is designed to test whether floral trait simplification is best understood as a universal island syndrome or as a regionally contingent process. The current workflow supports the latter interpretation: floral simplification emerges most clearly in northern-temperate islands, where *Bombus* environmental mismatch and self-compatibility jointly shape floral trait composition.

LOO-ELPD is used to compare predictive consistency among candidate pathways. It should not be interpreted as causal proof on its own; biological interpretation should combine scenario ranking, posterior uncertainty, path-effect estimates, and regional projection performance.

---

## Dependencies

Core R packages:

```r
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
library(usdm)
library(rnaturalearth)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
```

---
