############################################################
## Northern-temperate island floral syndrome
## Bayesian path-model analysis of island floral trait composition
##
## Standalone manuscript reproducibility script
##
## This script:
##   1. reads the original CSV directly
##   2. aggregates species-level records to island-level flora composition
##   3. fits Bayesian piecewise path models from scratch
##   4. compares candidate causal scenarios using summed LOO-ELPD
##   5. estimates direct and indirect effects
##   6. saves manuscript-ready tables and figures

############################################################

rm(list = ls())
gc()

############################################################
## 1. User settings
############################################################

INPUT_FILE <- "/Users/rachelzhang/Desktop/list_with_traits_island_distance_bom_with_PDI_reclassified_isolated_animal_trop_temp.csv"

OUT_ROOT <- "/Users/rachelzhang/Desktop/NorthernTemperate_IslandFloralSyndrome_BayesianPathModelAnalysis"

RUN_CHAINS <- 2
RUN_CORES  <- 1
RUN_ITER   <- 2000
RUN_WARMUP <- 1000
RUN_SEED   <- 20260525

USE_CMDSTANR_IF_AVAILABLE <- TRUE
LOO_CORES <- 1

SIG_Q_LOW  <- 0.025
SIG_Q_HIGH <- 0.975

DEBUG_SUBSET <- FALSE
DEBUG_N_ISLANDS <- 200

############################################################
## 2. Packages
############################################################

need_pkgs <- c(
  "data.table", "dplyr", "tidyr", "stringr", "ggplot2",
  "brms", "posterior", "loo", "purrr", "tibble", "scales",
  "DiagrammeR", "DiagrammeRsvg", "rsvg"
)

install_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    install.packages(miss, repos = "https://cran.rstudio.com")
  }
}

install_if_missing(need_pkgs)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(brms)
  library(posterior)
  library(loo)
  library(purrr)
  library(tibble)
  library(scales)
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
})

############################################################
## 3. Output folders
############################################################

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "01_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "02_bayesian_component_models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "03_model_comparison"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "04_path_effects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "05_figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "06_external_projection"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "07_metadata"), recursive = TRUE, showWarnings = FALSE)

options(mc.cores = RUN_CORES)
set.seed(RUN_SEED)

backend_use <- "rstan"
if (USE_CMDSTANR_IF_AVAILABLE && requireNamespace("cmdstanr", quietly = TRUE)) {
  ok_cmdstan <- FALSE
  try({
    ok_cmdstan <- !is.null(cmdstanr::cmdstan_path()) &&
      nzchar(cmdstanr::cmdstan_path())
  }, silent = TRUE)
  if (isTRUE(ok_cmdstan)) backend_use <- "cmdstanr"
}

message("============================================================")
message("Northern-temperate island floral syndrome analysis")
message("============================================================")
message("Input file: ", INPUT_FILE)
message("Output directory: ", OUT_ROOT)
message("Backend: ", backend_use)
message("chains=", RUN_CHAINS,
        " cores=", RUN_CORES,
        " iter=", RUN_ITER,
        " warmup=", RUN_WARMUP,
        " seed=", RUN_SEED)
message("============================================================")

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found: ", INPUT_FILE)
}

############################################################
## 4. Helper functions
############################################################

clean_chr <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_to_lower(x)
  x
}

zscale <- function(x) {
  as.numeric(scale(x))
}

safe_div <- function(a, b) {
  ifelse(is.na(b) | b <= 0, NA_real_, a / b)
}

pick_first_existing <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

invlogit <- function(x) {
  1 / (1 + exp(-x))
}

theme_manuscript <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Helvetica"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}

fit_brm_component <- function(formula, data, family, model_name) {
  message("\n--- fitting ", model_name, " | n islands = ", nrow(data), " ---")
  
  fit <- brm(
    formula = formula,
    data = data,
    family = family,
    chains = RUN_CHAINS,
    cores = RUN_CORES,
    iter = RUN_ITER,
    warmup = RUN_WARMUP,
    seed = RUN_SEED,
    backend = backend_use,
    refresh = max(100, floor(RUN_ITER / 6)),
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    save_pars = save_pars(all = TRUE),
    file = NULL
  )
  
  saveRDS(
    fit,
    file.path(
      OUT_ROOT,
      "02_bayesian_component_models",
      paste0(model_name, ".rds")
    )
  )
  
  fit
}

loo_component <- function(fit, model_name) {
  message("Computing LOO: ", model_name)
  
  out <- loo::loo(
    fit,
    cores = LOO_CORES,
    moment_match = FALSE
  )
  
  saveRDS(
    out,
    file.path(
      OUT_ROOT,
      "02_bayesian_component_models",
      paste0(model_name, "_loo.rds")
    )
  )
  
  out
}

extract_fixef <- function(fit, scenario_name, response_name) {
  fe <- as.data.frame(brms::fixef(fit, probs = c(SIG_Q_LOW, SIG_Q_HIGH)))
  fe$term <- rownames(fe)
  rownames(fe) <- NULL
  
  qlow_name <- paste0("Q", SIG_Q_LOW * 100)
  qhigh_name <- paste0("Q", SIG_Q_HIGH * 100)
  
  fe %>%
    transmute(
      scenario = scenario_name,
      response = response_name,
      term = term,
      mean = Estimate,
      se = Est.Error,
      q025 = .data[[qlow_name]],
      q975 = .data[[qhigh_name]],
      significant = q025 > 0 | q975 < 0
    )
}

posterior_beta <- function(fit, term) {
  dr <- posterior::as_draws_df(fit)
  nm <- paste0("b_", term)
  if (!nm %in% names(dr)) {
    return(rep(NA_real_, nrow(dr)))
  }
  as.numeric(dr[[nm]])
}

summ_draws <- function(x, label) {
  q025_val <- as.numeric(quantile(x, SIG_Q_LOW, na.rm = TRUE))
  q975_val <- as.numeric(quantile(x, SIG_Q_HIGH, na.rm = TRUE))
  
  tibble(
    effect = label,
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q025 = q025_val,
    q975 = q975_val,
    pd = max(mean(x > 0, na.rm = TRUE), mean(x < 0, na.rm = TRUE)),
    significant = q025_val > 0 | q975_val < 0
  )
}

############################################################
## 5. Load and clean input data
############################################################

raw <- data.table::fread(INPUT_FILE, encoding = "UTF-8", data.table = FALSE)

message("Rows loaded: ", nrow(raw))
message("Columns loaded: ", ncol(raw))

nms <- names(raw)

col_lat <- pick_first_existing(
  nms,
  c("centroid_lat", "latitude", "lat", "decimalLatitude")
)

col_lon <- pick_first_existing(
  nms,
  c("centroid_lon", "longitude", "lon", "decimalLongitude")
)

col_dist <- pick_first_existing(
  nms,
  c(
    "dist_to_continent_km",
    "distance_to_mainland",
    "dist_to_mainland",
    "min_distance_to_continent",
    "NEAR_DIST",
    "distance_km"
  )
)

col_area <- pick_first_existing(
  nms,
  c("area_km2", "IslandArea", "Area_Geode", "Shape_Area")
)

col_mismatch <- pick_first_existing(
  nms,
  c("mean_PDI", "PDI", "mean_pdi", "bumblebee_mismatch")
)

col_climate <- pick_first_existing(
  nms,
  c("climate_zone", "climate", "zone")
)

col_island <- pick_first_existing(
  nms,
  c("island_id", "USGS_ISID", "Name_USGSO", "NAME_wcmcI", "NAME_LOCAL")
)

col_species <- pick_first_existing(
  nms,
  c("species", "i.species", "taxon", "scientificName")
)

required_cols <- c(
  col_lat,
  col_dist,
  col_area,
  col_mismatch,
  col_climate,
  col_island,
  col_species
)

if (any(is.na(required_cols))) {
  stop(
    "Missing required columns. Detected: ",
    "latitude=", col_lat, "; ",
    "distance=", col_dist, "; ",
    "area=", col_area, "; ",
    "mismatch=", col_mismatch, "; ",
    "climate=", col_climate, "; ",
    "island=", col_island, "; ",
    "species=", col_species
  )
}

if (!"flower_color_simple" %in% nms) {
  stop("flower_color_simple column not found")
}

if (!"flower_shape_simple" %in% nms) {
  stop("flower_shape_simple column not found")
}

if (!"self_incompatibility" %in% nms) {
  stop("self_incompatibility column not found")
}

message("Detected columns:")
message("  latitude: ", col_lat)
message("  longitude: ", col_lon)
message("  distance: ", col_dist)
message("  area: ", col_area)
message("  bumblebee-environment mismatch: ", col_mismatch)
message("  climate zone: ", col_climate)
message("  island: ", col_island)
message("  species: ", col_species)

############################################################
## 6. Trait recoding
############################################################

plain_color_levels <- c(
  "white",
  "green_brown_inconspicuous",
  "other"
)

conspicuous_color_levels <- c(
  "blue_purple",
  "red_pink",
  "yellow_orange"
)

generalized_form_levels <- c(
  "actinomorphic_open",
  "brush_puff",
  "composite_head",
  "reduced_wind"
)

specialized_form_levels <- c(
  "tubular",
  "zygomorphic"
)

self_compatible_levels <- c(
  "sc",
  "likely_sc"
)

self_incompatible_levels <- c(
  "si",
  "likely_si",
  "obligate_si"
)

dat <- raw %>%
  mutate(
    island_id_use = as.character(.data[[col_island]]),
    species_use = as.character(.data[[col_species]]),
    
    lat = as.numeric(.data[[col_lat]]),
    lon = if (!is.na(col_lon)) as.numeric(.data[[col_lon]]) else NA_real_,
    
    dist_km = as.numeric(.data[[col_dist]]),
    area_km2 = as.numeric(.data[[col_area]]),
    bumblebee_mismatch = as.numeric(.data[[col_mismatch]]),
    
    climate_zone = clean_chr(.data[[col_climate]]),
    hemisphere = ifelse(lat >= 0, "northern", "southern"),
    
    analysis_group = case_when(
      hemisphere == "northern" & climate_zone == "temperate" ~ "temperate_NH",
      climate_zone == "tropical" ~ "tropical_all",
      hemisphere == "southern" & climate_zone != "tropical" ~ "southern_non_tropical",
      TRUE ~ NA_character_
    ),
    
    flower_color_simple_clean = clean_chr(flower_color_simple),
    flower_shape_simple_clean = clean_chr(flower_shape_simple),
    self_incompatibility_clean = clean_chr(self_incompatibility),
    
    color_plain = case_when(
      flower_color_simple_clean %in% plain_color_levels ~ 1L,
      flower_color_simple_clean %in% conspicuous_color_levels ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    form_generalized = case_when(
      flower_shape_simple_clean %in% generalized_form_levels ~ 1L,
      flower_shape_simple_clean %in% specialized_form_levels ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    self_compatible = case_when(
      self_incompatibility_clean %in% self_compatible_levels ~ 1L,
      self_incompatibility_clean %in% self_incompatible_levels ~ 0L,
      TRUE ~ NA_integer_
    )
  )

write.csv(
  dat %>%
    count(analysis_group, climate_zone, hemisphere),
  file.path(OUT_ROOT, "01_data", "diagnostic_analysis_group_counts.csv"),
  row.names = FALSE
)

write.csv(
  dat %>%
    count(flower_color_simple_clean, color_plain),
  file.path(OUT_ROOT, "01_data", "diagnostic_flower_color_recoding.csv"),
  row.names = FALSE
)

write.csv(
  dat %>%
    count(flower_shape_simple_clean, form_generalized),
  file.path(OUT_ROOT, "01_data", "diagnostic_flower_form_recoding.csv"),
  row.names = FALSE
)

write.csv(
  dat %>%
    count(self_incompatibility_clean, self_compatible),
  file.path(OUT_ROOT, "01_data", "diagnostic_self_compatibility_recoding.csv"),
  row.names = FALSE
)

############################################################
## 7. Island-level aggregation
############################################################

make_island_level_data <- function(x) {
  x %>%
    filter(
      !is.na(analysis_group),
      !is.na(island_id_use),
      is.finite(dist_km),
      is.finite(area_km2),
      is.finite(lat)
    ) %>%
    group_by(analysis_group, island_id_use) %>%
    summarise(
      lat = mean(lat, na.rm = TRUE),
      lon = mean(lon, na.rm = TRUE),
      dist_km = median(dist_km, na.rm = TRUE),
      area_km2 = median(area_km2, na.rm = TRUE),
      bumblebee_mismatch = mean(bumblebee_mismatch, na.rm = TRUE),
      
      n_species = n_distinct(species_use),
      
      n_color_total = sum(!is.na(color_plain)),
      n_plain = sum(color_plain == 1, na.rm = TRUE),
      
      n_form_total = sum(!is.na(form_generalized)),
      n_generalized = sum(form_generalized == 1, na.rm = TRUE),
      
      n_self_total = sum(!is.na(self_compatible)),
      n_self_compatible = sum(self_compatible == 1, na.rm = TRUE),
      
      prop_plain = safe_div(n_plain, n_color_total),
      prop_generalized = safe_div(n_generalized, n_form_total),
      prop_self_compatible = safe_div(n_self_compatible, n_self_total),
      
      .groups = "drop"
    ) %>%
    mutate(
      logdist = log1p(dist_km),
      logarea = log1p(area_km2),
      z_logdist = zscale(logdist),
      z_logarea = zscale(logarea),
      z_mismatch = zscale(bumblebee_mismatch)
    )
}

island_all <- make_island_level_data(dat)

write.csv(
  island_all,
  file.path(OUT_ROOT, "01_data", "island_level_all_regions_trait_counts.csv"),
  row.names = FALSE
)

region_summary <- island_all %>%
  group_by(analysis_group) %>%
  summarise(
    n_islands = n(),
    mean_distance_km = mean(dist_km, na.rm = TRUE),
    mean_area_km2 = mean(area_km2, na.rm = TRUE),
    mean_bumblebee_mismatch = mean(bumblebee_mismatch, na.rm = TRUE),
    mean_prop_self_compatible = mean(prop_self_compatible, na.rm = TRUE),
    mean_prop_plain = mean(prop_plain, na.rm = TRUE),
    mean_prop_generalized = mean(prop_generalized, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  region_summary,
  file.path(OUT_ROOT, "01_data", "regional_island_level_summary.csv"),
  row.names = FALSE
)

print(region_summary)

island_tnh <- island_all %>%
  filter(analysis_group == "temperate_NH") %>%
  filter(
    is.finite(z_logdist),
    is.finite(z_logarea),
    is.finite(z_mismatch),
    n_color_total > 0,
    n_form_total > 0,
    n_self_total > 0
  )

if (DEBUG_SUBSET) {
  island_tnh <- island_tnh %>%
    slice_sample(n = min(DEBUG_N_ISLANDS, nrow(island_tnh)))
}

message("Northern-temperate model-ready islands: ", nrow(island_tnh))

if (nrow(island_tnh) < 30) {
  stop("Too few northern-temperate islands after aggregation.")
}

write.csv(
  island_tnh,
  file.path(OUT_ROOT, "01_data", "island_level_northern_temperate_model_ready.csv"),
  row.names = FALSE
)

############################################################
## 8. Scenario definitions
############################################################
##
## M0: geographic filter
##   distance and island area directly predict final floral traits.
##
## M2: selfing syndrome
##   distance predicts self-compatible flora.
##   self-compatible flora predicts final floral traits.
##   island area is included only in final floral trait models.
##
## M3: full mediation
##   distance predicts bumblebee-environment mismatch.
##   mismatch predicts self-compatible flora.
##   self-compatible flora predicts final floral traits.
##   island area is included only in final floral trait models.
##
## M4: partial mediation
##   same as M3, but mismatch can also directly predict final floral traits.
##   island area is included only in final floral trait models.
############################################################

scenario_formulas <- list(
  M0_geographic_filter = list(
    mismatch = bf(z_mismatch ~ 1 + z_logdist),
    self = bf(n_self_compatible | trials(n_self_total) ~ 1 + z_logdist),
    color = bf(n_plain | trials(n_color_total) ~ 1 + z_logdist + z_logarea),
    form = bf(n_generalized | trials(n_form_total) ~ 1 + z_logdist + z_logarea)
  ),
  
  M2_selfing_syndrome = list(
    mismatch = bf(z_mismatch ~ 1 + z_logdist),
    self = bf(n_self_compatible | trials(n_self_total) ~ 1 + z_logdist),
    color = bf(n_plain | trials(n_color_total) ~ 1 + prop_self_compatible + z_logarea),
    form = bf(n_generalized | trials(n_form_total) ~ 1 + prop_self_compatible + z_logarea)
  ),
  
  M3_full_mediation = list(
    mismatch = bf(z_mismatch ~ 1 + z_logdist),
    self = bf(n_self_compatible | trials(n_self_total) ~ 1 + z_mismatch),
    color = bf(n_plain | trials(n_color_total) ~ 1 + prop_self_compatible + z_logarea),
    form = bf(n_generalized | trials(n_form_total) ~ 1 + prop_self_compatible + z_logarea)
  ),
  
  M4_partial_mediation = list(
    mismatch = bf(z_mismatch ~ 1 + z_logdist),
    self = bf(n_self_compatible | trials(n_self_total) ~ 1 + z_mismatch),
    color = bf(n_plain | trials(n_color_total) ~ 1 + z_mismatch + prop_self_compatible + z_logarea),
    form = bf(n_generalized | trials(n_form_total) ~ 1 + z_mismatch + prop_self_compatible + z_logarea)
  )
)

response_families <- list(
  mismatch = gaussian(),
  self = binomial(),
  color = binomial(),
  form = binomial()
)

response_labels <- c(
  mismatch = "Bumblebee-environment mismatch",
  self = "Self-compatible flora",
  color = "Plain-colored flowers",
  form = "Generalized floral forms"
)

scenario_labels <- c(
  M0_geographic_filter = "M0: geographic filter",
  M2_selfing_syndrome = "M2: selfing syndrome",
  M3_full_mediation = "M3: full mediation",
  M4_partial_mediation = "M4: partial mediation"
)

############################################################
## 9. Fit all Bayesian component models from scratch
############################################################

models <- list()
loos <- list()
fixed_effects_all <- list()

for (scenario_name in names(scenario_formulas)) {
  models[[scenario_name]] <- list()
  loos[[scenario_name]] <- list()
  
  for (response_name in names(scenario_formulas[[scenario_name]])) {
    
    model_name <- paste(scenario_name, response_name, sep = "__")
    
    fit <- fit_brm_component(
      formula = scenario_formulas[[scenario_name]][[response_name]],
      data = island_tnh,
      family = response_families[[response_name]],
      model_name = model_name
    )
    
    models[[scenario_name]][[response_name]] <- fit
    loos[[scenario_name]][[response_name]] <- loo_component(fit, model_name)
    
    fixed_effects_all[[model_name]] <- extract_fixef(
      fit = fit,
      scenario_name = scenario_name,
      response_name = response_name
    )
  }
}

fixed_effects_df <- bind_rows(fixed_effects_all)

write.csv(
  fixed_effects_df,
  file.path(OUT_ROOT, "03_model_comparison", "component_model_fixed_effects_all_scenarios.csv"),
  row.names = FALSE
)

############################################################
## 10. Scenario comparison by summed LOO-ELPD
############################################################

loo_by_component <- map_dfr(names(loos), function(scenario_name) {
  map_dfr(names(loos[[scenario_name]]), function(response_name) {
    x <- loos[[scenario_name]][[response_name]]
    
    tibble(
      scenario = scenario_name,
      scenario_label = scenario_labels[[scenario_name]],
      response = response_name,
      response_label = response_labels[[response_name]],
      elpd_loo = x$estimates["elpd_loo", "Estimate"],
      se_elpd_loo = x$estimates["elpd_loo", "SE"],
      p_loo = x$estimates["p_loo", "Estimate"],
      looic = x$estimates["looic", "Estimate"],
      se_looic = x$estimates["looic", "SE"]
    )
  })
}) %>%
  group_by(response) %>%
  mutate(
    best_elpd_response = max(elpd_loo, na.rm = TRUE),
    delta_elpd_response = elpd_loo - best_elpd_response,
    rank_response = rank(-elpd_loo, ties.method = "min")
  ) %>%
  ungroup()

write.csv(
  loo_by_component,
  file.path(OUT_ROOT, "03_model_comparison", "response_specific_LOO_ELPD_by_scenario.csv"),
  row.names = FALSE
)

loo_total <- loo_by_component %>%
  group_by(scenario, scenario_label) %>%
  summarise(
    elpd_total = sum(elpd_loo, na.rm = TRUE),
    se_total = sqrt(sum(se_elpd_loo^2, na.rm = TRUE)),
    looic_total = sum(looic, na.rm = TRUE),
    n_components = n(),
    .groups = "drop"
  ) %>%
  mutate(
    best_elpd_total = max(elpd_total, na.rm = TRUE),
    delta_elpd_total = elpd_total - best_elpd_total,
    rank_total = rank(-elpd_total, ties.method = "min"),
    relative_likelihood = exp(elpd_total - max(elpd_total, na.rm = TRUE)),
    approximate_weight = relative_likelihood / sum(relative_likelihood, na.rm = TRUE),
    within_2se_of_best = abs(delta_elpd_total) <= 2 * se_total
  ) %>%
  arrange(rank_total)

write.csv(
  loo_total,
  file.path(OUT_ROOT, "03_model_comparison", "scenario_comparison_total_LOO_ELPD.csv"),
  row.names = FALSE
)

print(loo_total)

best_scenario <- loo_total$scenario[which.max(loo_total$elpd_total)]

representative_scenario <- best_scenario
if ("M4_partial_mediation" %in% loo_total$scenario[loo_total$within_2se_of_best]) {
  representative_scenario <- "M4_partial_mediation"
} else if ("M3_full_mediation" %in% loo_total$scenario[loo_total$within_2se_of_best]) {
  representative_scenario <- "M3_full_mediation"
} else if ("M2_selfing_syndrome" %in% loo_total$scenario[loo_total$within_2se_of_best]) {
  representative_scenario <- "M2_selfing_syndrome"
}

message("Best scenario by total ELPD: ", best_scenario)
message("Representative scenario for path diagram: ", representative_scenario)

############################################################
## 11. Direct and indirect effects
############################################################

effect_rows <- list()

if ("M3_full_mediation" %in% names(models)) {
  b_distance_mismatch <- posterior_beta(models$M3_full_mediation$mismatch, "z_logdist")
  b_mismatch_self <- posterior_beta(models$M3_full_mediation$self, "z_mismatch")
  b_self_color <- posterior_beta(models$M3_full_mediation$color, "prop_self_compatible")
  b_self_form <- posterior_beta(models$M3_full_mediation$form, "prop_self_compatible")
  
  effect_rows[["M3_mismatch_self"]] <- summ_draws(
    b_mismatch_self,
    "M3: bumblebee-environment mismatch -> self-compatible flora"
  )
  
  effect_rows[["M3_self_color"]] <- summ_draws(
    b_self_color,
    "M3: self-compatible flora -> plain-colored flowers"
  )
  
  effect_rows[["M3_self_form"]] <- summ_draws(
    b_self_form,
    "M3: self-compatible flora -> generalized floral forms"
  )
  
  effect_rows[["M3_indirect_mismatch_self_color"]] <- summ_draws(
    b_mismatch_self * b_self_color,
    "M3 indirect: mismatch -> self-compatible flora -> plain-colored flowers"
  )
  
  effect_rows[["M3_indirect_mismatch_self_form"]] <- summ_draws(
    b_mismatch_self * b_self_form,
    "M3 indirect: mismatch -> self-compatible flora -> generalized floral forms"
  )
  
  effect_rows[["M3_chain_distance_color"]] <- summ_draws(
    b_distance_mismatch * b_mismatch_self * b_self_color,
    "M3 chain: distance -> mismatch -> self-compatible flora -> plain-colored flowers"
  )
  
  effect_rows[["M3_chain_distance_form"]] <- summ_draws(
    b_distance_mismatch * b_mismatch_self * b_self_form,
    "M3 chain: distance -> mismatch -> self-compatible flora -> generalized floral forms"
  )
}

if ("M4_partial_mediation" %in% names(models)) {
  b_distance_mismatch <- posterior_beta(models$M4_partial_mediation$mismatch, "z_logdist")
  b_mismatch_self <- posterior_beta(models$M4_partial_mediation$self, "z_mismatch")
  b_mismatch_color <- posterior_beta(models$M4_partial_mediation$color, "z_mismatch")
  b_mismatch_form <- posterior_beta(models$M4_partial_mediation$form, "z_mismatch")
  b_self_color <- posterior_beta(models$M4_partial_mediation$color, "prop_self_compatible")
  b_self_form <- posterior_beta(models$M4_partial_mediation$form, "prop_self_compatible")
  
  effect_rows[["M4_mismatch_self"]] <- summ_draws(
    b_mismatch_self,
    "M4: bumblebee-environment mismatch -> self-compatible flora"
  )
  
  effect_rows[["M4_direct_mismatch_color"]] <- summ_draws(
    b_mismatch_color,
    "M4 direct: bumblebee-environment mismatch -> plain-colored flowers"
  )
  
  effect_rows[["M4_direct_mismatch_form"]] <- summ_draws(
    b_mismatch_form,
    "M4 direct: bumblebee-environment mismatch -> generalized floral forms"
  )
  
  effect_rows[["M4_self_color"]] <- summ_draws(
    b_self_color,
    "M4: self-compatible flora -> plain-colored flowers"
  )
  
  effect_rows[["M4_self_form"]] <- summ_draws(
    b_self_form,
    "M4: self-compatible flora -> generalized floral forms"
  )
  
  effect_rows[["M4_indirect_mismatch_self_color"]] <- summ_draws(
    b_mismatch_self * b_self_color,
    "M4 indirect: mismatch -> self-compatible flora -> plain-colored flowers"
  )
  
  effect_rows[["M4_indirect_mismatch_self_form"]] <- summ_draws(
    b_mismatch_self * b_self_form,
    "M4 indirect: mismatch -> self-compatible flora -> generalized floral forms"
  )
  
  effect_rows[["M4_total_mismatch_color"]] <- summ_draws(
    b_mismatch_color + b_mismatch_self * b_self_color,
    "M4 total: mismatch -> plain-colored flowers"
  )
  
  effect_rows[["M4_total_mismatch_form"]] <- summ_draws(
    b_mismatch_form + b_mismatch_self * b_self_form,
    "M4 total: mismatch -> generalized floral forms"
  )
  
  effect_rows[["M4_chain_distance_color"]] <- summ_draws(
    b_distance_mismatch * b_mismatch_self * b_self_color,
    "M4 chain: distance -> mismatch -> self-compatible flora -> plain-colored flowers"
  )
  
  effect_rows[["M4_chain_distance_form"]] <- summ_draws(
    b_distance_mismatch * b_mismatch_self * b_self_form,
    "M4 chain: distance -> mismatch -> self-compatible flora -> generalized floral forms"
  )
}

path_effects_df <- bind_rows(effect_rows)

write.csv(
  path_effects_df,
  file.path(OUT_ROOT, "04_path_effects", "direct_indirect_and_total_path_effects.csv"),
  row.names = FALSE
)

print(path_effects_df)

############################################################
## 12. Manuscript figures
############################################################

p_total <- loo_total %>%
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = scenario_label[order(delta_elpd_total)]
    )
  ) %>%
  ggplot(aes(x = scenario_label, y = delta_elpd_total)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_col(width = 0.65) +
  coord_flip() +
  labs(
    title = "Bayesian comparison of island floral syndrome scenarios",
    x = NULL,
    y = "Delta ELPD from best scenario"
  ) +
  theme_manuscript()

ggsave(
  file.path(OUT_ROOT, "05_figures", "Figure_1_scenario_comparison_total_ELPD.png"),
  p_total,
  width = 7,
  height = 4,
  dpi = 300
)

ggsave(
  file.path(OUT_ROOT, "05_figures", "Figure_1_scenario_comparison_total_ELPD.pdf"),
  p_total,
  width = 7,
  height = 4
)

p_response <- loo_by_component %>%
  mutate(
    scenario_label = scenario_labels[scenario],
    response_label = response_labels[response]
  ) %>%
  ggplot(aes(x = scenario_label, y = delta_elpd_response)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(width = 0.65) +
  facet_wrap(~ response_label, scales = "free_y") +
  coord_flip() +
  labs(
    title = "Response-specific predictive performance",
    x = NULL,
    y = "Delta ELPD within response"
  ) +
  theme_manuscript(10)

ggsave(
  file.path(OUT_ROOT, "05_figures", "Figure_2_response_specific_ELPD.png"),
  p_response,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(OUT_ROOT, "05_figures", "Figure_2_response_specific_ELPD.pdf"),
  p_response,
  width = 8,
  height = 5
)

if (nrow(path_effects_df) > 0) {
  p_effects <- path_effects_df %>%
    mutate(
      effect_wrapped = stringr::str_wrap(effect, 55),
      effect_wrapped = factor(effect_wrapped, levels = rev(effect_wrapped[order(mean)]))
    ) %>%
    ggplot(aes(x = mean, y = effect_wrapped)) +
    geom_vline(xintercept = 0, linewidth = 0.35) +
    geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.18) +
    geom_point(size = 2.2) +
    labs(
      title = "Posterior estimates of direct, indirect, and total effects",
      x = "Posterior mean with 95% credible interval",
      y = NULL
    ) +
    theme_manuscript(9)
  
  ggsave(
    file.path(OUT_ROOT, "05_figures", "Figure_3_path_effect_estimates.png"),
    p_effects,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  ggsave(
    file.path(OUT_ROOT, "05_figures", "Figure_3_path_effect_estimates.pdf"),
    p_effects,
    width = 9,
    height = 7
  )
}

############################################################
## 13. Significant-path diagram
############################################################

make_edges_from_scenario <- function(scenario_name) {
  fixed_effects_df %>%
    filter(
      scenario == scenario_name,
      significant,
      term != "Intercept"
    ) %>%
    mutate(
      from = case_when(
        term == "z_logdist" ~ "Distance to continent",
        term == "z_logarea" ~ "Island area",
        term == "z_mismatch" ~ "Bumblebee-environment mismatch",
        term == "prop_self_compatible" ~ "Self-compatible flora",
        TRUE ~ term
      ),
      to = recode(response, !!!response_labels),
      edge_color = ifelse(mean >= 0, "#2166AC", "#B2182B"),
      edge_width = ifelse(
        length(unique(abs(mean))) <= 1,
        3,
        scales::rescale(abs(mean), to = c(1.5, 5.0))
      ),
      label = sprintf("%.2f", mean)
    ) %>%
    filter(from != to) %>%
    select(from, to, mean, q025, q975, edge_color, edge_width, label)
}

edges <- make_edges_from_scenario(representative_scenario)

write.csv(
  edges,
  file.path(
    OUT_ROOT,
    "04_path_effects",
    paste0("significant_edges_", representative_scenario, ".csv")
  ),
  row.names = FALSE
)

nodes <- c(
  "Distance to continent",
  "Island area",
  "Bumblebee-environment mismatch",
  "Self-compatible flora",
  "Plain-colored flowers",
  "Generalized floral forms"
)

node_defs <- paste0(
  '  "', nodes, '" [label="', nodes,
  '", shape=box, style="rounded,filled", fillcolor="white", ',
  'color="gray35", fontname="Helvetica"];',
  collapse = "\n"
)

edge_defs <- if (nrow(edges) == 0) {
  ""
} else {
  paste0(
    '  "', edges$from, '" -> "', edges$to,
    '" [label="', edges$label,
    '", color="', edges$edge_color,
    '", fontcolor="', edges$edge_color,
    '", penwidth=', sprintf("%.2f", edges$edge_width),
    ', arrowsize=0.8, fontname="Helvetica"];',
    collapse = "\n"
  )
}

graph_code <- paste0(
  'digraph G {\n',
  '  graph [rankdir=LR, bgcolor="white", layout=dot, splines=true, overlap=false, nodesep=0.7, ranksep=0.9];\n',
  '  node [fontsize=16, margin=0.16];\n',
  '  edge [fontsize=13];\n',
  node_defs, '\n',
  '  {rank=same; "Distance to continent"; "Island area";}\n',
  '  {rank=same; "Plain-colored flowers"; "Generalized floral forms";}\n',
  edge_defs, '\n',
  '}\n'
)

path_widget <- DiagrammeR::grViz(graph_code)
path_svg <- DiagrammeRsvg::export_svg(path_widget)

path_svg_file <- file.path(
  OUT_ROOT,
  "05_figures",
  paste0("Figure_4_significant_path_diagram_", representative_scenario, ".svg")
)

path_png_file <- file.path(
  OUT_ROOT,
  "05_figures",
  paste0("Figure_4_significant_path_diagram_", representative_scenario, ".png")
)

path_pdf_file <- file.path(
  OUT_ROOT,
  "05_figures",
  paste0("Figure_4_significant_path_diagram_", representative_scenario, ".pdf")
)

writeLines(path_svg, path_svg_file)
rsvg::rsvg_png(charToRaw(path_svg), path_png_file, width = 2200, height = 1200)
rsvg::rsvg_pdf(charToRaw(path_svg), path_pdf_file, width = 11, height = 6)

############################################################
## 14. External projection test
############################################################

fit_glm_component <- function(scenario_name, response_name, train_df) {
  train_df <- as.data.frame(train_df)
  
  if (response_name == "mismatch") {
    return(lm(z_mismatch ~ z_logdist, data = train_df))
  }
  
  if (response_name == "self") {
    if (scenario_name %in% c("M3_full_mediation", "M4_partial_mediation")) {
      return(
        glm(
          cbind(n_self_compatible, n_self_total - n_self_compatible) ~ z_mismatch,
          data = train_df,
          family = binomial()
        )
      )
    } else {
      return(
        glm(
          cbind(n_self_compatible, n_self_total - n_self_compatible) ~ z_logdist,
          data = train_df,
          family = binomial()
        )
      )
    }
  }
  
  if (response_name == "color") {
    if (scenario_name == "M0_geographic_filter") {
      return(
        glm(
          cbind(n_plain, n_color_total - n_plain) ~ z_logdist + z_logarea,
          data = train_df,
          family = binomial()
        )
      )
    } else if (scenario_name %in% c("M2_selfing_syndrome", "M3_full_mediation")) {
      return(
        glm(
          cbind(n_plain, n_color_total - n_plain) ~ prop_self_compatible + z_logarea,
          data = train_df,
          family = binomial()
        )
      )
    } else {
      return(
        glm(
          cbind(n_plain, n_color_total - n_plain) ~ z_mismatch + prop_self_compatible + z_logarea,
          data = train_df,
          family = binomial()
        )
      )
    }
  }
  
  if (response_name == "form") {
    if (scenario_name == "M0_geographic_filter") {
      return(
        glm(
          cbind(n_generalized, n_form_total - n_generalized) ~ z_logdist + z_logarea,
          data = train_df,
          family = binomial()
        )
      )
    } else if (scenario_name %in% c("M2_selfing_syndrome", "M3_full_mediation")) {
      return(
        glm(
          cbind(n_generalized, n_form_total - n_generalized) ~ prop_self_compatible + z_logarea,
          data = train_df,
          family = binomial()
        )
      )
    } else {
      return(
        glm(
          cbind(n_generalized, n_form_total - n_generalized) ~ z_mismatch + prop_self_compatible + z_logarea,
          data = train_df,
          family = binomial()
        )
      )
    }
  }
  
  stop("Unknown response: ", response_name)
}

nll_binom <- function(y, n, p) {
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  -sum(dbinom(y, size = n, prob = p, log = TRUE), na.rm = TRUE)
}

nll_gauss <- function(y, mu, sigma) {
  sigma <- max(sigma, 1e-8)
  -sum(dnorm(y, mean = mu, sd = sigma, log = TRUE), na.rm = TRUE)
}

external_eval <- function(train_df, test_df, scenario_name, test_label) {
  train_df <- as.data.frame(train_df)
  test_df <- as.data.frame(test_df)
  
  if (nrow(test_df) < 5) {
    return(NULL)
  }
  
  fit_mismatch <- fit_glm_component(scenario_name, "mismatch", train_df)
  fit_self <- fit_glm_component(scenario_name, "self", train_df)
  fit_color <- fit_glm_component(scenario_name, "color", train_df)
  fit_form <- fit_glm_component(scenario_name, "form", train_df)
  
  pred_mismatch <- as.numeric(predict(fit_mismatch, newdata = test_df, type = "response"))
  pred_self <- as.numeric(predict(fit_self, newdata = test_df, type = "response"))
  pred_color <- as.numeric(predict(fit_color, newdata = test_df, type = "response"))
  pred_form <- as.numeric(predict(fit_form, newdata = test_df, type = "response"))
  
  sigma_mismatch <- sd(resid(fit_mismatch), na.rm = TRUE)
  if (!is.finite(sigma_mismatch) || sigma_mismatch <= 0) {
    sigma_mismatch <- 1e-6
  }
  
  nll_mismatch <- nll_gauss(test_df$z_mismatch, pred_mismatch, sigma_mismatch)
  nll_self <- nll_binom(test_df$n_self_compatible, test_df$n_self_total, pred_self)
  nll_color <- nll_binom(test_df$n_plain, test_df$n_color_total, pred_color)
  nll_form <- nll_binom(test_df$n_generalized, test_df$n_form_total, pred_form)
  
  tibble(
    train_region = "temperate_NH",
    test_region = test_label,
    scenario = scenario_name,
    scenario_label = scenario_labels[[scenario_name]],
    n_islands = nrow(test_df),
    nll_mismatch = nll_mismatch,
    nll_self = nll_self,
    nll_color = nll_color,
    nll_form = nll_form,
    nll_traits = nll_color + nll_form,
    nll_full = nll_mismatch + nll_self + nll_color + nll_form,
    nll_traits_per_island = nll_traits / nrow(test_df),
    nll_full_per_island = nll_full / nrow(test_df)
  )
}

model_ready_all <- island_all %>%
  filter(
    is.finite(z_logdist),
    is.finite(z_logarea),
    is.finite(z_mismatch),
    n_color_total > 0,
    n_form_total > 0,
    n_self_total > 0
  )

train_tnh <- model_ready_all %>%
  filter(analysis_group == "temperate_NH")

external_tropical <- model_ready_all %>%
  filter(analysis_group == "tropical_all")

external_southern <- model_ready_all %>%
  filter(analysis_group == "southern_non_tropical")

external_projection <- bind_rows(
  external_eval(train_tnh, train_tnh, representative_scenario, "temperate_NH_in_sample"),
  external_eval(train_tnh, external_tropical, representative_scenario, "tropical_external"),
  external_eval(train_tnh, external_southern, representative_scenario, "southern_non_tropical_external")
)

write.csv(
  external_projection,
  file.path(
    OUT_ROOT,
    "06_external_projection",
    "external_projection_negative_log_likelihood.csv"
  ),
  row.names = FALSE
)

print(external_projection)

if (!is.null(external_projection) && nrow(external_projection) > 0) {
  p_external <- external_projection %>%
    ggplot(aes(x = test_region, y = nll_traits_per_island)) +
    geom_col(width = 0.65) +
    coord_flip() +
    labs(
      title = "External projection of the northern-temperate path scenario",
      x = NULL,
      y = "Trait prediction negative log-likelihood per island"
    ) +
    theme_manuscript()
  
  ggsave(
    file.path(OUT_ROOT, "05_figures", "Figure_S1_external_projection_NLL.png"),
    p_external,
    width = 7,
    height = 4,
    dpi = 300
  )
  
  ggsave(
    file.path(OUT_ROOT, "05_figures", "Figure_S1_external_projection_NLL.pdf"),
    p_external,
    width = 7,
    height = 4
  )
}

############################################################
## 15. Save metadata and README
############################################################

metadata_file <- file.path(
  OUT_ROOT,
  "07_metadata",
  "analysis_metadata.txt"
)

sink(metadata_file)
cat("Northern-temperate island floral syndrome Bayesian path-model analysis\n")
cat("=====================================================================\n\n")
cat("Run time: ", as.character(Sys.time()), "\n", sep = "")
cat("Input file: ", INPUT_FILE, "\n", sep = "")
cat("Output directory: ", OUT_ROOT, "\n", sep = "")
cat("Backend: ", backend_use, "\n", sep = "")
cat("Chains: ", RUN_CHAINS, "\n", sep = "")
cat("Cores: ", RUN_CORES, "\n", sep = "")
cat("Iterations: ", RUN_ITER, "\n", sep = "")
cat("Warmup: ", RUN_WARMUP, "\n", sep = "")
cat("Seed: ", RUN_SEED, "\n\n", sep = "")
cat("Best scenario by total ELPD: ", best_scenario, "\n", sep = "")
cat("Representative scenario: ", representative_scenario, "\n\n", sep = "")
cat("Scenario definitions:\n")
cat("M0_geographic_filter: geographic filter\n")
cat("M2_selfing_syndrome: selfing syndrome\n")
cat("M3_full_mediation: full mediation through mismatch and self-compatibility\n")
cat("M4_partial_mediation: M3 plus direct mismatch-to-trait effects\n\n")
cat("Session information:\n")
print(sessionInfo())
sink()

readme_file <- file.path(
  OUT_ROOT,
  "README_NorthernTemperate_IslandFloralSyndrome_BayesianPathModelAnalysis.txt"
)

readme <- c(
  "Northern-temperate island floral syndrome Bayesian path-model analysis",
  "",
  "This is a standalone manuscript reproducibility run.",
  "The script reads the original CSV directly and refits all Bayesian component models from scratch.",
  "",
  "Main outputs:",
  "01_data/island_level_northern_temperate_model_ready.csv",
  "03_model_comparison/scenario_comparison_total_LOO_ELPD.csv",
  "03_model_comparison/response_specific_LOO_ELPD_by_scenario.csv",
  "04_path_effects/direct_indirect_and_total_path_effects.csv",
  "04_path_effects/significant_edges_<representative_scenario>.csv",
  "05_figures/Figure_1_scenario_comparison_total_ELPD.png/pdf",
  "05_figures/Figure_2_response_specific_ELPD.png/pdf",
  "05_figures/Figure_3_path_effect_estimates.png/pdf",
  "05_figures/Figure_4_significant_path_diagram_<representative_scenario>.png/pdf/svg",
  "06_external_projection/external_projection_negative_log_likelihood.csv",
  "07_metadata/analysis_metadata.txt",
  "",
  "Interpretation note:",
  "LOO-ELPD evaluates predictive consistency among candidate path scenarios.",
  "It does not by itself prove causality.",
  "Biological interpretation should combine scenario ranking, uncertainty, and posterior path effects."
)

writeLines(readme, readme_file)

############################################################
## 16. Final message
############################################################

message("\nDONE.")
message("Results saved to: ", OUT_ROOT)
message("Best scenario by total ELPD: ", best_scenario)
message("Representative scenario: ", representative_scenario)
message("Key outputs:")
message(" - 03_model_comparison/scenario_comparison_total_LOO_ELPD.csv")
message(" - 03_model_comparison/response_specific_LOO_ELPD_by_scenario.csv")
message(" - 04_path_effects/direct_indirect_and_total_path_effects.csv")
message(" - 05_figures/")
message(" - 06_external_projection/external_projection_negative_log_likelihood.csv")
message(" - 07_metadata/analysis_metadata.txt")
