############################################################
## Northern_island.R
##
## Purpose
##   Northern-temperate islands only:
##   Compare island-syndrome causal scenarios using Bayesian
##   piecewise path models, but keep computation light.
##
## Core design
##   Unit: island-level flora composition
##   Responses:
##     1) bumblebee-environment mismatch index (continuous)
##     2) self-compatible flora proportion (binomial counts)
##     3) plain-colored flowers proportion (binomial counts)
##     4) generalized floral forms proportion (binomial counts)
##
## Scenarios
##   M1: geographic filter
##       distance + area -> floral traits
##
##   M2: selfing syndrome
##       distance -> self-compatible -> floral traits
##       area -> floral traits
##
##   M3: full mediation
##       distance -> bumblebee mismatch -> self-compatible -> floral traits
##       area -> floral traits
##
##   M4: partial mediation
##       distance -> bumblebee mismatch -> self-compatible -> floral traits
##       bumblebee mismatch -> floral traits
##       area -> floral traits
##
## Notes
##   - Island area is NOT used as a cause of mismatch or selfing.
##     It only enters the final floral trait models.
##   - "PDI" is treated internally as mean_PDI but output labels use
##     "bumblebee-environment mismatch".
##   - This script keeps computation lighter by:
##       * fitting separate brms component models (piecewise)
##       * using island-level data, not species x island rows
##       * comparing scenarios using the SUM of component model LOO values
##       * using modest MCMC defaults
##       * optionally using cmdstanr if available
############################################################

rm(list = ls())
gc()

## =========================
## User settings
## =========================
INPUT_FILE <- "/Users/rachelzhang/Desktop/list_with_traits_island_distance_bom_with_PDI_reclassified_isolated_animal_trop_temp.csv"
OUT_ROOT <- "/Users/rachelzhang/Desktop/northern_island"

## Light but usable MCMC settings. Increase for final paper if needed.
RUN_CHAINS <- 2
RUN_CORES  <- 1
RUN_ITER   <- 1200
RUN_WARMUP <- 600
RUN_SEED   <- 20260522

## Reproducibility/cache behavior
## FALSE = if an .rds model exists, load it; if not, fit from scratch and save it.
## TRUE  = ignore existing .rds files and refit all component models.
FORCE_REFIT_MODELS <- FALSE

## FALSE = if a *_loo.rds exists, load it; if not, compute and save it.
## TRUE  = recompute all LOO objects from fitted models.
FORCE_RECOMPUTE_LOO <- FALSE

## Optional: sample islands for quick debugging. For production keep FALSE.
DEBUG_SUBSET <- FALSE
DEBUG_N_ISLANDS <- 200

## If TRUE, fit with cmdstanr backend when installed and configured.
USE_CMDSTANR_IF_AVAILABLE <- TRUE

## LOO settings
LOO_CORES <- 1

## Significance threshold for path diagram: 95% credible interval excludes zero.
SIG_Q_LOW <- 0.025
SIG_Q_HIGH <- 0.975

## =========================
## Packages
## =========================
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

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "01_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "02_bayes_component_models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "03_scenario_comparison"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "04_effects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "05_figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "06_supp_external_projection"), recursive = TRUE, showWarnings = FALSE)

options(mc.cores = RUN_CORES)
set.seed(RUN_SEED)

backend_use <- "rstan"
if (USE_CMDSTANR_IF_AVAILABLE && requireNamespace("cmdstanr", quietly = TRUE)) {
  ok_cmdstan <- FALSE
  try({ ok_cmdstan <- !is.null(cmdstanr::cmdstan_path()) && nzchar(cmdstanr::cmdstan_path()) }, silent = TRUE)
  if (isTRUE(ok_cmdstan)) backend_use <- "cmdstanr"
}
message("[Backend] ", backend_use)
message("[Input] ", INPUT_FILE)
message("[Output] ", OUT_ROOT)

## =========================
## Helper functions
## =========================
clean_chr <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_to_lower(x)
  x
}

zscale <- function(x) as.numeric(scale(x))

safe_div <- function(a, b) ifelse(is.na(b) | b <= 0, NA_real_, a / b)

pick_first_existing <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

make_brms <- function(formula, data, family, name) {
  message("
--- Fitting ", name, " | n = ", nrow(data), " ---")
  fit_file <- file.path(OUT_ROOT, "02_bayes_component_models", paste0(name, ".rds"))

  if (file.exists(fit_file) && !isTRUE(FORCE_REFIT_MODELS)) {
    message("Loading cached model: ", fit_file)
    cached <- tryCatch(readRDS(fit_file), error = function(e) e)
    if (!inherits(cached, "error")) {
      return(cached)
    }
    message("Cached model could not be read; refitting: ", conditionMessage(cached))
  }

  message("No usable cached .rds found, so fitting from scratch and saving to: ", fit_file)
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
  saveRDS(fit, fit_file)
  fit
}

loo_one <- function(fit, name) {
  message("Computing LOO: ", name)
  loo_file <- file.path(OUT_ROOT, "02_bayes_component_models", paste0(name, "_loo.rds"))

  if (file.exists(loo_file) && !isTRUE(FORCE_RECOMPUTE_LOO)) {
    cached <- tryCatch(readRDS(loo_file), error = function(e) e)
    if (!inherits(cached, "error")) {
      return(cached)
    }
    message("Cached LOO could not be read; recomputing: ", conditionMessage(cached))
  }

  out <- loo::loo(fit, cores = LOO_CORES, moment_match = FALSE)
  saveRDS(out, loo_file)
  out
}

extract_fixef <- function(fit, model_name, response_name) {
  fe <- as.data.frame(brms::fixef(fit, probs = c(SIG_Q_LOW, SIG_Q_HIGH)))
  fe$term <- rownames(fe)
  rownames(fe) <- NULL
  fe %>%
    transmute(
      model = model_name,
      response = response_name,
      term = term,
      mean = Estimate,
      se = Est.Error,
      q025 = .data[[paste0("Q", SIG_Q_LOW * 100)]],
      q975 = .data[[paste0("Q", SIG_Q_HIGH * 100)]],
      significant = q025 > 0 | q975 < 0
    )
}

posterior_beta <- function(fit, term) {
  dr <- as_draws_df(fit)
  nm <- paste0("b_", term)
  if (!nm %in% names(dr)) return(rep(NA_real_, nrow(dr)))
  as.numeric(dr[[nm]])
}

summ_draws <- function(x, label) {
  tibble(
    effect = label,
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q025 = quantile(x, SIG_Q_LOW, na.rm = TRUE),
    q975 = quantile(x, SIG_Q_HIGH, na.rm = TRUE),
    pd = max(mean(x > 0, na.rm = TRUE), mean(x < 0, na.rm = TRUE)),
    significant = q025 > 0 | q975 < 0
  )
}

## Extract y-level log-likelihood for binomial brms model works through loo.
## For external projection, use simple fixed-effect posterior means from best scenario.
invlogit <- function(x) 1 / (1 + exp(-x))

## =========================
## Load and clean data
## =========================
stopifnot(file.exists(INPUT_FILE))
raw <- data.table::fread(INPUT_FILE, encoding = "UTF-8", data.table = FALSE)
message("Rows loaded: ", nrow(raw))

nms <- names(raw)
col_lat  <- pick_first_existing(nms, c("centroid_lat", "latitude", "lat", "decimalLatitude"))
col_lon  <- pick_first_existing(nms, c("centroid_lon", "longitude", "lon", "decimalLongitude"))
col_dist <- pick_first_existing(nms, c("dist_to_continent_km", "distance_to_mainland", "dist_to_mainland", "min_distance_to_continent", "NEAR_DIST", "distance_km"))
col_area <- pick_first_existing(nms, c("area_km2", "IslandArea", "Area_Geode", "Shape_Area"))
col_pdi  <- pick_first_existing(nms, c("mean_PDI", "PDI", "mean_pdi", "bumblebee_mismatch"))
col_clim <- pick_first_existing(nms, c("climate_zone", "climate", "zone"))
col_island <- pick_first_existing(nms, c("island_id", "USGS_ISID", "Name_USGSO", "NAME_wcmcI", "NAME_LOCAL"))
col_species <- pick_first_existing(nms, c("species", "i.species", "taxon", "scientificName"))

required <- c(col_lat, col_dist, col_area, col_pdi, col_clim, col_island, col_species)
if (any(is.na(required))) {
  stop("Missing required columns. Detected: lat=", col_lat, ", dist=", col_dist,
       ", area=", col_area, ", pdi=", col_pdi, ", climate=", col_clim,
       ", island=", col_island, ", species=", col_species)
}

## Trait columns from the user's current dataset
if (!"flower_color_simple" %in% nms) stop("flower_color_simple column not found")
if (!"flower_shape_simple" %in% nms) stop("flower_shape_simple column not found")
if (!"self_incompatibility" %in% nms) stop("self_incompatibility column not found")

## Map traits using previous validated rules
plain_levels <- c("white", "green_brown_inconspicuous", "other")
consp_levels <- c("blue_purple", "red_pink", "yellow_orange")
generalized_levels <- c("actinomorphic_open", "brush_puff", "composite_head", "reduced_wind")
specialized_levels <- c("tubular", "zygomorphic")
self_yes <- c("sc", "likely_sc")
self_no  <- c("si", "likely_si", "obligate_si")

 dat <- raw %>%
  mutate(
    island_id_use = as.character(.data[[col_island]]),
    species_use = as.character(.data[[col_species]]),
    lat = as.numeric(.data[[col_lat]]),
    lon = if (!is.na(col_lon)) as.numeric(.data[[col_lon]]) else NA_real_,
    dist_km = as.numeric(.data[[col_dist]]),
    area_km2 = as.numeric(.data[[col_area]]),
    bumblebee_mismatch = as.numeric(.data[[col_pdi]]),
    climate_zone = clean_chr(.data[[col_clim]]),
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
      flower_color_simple_clean %in% plain_levels ~ 1L,
      flower_color_simple_clean %in% consp_levels ~ 0L,
      TRUE ~ NA_integer_
    ),
    form_generalized = case_when(
      flower_shape_simple_clean %in% generalized_levels ~ 1L,
      flower_shape_simple_clean %in% specialized_levels ~ 0L,
      TRUE ~ NA_integer_
    ),
    self_compatible = case_when(
      self_incompatibility_clean %in% self_yes ~ 1L,
      self_incompatibility_clean %in% self_no ~ 0L,
      TRUE ~ NA_integer_
    )
  )

## Diagnostic tables
write.csv(dat %>% count(analysis_group, climate_zone, hemisphere),
          file.path(OUT_ROOT, "01_data", "diagnostic_group_counts_raw.csv"), row.names = FALSE)
write.csv(dat %>% count(flower_color_simple_clean, color_plain),
          file.path(OUT_ROOT, "01_data", "diagnostic_color_mapping.csv"), row.names = FALSE)
write.csv(dat %>% count(flower_shape_simple_clean, form_generalized),
          file.path(OUT_ROOT, "01_data", "diagnostic_form_mapping.csv"), row.names = FALSE)
write.csv(dat %>% count(self_incompatibility_clean, self_compatible),
          file.path(OUT_ROOT, "01_data", "diagnostic_self_mapping.csv"), row.names = FALSE)

## =========================
## Island-level aggregation
## =========================
make_island_counts <- function(x) {
  x %>%
    filter(!is.na(analysis_group), !is.na(island_id_use), is.finite(dist_km), is.finite(area_km2), is.finite(lat)) %>%
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

island_all <- make_island_counts(dat)
write.csv(island_all, file.path(OUT_ROOT, "01_data", "island_level_all_regions_counts.csv"), row.names = FALSE)

island_tnh <- island_all %>%
  filter(analysis_group == "temperate_NH") %>%
  filter(
    is.finite(z_logdist), is.finite(z_logarea), is.finite(z_mismatch),
    n_color_total > 0, n_form_total > 0, n_self_total > 0
  )

if (DEBUG_SUBSET) {
  island_tnh <- island_tnh %>% slice_sample(n = min(DEBUG_N_ISLANDS, nrow(island_tnh)))
}

message("TNH island-level rows: ", nrow(island_tnh))
write.csv(island_tnh, file.path(OUT_ROOT, "01_data", "island_level_temperate_NH_counts.csv"), row.names = FALSE)

if (nrow(island_tnh) < 30) stop("Too few TNH islands after aggregation.")

## =========================
## Model formulas
## =========================
## Component models use common responses where meaningful.
## For fair scenario comparison, each scenario contributes the same response components:
##   mismatch, self, color, form.
## Scenario structure differs in predictors.

## M0: geography only
scenario_formulas <- list(
  M0 = list(
    mismatch = bf(z_mismatch ~ 1 + z_logdist),  ## distance describes mismatch; included for same response set
    self     = bf(n_self_compatible | trials(n_self_total) ~ 1 + z_logdist),
    color    = bf(n_plain | trials(n_color_total) ~ 1 + z_logdist + z_logarea),
    form     = bf(n_generalized | trials(n_form_total) ~ 1 + z_logdist + z_logarea)
  ),
  M2 = list(
    mismatch = bf(z_mismatch ~ 1 + z_logdist),
    self     = bf(n_self_compatible | trials(n_self_total) ~ 1 + z_logdist),
    color    = bf(n_plain | trials(n_color_total) ~ 1 + prop_self_compatible + z_logarea),
    form     = bf(n_generalized | trials(n_form_total) ~ 1 + prop_self_compatible + z_logarea)
  ),
  M3 = list(
    mismatch = bf(z_mismatch ~ 1 + z_logdist),
    self     = bf(n_self_compatible | trials(n_self_total) ~ 1 + z_mismatch),
    color    = bf(n_plain | trials(n_color_total) ~ 1 + prop_self_compatible + z_logarea),
    form     = bf(n_generalized | trials(n_form_total) ~ 1 + prop_self_compatible + z_logarea)
  ),
  M4 = list(
    mismatch = bf(z_mismatch ~ 1 + z_logdist),
    self     = bf(n_self_compatible | trials(n_self_total) ~ 1 + z_mismatch),
    color    = bf(n_plain | trials(n_color_total) ~ 1 + z_mismatch + prop_self_compatible + z_logarea),
    form     = bf(n_generalized | trials(n_form_total) ~ 1 + z_mismatch + prop_self_compatible + z_logarea)
  )
)

families <- list(
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

## =========================
## Fit all Bayesian piecewise scenario models
## =========================
models <- list()
loos <- list()
fixefs_all <- list()

for (sc in names(scenario_formulas)) {
  models[[sc]] <- list()
  loos[[sc]] <- list()
  for (resp in names(scenario_formulas[[sc]])) {
    nm <- paste(sc, resp, sep = "_")
    fit <- make_brms(
      formula = scenario_formulas[[sc]][[resp]],
      data = island_tnh,
      family = families[[resp]],
      name = nm
    )
    models[[sc]][[resp]] <- fit
    loos[[sc]][[resp]] <- loo_one(fit, nm)
    fixefs_all[[nm]] <- extract_fixef(fit, sc, resp)
  }
}

fixefs_df <- bind_rows(fixefs_all)
write.csv(fixefs_df, file.path(OUT_ROOT, "03_scenario_comparison", "TNH_all_scenario_component_fixed_effects.csv"), row.names = FALSE)

## =========================
## Scenario comparison: total and by response
## =========================
loo_component <- map_dfr(names(loos), function(sc) {
  map_dfr(names(loos[[sc]]), function(resp) {
    x <- loos[[sc]][[resp]]
    tibble(
      scenario = sc,
      response = resp,
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

write.csv(loo_component, file.path(OUT_ROOT, "03_scenario_comparison", "TNH_LOO_components_by_response.csv"), row.names = FALSE)

loo_total <- loo_component %>%
  group_by(scenario) %>%
  summarise(
    elpd_total = sum(elpd_loo, na.rm = TRUE),
    ## Approx SE of sum assuming independence among component elpd estimates.
    se_total = sqrt(sum(se_elpd_loo^2, na.rm = TRUE)),
    looic_total = sum(looic, na.rm = TRUE),
    n_components = n(),
    .groups = "drop"
  ) %>%
  mutate(
    best_elpd_total = max(elpd_total, na.rm = TRUE),
    delta_elpd_total = elpd_total - best_elpd_total,
    rank_total = rank(-elpd_total, ties.method = "min")
  ) %>%
  arrange(rank_total)

write.csv(loo_total, file.path(OUT_ROOT, "03_scenario_comparison", "TNH_M0_M2_M3_M4_TOTAL_LOO.csv"), row.names = FALSE)
print(loo_total)

## Use model weights from total elpd as an approximate softmax weight.
## This is not formal stacking across piecewise bundles, but useful for visual summary.
loo_total <- loo_total %>%
  mutate(
    rel_weight = exp(elpd_total - max(elpd_total)),
    approx_weight = rel_weight / sum(rel_weight)
  )
write.csv(loo_total, file.path(OUT_ROOT, "03_scenario_comparison", "TNH_M0_M2_M3_M4_TOTAL_LOO_with_approx_weights.csv"), row.names = FALSE)

## =========================
## Direct and indirect effects for M3 and M4
## =========================
effect_rows <- list()

if ("M3" %in% names(models)) {
  b_dist_mis <- posterior_beta(models$M3$mismatch, "z_logdist")
  b_mis_self <- posterior_beta(models$M3$self, "z_mismatch")
  b_self_col <- posterior_beta(models$M3$color, "prop_self_compatible")
  b_self_for <- posterior_beta(models$M3$form, "prop_self_compatible")
  effect_rows[["M3_mis_self"]] <- summ_draws(b_mis_self, "M3: mismatch -> self-compatible")
  effect_rows[["M3_self_color"]] <- summ_draws(b_self_col, "M3: self-compatible -> plain-colored flowers")
  effect_rows[["M3_self_form"]] <- summ_draws(b_self_for, "M3: self-compatible -> generalized floral forms")
  effect_rows[["M3_mis_self_color"]] <- summ_draws(b_mis_self * b_self_col, "M3 indirect: mismatch -> self -> plain-colored flowers")
  effect_rows[["M3_mis_self_form"]] <- summ_draws(b_mis_self * b_self_for, "M3 indirect: mismatch -> self -> generalized floral forms")
  effect_rows[["M3_dist_mis_self_form"]] <- summ_draws(b_dist_mis * b_mis_self * b_self_for, "M3 chain: distance -> mismatch -> self -> generalized floral forms")
  effect_rows[["M3_dist_mis_self_color"]] <- summ_draws(b_dist_mis * b_mis_self * b_self_col, "M3 chain: distance -> mismatch -> self -> plain-colored flowers")
}

if ("M4" %in% names(models)) {
  b_dist_mis <- posterior_beta(models$M4$mismatch, "z_logdist")
  b_mis_self <- posterior_beta(models$M4$self, "z_mismatch")
  b_mis_col <- posterior_beta(models$M4$color, "z_mismatch")
  b_mis_for <- posterior_beta(models$M4$form, "z_mismatch")
  b_self_col <- posterior_beta(models$M4$color, "prop_self_compatible")
  b_self_for <- posterior_beta(models$M4$form, "prop_self_compatible")
  effect_rows[["M4_mis_self"]] <- summ_draws(b_mis_self, "M4: mismatch -> self-compatible")
  effect_rows[["M4_direct_mis_color"]] <- summ_draws(b_mis_col, "M4 direct: mismatch -> plain-colored flowers")
  effect_rows[["M4_direct_mis_form"]] <- summ_draws(b_mis_for, "M4 direct: mismatch -> generalized floral forms")
  effect_rows[["M4_self_color"]] <- summ_draws(b_self_col, "M4: self-compatible -> plain-colored flowers")
  effect_rows[["M4_self_form"]] <- summ_draws(b_self_for, "M4: self-compatible -> generalized floral forms")
  effect_rows[["M4_ind_mis_self_color"]] <- summ_draws(b_mis_self * b_self_col, "M4 indirect: mismatch -> self -> plain-colored flowers")
  effect_rows[["M4_ind_mis_self_form"]] <- summ_draws(b_mis_self * b_self_for, "M4 indirect: mismatch -> self -> generalized floral forms")
  effect_rows[["M4_total_mis_color"]] <- summ_draws(b_mis_col + b_mis_self * b_self_col, "M4 total: mismatch -> plain-colored flowers")
  effect_rows[["M4_total_mis_form"]] <- summ_draws(b_mis_for + b_mis_self * b_self_for, "M4 total: mismatch -> generalized floral forms")
  effect_rows[["M4_dist_mis_self_form"]] <- summ_draws(b_dist_mis * b_mis_self * b_self_for, "M4 chain: distance -> mismatch -> self -> generalized floral forms")
  effect_rows[["M4_dist_mis_self_color"]] <- summ_draws(b_dist_mis * b_mis_self * b_self_col, "M4 chain: distance -> mismatch -> self -> plain-colored flowers")
}

effects_df <- bind_rows(effect_rows)
write.csv(effects_df, file.path(OUT_ROOT, "04_effects", "TNH_M3_M4_direct_indirect_effects.csv"), row.names = FALSE)
print(effects_df)

## =========================
## Select representative scenario for path diagram and external projection
## =========================
## Because scenario comparison is the main goal, we do not blindly force a winner.
## We pick top scenario by total ELPD but mark models with delta within 2 SE as comparable.

best_scenario <- loo_total$scenario[which.max(loo_total$elpd_total)]
loo_total <- loo_total %>%
  mutate(
    within_2se_of_best = abs(delta_elpd_total) <= 2 * se_total
  )
write.csv(loo_total, file.path(OUT_ROOT, "03_scenario_comparison", "TNH_TOTAL_LOO_interpretation_with_2SE.csv"), row.names = FALSE)
message("Top scenario by total ELPD: ", best_scenario)
message("Scenarios within approx. 2 SE of best: ", paste(loo_total$scenario[loo_total$within_2se_of_best], collapse = ", "))

## If M3 or M4 is statistically comparable, use the more mechanistic one for pathway visualization.
representative_scenario <- best_scenario
if ("M4" %in% loo_total$scenario[loo_total$within_2se_of_best]) {
  representative_scenario <- "M4"
} else if ("M3" %in% loo_total$scenario[loo_total$within_2se_of_best]) {
  representative_scenario <- "M3"
} else if ("M2" %in% loo_total$scenario[loo_total$within_2se_of_best]) {
  representative_scenario <- "M2"
}
message("Representative scenario for path diagram: ", representative_scenario)

## =========================
## Figures: model comparison and effects
## =========================
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Helvetica"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

p_total <- loo_total %>%
  mutate(scenario = factor(scenario, levels = scenario[order(delta_elpd_total)])) %>%
  ggplot(aes(x = scenario, y = delta_elpd_total)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_col(width = 0.65) +
  coord_flip() +
  labs(
    title = "Bayesian scenario comparison in northern-temperate islands",
    x = NULL,
    y = "Delta ELPD from best model (larger is better)"
  ) +
  theme_pub()

ggsave(file.path(OUT_ROOT, "05_figures", "Fig_model_comparison_TOTAL_ELPD.png"), p_total, width = 7, height = 4, dpi = 300)
ggsave(file.path(OUT_ROOT, "05_figures", "Fig_model_comparison_TOTAL_ELPD.pdf"), p_total, width = 7, height = 4)

p_resp <- loo_component %>%
  mutate(response = recode(response, !!!response_labels)) %>%
  ggplot(aes(x = scenario, y = delta_elpd_response, fill = scenario)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(width = 0.65, show.legend = FALSE) +
  facet_wrap(~ response, scales = "free_y") +
  coord_flip() +
  labs(
    title = "Response-specific scenario comparison",
    x = NULL,
    y = "Delta ELPD within response"
  ) +
  theme_pub(10)

ggsave(file.path(OUT_ROOT, "05_figures", "Fig_response_specific_ELPD.png"), p_resp, width = 8, height = 5, dpi = 300)
ggsave(file.path(OUT_ROOT, "05_figures", "Fig_response_specific_ELPD.pdf"), p_resp, width = 8, height = 5)

if (nrow(effects_df) > 0) {
  p_eff <- effects_df %>%
    filter(str_detect(effect, "indirect|direct|total|chain")) %>%
    mutate(effect = stringr::str_wrap(effect, 45),
           effect = factor(effect, levels = rev(effect[order(mean)]))) %>%
    ggplot(aes(x = mean, y = effect)) +
    geom_vline(xintercept = 0, linewidth = 0.35) +
    geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.18) +
    geom_point(size = 2.2) +
    labs(
      title = "Direct and indirect effects in mechanistic scenarios",
      x = "Posterior mean with 95% credible interval",
      y = NULL
    ) +
    theme_pub(10)
  ggsave(file.path(OUT_ROOT, "05_figures", "Fig_direct_indirect_effects.png"), p_eff, width = 8, height = 6, dpi = 300)
  ggsave(file.path(OUT_ROOT, "05_figures", "Fig_direct_indirect_effects.pdf"), p_eff, width = 8, height = 6)
}

## =========================
## Path diagram for representative scenario: significant paths only
## =========================
make_edges_from_scenario <- function(sc) {
  fe <- fixefs_df %>% filter(model == sc, significant, term != "Intercept")
  edge <- fe %>%
    mutate(
      from = case_when(
        term == "z_logdist" ~ "Distance to continent",
        term == "z_logarea" ~ "Island area",
        term == "z_mismatch" ~ "Bumblebee-environment mismatch",
        term == "prop_self_compatible" ~ "Self-compatible flora",
        TRUE ~ term
      ),
      to = recode(response, !!!response_labels),
      color = ifelse(mean >= 0, "#2166AC", "#B2182B"),
      penwidth = scales::rescale(abs(mean), to = c(1.4, 5.0), from = range(abs(mean), na.rm = TRUE)),
      label = sprintf("%.2f", mean)
    ) %>%
    filter(from != to) %>%
    select(from, to, mean, q025, q975, color, penwidth, label)
  edge
}

edges <- make_edges_from_scenario(representative_scenario)
write.csv(edges, file.path(OUT_ROOT, "04_effects", paste0("TNH_", representative_scenario, "_significant_edges.csv")), row.names = FALSE)

nodes <- c("Distance to continent", "Island area", "Bumblebee-environment mismatch", "Self-compatible flora", "Plain-colored flowers", "Generalized floral forms")
node_defs <- paste0(
  '  "', nodes, '" [label="', nodes, '", shape=box, style="rounded,filled", fillcolor="white", color="gray35", fontname="Helvetica"];',
  collapse = "\n"
)

edge_defs <- if (nrow(edges) == 0) {
  ""
} else {
  paste0(
    '  "', edges$from, '" -> "', edges$to, '" [label="', edges$label,
    '", color="', edges$color, '", fontcolor="', edges$color,
    '", penwidth=', sprintf("%.2f", edges$penwidth), ', arrowsize=0.8, fontname="Helvetica"];',
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

htmlwidget <- DiagrammeR::grViz(graph_code)
svg_txt <- DiagrammeRsvg::export_svg(htmlwidget)
svg_file <- file.path(OUT_ROOT, "05_figures", paste0("Fig_path_diagram_", representative_scenario, "_sig_only.svg"))
png_file <- file.path(OUT_ROOT, "05_figures", paste0("Fig_path_diagram_", representative_scenario, "_sig_only.png"))
pdf_file <- file.path(OUT_ROOT, "05_figures", paste0("Fig_path_diagram_", representative_scenario, "_sig_only.pdf"))
writeLines(svg_txt, svg_file)
rsvg::rsvg_png(charToRaw(svg_txt), png_file, width = 2200, height = 1200)
rsvg::rsvg_pdf(charToRaw(svg_txt), pdf_file, width = 11, height = 6)

## =========================
## Supplement: external projection to other regions
## =========================
## Lightweight external projection. Fit fixed-effect GLMs on TNH for the representative scenario,
## then evaluate negative log likelihood on tropical/southern. This is supplemental and fast.

fit_glm_component <- function(sc, resp, train) {
  ## Use same formulas conceptually, but base glm.
  if (resp == "mismatch") {
    return(lm(z_mismatch ~ z_logdist, data = train))
  }
  if (resp == "self") {
    if (sc %in% c("M3", "M4")) {
      return(glm(cbind(n_self_compatible, n_self_total - n_self_compatible) ~ z_mismatch, data = train, family = binomial()))
    } else {
      return(glm(cbind(n_self_compatible, n_self_total - n_self_compatible) ~ z_logdist, data = train, family = binomial()))
    }
  }
  if (resp == "color") {
    if (sc == "M0") {
      return(glm(cbind(n_plain, n_color_total - n_plain) ~ z_logdist + z_logarea, data = train, family = binomial()))
    } else if (sc %in% c("M2", "M3")) {
      return(glm(cbind(n_plain, n_color_total - n_plain) ~ prop_self_compatible + z_logarea, data = train, family = binomial()))
    } else {
      return(glm(cbind(n_plain, n_color_total - n_plain) ~ z_mismatch + prop_self_compatible + z_logarea, data = train, family = binomial()))
    }
  }
  if (resp == "form") {
    if (sc == "M0") {
      return(glm(cbind(n_generalized, n_form_total - n_generalized) ~ z_logdist + z_logarea, data = train, family = binomial()))
    } else if (sc %in% c("M2", "M3")) {
      return(glm(cbind(n_generalized, n_form_total - n_generalized) ~ prop_self_compatible + z_logarea, data = train, family = binomial()))
    } else {
      return(glm(cbind(n_generalized, n_form_total - n_generalized) ~ z_mismatch + prop_self_compatible + z_logarea, data = train, family = binomial()))
    }
  }
}

nll_binom <- function(y, n, p) {
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  -sum(dbinom(y, size = n, prob = p, log = TRUE), na.rm = TRUE)
}

nll_gauss <- function(y, mu, sigma) {
  sigma <- max(sigma, 1e-8)
  -sum(dnorm(y, mean = mu, sd = sigma, log = TRUE), na.rm = TRUE)
}

external_eval <- function(train_df, test_df, sc, label) {
  train_df <- as.data.frame(train_df)
  test_df  <- as.data.frame(test_df)

  if (nrow(test_df) < 5) return(NULL)

  fm <- fit_glm_component(sc, "mismatch", train_df)
  fs <- fit_glm_component(sc, "self", train_df)
  fc <- fit_glm_component(sc, "color", train_df)
  ff <- fit_glm_component(sc, "form", train_df)

  pred_m <- as.numeric(predict(fm, newdata = test_df, type = "response"))
  pred_s <- as.numeric(predict(fs, newdata = test_df, type = "response"))
  pred_c <- as.numeric(predict(fc, newdata = test_df, type = "response"))
  pred_f <- as.numeric(predict(ff, newdata = test_df, type = "response"))

  sigma_m <- stats::sd(stats::resid(fm), na.rm = TRUE)
  if (!is.finite(sigma_m) || sigma_m <= 0) sigma_m <- 1e-6

  nll_mismatch_val <- nll_gauss(test_df$z_mismatch, pred_m, sigma_m)
  nll_self_val     <- nll_binom(test_df$n_self_compatible, test_df$n_self_total, pred_s)
  nll_color_val    <- nll_binom(test_df$n_plain, test_df$n_color_total, pred_c)
  nll_form_val     <- nll_binom(test_df$n_generalized, test_df$n_form_total, pred_f)

  tibble(
    train = "temperate_NH",
    test = label,
    scenario = sc,
    n_islands = nrow(test_df),
    nll_mismatch = nll_mismatch_val,
    nll_self = nll_self_val,
    nll_color = nll_color_val,
    nll_form = nll_form_val,
    nll_traits = nll_color_val + nll_form_val,
    nll_full = nll_mismatch_val + nll_self_val + nll_color_val + nll_form_val,
    nll_traits_per_island = nll_traits / n_islands,
    nll_full_per_island = nll_full / n_islands
  )
}

train_tnh <- island_all %>%
  filter(analysis_group == "temperate_NH") %>%
  filter(is.finite(z_logdist), is.finite(z_logarea), is.finite(z_mismatch),
         n_color_total > 0, n_form_total > 0, n_self_total > 0)

trop <- island_all %>%
  filter(analysis_group == "tropical_all") %>%
  filter(is.finite(z_logdist), is.finite(z_logarea), is.finite(z_mismatch),
         n_color_total > 0, n_form_total > 0, n_self_total > 0)

south <- island_all %>%
  filter(analysis_group == "southern_non_tropical") %>%
  filter(is.finite(z_logdist), is.finite(z_logarea), is.finite(z_mismatch),
         n_color_total > 0, n_form_total > 0, n_self_total > 0)

ext <- bind_rows(
  external_eval(train_tnh, train_tnh, representative_scenario, "temperate_NH_in_sample_glm"),
  external_eval(train_tnh, trop, representative_scenario, "tropical_external"),
  external_eval(train_tnh, south, representative_scenario, "southern_non_tropical_external")
)
write.csv(ext, file.path(OUT_ROOT, "06_supp_external_projection", "TNH_representative_scenario_external_projection_NLL.csv"), row.names = FALSE)
print(ext)

if (!is.null(ext) && nrow(ext) > 0) {
  p_ext <- ext %>%
    ggplot(aes(x = test, y = nll_traits_per_island)) +
    geom_col(width = 0.65) +
    coord_flip() +
    labs(
      title = "External projection of the northern-temperate scenario",
      x = NULL,
      y = "Trait prediction negative log-likelihood per island"
    ) +
    theme_pub()
  ggsave(file.path(OUT_ROOT, "05_figures", "Fig_supp_external_projection_NLL.png"), p_ext, width = 7, height = 4, dpi = 300)
  ggsave(file.path(OUT_ROOT, "05_figures", "Fig_supp_external_projection_NLL.pdf"), p_ext, width = 7, height = 4)
}

## =========================
## README
## =========================
readme <- c(
  "Final TNH Bayesian scenario comparison, light computation version",
  "", 
  "Main purpose:",
  "Compare M1, M2, M3, M4 as biologically motivated island-syndrome scenarios using Bayesian piecewise path models.",
  "", 
  "Scenario meanings:",
  "M1: geographic filter; distance and island area directly predict floral traits.",
  "M2: selfing syndrome; distance predicts self-compatible flora, and self-compatible flora predicts floral traits; island area only predicts final floral traits.",
  "M3: full mediation; distance predicts bumblebee-environment mismatch, mismatch predicts self-compatible flora, self-compatible flora predicts floral traits; island area only predicts final floral traits.",
  "M4: partial mediation; M3 plus direct mismatch-to-floral-trait paths; island area only predicts final floral traits.",
  "", 
  "Important interpretation:",
  "LOO-CV/ELPD ranks predictive consistency, not definitive causality. If models are within uncertainty, interpret them as comparable and use path effects to discuss biological plausibility.",
  "", 
  "Cache/reproducibility:",
  "If component-model .rds files are missing, this script fits them from scratch and saves them under 02_bayes_component_models/.",
  "Set FORCE_REFIT_MODELS <- TRUE near the top if you want to ignore old .rds files and refit everything.",
  "",
  "Key outputs:",
  "03_scenario_comparison/TNH_M0_M2_M3_M4_TOTAL_LOO.csv",
  "03_scenario_comparison/TNH_LOO_components_by_response.csv",
  "04_effects/TNH_M3_M4_direct_indirect_effects.csv",
  "05_figures/Fig_model_comparison_TOTAL_ELPD.png",
  "05_figures/Fig_response_specific_ELPD.png",
  paste0("05_figures/Fig_path_diagram_", representative_scenario, "_sig_only.png"),
  "06_supp_external_projection/TNH_representative_scenario_external_projection_NLL.csv"
)
writeLines(readme, file.path(OUT_ROOT, "README_analysis_design.txt"))

message("\nDONE.")
message("Results saved to: ", OUT_ROOT)
message("Representative scenario: ", representative_scenario)
message("Key files:")
message(" - 03_scenario_comparison/TNH_M0_M2_M3_M4_TOTAL_LOO.csv")
message(" - 03_scenario_comparison/TNH_LOO_components_by_response.csv")
message(" - 04_effects/TNH_M3_M4_direct_indirect_effects.csv")
message(" - 05_figures/")
message(" - 06_supp_external_projection/TNH_representative_scenario_external_projection_NLL.csv")
