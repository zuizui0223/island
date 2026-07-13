#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: run_inla_final_flora_composition.R <input.csv> <output_dir>")
infile <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260714)

d <- fread(infile)
colors <- c("white", "blue", "yellow", "red", "green", "cream", "rare_other")
forms <- c("tubular", "composite_head", "star", "bell_campanulate", "open_radial", "trumpet_funnel", "rare_other")
endpoints <- list(
  color = list(
    all_prefix = "all_fine_color_", all_trials = "all_fine_color_trials",
    animal_prefix = "fine_color_", animal_trials = "fine_color_trials", categories = colors
  ),
  form = list(
    all_prefix = "all_fine_form_", all_trials = "all_fine_form_trials",
    animal_prefix = "fine_form_", animal_trials = "fine_form_trials", categories = forms
  )
)

regime_levels <- c(
  "northern_midlatitude", "northern_high_latitude",
  "tropical", "southern_extratropical"
)
d[, analysis_regime := factor(analysis_regime, levels = regime_levels)]
d[, wind_share := fifelse(animal_status_observed > 0, n_wind_species / animal_status_observed, NA_real_)]
d[, mixed_share := fifelse(animal_status_observed > 0, n_mixed_species / animal_status_observed, NA_real_)]

base_required <- c(
  "island_id", "analysis_regime", "spatial_block",
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
  "bombus_deficit", "sc_share", "animal_status_observed", "n_wind_species", "n_mixed_species"
)
required <- unique(c(base_required, unlist(lapply(endpoints, function(x) {
  c(x$all_trials, x$animal_trials,
    paste0(x$all_prefix, x$categories), paste0(x$animal_prefix, x$categories))
}))))
missing <- setdiff(required, names(d))
if (length(missing)) stop("missing final-model columns: ", paste(missing, collapse = ", "))

for (ep in names(endpoints)) {
  x <- endpoints[[ep]]
  all_cols <- paste0(x$all_prefix, x$categories)
  animal_cols <- paste0(x$animal_prefix, x$categories)
  if (any(rowSums(d[, ..all_cols], na.rm = TRUE) != fifelse(is.na(d[[x$all_trials]]), 0, d[[x$all_trials]]))) {
    stop(ep, " all-flora categories do not sum to trials")
  }
  if (any(rowSums(d[, ..animal_cols], na.rm = TRUE) != fifelse(is.na(d[[x$animal_trials]]), 0, d[[x$animal_trials]]))) {
    stop(ep, " animal-flora categories do not sum to trials")
  }
}

# Scale isolation/geography on all positive-distance islands.
global_scale_vars <- c(
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"
)
for (nm in global_scale_vars) {
  mu <- mean(d[[nm]], na.rm = TRUE)
  sig <- sd(d[[nm]], na.rm = TRUE)
  if (!is.finite(sig) || sig <= 0) stop("invalid scale for ", nm)
  d[[paste0("z_", nm)]] <- (d[[nm]] - mu) / sig
}
set(d, j = "z_isolation_sq", value = d$z_log_distance_to_continent_km^2)

# Mechanism variables are scaled on the northern reference population. A
# zero-variance variable is not statistically identifiable and is omitted from
# fitted equations, but its availability and reason are written to an audit.
north <- d[analysis_regime == "northern_midlatitude"]
mechanism_candidates <- c("bombus_deficit", "sc_share", "wind_share", "mixed_share")
mechanism_audit <- rbindlist(lapply(mechanism_candidates, function(nm) {
  x <- north[[nm]]
  n_nonmissing <- sum(is.finite(x))
  n_unique <- uniqueN(x[is.finite(x)])
  sig <- sd(x, na.rm = TRUE)
  estimable <- n_nonmissing >= 2 && n_unique >= 2 && is.finite(sig) && sig > 0
  data.table(
    variable = nm,
    n_nonmissing = n_nonmissing,
    n_unique = n_unique,
    sd = sig,
    estimable = estimable,
    reason = if (estimable) "estimable" else "zero_or_insufficient_variation"
  )
}))
fwrite(mechanism_audit, file.path(outdir, "final_mechanism_variable_audit.csv"))

for (nm in mechanism_audit[estimable == TRUE, variable]) {
  mu <- mean(north[[nm]], na.rm = TRUE)
  sig <- sd(north[[nm]], na.rm = TRUE)
  d[[paste0("z_", nm)]] <- (d[[nm]] - mu) / sig
}

if (!"z_bombus_deficit" %in% names(d)) stop("bombus_deficit is not estimable in northern mechanism population")
if (!"z_sc_share" %in% names(d)) stop("sc_share is not estimable in northern mechanism population")

pollination_controls <- mechanism_audit[
  variable %in% c("wind_share", "mixed_share") & estimable == TRUE,
  variable
]
pollination_control_terms <- paste0("category:z_", pollination_controls)
pollination_control_sc_terms <- paste0("z_", pollination_controls)

compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
fixed_rows <- list()
score_rows <- list()
support_rows <- list()

fit_poisson_composition <- function(data, endpoint, model_name, prefix, trials, categories, rhs) {
  cols <- paste0(prefix, categories)
  candidate_id_vars <- c(
    "island_id", "spatial_block", "analysis_regime", trials,
    "z_log_distance_to_continent_km", "z_isolation_sq", "z_log_island_area_km2",
    "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
    "z_bombus_deficit", "z_sc_share", "z_wind_share", "z_mixed_share"
  )
  long <- melt(
    data,
    id.vars = intersect(candidate_id_vars, names(data)),
    measure.vars = cols,
    variable.name = "category", value.name = "count"
  )
  long[, category := factor(sub(paste0("^", prefix), "", category), levels = categories)]
  long[, island_id_re := as.integer(factor(island_id))]
  long[, spatial_id := as.integer(factor(spatial_block))]
  long[, obs_id := seq_len(.N)]
  long[, log_trials := log(get(trials))]
  fit <- inla(
    as.formula(paste("count ~", rhs, "+ offset(log_trials)")),
    family = "poisson", data = long,
    control.compute = compute,
    control.predictor = list(compute = TRUE), verbose = FALSE
  )
  fx <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
  fx[, `:=`(endpoint = endpoint, model = model_name)]
  fixed_rows[[length(fixed_rows) + 1L]] <<- fx
  cpo <- fit$cpo$cpo
  cpo <- cpo[is.finite(cpo) & cpo > 0]
  score_rows[[length(score_rows) + 1L]] <<- data.table(
    endpoint = endpoint, model = model_name,
    n_islands = uniqueN(data$island_id), n_rows = nrow(long),
    log_cpo_sum = if (length(cpo)) sum(log(cpo)) else NA_real_,
    waic = fit$waic$waic, dic = fit$dic$dic
  )
}

geo_terms <- paste(
  "category:z_log_island_area_km2 +",
  "category:z_climate_pc1 + category:z_climate_pc2 +",
  "category:z_climate_pc3 + category:z_climate_pc4"
)
random_terms <- "+ f(island_id_re, model='iid') + f(spatial_id, model='iid') + f(obs_id, model='iid')"

for (ep in names(endpoints)) {
  x <- endpoints[[ep]]

  # GLOBAL PRIMARY: all trait-resolved flora, irrespective of pollination mode.
  global_support <- d[
    !is.na(analysis_regime) & get(x$all_trials) > 0 &
      complete.cases(d[, .(
        z_log_distance_to_continent_km, z_isolation_sq, z_log_island_area_km2,
        z_climate_pc1, z_climate_pc2, z_climate_pc3, z_climate_pc4
      )])
  ]
  global_rhs <- paste(
    "0 + category +",
    "category:analysis_regime +",
    "category:analysis_regime:z_log_distance_to_continent_km +",
    "category:analysis_regime:z_isolation_sq +",
    geo_terms, random_terms
  )
  fit_poisson_composition(
    global_support, paste0("all_flora_", ep), "G_global_isolation_by_regime",
    x$all_prefix, x$all_trials, x$categories, global_rhs
  )
  support_rows[[length(support_rows) + 1L]] <- data.table(
    endpoint = paste0("all_flora_", ep), layer = "global_primary_all_flora",
    n_islands = nrow(global_support)
  )

  # NORTHERN MECHANISM: animal-pollinated trait-resolved flora only.
  mech_vars <- c(
    "z_log_distance_to_continent_km", "z_isolation_sq", "z_log_island_area_km2",
    "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
    "z_bombus_deficit", "z_sc_share", paste0("z_", pollination_controls)
  )
  north_support <- d[
    analysis_regime == "northern_midlatitude" & get(x$animal_trials) > 0 &
      complete.cases(d[, ..mech_vars])
  ]
  base_rhs <- paste(
    "0 + category + category:z_log_distance_to_continent_km + category:z_isolation_sq +",
    geo_terms, random_terms
  )
  n3_extra <- c("category:z_bombus_deficit", "category:z_sc_share", pollination_control_terms)
  north_models <- list(
    N0_isolation = base_rhs,
    N1_bombus = paste(base_rhs, "+ category:z_bombus_deficit"),
    N2_bombus_sc = paste(base_rhs, "+ category:z_bombus_deficit + category:z_sc_share"),
    N3_full = paste(base_rhs, "+", paste(n3_extra, collapse = " + "))
  )
  for (nm in names(north_models)) {
    fit_poisson_composition(
      north_support, paste0("animal_flora_", ep), nm,
      x$animal_prefix, x$animal_trials, x$categories, north_models[[nm]]
    )
  }
  support_rows[[length(support_rows) + 1L]] <- data.table(
    endpoint = paste0("animal_flora_", ep), layer = "northern_bombus_mechanism",
    n_islands = nrow(north_support)
  )
}

# Explicit reproductive-assurance equation for the proposed indirect pathway.
sc_mech_vars <- c(
  "z_log_distance_to_continent_km", "z_isolation_sq", "z_log_island_area_km2",
  "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
  "z_bombus_deficit", pollination_control_sc_terms
)
sc_support <- d[
  analysis_regime == "northern_midlatitude" & !is.na(sc_share) & sc_trials > 0 &
    complete.cases(d[, ..sc_mech_vars]) & !is.na(sc_successes)
]
sc_rhs <- c(
  "z_log_distance_to_continent_km", "z_isolation_sq", "z_bombus_deficit",
  pollination_control_sc_terms,
  "z_log_island_area_km2", "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
  "f(spatial_id_sc, model='iid')"
)
sc_support[, spatial_id_sc := as.integer(factor(spatial_block))]
sc_formula <- as.formula(paste("sc_successes ~", paste(sc_rhs, collapse = " + ")))
sc_fit <- inla(
  sc_formula, family = "betabinomial", data = sc_support, Ntrials = sc_trials,
  control.compute = compute, control.predictor = list(compute = TRUE), verbose = FALSE
)
sc_fx <- as.data.table(sc_fit$summary.fixed, keep.rownames = "parameter")
sc_fx[, `:=`(endpoint = "self_compatibility", model = "SC_bombus_reproductive_assurance")]
fixed_rows[[length(fixed_rows) + 1L]] <- sc_fx

fixed <- rbindlist(fixed_rows, fill = TRUE)
fixed[, excludes_zero := (`0.025quant` > 0 & `0.975quant` > 0) |
  (`0.025quant` < 0 & `0.975quant` < 0)]
scores <- rbindlist(score_rows, fill = TRUE)
support <- rbindlist(support_rows, fill = TRUE)

fwrite(fixed, file.path(outdir, "final_flora_fixed_effects.csv"))
fwrite(scores, file.path(outdir, "final_flora_model_comparisons.csv"))
fwrite(support, file.path(outdir, "final_flora_support.csv"))
fwrite(sc_fx, file.path(outdir, "final_sc_pathway.csv"))

key <- fixed[
  grepl("z_log_distance_to_continent_km|z_isolation_sq|z_bombus_deficit|z_sc_share|z_wind_share|z_mixed_share", parameter)
]
fwrite(key, file.path(outdir, "final_flora_key_effects.csv"))

meta <- list(
  primary_response_definition = paste(
    "Conditional relative category probability within the observed, trait-resolved island flora;",
    "not raw detection probability and not latent true occupancy."
  ),
  global_primary_population = "all trait-resolved flora, irrespective of animal, wind, mixed, or unresolved pollination mode",
  primary_endpoints = c("fine colour composition", "fine floral-form composition"),
  isolation = "continuous standardized log distance plus category-specific quadratic term",
  global_test = "category-specific isolation and isolation-squared effects by four analysis regimes",
  northern_mechanism_population = "animal-pollinated trait-resolved flora in northern_midlatitude islands",
  northern_mechanisms = c(
    "Bombus deficit: reduced/relaxed Bombus-mediated selection",
    "self-compatibility: reproductive-assurance pathway consistent with selfing-syndrome evolution",
    paste0(
      "pollination-mode composition controls included only when estimable: ",
      if (length(pollination_controls)) paste(pollination_controls, collapse = ", ") else "none"
    )
  ),
  nonestimable_mechanism_controls = mechanism_audit[estimable == FALSE, variable],
  joint_colour_form = "sensitivity only, not part of the final primary model"
)
write_json(meta, file.path(outdir, "final_flora_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
