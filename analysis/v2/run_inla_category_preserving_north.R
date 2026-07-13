#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: run_inla_category_preserving_north.R <input.csv> <output_dir>")
infile <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260713)

d <- fread(infile)
colour_outcomes <- list(
  plain = c("color_plain", "color_trials"),
  yellow_orange = c("color_yellow_orange", "color_trials"),
  red_pink = c("color_red_pink", "color_trials"),
  blue_purple = c("color_blue_purple", "color_trials")
)
form_outcomes <- list(
  open_generalized = c("form_open_generalized", "form_trials"),
  tubular_trumpet = c("form_tubular_trumpet", "form_trials"),
  zygomorphic_specialized = c("form_zygomorphic_specialized", "form_trials"),
  composite_brush = c("form_composite_brush", "form_trials")
)
trait_outcomes <- c(colour_outcomes, form_outcomes)

required <- c(
  "island_id", "analysis_tier", "analysis_regime",
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
  "bombus_deficit", "sc_trials", "sc_successes", "sc_share", "spatial_block",
  "color_trials", unname(vapply(colour_outcomes, `[[`, character(1), 1)),
  "form_trials", unname(vapply(form_outcomes, `[[`, character(1), 1))
)
missing <- setdiff(required, names(d))
if (length(missing)) stop("missing required columns: ", paste(missing, collapse = ", "))

# Verify that the retained categories exhaust the recorded totals.
colour_sum <- Reduce(`+`, lapply(colour_outcomes, function(x) fifelse(is.na(d[[x[[1]]]]), 0, d[[x[[1]]]])))
form_sum <- Reduce(`+`, lapply(form_outcomes, function(x) fifelse(is.na(d[[x[[1]]]]), 0, d[[x[[1]]]])))
if (any(colour_sum != fifelse(is.na(d$color_trials), 0, d$color_trials))) {
  stop("retained colour categories do not sum to color_trials")
}
if (any(form_sum != fifelse(is.na(d$form_trials), 0, d$form_trials))) {
  stop("retained floral-form categories do not sum to form_trials")
}

scale_vars <- c(
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
  "bombus_deficit", "sc_share"
)
north_ref <- d[analysis_regime == "northern_midlatitude"]
scaling <- rbindlist(lapply(scale_vars, function(nm) {
  x <- north_ref[[nm]]
  data.table(variable = nm, center = mean(x, na.rm = TRUE), scale = sd(x, na.rm = TRUE))
}))
if (any(!is.finite(scaling$scale) | scaling$scale <= 0)) stop("invalid northern scaling")
for (i in seq_len(nrow(scaling))) {
  nm <- scaling$variable[[i]]
  d[[paste0("z_", nm)]] <- (d[[nm]] - scaling$center[[i]]) / scaling$scale[[i]]
}
fwrite(scaling, file.path(outdir, "northern_training_scaling.csv"))

d[, spatial_id := as.integer(factor(spatial_block))]
d[, regime := factor(analysis_regime, levels = c(
  "northern_midlatitude", "tropical", "southern_extratropical"
))]

geo_vars <- c(
  "z_log_distance_to_continent_km", "z_log_island_area_km2",
  "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4"
)
base_vars <- c(geo_vars, "z_bombus_deficit")
geo_fixed <- paste(
  "z_log_distance_to_continent_km * z_log_island_area_km2 +",
  "z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4"
)
geo_fit <- paste(geo_fixed, "+ f(spatial_id, model='iid')")
compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)

fit_grouped <- function(success_col, trials_col, rhs, dd) {
  x <- copy(dd)
  if (nrow(x) < 30) stop("too few rows for ", success_col, ": ", nrow(x))
  x[, obs_id := seq_len(.N)]
  inla(
    as.formula(paste(success_col, "~", rhs, "+ f(obs_id, model='iid')")),
    family = "binomial", data = x, Ntrials = x[[trials_col]],
    control.compute = compute,
    control.predictor = list(compute = TRUE),
    verbose = FALSE
  )
}

score_fit <- function(fit) {
  cpo <- fit$cpo$cpo
  cpo <- cpo[is.finite(cpo) & cpo > 0]
  data.table(
    log_cpo_sum = if (length(cpo)) sum(log(cpo)) else NA_real_,
    n_cpo = length(cpo), waic = fit$waic$waic, dic = fit$dic$dic
  )
}

fixed_table <- function(fit, outcome_name, model_name, support_name, n_support) {
  x <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
  x[, `:=`(
    outcome = outcome_name,
    model = model_name,
    support = support_name,
    n_support = n_support
  )]
  x
}

sample_fixed <- function(fit, parameter, n = 50000L) {
  if (!parameter %in% names(fit$marginals.fixed)) stop("missing parameter: ", parameter)
  inla.rmarginal(n, fit$marginals.fixed[[parameter]])
}

summarize_draws <- function(x, outcome_name, effect_name, n_support) {
  q <- quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)
  data.table(
    outcome = outcome_name,
    effect = effect_name,
    n_support = n_support,
    mean = mean(x, na.rm = TRUE),
    median = unname(q[[2]]),
    `0.025quant` = unname(q[[1]]),
    `0.975quant` = unname(q[[3]]),
    excludes_zero = unname(q[[1]] > 0 | q[[3]] < 0)
  )
}

base_complete <- d[complete.cases(d[, ..base_vars]) & !is.na(regime)]
north_base <- base_complete[regime == "northern_midlatitude"]
if (nrow(north_base) < 30) stop("too few northern-midlatitude islands")

fits <- list()
score_rows <- list()
fixed_rows <- list()
support_rows <- list()

sc_support <- north_base[sc_trials > 0]
support_rows[["self_compatibility"]] <- data.table(
  outcome = "self_compatibility", support = "maximum_equation_support", n_islands = nrow(sc_support)
)
sc_models <- list(
  N0_geo = geo_fit,
  N1_bombus = paste("z_bombus_deficit +", geo_fit)
)
for (model_name in names(sc_models)) {
  key <- paste("self_compatibility", model_name, sep = "__")
  fit <- fit_grouped("sc_successes", "sc_trials", sc_models[[model_name]], sc_support)
  fits[[key]] <- fit
  s <- score_fit(fit)
  s[, `:=`(
    outcome = "self_compatibility", model = model_name,
    support = "maximum_equation_support", n_support = nrow(sc_support)
  )]
  score_rows[[key]] <- s
  fixed_rows[[key]] <- fixed_table(
    fit, "self_compatibility", model_name, "maximum_equation_support", nrow(sc_support)
  )
}

for (outcome_name in names(trait_outcomes)) {
  success_col <- trait_outcomes[[outcome_name]][[1]]
  trials_col <- trait_outcomes[[outcome_name]][[2]]

  max_support <- north_base[get(trials_col) > 0]
  path_support <- max_support[is.finite(z_sc_share)]
  if (nrow(max_support) < 30) stop("too few maximum-support rows for ", outcome_name)
  if (nrow(path_support) < 30) stop("too few path-support rows for ", outcome_name)

  support_rows[[paste0(outcome_name, "__max")]] <- data.table(
    outcome = outcome_name, support = "maximum_category_support", n_islands = nrow(max_support)
  )
  support_rows[[paste0(outcome_name, "__path")]] <- data.table(
    outcome = outcome_name, support = "sc_path_support", n_islands = nrow(path_support)
  )

  ladders <- list(
    maximum_category_support = list(
      data = max_support,
      models = list(
        N0_geo_max = geo_fit,
        N1_bombus_max = paste("z_bombus_deficit +", geo_fit)
      )
    ),
    sc_path_support = list(
      data = path_support,
      models = list(
        N0_geo_path = geo_fit,
        N1_bombus_path = paste("z_bombus_deficit +", geo_fit),
        N2_sc_path = paste("z_sc_share +", geo_fit),
        N3_direct_indirect_path = paste("z_bombus_deficit + z_sc_share +", geo_fit)
      )
    )
  )

  for (support_name in names(ladders)) {
    dd <- ladders[[support_name]]$data
    for (model_name in names(ladders[[support_name]]$models)) {
      key <- paste(outcome_name, model_name, sep = "__")
      fit <- fit_grouped(success_col, trials_col, ladders[[support_name]]$models[[model_name]], dd)
      fits[[key]] <- fit
      s <- score_fit(fit)
      s[, `:=`(
        outcome = outcome_name, model = model_name,
        support = support_name, n_support = nrow(dd)
      )]
      score_rows[[key]] <- s
      fixed_rows[[key]] <- fixed_table(fit, outcome_name, model_name, support_name, nrow(dd))
    }
  }
}

scores <- rbindlist(score_rows, fill = TRUE)
scores[, delta_log_cpo_from_best := log_cpo_sum - max(log_cpo_sum, na.rm = TRUE), by = .(outcome, support)]
setorder(scores, outcome, support, -log_cpo_sum)
fwrite(scores, file.path(outdir, "northern_model_comparisons.csv"))

fixed <- rbindlist(fixed_rows, fill = TRUE)
fixed[, excludes_zero := (`0.025quant` > 0 & `0.975quant` > 0) |
  (`0.025quant` < 0 & `0.975quant` < 0)]
fwrite(fixed, file.path(outdir, "northern_fixed_effects.csv"))
fwrite(rbindlist(support_rows), file.path(outdir, "northern_equation_support.csv"))

a_fit <- fits[["self_compatibility__N1_bombus"]]
a_draw <- sample_fixed(a_fit, "z_bombus_deficit")
mediation_rows <- list()
for (outcome_name in names(trait_outcomes)) {
  trait_fit <- fits[[paste0(outcome_name, "__N3_direct_indirect_path")]]
  direct_draw <- sample_fixed(trait_fit, "z_bombus_deficit")
  b_draw <- sample_fixed(trait_fit, "z_sc_share")
  indirect_draw <- a_draw * b_draw
  total_draw <- direct_draw + indirect_draw
  n_support <- unique(fixed[
    outcome == outcome_name & model == "N3_direct_indirect_path", n_support
  ])

  mediation_rows[[paste0(outcome_name, "_a")]] <- summarize_draws(
    a_draw, outcome_name, "a_bombus_to_self_compatibility", nrow(sc_support)
  )
  mediation_rows[[paste0(outcome_name, "_b")]] <- summarize_draws(
    b_draw, outcome_name, "b_self_compatibility_to_category", n_support
  )
  mediation_rows[[paste0(outcome_name, "_direct")]] <- summarize_draws(
    direct_draw, outcome_name, "direct_bombus_to_category", n_support
  )
  mediation_rows[[paste0(outcome_name, "_indirect")]] <- summarize_draws(
    indirect_draw, outcome_name, "indirect_bombus_via_self_compatibility", n_support
  )
  mediation_rows[[paste0(outcome_name, "_total")]] <- summarize_draws(
    total_draw, outcome_name, "total_bombus_effect", n_support
  )
}
mediation_effects <- rbindlist(mediation_rows)
fwrite(mediation_effects, file.path(outdir, "northern_direct_indirect_total_effects.csv"))

key_effects <- fixed[
  (outcome == "self_compatibility" & model == "N1_bombus" & parameter == "z_bombus_deficit") |
  (outcome %in% names(trait_outcomes) & model == "N3_direct_indirect_path" &
     parameter %in% c("z_bombus_deficit", "z_sc_share"))
]
fwrite(key_effects, file.path(outdir, "northern_bombus_syndrome_effects.csv"))

falsification_outcomes <- c(
  list(self_compatibility = c("sc_successes", "sc_trials")),
  trait_outcomes
)
regime_rows <- list()
regime_fits <- list()
regime_support_rows <- list()
for (outcome_name in names(falsification_outcomes)) {
  success_col <- falsification_outcomes[[outcome_name]][[1]]
  trials_col <- falsification_outcomes[[outcome_name]][[2]]
  dd <- copy(base_complete[get(trials_col) > 0])
  support <- dd[, .N, by = analysis_regime]
  support[, outcome := outcome_name]
  setcolorder(support, c("outcome", "analysis_regime", "N"))
  regime_support_rows[[outcome_name]] <- support
  dd[, obs_id := seq_len(.N)]
  rhs <- paste(
    "regime * z_log_distance_to_continent_km + z_log_island_area_km2 +",
    "z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4 +",
    "f(spatial_id, model='iid') + f(obs_id, model='iid')"
  )
  fit <- inla(
    as.formula(paste(success_col, "~", rhs)),
    family = "binomial", data = dd, Ntrials = dd[[trials_col]],
    control.compute = compute,
    control.predictor = list(compute = TRUE),
    verbose = FALSE
  )
  regime_fits[[outcome_name]] <- fit
  f <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
  f[, outcome := outcome_name]
  regime_rows[[outcome_name]] <- f
}

regime_effects <- rbindlist(regime_rows, fill = TRUE)
regime_effects[, excludes_zero := (`0.025quant` > 0 & `0.975quant` > 0) |
  (`0.025quant` < 0 & `0.975quant` < 0)]
fwrite(regime_effects, file.path(outdir, "cross_regime_isolation_effects.csv"))
fwrite(
  regime_effects[grepl("z_log_distance_to_continent_km", parameter),
    .(outcome, parameter, mean, `0.025quant`, `0.975quant`, excludes_zero)],
  file.path(outdir, "cross_regime_isolation_summary.csv")
)
fwrite(rbindlist(regime_support_rows), file.path(outdir, "cross_regime_support.csv"))

meta <- list(
  contract = "v2_all_retained_categories_northern_paths_global_falsification_v4",
  analysis_tier = unique(d$analysis_tier),
  primary_engine = "INLA",
  primary_domain = "northern_midlatitude",
  retained_colour_categories = names(colour_outcomes),
  retained_floral_form_categories = names(form_outcomes),
  n_retained_colour_categories = length(colour_outcomes),
  n_retained_floral_form_categories = length(form_outcomes),
  support_strategy = paste(
    "maximum equation support for SC and N0/N1 category models;",
    "matched SC-path support for N0-N3 and mediation"
  ),
  category_totals_verified = TRUE,
  zero_distance_rule = "excluded before this script by workflow",
  mediation_uncertainty = "product-of-INLA-posterior-marginals approximation",
  falsification = "isolation-by-regime slopes for SC and every retained colour/form category",
  alternative_pollinator_covariates = "excluded from primary models"
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(
  c(lapply(fits, summary), lapply(regime_fits, summary)),
  file = file.path(outdir, "model_summaries.txt")
)
