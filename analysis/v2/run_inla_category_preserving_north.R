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
set.seed(20260714)

d <- fread(infile)
required <- c(
  "island_id", "analysis_tier", "analysis_regime",
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4", "spatial_block",
  "bombus_deficit", "bombus_environmental_compatibility",
  "environmental_compatibility_max", "environmental_compatibility_mean",
  "environmental_extrapolation_share", "n_bombus_species_total",
  "n_bombus_species_scored", "n_bombus_species_unscored",
  "n_environmentally_compatible_species",
  "sc_trials", "sc_successes", "sc_share",
  "wind_mixed_trials", "wind_mixed_successes", "wind_mixed_share",
  "animal_status_observed", "n_animal_species", "n_wind_species", "n_species_total",
  "showy_alternative_guild_trials", "showy_alternative_guild_successes",
  "other_bee_guild_trials", "other_bee_guild_successes",
  "generalist_insect_guild_trials", "generalist_insect_guild_successes",
  "color_trials", "color_plain", "color_yellow_orange", "color_red_pink", "color_blue_purple",
  "form_trials", "form_open_generalized", "form_tubular_trumpet",
  "form_zygomorphic_specialized", "form_composite_brush"
)
missing <- setdiff(required, names(d))
if (length(missing)) stop("missing required columns: ", paste(missing, collapse = ", "))

if (any(rowSums(d[, .(color_plain, color_yellow_orange, color_red_pink, color_blue_purple)], na.rm = TRUE) !=
        fifelse(is.na(d$color_trials), 0, d$color_trials))) stop("colour categories do not sum to color_trials")
if (any(rowSums(d[, .(form_open_generalized, form_tubular_trumpet,
                      form_zygomorphic_specialized, form_composite_brush)], na.rm = TRUE) !=
        fifelse(is.na(d$form_trials), 0, d$form_trials))) stop("form categories do not sum to form_trials")
if (any(d$n_animal_species + d$n_wind_species != d$animal_status_observed, na.rm = TRUE))
  stop("animal + wind species do not sum to observed pollen-vector trials")
if (any(d$wind_mixed_trials != d$animal_status_observed, na.rm = TRUE))
  stop("wind_mixed_trials and animal_status_observed disagree")
if (any(d$wind_mixed_successes != d$n_wind_species, na.rm = TRUE))
  stop("wind_mixed_successes and n_wind_species disagree")

d[, animal_pollinated_trials := animal_status_observed]
d[, animal_pollinated_successes := n_animal_species]
clip01 <- function(x, eps = 1e-5) pmin(pmax(as.numeric(x), eps), 1 - eps)
d[, logit_bombus_deficit := qlogis(clip01(bombus_deficit))]

regime_levels <- c(
  "northern_midlatitude", "northern_high_latitude",
  "tropical", "southern_extratropical"
)
d[, regime := factor(analysis_regime, levels = regime_levels)]
d[, spatial_id := as.integer(factor(spatial_block))]

scale_vars <- c(
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
  "bombus_deficit", "sc_share", "wind_mixed_share"
)
north_ref <- d[analysis_regime == "northern_midlatitude"]
scaling <- rbindlist(lapply(scale_vars, function(nm) {
  x <- north_ref[[nm]]
  data.table(variable = nm, center = mean(x, na.rm = TRUE), scale = sd(x, na.rm = TRUE))
}))
if (any(!is.finite(scaling$scale) | scaling$scale <= 0)) stop("invalid northern scaling")
for (i in seq_len(nrow(scaling))) {
  nm <- scaling$variable[i]
  d[[paste0("z_", nm)]] <- (d[[nm]] - scaling$center[i]) / scaling$scale[i]
}
fwrite(scaling, file.path(outdir, "northern_training_scaling.csv"))

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
    control.compute = compute, control.predictor = list(compute = TRUE), verbose = FALSE
  )
}

fit_gaussian <- function(response_col, rhs, dd) {
  x <- copy(dd)
  if (nrow(x) < 30) stop("too few rows for ", response_col, ": ", nrow(x))
  inla(
    as.formula(paste(response_col, "~", rhs)), family = "gaussian", data = x,
    control.compute = compute, control.predictor = list(compute = TRUE), verbose = FALSE
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

fixed_table <- function(fit, outcome_name, model_name, support_name, n_support, layer) {
  x <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
  x[, `:=`(
    outcome = outcome_name, model = model_name, support = support_name,
    n_support = n_support, analysis_layer = layer
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
    outcome = outcome_name, effect = effect_name, n_support = n_support,
    mean = mean(x, na.rm = TRUE), median = unname(q[[2]]),
    `0.025quant` = unname(q[[1]]), `0.975quant` = unname(q[[3]]),
    excludes_zero = unname(q[[1]] > 0 | q[[3]] < 0)
  )
}

complete_geo <- d[complete.cases(d[, ..geo_vars]) & !is.na(regime)]
base_complete <- d[complete.cases(d[, ..base_vars]) & !is.na(regime)]
north_geo <- complete_geo[regime == "northern_midlatitude"]
north_base <- base_complete[regime == "northern_midlatitude"]
if (nrow(north_base) < 30) stop("too few northern-midlatitude islands")

fits <- list()
score_rows <- list()
fixed_rows <- list()
support_rows <- list()

full_flora_outcomes <- list(
  animal_pollinated = c("animal_pollinated_successes", "animal_pollinated_trials"),
  wind_mixed = c("wind_mixed_successes", "wind_mixed_trials"),
  self_compatibility = c("sc_successes", "sc_trials")
)
for (outcome_name in names(full_flora_outcomes)) {
  success_col <- full_flora_outcomes[[outcome_name]][1]
  trials_col <- full_flora_outcomes[[outcome_name]][2]
  dd <- north_base[get(trials_col) > 0]
  support_rows[[paste0("full_flora__", outcome_name)]] <- data.table(
    analysis_layer = "full_flora", outcome = outcome_name,
    support = "maximum_equation_support", n_islands = nrow(dd)
  )
  models <- if (outcome_name == "self_compatibility") {
    list(
      F0_geo = geo_fit,
      F1_bombus = paste("z_bombus_deficit +", geo_fit),
      F2_wind = paste("z_wind_mixed_share +", geo_fit),
      F3_bombus_wind = paste("z_bombus_deficit + z_wind_mixed_share +", geo_fit)
    )
  } else {
    list(F0_geo = geo_fit, F1_bombus = paste("z_bombus_deficit +", geo_fit))
  }
  if (outcome_name == "self_compatibility") dd <- dd[is.finite(z_wind_mixed_share)]
  for (model_name in names(models)) {
    key <- paste("full_flora", outcome_name, model_name, sep = "__")
    fit <- fit_grouped(success_col, trials_col, models[[model_name]], dd)
    fits[[key]] <- fit
    s <- score_fit(fit)
    s[, `:=`(outcome = outcome_name, model = model_name,
             support = "matched_full_flora_support", n_support = nrow(dd),
             analysis_layer = "full_flora")]
    score_rows[[key]] <- s
    fixed_rows[[key]] <- fixed_table(
      fit, outcome_name, model_name, "matched_full_flora_support", nrow(dd), "full_flora"
    )
  }
}

bombus_support <- north_geo[is.finite(logit_bombus_deficit)]
bombus_fit <- fit_gaussian("logit_bombus_deficit", geo_fit, bombus_support)
fits[["bombus_deficit__B0_geo"]] <- bombus_fit
fixed_rows[["bombus_deficit__B0_geo"]] <- fixed_table(
  bombus_fit, "bombus_deficit", "B0_geo", "maximum_bombus_support",
  nrow(bombus_support), "bombus_channel"
)
score_rows[["bombus_deficit__B0_geo"]] <- cbind(
  score_fit(bombus_fit),
  data.table(outcome = "bombus_deficit", model = "B0_geo",
             support = "maximum_bombus_support", n_support = nrow(bombus_support),
             analysis_layer = "bombus_channel")
)
support_rows[["bombus_deficit"]] <- data.table(
  analysis_layer = "bombus_channel", outcome = "bombus_deficit",
  support = "maximum_bombus_support", n_islands = nrow(bombus_support)
)

trait_outcomes <- list(
  plain = c("color_plain", "color_trials"),
  yellow_orange = c("color_yellow_orange", "color_trials"),
  red_pink = c("color_red_pink", "color_trials"),
  blue_purple = c("color_blue_purple", "color_trials"),
  open_generalized = c("form_open_generalized", "form_trials"),
  tubular_trumpet = c("form_tubular_trumpet", "form_trials"),
  zygomorphic_specialized = c("form_zygomorphic_specialized", "form_trials"),
  composite_brush = c("form_composite_brush", "form_trials")
)

sc_path_support <- north_base[sc_trials > 0]
a_fit <- fits[["full_flora__self_compatibility__F1_bombus"]]
a_adjusted_fit <- fits[["full_flora__self_compatibility__F3_bombus_wind"]]
a_draw <- sample_fixed(a_fit, "z_bombus_deficit")
a_adjusted_draw <- sample_fixed(a_adjusted_fit, "z_bombus_deficit")
mediation_rows <- list()

for (outcome_name in names(trait_outcomes)) {
  success_col <- trait_outcomes[[outcome_name]][1]
  trials_col <- trait_outcomes[[outcome_name]][2]
  max_support <- north_base[get(trials_col) > 0]
  path_support <- max_support[is.finite(z_sc_share)]
  adjusted_support <- path_support[is.finite(z_wind_mixed_share)]
  if (min(nrow(max_support), nrow(path_support), nrow(adjusted_support)) < 30)
    stop("too few support rows for ", outcome_name)

  support_rows[[paste0(outcome_name, "__max")]] <- data.table(
    analysis_layer = "animal_pollinated_flower", outcome = outcome_name,
    support = "maximum_category_support", n_islands = nrow(max_support)
  )
  support_rows[[paste0(outcome_name, "__path")]] <- data.table(
    analysis_layer = "animal_pollinated_flower", outcome = outcome_name,
    support = "sc_path_support", n_islands = nrow(path_support)
  )
  support_rows[[paste0(outcome_name, "__adjusted")]] <- data.table(
    analysis_layer = "animal_pollinated_flower", outcome = outcome_name,
    support = "sc_wind_adjusted_support", n_islands = nrow(adjusted_support)
  )

  model_sets <- list(
    maximum_category_support = list(
      data = max_support,
      models = list(N0_geo_max = geo_fit, N1_bombus_max = paste("z_bombus_deficit +", geo_fit))
    ),
    sc_path_support = list(
      data = path_support,
      models = list(
        N0_geo_path = geo_fit,
        N1_bombus_path = paste("z_bombus_deficit +", geo_fit),
        N2_sc_path = paste("z_sc_share +", geo_fit),
        N3_direct_indirect_path = paste("z_bombus_deficit + z_sc_share +", geo_fit)
      )
    ),
    sc_wind_adjusted_support = list(
      data = adjusted_support,
      models = list(
        N0_geo_adjusted = geo_fit,
        N1_bombus_adjusted = paste("z_bombus_deficit +", geo_fit),
        N2_sc_adjusted = paste("z_sc_share +", geo_fit),
        N3_direct_indirect_adjusted = paste("z_bombus_deficit + z_sc_share +", geo_fit),
        N4_wind_adjusted = paste("z_bombus_deficit + z_sc_share + z_wind_mixed_share +", geo_fit)
      )
    )
  )

  for (support_name in names(model_sets)) {
    dd <- model_sets[[support_name]]$data
    models <- model_sets[[support_name]]$models
    for (model_name in names(models)) {
      key <- paste(outcome_name, model_name, sep = "__")
      fit <- fit_grouped(success_col, trials_col, models[[model_name]], dd)
      fits[[key]] <- fit
      s <- score_fit(fit)
      s[, `:=`(outcome = outcome_name, model = model_name, support = support_name,
               n_support = nrow(dd), analysis_layer = "animal_pollinated_flower")]
      score_rows[[key]] <- s
      fixed_rows[[key]] <- fixed_table(
        fit, outcome_name, model_name, support_name, nrow(dd), "animal_pollinated_flower"
      )
    }
  }

  trait_fit <- fits[[paste0(outcome_name, "__N3_direct_indirect_path")]]
  direct_draw <- sample_fixed(trait_fit, "z_bombus_deficit")
  b_draw <- sample_fixed(trait_fit, "z_sc_share")
  indirect_draw <- a_draw * b_draw
  total_draw <- direct_draw + indirect_draw
  mediation_rows[[paste0(outcome_name, "__unadjusted_a")]] <- summarize_draws(
    a_draw, outcome_name, "a_bombus_to_sc", nrow(sc_path_support)
  )
  mediation_rows[[paste0(outcome_name, "__unadjusted_b")]] <- summarize_draws(
    b_draw, outcome_name, "b_sc_to_category", nrow(path_support)
  )
  mediation_rows[[paste0(outcome_name, "__unadjusted_direct")]] <- summarize_draws(
    direct_draw, outcome_name, "direct_bombus_to_category", nrow(path_support)
  )
  mediation_rows[[paste0(outcome_name, "__unadjusted_indirect")]] <- summarize_draws(
    indirect_draw, outcome_name, "indirect_bombus_via_sc", nrow(path_support)
  )
  mediation_rows[[paste0(outcome_name, "__unadjusted_total")]] <- summarize_draws(
    total_draw, outcome_name, "total_bombus_effect", nrow(path_support)
  )

  adjusted_trait_fit <- fits[[paste0(outcome_name, "__N4_wind_adjusted")]]
  direct_adjusted_draw <- sample_fixed(adjusted_trait_fit, "z_bombus_deficit")
  b_adjusted_draw <- sample_fixed(adjusted_trait_fit, "z_sc_share")
  indirect_adjusted_draw <- a_adjusted_draw * b_adjusted_draw
  total_adjusted_draw <- direct_adjusted_draw + indirect_adjusted_draw
  mediation_rows[[paste0(outcome_name, "__adjusted_a")]] <- summarize_draws(
    a_adjusted_draw, outcome_name, "a_bombus_to_sc_adjusted_for_wind_mixed", nrow(sc_path_support)
  )
  mediation_rows[[paste0(outcome_name, "__adjusted_b")]] <- summarize_draws(
    b_adjusted_draw, outcome_name, "b_sc_to_category_adjusted_for_wind_mixed", nrow(adjusted_support)
  )
  mediation_rows[[paste0(outcome_name, "__adjusted_direct")]] <- summarize_draws(
    direct_adjusted_draw, outcome_name, "direct_bombus_to_category_adjusted_for_wind_mixed", nrow(adjusted_support)
  )
  mediation_rows[[paste0(outcome_name, "__adjusted_indirect")]] <- summarize_draws(
    indirect_adjusted_draw, outcome_name, "indirect_bombus_via_sc_adjusted_for_wind_mixed", nrow(adjusted_support)
  )
  mediation_rows[[paste0(outcome_name, "__adjusted_total")]] <- summarize_draws(
    total_adjusted_draw, outcome_name, "total_bombus_effect_adjusted_for_wind_mixed", nrow(adjusted_support)
  )
}

scores <- rbindlist(score_rows, fill = TRUE)
scores[, delta_log_cpo_from_best := log_cpo_sum - max(log_cpo_sum, na.rm = TRUE),
       by = .(analysis_layer, outcome, support)]
setorder(scores, analysis_layer, outcome, support, -log_cpo_sum)
fwrite(scores, file.path(outdir, "model_comparisons_all_layers.csv"))

fixed <- rbindlist(fixed_rows, fill = TRUE)
fixed[, excludes_zero := (`0.025quant` > 0 & `0.975quant` > 0) |
  (`0.025quant` < 0 & `0.975quant` < 0)]
fwrite(fixed, file.path(outdir, "fixed_effects_all_layers.csv"))
fwrite(rbindlist(support_rows, fill = TRUE), file.path(outdir, "equation_support_all_layers.csv"))
fwrite(rbindlist(mediation_rows, fill = TRUE),
       file.path(outdir, "direct_indirect_total_effects_all_categories.csv"))

bombus_audit_cols <- c(
  "island_id", "analysis_regime", "n_bombus_species_total", "n_bombus_species_scored",
  "n_bombus_species_unscored", "n_environmentally_compatible_species",
  "bombus_environmental_compatibility", "environmental_compatibility_max",
  "environmental_compatibility_mean", "environmental_extrapolation_share", "bombus_deficit"
)
fwrite(d[, ..bombus_audit_cols], file.path(outdir, "bombus_channel_component_audit.csv"))

guild_outcomes <- list(
  showy_alternative_guild = c("showy_alternative_guild_successes", "showy_alternative_guild_trials"),
  other_bee_guild = c("other_bee_guild_successes", "other_bee_guild_trials"),
  generalist_insect_guild = c("generalist_insect_guild_successes", "generalist_insect_guild_trials")
)
global_outcomes <- c(full_flora_outcomes, trait_outcomes, guild_outcomes)
regime_rows <- list()
regime_support_rows <- list()
regime_fits <- list()
for (outcome_name in names(global_outcomes)) {
  success_col <- global_outcomes[[outcome_name]][1]
  trials_col <- global_outcomes[[outcome_name]][2]
  dd <- copy(complete_geo[get(trials_col) > 0])
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
    as.formula(paste(success_col, "~", rhs)), family = "binomial", data = dd,
    Ntrials = dd[[trials_col]], control.compute = compute,
    control.predictor = list(compute = TRUE), verbose = FALSE
  )
  regime_fits[[outcome_name]] <- fit
  f <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
  f[, outcome := outcome_name]
  regime_rows[[outcome_name]] <- f
}
regime_effects <- rbindlist(regime_rows, fill = TRUE)
regime_effects[, excludes_zero := (`0.025quant` > 0 & `0.975quant` > 0) |
  (`0.025quant` < 0 & `0.975quant` < 0)]
fwrite(regime_effects, file.path(outdir, "cross_regime_effects_all_outcomes.csv"))
fwrite(
  regime_effects[grepl("z_log_distance_to_continent_km", parameter),
                 .(outcome, parameter, mean, `0.025quant`, `0.975quant`, excludes_zero)],
  file.path(outdir, "cross_regime_isolation_summary_all_outcomes.csv")
)
fwrite(rbindlist(regime_support_rows, fill = TRUE),
       file.path(outdir, "cross_regime_support_all_outcomes.csv"))

meta <- list(
  contract = "v2_all_data_layered_inla_v4",
  analysis_tier = unique(d$analysis_tier),
  engine = "INLA",
  regime_levels = regime_levels,
  layers = list(
    full_flora = names(full_flora_outcomes),
    bombus_channel = c("isolation_to_bombus_deficit", bombus_audit_cols[-1]),
    animal_pollinated_flower = names(trait_outcomes),
    global_descriptive_falsification = c(names(global_outcomes))
  ),
  flower_path_models = c(
    "unadjusted Bombus -> SC -> category",
    "wind/mixed-adjusted Bombus -> SC -> category"
  ),
  alternative_guild_role = "descriptive global outcomes only; not primary causal predictors",
  bombus_component_role = "fully exported and audited; collinear construction components are not simultaneously entered as causal covariates",
  support_rule = "maximum support per equation; model comparisons only within matched support",
  distance_rule = "distance_to_continent_km > 0 enforced by workflow",
  n_model_islands = nrow(d)
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(
  c(lapply(fits, summary), lapply(regime_fits, summary)),
  file = file.path(outdir, "model_summaries.txt")
)
