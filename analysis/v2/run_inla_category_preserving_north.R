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
required <- c(
  "analysis_regime", "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
  "bombus_deficit", "sc_trials", "sc_successes", "sc_share", "spatial_block",
  "color_trials", "color_plain", "color_yellow_orange", "color_red_pink",
  "form_trials", "form_open_generalized", "form_tubular_trumpet"
)
missing <- setdiff(required, names(d))
if (length(missing)) stop("missing required columns: ", paste(missing, collapse = ", "))

# Engine split:
#   INLA owns all category-preserving flower colour/form analyses.
#   brms is not used for category models; it may only be used later for compact
#   scalar/aggregate replication of simple Bombus/SC paths.
# Primary northern mechanism:
#   Bombus deficit -> dull/generalized flowers                 (direct)
#   Bombus deficit -> self-compatibility -> dull/generalized   (indirect)
# Tropical and southern-extratropical regimes are falsification domains.

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
  nm <- scaling$variable[i]
  d[[paste0("z_", nm)]] <- (d[[nm]] - scaling$center[i]) / scaling$scale[i]
}
fwrite(scaling, file.path(outdir, "northern_training_scaling.csv"))
d[, spatial_id := as.integer(factor(spatial_block))]
d[, regime := factor(
  analysis_regime,
  levels = c("northern_midlatitude", "tropical", "southern_extratropical")
)]

geo_fixed <- paste(
  "z_log_distance_to_continent_km * z_log_island_area_km2 +",
  "z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4"
)
geo_fit <- paste(geo_fixed, "+ f(spatial_id, model='iid')")
compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)

fit_grouped <- function(success_col, trials_col, rhs, dd) {
  x <- copy(dd)
  x[, obs_id := seq_len(.N)]
  inla(
    as.formula(paste(success_col, "~", rhs, "+ f(obs_id, model='iid')")),
    family = "binomial", data = x, Ntrials = x[[trials_col]],
    control.compute = compute, control.predictor = list(compute = TRUE),
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

fixed_table <- function(fit, outcome_name, model_name) {
  x <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
  x[, `:=`(outcome = outcome_name, model = model_name)]
  x
}

core_z <- c(
  "z_log_distance_to_continent_km", "z_log_island_area_km2",
  "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
  "z_bombus_deficit", "z_sc_share"
)
complete <- d[complete.cases(d[, ..core_z]) & !is.na(regime)]
north <- complete[regime == "northern_midlatitude"]
if (nrow(north) < 30) stop("Too few complete northern-midlatitude islands")

# Common support keeps direct and indirect paths comparable.
mediation_support <- north[sc_trials > 0 & color_trials > 0 & form_trials > 0]
if (nrow(mediation_support) < 30) stop("Too few islands on common mediation support")
fwrite(
  mediation_support[, .(island_id, analysis_regime, sc_trials, color_trials, form_trials)],
  file.path(outdir, "northern_mediation_support.csv")
)

sc_models <- list(
  N0_geo = geo_fit,
  N1_bombus = paste("z_bombus_deficit +", geo_fit)
)
trait_models <- list(
  N0_geo = geo_fit,
  N1_bombus = paste("z_bombus_deficit +", geo_fit),
  N2_sc = paste("z_sc_share +", geo_fit),
  N3_direct_indirect = paste("z_bombus_deficit + z_sc_share +", geo_fit)
)
trait_outcomes <- list(
  plain_flower = c("color_plain", "color_trials"),
  generalized_flower = c("form_open_generalized", "form_trials")
)

fits <- list(); score_rows <- list(); fixed_rows <- list()
for (model_name in names(sc_models)) {
  key <- paste("self_compatibility", model_name, sep = "__")
  fit <- fit_grouped("sc_successes", "sc_trials", sc_models[[model_name]], mediation_support)
  fits[[key]] <- fit
  s <- score_fit(fit)
  s[, `:=`(outcome = "self_compatibility", model = model_name)]
  score_rows[[key]] <- s
  fixed_rows[[key]] <- fixed_table(fit, "self_compatibility", model_name)
}
for (outcome_name in names(trait_outcomes)) {
  success_col <- trait_outcomes[[outcome_name]][1]
  trials_col <- trait_outcomes[[outcome_name]][2]
  for (model_name in names(trait_models)) {
    key <- paste(outcome_name, model_name, sep = "__")
    fit <- fit_grouped(success_col, trials_col, trait_models[[model_name]], mediation_support)
    fits[[key]] <- fit
    s <- score_fit(fit)
    s[, `:=`(outcome = outcome_name, model = model_name)]
    score_rows[[key]] <- s
    fixed_rows[[key]] <- fixed_table(fit, outcome_name, model_name)
  }
}

scores <- rbindlist(score_rows, fill = TRUE)
scores[, delta_log_cpo_from_best := log_cpo_sum - max(log_cpo_sum, na.rm = TRUE), by = outcome]
setorder(scores, outcome, -log_cpo_sum)
fwrite(scores, file.path(outdir, "northern_model_comparisons.csv"))

fixed <- rbindlist(fixed_rows, fill = TRUE)
fixed[, excludes_zero := (`0.025quant` > 0 & `0.975quant` > 0) |
  (`0.025quant` < 0 & `0.975quant` < 0)]
fwrite(fixed, file.path(outdir, "northern_fixed_effects.csv"))

sample_fixed <- function(fit, parameter, n = 50000L) {
  if (!parameter %in% names(fit$marginals.fixed)) stop("missing parameter: ", parameter)
  inla.rmarginal(n, fit$marginals.fixed[[parameter]])
}
summarize_draws <- function(x, outcome_name, effect_name) {
  q <- quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)
  data.table(
    outcome = outcome_name, effect = effect_name, mean = mean(x, na.rm = TRUE),
    median = unname(q[[2]]), `0.025quant` = unname(q[[1]]),
    `0.975quant` = unname(q[[3]]),
    excludes_zero = unname(q[[1]] > 0 | q[[3]] < 0)
  )
}

a_fit <- fits[["self_compatibility__N1_bombus"]]
a_draw <- sample_fixed(a_fit, "z_bombus_deficit")
mediation_rows <- list()
for (outcome_name in names(trait_outcomes)) {
  trait_fit <- fits[[paste0(outcome_name, "__N3_direct_indirect")]]
  direct_draw <- sample_fixed(trait_fit, "z_bombus_deficit")
  b_draw <- sample_fixed(trait_fit, "z_sc_share")
  indirect_draw <- a_draw * b_draw
  total_draw <- direct_draw + indirect_draw
  mediation_rows[[paste0(outcome_name, "_a")]] <- summarize_draws(a_draw, outcome_name, "a_bombus_to_self_compatibility")
  mediation_rows[[paste0(outcome_name, "_b")]] <- summarize_draws(b_draw, outcome_name, "b_self_compatibility_to_trait")
  mediation_rows[[paste0(outcome_name, "_direct")]] <- summarize_draws(direct_draw, outcome_name, "direct_bombus_to_trait")
  mediation_rows[[paste0(outcome_name, "_indirect")]] <- summarize_draws(indirect_draw, outcome_name, "indirect_bombus_via_self_compatibility")
  mediation_rows[[paste0(outcome_name, "_total")]] <- summarize_draws(total_draw, outcome_name, "total_bombus_effect")
}
mediation_effects <- rbindlist(mediation_rows)
fwrite(mediation_effects, file.path(outdir, "northern_direct_indirect_total_effects.csv"))

key_effects <- fixed[
  (outcome == "self_compatibility" & model == "N1_bombus" & parameter == "z_bombus_deficit") |
  (outcome %in% names(trait_outcomes) & model == "N3_direct_indirect" &
    parameter %in% c("z_bombus_deficit", "z_sc_share"))
]
fwrite(key_effects, file.path(outdir, "northern_bombus_syndrome_effects.csv"))

# INLA category-preserving falsification analysis. Each retained category is modelled
# against its category total on the same engine; no brms multinomial/category fit.
falsification_outcomes <- list(
  self_compatibility = c("sc_successes", "sc_trials"),
  plain_flower = c("color_plain", "color_trials"),
  generalized_flower = c("form_open_generalized", "form_trials"),
  yellow_orange = c("color_yellow_orange", "color_trials"),
  red_pink = c("color_red_pink", "color_trials"),
  tubular_trumpet = c("form_tubular_trumpet", "form_trials")
)
regime_rows <- list(); regime_fits <- list(); regime_support_rows <- list()
for (outcome_name in names(falsification_outcomes)) {
  success_col <- falsification_outcomes[[outcome_name]][1]
  trials_col <- falsification_outcomes[[outcome_name]][2]
  dd <- copy(complete[get(trials_col) > 0])
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
fwrite(regime_effects, file.path(outdir, "cross_regime_isolation_effects.csv"))
isolation_terms <- regime_effects[
  grepl("z_log_distance_to_continent_km", parameter),
  .(outcome, parameter, mean, `0.025quant`, `0.975quant`, excludes_zero)
]
fwrite(isolation_terms, file.path(outdir, "cross_regime_isolation_summary.csv"))
fwrite(rbindlist(regime_support_rows), file.path(outdir, "cross_regime_support.csv"))

meta <- list(
  contract = "v2_northern_bombus_direct_indirect_with_global_falsification_v2",
  analysis_tier = unique(d$analysis_tier),
  primary_engine = "INLA",
  category_engine = "INLA",
  bayesian_role = "optional scalar/aggregate replication only; no category models",
  primary_domain = "northern_midlatitude",
  primary_claim = "isolation-associated self-compatibility and dull/generalized floral composition are tested as a northern Bombus-linked syndrome",
  causal_decomposition = list(
    direct = "Bombus deficit -> plain/generalized floral composition",
    indirect = "Bombus deficit -> self-compatibility -> plain/generalized floral composition",
    total = "direct + indirect"
  ),
  mediation_support = "common northern island support with SC, colour, and floral-form trials",
  mediation_uncertainty = "product-of-INLA-posterior-marginals approximation across separately fitted component models",
  falsification = "tropical and southern-extratropical regimes test whether northern isolation slopes weaken, disappear, or reverse; conspicuous categories are retained in INLA as counter-patterns",
  excluded_from_primary_mechanism = "alternative-pollinator covariates and Bombus-by-alternative interactions",
  n_common_northern_mediation_islands = nrow(mediation_support),
  n_complete_all_regimes = nrow(complete)
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(
  c(lapply(fits, summary), lapply(regime_fits, summary)),
  file = file.path(outdir, "model_summaries.txt")
)
