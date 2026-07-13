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
  "island_id", "analysis_regime", "log_distance_to_continent_km",
  "log_island_area_km2", "climate_pc1", "climate_pc2", "climate_pc3",
  "climate_pc4", "bombus_deficit", "sc_share", "spatial_block",
  "showy_alt_guild_share", "other_bee_guild_share_1",
  "generalist_insect_guild_share_1", "color_trials", "color_plain",
  "color_yellow_orange", "color_red_pink", "color_blue_purple",
  "form_trials", "form_open_generalized", "form_tubular_trumpet",
  "form_zygomorphic_specialized", "form_composite_brush"
)
missing <- setdiff(required, names(d))
if (length(missing)) stop("missing required columns: ", paste(missing, collapse = ", "))

# Freeze scaling on the northern-midlatitude training domain, then apply those
# same transformations to every regime for honest out-of-domain projection.
raw_predictors <- c(
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
  "bombus_deficit", "sc_share", "showy_alt_guild_share",
  "other_bee_guild_share_1", "generalist_insect_guild_share_1"
)
north_ref <- d[analysis_regime == "northern_midlatitude"]
scaling <- rbindlist(lapply(raw_predictors, function(nm) {
  x <- north_ref[[nm]]
  data.table(variable = nm, center = mean(x, na.rm = TRUE), scale = sd(x, na.rm = TRUE))
}))
if (any(!is.finite(scaling$scale) | scaling$scale <= 0)) {
  stop("non-finite or zero northern training scale for: ", paste(scaling[!is.finite(scale) | scale <= 0, variable], collapse = ", "))
}
for (i in seq_len(nrow(scaling))) {
  nm <- scaling$variable[i]
  d[[paste0("z_", nm)]] <- (d[[nm]] - scaling$center[i]) / scaling$scale[i]
}
fwrite(scaling, file.path(outdir, "northern_training_scaling.csv"))
d[, spatial_id := as.integer(factor(spatial_block))]

generic_geo <- paste(
  "z_log_distance_to_continent_km * z_log_island_area_km2 +",
  "z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4"
)
alt_terms <- paste(
  "z_showy_alt_guild_share + z_other_bee_guild_share_1 +",
  "z_generalist_insect_guild_share_1"
)
replacement_interactions <- paste(
  "z_bombus_deficit:z_showy_alt_guild_share +",
  "z_bombus_deficit:z_other_bee_guild_share_1 +",
  "z_bombus_deficit:z_generalist_insect_guild_share_1"
)
fixed_rhs <- list(
  M0_geo = generic_geo,
  M1_bombus = paste("z_bombus_deficit +", generic_geo),
  M2_sc = paste("z_sc_share +", generic_geo),
  M3_bombus_sc = paste("z_bombus_deficit + z_sc_share +", generic_geo),
  M4_alt = paste(alt_terms, "+", generic_geo),
  M5_bombus_alt = paste("z_bombus_deficit +", alt_terms, "+", generic_geo),
  M6_sc_alt = paste("z_sc_share +", alt_terms, "+", generic_geo),
  M7_full = paste("z_bombus_deficit + z_sc_share +", alt_terms, "+", generic_geo),
  M8_replacement = paste(
    "z_bombus_deficit + z_sc_share +", alt_terms, "+",
    replacement_interactions, "+", generic_geo
  )
)
fit_rhs <- lapply(fixed_rhs, function(x) paste(x, "+ f(spatial_id, model='iid')"))

specs <- list(
  color_plain = c("color_plain", "color_trials"),
  color_yellow_orange = c("color_yellow_orange", "color_trials"),
  color_red_pink = c("color_red_pink", "color_trials"),
  color_blue_purple = c("color_blue_purple", "color_trials"),
  form_open_generalized = c("form_open_generalized", "form_trials"),
  form_tubular_trumpet = c("form_tubular_trumpet", "form_trials"),
  form_zygomorphic_specialized = c("form_zygomorphic_specialized", "form_trials"),
  form_composite_brush = c("form_composite_brush", "form_trials")
)

all_z <- c(
  "z_log_distance_to_continent_km", "z_log_island_area_km2",
  "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
  "z_bombus_deficit", "z_sc_share", "z_showy_alt_guild_share",
  "z_other_bee_guild_share_1", "z_generalist_insect_guild_share_1"
)
complete_full <- d[complete.cases(d[, ..all_z])]
north <- complete_full[analysis_regime == "northern_midlatitude"]
if (nrow(north) < 30) stop("Too few complete northern-midlatitude islands")

compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
fit_grouped <- function(success_col, trials_col, rhs, dd) {
  x <- copy(dd)
  success <- x[[success_col]]
  trials <- x[[trials_col]]
  bad <- !is.finite(success) | !is.finite(trials) | success < 0 | trials <= 0 |
    success > trials | success != floor(success) | trials != floor(trials)
  if (any(bad)) stop("Invalid grouped counts for ", success_col, "/", trials_col)
  x[, obs_id := seq_len(.N)]
  inla(
    as.formula(paste(success_col, "~", rhs, "+ f(obs_id, model='iid')")),
    family = "binomial", data = x, Ntrials = trials,
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

# Population-level projection deliberately excludes spatial and observation-level
# random effects. This tests transfer of the northern fixed-effect mechanism itself.
project_fixed <- function(fit, rhs, dd) {
  mm <- model.matrix(as.formula(paste("~", rhs)), data = dd)
  beta <- fit$summary.fixed$mean
  names(beta) <- rownames(fit$summary.fixed)
  common <- intersect(colnames(mm), names(beta))
  missing_coef <- setdiff(colnames(mm), names(beta))
  if (length(missing_coef)) stop("projection coefficients missing: ", paste(missing_coef, collapse = ", "))
  eta <- as.numeric(mm[, common, drop = FALSE] %*% beta[common])
  plogis(eta)
}
projection_metrics <- function(success, trials, p) {
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  observed <- success / trials
  ll <- dbinom(success, size = trials, prob = p, log = TRUE)
  data.table(
    n_islands = length(success),
    n_trials = sum(trials),
    log_score_sum = sum(ll),
    log_score_per_trial = sum(ll) / sum(trials),
    brier_share = mean((observed - p)^2),
    mean_observed_share = mean(observed),
    mean_predicted_share = mean(p)
  )
}

fits <- list()
score_rows <- list()
fixed_rows <- list()
support_rows <- list()
projection_rows <- list()

for (outcome in names(specs)) {
  success_col <- specs[[outcome]][1]
  trials_col <- specs[[outcome]][2]
  train <- north[get(trials_col) > 0]
  global_eval <- complete_full[
    analysis_regime %in% c("northern_midlatitude", "tropical", "southern_extratropical") &
      get(trials_col) > 0
  ]
  support_rows[[outcome]] <- data.table(
    outcome = outcome,
    n_train_northern = nrow(train),
    n_eval_northern = nrow(global_eval[analysis_regime == "northern_midlatitude"]),
    n_eval_tropical = nrow(global_eval[analysis_regime == "tropical"]),
    n_eval_southern_extratropical = nrow(global_eval[analysis_regime == "southern_extratropical"])
  )

  for (model in names(fit_rhs)) {
    key <- paste(outcome, model, sep = "__")
    fit <- fit_grouped(success_col, trials_col, fit_rhs[[model]], train)
    fits[[key]] <- fit

    s <- score_fit(fit)
    s[, `:=`(outcome = outcome, model = model)]
    score_rows[[key]] <- s

    f <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
    f[, `:=`(outcome = outcome, model = model)]
    fixed_rows[[key]] <- f

    for (regime in c("northern_midlatitude", "tropical", "southern_extratropical")) {
      eval <- global_eval[analysis_regime == regime]
      if (!nrow(eval)) next
      p <- project_fixed(fit, fixed_rhs[[model]], eval)
      pm <- projection_metrics(eval[[success_col]], eval[[trials_col]], p)
      pm[, `:=`(outcome = outcome, model = model, regime = regime)]
      projection_rows[[paste(key, regime, sep = "__")]] <- pm
    }
  }
}

scores <- rbindlist(score_rows, fill = TRUE)
scores[, delta_log_cpo_from_best := log_cpo_sum - max(log_cpo_sum, na.rm = TRUE), by = outcome]
scores[, best_model := model[which.max(log_cpo_sum)], by = outcome]
setorder(scores, outcome, -log_cpo_sum)
fwrite(scores, file.path(outdir, "category_model_comparisons.csv"))

fixed <- rbindlist(fixed_rows, fill = TRUE)
setcolorder(fixed, c("outcome", "model", "parameter"))
fwrite(fixed, file.path(outdir, "category_fixed_effects.csv"))

replacement_effects <- fixed[
  model == "M8_replacement" & parameter %in% c(
    "z_bombus_deficit", "z_sc_share", "z_showy_alt_guild_share",
    "z_other_bee_guild_share_1", "z_generalist_insect_guild_share_1",
    "z_bombus_deficit:z_showy_alt_guild_share",
    "z_bombus_deficit:z_other_bee_guild_share_1",
    "z_bombus_deficit:z_generalist_insect_guild_share_1"
  )
]
replacement_effects[, excludes_zero :=
  (`0.025quant` > 0 & `0.975quant` > 0) |
    (`0.025quant` < 0 & `0.975quant` < 0)
]
replacement_effects[, direction := fifelse(mean > 0, "positive", "negative")]
setorder(replacement_effects, outcome, parameter)
fwrite(replacement_effects, file.path(outdir, "replacement_category_channel_effects.csv"))

best <- scores[, .SD[which.max(log_cpo_sum)], by = outcome]
best <- best[, .(outcome, best_model = model, log_cpo_sum, waic, dic, n_cpo)]
fwrite(best, file.path(outdir, "category_best_models.csv"))

projection <- rbindlist(projection_rows, fill = TRUE)
# Projection failure is defined relative to each northern-trained model's score
# in its own northern training regime. Negative values mean worse transfer.
projection[, northern_reference := log_score_per_trial[regime == "northern_midlatitude"][1], by = .(outcome, model)]
projection[, transfer_delta_log_score_per_trial := log_score_per_trial - northern_reference]
projection[, projection_failure := regime != "northern_midlatitude" & transfer_delta_log_score_per_trial < 0]
setorder(projection, outcome, model, regime)
fwrite(projection, file.path(outdir, "cross_regime_projection_performance.csv"))

projection_summary <- projection[, .(
  n_outcome_models = .N,
  mean_transfer_delta = mean(transfer_delta_log_score_per_trial, na.rm = TRUE),
  median_transfer_delta = median(transfer_delta_log_score_per_trial, na.rm = TRUE),
  share_projection_failure = mean(projection_failure, na.rm = TRUE)
), by = regime]
fwrite(projection_summary, file.path(outdir, "cross_regime_projection_summary.csv"))

support <- rbindlist(support_rows)
fwrite(support, file.path(outdir, "category_support.csv"))

meta <- list(
  contract = "v2_inla_category_preserving_replacement_projection_v3",
  analysis_tier = unique(d$analysis_tier),
  training_regime = "northern_midlatitude",
  evaluation_regimes = c("northern_midlatitude", "tropical", "southern_extratropical"),
  support_policy = "for each category, M0-M8 use the same maximum complete-case northern training support; projection uses all complete islands with category trials in each evaluation regime",
  category_policy = "one-vs-rest grouped binomial-logit-normal; no color/form binary collapse",
  channel_policy = "compare geography, Bombus deficit, reproductive assurance, alternative-pollinator replacement, additive combinations, and Bombus-deficit by replacement interactions",
  projection_policy = "northern scaling and fixed effects are frozen before population-level projection; spatial and observation-level random effects are excluded from out-of-domain prediction",
  alternative_pollinator_terms = c(
    "showy_alt_guild_share", "other_bee_guild_share_1", "generalist_insect_guild_share_1"
  ),
  note = "showy_alt_guild_share is the currently available combined showy-alternative functional channel; this artifact does not yet separate bird from butterfly/moth evidence",
  n_complete_northern_midlatitude = nrow(north),
  n_complete_global_evaluation = nrow(complete_full),
  n_spatial_blocks_north = uniqueN(north$spatial_id)
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)

capture.output(lapply(fits, summary), file = file.path(outdir, "model_summaries.txt"))
