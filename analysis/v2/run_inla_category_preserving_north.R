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

# Simple v2 question:
# 1) In the northern-midlatitude Bombus domain, does isolation predict more
#    self-compatibility and plainer/generalized flowers?
# 2) Does Bombus deficit improve that northern explanation?
# 3) Do those isolation relationships fail outside the northern domain?
# Alternative-pollinator covariates are deliberately excluded from the
# northern mechanism; tropical/southern regimes are falsification domains.

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
d[, regime := factor(analysis_regime,
  levels = c("northern_midlatitude", "tropical", "southern_extratropical"))]

geo_fixed <- paste(
  "z_log_distance_to_continent_km * z_log_island_area_km2 +",
  "z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4"
)
geo_fit <- paste(geo_fixed, "+ f(spatial_id, model='iid')")

outcomes <- list(
  self_compatibility = c("sc_successes", "sc_trials"),
  plain_flower = c("color_plain", "color_trials"),
  generalized_flower = c("form_open_generalized", "form_trials")
)
model_rhs <- list(
  self_compatibility = list(
    N0_geo = geo_fit,
    N1_bombus = paste("z_bombus_deficit +", geo_fit)
  ),
  plain_flower = list(
    N0_geo = geo_fit,
    N1_bombus = paste("z_bombus_deficit +", geo_fit),
    N2_sc = paste("z_sc_share +", geo_fit),
    N3_bombus_sc = paste("z_bombus_deficit + z_sc_share +", geo_fit)
  ),
  generalized_flower = list(
    N0_geo = geo_fit,
    N1_bombus = paste("z_bombus_deficit +", geo_fit),
    N2_sc = paste("z_sc_share +", geo_fit),
    N3_bombus_sc = paste("z_bombus_deficit + z_sc_share +", geo_fit)
  )
)

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

core_z <- c(
  "z_log_distance_to_continent_km", "z_log_island_area_km2",
  "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
  "z_bombus_deficit", "z_sc_share"
)
complete <- d[complete.cases(d[, ..core_z]) & !is.na(regime)]
north <- complete[regime == "northern_midlatitude"]
if (nrow(north) < 30) stop("Too few complete northern-midlatitude islands")

fits <- list(); score_rows <- list(); fixed_rows <- list(); support_rows <- list()
for (outcome in names(outcomes)) {
  success_col <- outcomes[[outcome]][1]
  trials_col <- outcomes[[outcome]][2]
  dd <- north[get(trials_col) > 0]
  bad <- !is.finite(dd[[success_col]]) | !is.finite(dd[[trials_col]]) |
    dd[[success_col]] < 0 | dd[[trials_col]] <= 0 | dd[[success_col]] > dd[[trials_col]]
  if (any(bad)) stop("invalid grouped counts for ", outcome)
  support_rows[[outcome]] <- data.table(outcome = outcome, n_northern_islands = nrow(dd))
  for (model in names(model_rhs[[outcome]])) {
    key <- paste(outcome, model, sep = "__")
    fit <- fit_grouped(success_col, trials_col, model_rhs[[outcome]][[model]], dd)
    fits[[key]] <- fit
    s <- score_fit(fit); s[, `:=`(outcome = outcome, model = model)]; score_rows[[key]] <- s
    f <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
    f[, `:=`(outcome = outcome, model = model)]; fixed_rows[[key]] <- f
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
key_effects <- fixed[
  (outcome == "self_compatibility" & model == "N1_bombus" & parameter == "z_bombus_deficit") |
  (outcome %in% c("plain_flower", "generalized_flower") & model == "N3_bombus_sc" &
    parameter %in% c("z_bombus_deficit", "z_sc_share"))
]
fwrite(key_effects, file.path(outdir, "northern_bombus_syndrome_effects.csv"))
fwrite(rbindlist(support_rows), file.path(outdir, "northern_support.csv"))

# Falsification across regimes. The first three outcomes test whether the
# northern isolation syndrome is geographically restricted. The last three
# are descriptive counter-patterns for conspicuous/specialized flowers only.
falsification_outcomes <- list(
  self_compatibility = c("sc_successes", "sc_trials"),
  plain_flower = c("color_plain", "color_trials"),
  generalized_flower = c("form_open_generalized", "form_trials"),
  yellow_orange = c("color_yellow_orange", "color_trials"),
  red_pink = c("color_red_pink", "color_trials"),
  tubular_trumpet = c("form_tubular_trumpet", "form_trials")
)
regime_rows <- list(); regime_fits <- list()
for (outcome in names(falsification_outcomes)) {
  success_col <- falsification_outcomes[[outcome]][1]
  trials_col <- falsification_outcomes[[outcome]][2]
  dd <- complete[get(trials_col) > 0]
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
  regime_fits[[outcome]] <- fit
  f <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
  f[, outcome := outcome]; regime_rows[[outcome]] <- f
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

meta <- list(
  contract = "v2_simple_northern_bombus_syndrome_falsification_v1",
  analysis_tier = unique(d$analysis_tier),
  primary_claim = "isolation-associated increases in self-compatibility and plain/generalized floral composition are tested as a northern-midlatitude Bombus-linked syndrome",
  northern_mechanism = "geography -> Bombus deficit -> self-compatibility and plain/generalized floral composition",
  falsification = "the same isolation slopes are estimated in tropical and southern-extratropical regimes; yellow/orange, red/pink and tubular/trumpet categories are descriptive counter-patterns only",
  excluded_from_primary_mechanism = "alternative-pollinator covariates and Bombus-by-alternative interactions",
  n_complete_northern_midlatitude = nrow(north),
  n_complete_all_regimes = nrow(complete)
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(c(lapply(fits, summary), lapply(regime_fits, summary)),
  file = file.path(outdir, "model_summaries.txt"))
