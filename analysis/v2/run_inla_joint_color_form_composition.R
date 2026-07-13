#!/usr/bin/env Rscript

# Jointly analyse the 4 x 4 colour-form cross-classification.
# The Poisson log-linear representation treats the 16 counts from each island
# as one composition through a shared island intercept. This avoids interpreting
# eight separate one-vs-rest binomial fits as independent responses.

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: run_inla_joint_color_form_composition.R <input.csv> <output_dir>")
}
infile <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260714)

d <- fread(infile)
colors <- c("plain", "yellow_orange", "red_pink", "blue_purple")
forms <- c(
  "open_generalized", "tubular_trumpet",
  "zygomorphic_specialized", "composite_brush"
)
joint_categories <- as.vector(outer(colors, forms, paste, sep = "__"))
joint_cols <- paste0("joint_", joint_categories)
required <- c(
  "island_id", "analysis_regime", "spatial_block", "joint_color_form_trials",
  joint_cols,
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
  "bombus_deficit", "sc_share", "wind_mixed_share"
)
missing <- setdiff(required, names(d))
if (length(missing)) stop("missing joint-composition columns: ", paste(missing, collapse = ", "))

joint_sum <- rowSums(d[, ..joint_cols], na.rm = TRUE)
expected_trials <- fifelse(is.na(d$joint_color_form_trials), 0, d$joint_color_form_trials)
if (any(joint_sum != expected_trials)) {
  stop("joint colour-form cells do not sum to joint_color_form_trials")
}

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
if (any(!is.finite(scaling$scale) | scaling$scale <= 0)) stop("invalid joint-model scaling")
for (i in seq_len(nrow(scaling))) {
  nm <- scaling$variable[i]
  d[[paste0("z_", nm)]] <- (d[[nm]] - scaling$center[i]) / scaling$scale[i]
}
fwrite(scaling, file.path(outdir, "joint_color_form_scaling.csv"))

predictors <- c(
  "z_log_distance_to_continent_km", "z_log_island_area_km2",
  "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
  "z_bombus_deficit", "z_sc_share", "z_wind_mixed_share"
)
support <- d[
  analysis_regime == "northern_midlatitude" &
    joint_color_form_trials > 0 &
    complete.cases(d[, ..predictors])
]
if (nrow(support) < 30) stop("too few northern islands for joint colour-form model")

long <- melt(
  support,
  id.vars = c("island_id", "spatial_block", "joint_color_form_trials", predictors),
  measure.vars = joint_cols,
  variable.name = "joint_category",
  value.name = "count"
)
long[, joint_category := factor(
  sub("^joint_", "", joint_category),
  levels = joint_categories
)]
long[, island_comp_id := as.integer(factor(island_id))]
long[, spatial_id := as.integer(factor(spatial_block))]
long[, obs_id := seq_len(.N)]
long[, log_trials := log(joint_color_form_trials)]

geo_terms <- paste(
  "joint_category:z_log_distance_to_continent_km +",
  "joint_category:z_log_island_area_km2 +",
  "joint_category:z_climate_pc1 + joint_category:z_climate_pc2 +",
  "joint_category:z_climate_pc3 + joint_category:z_climate_pc4"
)
base_terms <- paste(
  "0 + joint_category +", geo_terms,
  "+ f(island_comp_id, model='iid') + f(spatial_id, model='iid') + f(obs_id, model='iid')"
)
models <- list(
  J0_geo = base_terms,
  J1_bombus = paste(base_terms, "+ joint_category:z_bombus_deficit"),
  J2_sc_wind = paste(base_terms, "+ joint_category:z_sc_share + joint_category:z_wind_mixed_share"),
  J3_full = paste(
    base_terms,
    "+ joint_category:z_bombus_deficit + joint_category:z_sc_share +",
    "joint_category:z_wind_mixed_share"
  )
)
compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
fits <- list()
score_rows <- list()
fixed_rows <- list()
for (model_name in names(models)) {
  fit <- inla(
    as.formula(paste("count ~", models[[model_name]], "+ offset(log_trials)")),
    family = "poisson",
    data = long,
    control.compute = compute,
    control.predictor = list(compute = TRUE),
    verbose = FALSE
  )
  fits[[model_name]] <- fit
  cpo <- fit$cpo$cpo
  cpo <- cpo[is.finite(cpo) & cpo > 0]
  score_rows[[model_name]] <- data.table(
    model = model_name,
    n_islands = nrow(support),
    n_joint_rows = nrow(long),
    log_cpo_sum = if (length(cpo)) sum(log(cpo)) else NA_real_,
    waic = fit$waic$waic,
    dic = fit$dic$dic
  )
  fx <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
  fx[, model := model_name]
  fixed_rows[[model_name]] <- fx
}

scores <- rbindlist(score_rows)
scores[, delta_log_cpo_from_best := log_cpo_sum - max(log_cpo_sum, na.rm = TRUE)]
setorder(scores, -log_cpo_sum)
fwrite(scores, file.path(outdir, "joint_color_form_model_comparisons.csv"))

fixed <- rbindlist(fixed_rows, fill = TRUE)
fixed[, excludes_zero := (`0.025quant` > 0 & `0.975quant` > 0) |
  (`0.025quant` < 0 & `0.975quant` < 0)]
fwrite(fixed, file.path(outdir, "joint_color_form_fixed_effects.csv"))

joint_effects <- fixed[
  model == "J3_full" & grepl("joint_category.*:(z_bombus_deficit|z_sc_share|z_wind_mixed_share)$", parameter)
]
joint_effects[, predictor := fifelse(
  grepl("z_bombus_deficit$", parameter), "bombus_deficit",
  fifelse(grepl("z_sc_share$", parameter), "self_compatibility", "wind_mixed")
)]
joint_effects[, joint_category := sub(
  "^joint_category([^:]+):.*$", "\\1", parameter
)]
fwrite(joint_effects, file.path(outdir, "joint_color_form_key_effects.csv"))

meta <- list(
  model = "joint_4x4_colour_form_poisson_loglinear",
  interpretation = paste(
    "Sixteen colour-form cells are fitted together with a shared island effect;",
    "coefficients describe compositional reallocation, not independent binary outcomes."
  ),
  n_islands = nrow(support),
  n_joint_rows = nrow(long),
  categories = joint_categories,
  support = "northern_midlatitude islands with both colour and form resolved and complete predictors"
)
write_json(meta, file.path(outdir, "joint_color_form_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(lapply(fits, summary), file = file.path(outdir, "joint_color_form_model_summaries.txt"))
