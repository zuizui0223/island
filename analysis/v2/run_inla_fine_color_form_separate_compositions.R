#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: run_inla_fine_color_form_separate_compositions.R <input.csv> <output_dir>")
infile <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260714)

d <- fread(infile)
colors <- c("white", "blue", "yellow", "red", "green", "cream", "rare_other")
forms <- c("tubular", "composite_head", "star", "bell_campanulate", "open_radial", "trumpet_funnel", "rare_other")

endpoints <- list(
  fine_color = list(categories = colors, prefix = "fine_color_", trials = "fine_color_trials"),
  fine_form = list(categories = forms, prefix = "fine_form_", trials = "fine_form_trials")
)

base_required <- c(
  "island_id", "analysis_regime", "spatial_block",
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
  "bombus_deficit", "sc_share", "wind_mixed_share"
)
required <- unique(c(
  base_required,
  unlist(lapply(endpoints, function(x) c(x$trials, paste0(x$prefix, x$categories))))
))
missing <- setdiff(required, names(d))
if (length(missing)) stop("missing fine composition columns: ", paste(missing, collapse = ", "))

for (endpoint_name in names(endpoints)) {
  spec <- endpoints[[endpoint_name]]
  cols <- paste0(spec$prefix, spec$categories)
  totals <- rowSums(d[, ..cols], na.rm = TRUE)
  expected <- fifelse(is.na(d[[spec$trials]]), 0, d[[spec$trials]])
  if (any(totals != expected)) stop(endpoint_name, " categories do not sum to trials")
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
if (any(!is.finite(scaling$scale) | scaling$scale <= 0)) stop("invalid fine-composition scaling")
for (i in seq_len(nrow(scaling))) {
  nm <- scaling$variable[i]
  d[[paste0("z_", nm)]] <- (d[[nm]] - scaling$center[i]) / scaling$scale[i]
}
d[, z_isolation_sq := z_log_distance_to_continent_km^2]
fwrite(scaling, file.path(outdir, "fine_separate_scaling.csv"))

complete_predictors <- c(
  "z_log_distance_to_continent_km", "z_isolation_sq", "z_log_island_area_km2",
  "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
  "z_bombus_deficit", "z_sc_share", "z_wind_mixed_share"
)
compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
all_scores <- list()
all_fixed <- list()
all_support <- list()

fit_endpoint <- function(endpoint_name, spec) {
  count_cols <- paste0(spec$prefix, spec$categories)
  support <- d[
    analysis_regime == "northern_midlatitude" &
      get(spec$trials) > 0 &
      complete.cases(d[, ..complete_predictors])
  ]
  if (nrow(support) < 30) stop("too few islands for ", endpoint_name)

  long <- melt(
    support,
    id.vars = c("island_id", "spatial_block", spec$trials, complete_predictors),
    measure.vars = count_cols,
    variable.name = "category",
    value.name = "count"
  )
  long[, category := factor(sub(paste0("^", spec$prefix), "", category), levels = spec$categories)]
  long[, island_comp_id := as.integer(factor(island_id))]
  long[, spatial_id := as.integer(factor(spatial_block))]
  long[, obs_id := seq_len(.N)]
  long[, log_trials := log(get(spec$trials))]

  geo_controls <- paste(
    "category:z_log_island_area_km2 +",
    "category:z_climate_pc1 + category:z_climate_pc2 +",
    "category:z_climate_pc3 + category:z_climate_pc4"
  )
  random_terms <- "+ f(island_comp_id, model='iid') + f(spatial_id, model='iid') + f(obs_id, model='iid')"
  linear_iso <- paste("0 + category + category:z_log_distance_to_continent_km +", geo_controls, random_terms)
  accel_iso <- paste(linear_iso, "+ category:z_isolation_sq")
  models <- list(
    I0_linear_isolation = linear_iso,
    I1_accelerating_isolation = accel_iso,
    I2_bombus = paste(accel_iso, "+ category:z_bombus_deficit"),
    I3_sc_wind = paste(accel_iso, "+ category:z_sc_share + category:z_wind_mixed_share"),
    I4_full = paste(accel_iso, "+ category:z_bombus_deficit + category:z_sc_share + category:z_wind_mixed_share")
  )

  scores <- list()
  fixed_rows <- list()
  for (model_name in names(models)) {
    fit <- inla(
      as.formula(paste("count ~", models[[model_name]], "+ offset(log_trials)")),
      family = "poisson", data = long,
      control.compute = compute,
      control.predictor = list(compute = TRUE),
      verbose = FALSE
    )
    cpo <- fit$cpo$cpo
    cpo <- cpo[is.finite(cpo) & cpo > 0]
    scores[[model_name]] <- data.table(
      endpoint = endpoint_name,
      model = model_name,
      n_islands = nrow(support),
      n_rows = nrow(long),
      log_cpo_sum = if (length(cpo)) sum(log(cpo)) else NA_real_,
      waic = fit$waic$waic,
      dic = fit$dic$dic
    )
    fx <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
    fx[, `:=`(endpoint = endpoint_name, model = model_name)]
    fixed_rows[[model_name]] <- fx
  }
  score_dt <- rbindlist(scores)
  score_dt[, delta_log_cpo_from_best := log_cpo_sum - max(log_cpo_sum, na.rm = TRUE)]
  setorder(score_dt, -log_cpo_sum)
  fixed_dt <- rbindlist(fixed_rows, fill = TRUE)
  fixed_dt[, excludes_zero := (`0.025quant` > 0 & `0.975quant` > 0) |
    (`0.025quant` < 0 & `0.975quant` < 0)]

  all_scores[[endpoint_name]] <<- score_dt
  all_fixed[[endpoint_name]] <<- fixed_dt
  all_support[[endpoint_name]] <<- data.table(
    endpoint = endpoint_name,
    n_islands = nrow(support),
    n_categories = length(spec$categories),
    categories = paste(spec$categories, collapse = ";")
  )
}

for (endpoint_name in names(endpoints)) fit_endpoint(endpoint_name, endpoints[[endpoint_name]])

scores <- rbindlist(all_scores, fill = TRUE)
fixed <- rbindlist(all_fixed, fill = TRUE)
support <- rbindlist(all_support, fill = TRUE)
fwrite(scores, file.path(outdir, "fine_separate_model_comparisons.csv"))
fwrite(fixed, file.path(outdir, "fine_separate_fixed_effects.csv"))
fwrite(support, file.path(outdir, "fine_separate_support.csv"))

key <- fixed[
  model == "I4_full" & grepl("category.*:(z_log_distance_to_continent_km|z_isolation_sq|z_bombus_deficit|z_sc_share|z_wind_mixed_share)$", parameter)
]
key[, predictor := fcase(
  grepl("z_log_distance_to_continent_km$", parameter), "isolation_linear",
  grepl("z_isolation_sq$", parameter), "isolation_acceleration",
  grepl("z_bombus_deficit$", parameter), "bombus_deficit",
  grepl("z_sc_share$", parameter), "self_compatibility",
  default = "wind_mixed"
)]
key[, category := sub("^category([^:]+):.*$", "\\1", parameter)]
fwrite(key, file.path(outdir, "fine_separate_key_effects.csv"))

meta <- list(
  model = "separate_fine_colour_and_form_compositions",
  primary_endpoints = c("fine_color", "fine_form"),
  joint_colour_form_primary = false,
  isolation = list(
    scale = "continuous standardized log distance to continent",
    acceleration_test = "category-specific quadratic term z_isolation_sq",
    interpretation = "linear term gives the local isolation gradient; quadratic term tests whether compositional change strengthens or bends at greater isolation"
  ),
  model_sequence = c(
    "I0_linear_isolation",
    "I1_accelerating_isolation",
    "I2_bombus",
    "I3_sc_wind",
    "I4_full"
  ),
  sensitivity_only = "colour-by-form joint composition"
)
write_json(meta, file.path(outdir, "fine_separate_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
