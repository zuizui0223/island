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
z <- function(x) as.numeric(scale(x))
for (nm in c("log_distance_to_continent_km", "log_island_area_km2", "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4")) {
  d[[paste0("z_", nm)]] <- z(d[[nm]])
}
d[, z_bombus_deficit := z(bombus_deficit)]
d[, z_sc_share := z(sc_share)]
d[, spatial_id := as.integer(factor(spatial_block))]

geo_terms <- paste(
  "z_log_distance_to_continent_km * z_log_island_area_km2 +",
  "z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4 +",
  "f(spatial_id, model='iid')"
)
model_rhs <- list(
  M0_geo = geo_terms,
  M1_bombus = paste("z_bombus_deficit +", geo_terms),
  M2_sc = paste("z_sc_share +", geo_terms),
  M3_bombus_sc = paste("z_bombus_deficit + z_sc_share +", geo_terms)
)

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

compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
fit_grouped <- function(success_col, trials_col, rhs, dd) {
  x <- copy(dd)
  x[, obs_id := seq_len(.N)]
  fit <- inla(
    as.formula(paste(success_col, "~", rhs, "+ f(obs_id, model='iid')")),
    family = "binomial",
    data = x,
    Ntrials = x[[trials_col]],
    control.compute = compute,
    control.predictor = list(compute = TRUE),
    verbose = FALSE
  )
  fit
}
score_fit <- function(fit) {
  cpo <- fit$cpo$cpo
  cpo <- cpo[is.finite(cpo) & cpo > 0]
  data.table(
    log_cpo_sum = if (length(cpo)) sum(log(cpo)) else NA_real_,
    n_cpo = length(cpo),
    waic = fit$waic$waic,
    dic = fit$dic$dic
  )
}

north <- d[
  analysis_regime == "northern_midlatitude" &
    complete.cases(
      z_log_distance_to_continent_km, z_log_island_area_km2,
      z_climate_pc1, z_climate_pc2, z_climate_pc3, z_climate_pc4,
      z_bombus_deficit, z_sc_share, spatial_id
    )
]
if (nrow(north) < 30) stop("Too few complete northern-midlatitude islands")

fits <- list()
score_rows <- list()
fixed_rows <- list()
support_rows <- list()

for (outcome in names(specs)) {
  success_col <- specs[[outcome]][1]
  trials_col <- specs[[outcome]][2]
  dd <- north[get(trials_col) > 0]
  bad <- !is.finite(dd[[success_col]]) | !is.finite(dd[[trials_col]]) |
    dd[[success_col]] < 0 | dd[[trials_col]] <= 0 |
    dd[[success_col]] > dd[[trials_col]]
  if (any(bad)) stop("Invalid grouped counts for ", outcome)
  support_rows[[outcome]] <- data.table(outcome = outcome, n_islands = nrow(dd))

  for (model in names(model_rhs)) {
    key <- paste(outcome, model, sep = "__")
    fit <- fit_grouped(success_col, trials_col, model_rhs[[model]], dd)
    fits[[key]] <- fit

    s <- score_fit(fit)
    s[, `:=`(outcome = outcome, model = model)]
    score_rows[[key]] <- s

    f <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
    f[, `:=`(outcome = outcome, model = model)]
    fixed_rows[[key]] <- f
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

m3 <- fixed[model == "M3_bombus_sc" & parameter %in% c("z_bombus_deficit", "z_sc_share")]
m3[, excludes_zero := (`0.025quant` > 0 & `0.975quant` > 0) | (`0.025quant` < 0 & `0.975quant` < 0)]
m3[, direction := fifelse(mean > 0, "positive", "negative")]
setorder(m3, outcome, parameter)
fwrite(m3, file.path(outdir, "m3_category_channel_effects.csv"))

best <- scores[, .SD[which.max(log_cpo_sum)], by = outcome]
best <- best[, .(outcome, best_model = model, log_cpo_sum, waic, dic, n_cpo)]
fwrite(best, file.path(outdir, "category_best_models.csv"))

support <- rbindlist(support_rows)
fwrite(support, file.path(outdir, "category_support.csv"))

meta <- list(
  contract = "v2_inla_category_preserving_northern_midlatitude_max_data_v1",
  analysis_tier = unique(d$analysis_tier),
  regime = "northern_midlatitude",
  support_policy = "for each category, M0-M3 are fit on the same maximum complete-case support required by M3",
  category_policy = "one-vs-rest grouped binomial-logit-normal; no color/form binary collapse",
  n_complete_northern_midlatitude = nrow(north),
  n_spatial_blocks = uniqueN(north$spatial_id)
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)

capture.output(lapply(fits, summary), file = file.path(outdir, "model_summaries.txt"))
