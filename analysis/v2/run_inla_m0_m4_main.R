#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: run_inla_m0_m4_main.R <input.csv> <output_dir>")
infile <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260713)
d <- fread(infile)

z <- function(x) as.numeric(scale(x))
clamp01 <- function(x, eps = 1e-5) pmin(pmax(x, eps), 1 - eps)

required <- c(
  "island_id", "analysis_regime", "log_distance_to_continent_km",
  "log_island_area_km2", "climate_pc1", "climate_pc2", "climate_pc3",
  "climate_pc4", "bombus_deficit", "spatial_block", "color_trials",
  "color_plain", "form_trials", "form_open_generalized", "sc_trials",
  "sc_successes", "sc_share", "showy_alt_guild_share",
  "other_bee_guild_share_1", "generalist_insect_guild_share_1"
)
missing <- setdiff(required, names(d))
if (length(missing)) stop("missing required columns: ", paste(missing, collapse = ", "))

# Standardize on the full model-input universe (distance > 0 is enforced upstream).
for (nm in c(
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"
)) {
  d[[paste0("z_", nm)]] <- z(d[[nm]])
}
d[, z_bombus_deficit := z(bombus_deficit)]
d[, z_sc_share := z(sc_share)]
d[, z_showy_alt_guild_share := z(showy_alt_guild_share)]
d[, z_other_bee_guild_share := z(other_bee_guild_share_1)]
d[, z_generalist_insect_guild_share := z(generalist_insect_guild_share_1)]
d[, spatial_id := as.integer(factor(spatial_block))]
d[, bombus_deficit_logit := qlogis(clamp01(bombus_deficit))]

geo_terms <- paste(
  "z_log_distance_to_continent_km * z_log_island_area_km2 +",
  "z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4 +",
  "f(spatial_id, model='iid')"
)

base_north <- d[
  analysis_regime == "northern_midlatitude" &
    complete.cases(
      z_log_distance_to_continent_km, z_log_island_area_km2,
      z_climate_pc1, z_climate_pc2, z_climate_pc3, z_climate_pc4,
      spatial_id
    )
]
if (nrow(base_north) < 30) stop("Too few northern-midlatitude islands after geographic filtering")

bombus_data <- base_north[complete.cases(z_bombus_deficit, bombus_deficit_logit)]
sc_geo_data <- base_north[sc_trials > 0]
sc_bombus_data <- base_north[sc_trials > 0 & complete.cases(z_bombus_deficit)]
color_geo_data <- base_north[color_trials > 0]
color_bombus_data <- base_north[color_trials > 0 & complete.cases(z_bombus_deficit)]
color_sc_data <- base_north[color_trials > 0 & complete.cases(z_sc_share)]
color_bombus_sc_data <- base_north[color_trials > 0 & complete.cases(z_bombus_deficit, z_sc_share)]
form_geo_data <- base_north[form_trials > 0]
form_bombus_data <- base_north[form_trials > 0 & complete.cases(z_bombus_deficit)]
form_sc_data <- base_north[form_trials > 0 & complete.cases(z_sc_share)]
form_bombus_sc_data <- base_north[form_trials > 0 & complete.cases(z_bombus_deficit, z_sc_share)]

compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)

fit_bb <- function(lhs_success, lhs_trials, rhs, data) {
  dd <- copy(data)
  success <- dd[[lhs_success]]
  trials <- dd[[lhs_trials]]
  bad <- !is.finite(success) | !is.finite(trials) |
    success < 0 | trials <= 0 | success > trials |
    success != floor(success) | trials != floor(trials)
  if (any(bad)) {
    stop(
      "Invalid grouped-binomial counts for ", lhs_success, "/", lhs_trials,
      ": ", sum(bad), " invalid rows"
    )
  }
  dd[, obs_id := seq_len(.N)]
  formula <- as.formula(paste(lhs_success, "~", rhs, "+ f(obs_id, model='iid')"))
  inla(
    formula,
    family = "binomial",
    data = dd,
    Ntrials = trials,
    control.compute = compute,
    control.predictor = list(compute = TRUE),
    verbose = FALSE
  )
}

fit_gaussian <- function(lhs, rhs, data) {
  inla(
    as.formula(paste(lhs, "~", rhs)),
    family = "gaussian",
    data = data,
    control.compute = compute,
    control.predictor = list(compute = TRUE),
    verbose = FALSE
  )
}

components <- list()
components$bombus_geo <- fit_gaussian("bombus_deficit_logit", geo_terms, bombus_data)
components$sc_geo <- fit_bb("sc_successes", "sc_trials", geo_terms, sc_geo_data)
components$sc_bombus_geo <- fit_bb(
  "sc_successes", "sc_trials", paste("z_bombus_deficit +", geo_terms), sc_bombus_data
)
components$plain_geo <- fit_bb("color_plain", "color_trials", geo_terms, color_geo_data)
components$plain_bombus_geo <- fit_bb(
  "color_plain", "color_trials", paste("z_bombus_deficit +", geo_terms), color_bombus_data
)
components$plain_sc_geo <- fit_bb(
  "color_plain", "color_trials", paste("z_sc_share +", geo_terms), color_sc_data
)
components$plain_bombus_sc_geo <- fit_bb(
  "color_plain", "color_trials", paste("z_bombus_deficit + z_sc_share +", geo_terms), color_bombus_sc_data
)
components$form_geo <- fit_bb("form_open_generalized", "form_trials", geo_terms, form_geo_data)
components$form_bombus_geo <- fit_bb(
  "form_open_generalized", "form_trials", paste("z_bombus_deficit +", geo_terms), form_bombus_data
)
components$form_sc_geo <- fit_bb(
  "form_open_generalized", "form_trials", paste("z_sc_share +", geo_terms), form_sc_data
)
components$form_bombus_sc_geo <- fit_bb(
  "form_open_generalized", "form_trials", paste("z_bombus_deficit + z_sc_share +", geo_terms), form_bombus_sc_data
)

fixed_rows <- rbindlist(lapply(names(components), function(nm) {
  x <- as.data.table(components[[nm]]$summary.fixed, keep.rownames = "parameter")
  x[, component := nm]
  setcolorder(x, c("component", "parameter"))
  x
}), fill = TRUE)
fwrite(fixed_rows, file.path(outdir, "inla_component_fixed_effects.csv"))

component_score <- function(fit) {
  cpo <- fit$cpo$cpo
  cpo <- cpo[is.finite(cpo) & cpo > 0]
  data.table(
    log_cpo_sum = if (length(cpo)) sum(log(cpo)) else NA_real_,
    n_cpo = length(cpo),
    waic = fit$waic$waic,
    dic = fit$dic$dic
  )
}

score_rows <- rbindlist(lapply(names(components), function(nm) {
  x <- component_score(components[[nm]])
  x[, component := nm]
  setcolorder(x, c("component", "log_cpo_sum", "n_cpo", "waic", "dic"))
  x
}))
fwrite(score_rows, file.path(outdir, "inla_component_scores.csv"))

# Compare alternatives only within the same biological endpoint. Different endpoints
# deliberately use all available observations and therefore need not share one sample size.
endpoint_comparisons <- rbindlist(list(
  data.table(endpoint = "self_compatibility", model = c("M0_geo", "M1_bombus"), component = c("sc_geo", "sc_bombus_geo")),
  data.table(endpoint = "flower_color", model = c("M0_geo", "M1_bombus", "M2_sc", "M3_bombus_sc"), component = c("plain_geo", "plain_bombus_geo", "plain_sc_geo", "plain_bombus_sc_geo")),
  data.table(endpoint = "floral_form", model = c("M0_geo", "M1_bombus", "M2_sc", "M3_bombus_sc"), component = c("form_geo", "form_bombus_geo", "form_sc_geo", "form_bombus_sc_geo"))
))
endpoint_comparisons <- merge(endpoint_comparisons, score_rows, by = "component", all.x = TRUE)
endpoint_comparisons[, delta_log_cpo_from_best := log_cpo_sum - max(log_cpo_sum, na.rm = TRUE), by = endpoint]
setorder(endpoint_comparisons, endpoint, -log_cpo_sum)
fwrite(endpoint_comparisons, file.path(outdir, "inla_endpoint_model_comparisons.csv"))

# M3 mediation remains estimated on the observations required by each path.
draw_coef <- function(fit, parameter, n = 50000L) {
  m <- fit$marginals.fixed[[parameter]]
  if (is.null(m)) stop("fixed-effect marginal not found: ", parameter)
  inla.rmarginal(n, m)
}
b_bombus_sc <- draw_coef(components$sc_bombus_geo, "z_bombus_deficit")
b_sc_plain <- draw_coef(components$plain_bombus_sc_geo, "z_sc_share")
b_sc_form <- draw_coef(components$form_bombus_sc_geo, "z_sc_share")
plain_indirect <- b_bombus_sc * b_sc_plain
form_indirect <- b_bombus_sc * b_sc_form
indirect <- data.table(
  path = c("Bombus_deficit_to_SC_to_plain", "Bombus_deficit_to_SC_to_generalized"),
  mean = c(mean(plain_indirect), mean(form_indirect)),
  q025 = c(quantile(plain_indirect, 0.025), quantile(form_indirect, 0.025)),
  q975 = c(quantile(plain_indirect, 0.975), quantile(form_indirect, 0.975))
)
fwrite(indirect, file.path(outdir, "inla_m3_indirect_effects.csv"))

# M4: category-preserving one-vs-rest models, each using every island with trials > 0.
m4 <- d[
  analysis_regime %in% c("northern_midlatitude", "tropical", "southern_extratropical") &
    complete.cases(
      z_log_distance_to_continent_km, z_log_island_area_km2,
      z_showy_alt_guild_share, z_other_bee_guild_share,
      z_generalist_insect_guild_share, spatial_id
    )
]
m4[, analysis_regime := factor(analysis_regime)]
m4_rhs <- paste(
  "z_log_distance_to_continent_km * z_log_island_area_km2 * analysis_regime +",
  "z_showy_alt_guild_share + z_other_bee_guild_share +",
  "z_generalist_insect_guild_share + f(spatial_id, model='iid')"
)
m4_specs <- list(
  color_yellow_orange = c("color_yellow_orange", "color_trials"),
  color_red_pink = c("color_red_pink", "color_trials"),
  color_blue_purple = c("color_blue_purple", "color_trials"),
  form_tubular_trumpet = c("form_tubular_trumpet", "form_trials"),
  form_zygomorphic_specialized = c("form_zygomorphic_specialized", "form_trials"),
  form_composite_brush = c("form_composite_brush", "form_trials")
)
m4_fits <- list()
for (nm in names(m4_specs)) {
  success <- m4_specs[[nm]][1]
  trials <- m4_specs[[nm]][2]
  dd <- m4[get(trials) > 0]
  m4_fits[[nm]] <- fit_bb(success, trials, m4_rhs, dd)
}
m4_fixed <- rbindlist(lapply(names(m4_fits), function(nm) {
  x <- as.data.table(m4_fits[[nm]]$summary.fixed, keep.rownames = "parameter")
  x[, outcome := nm]
  setcolorder(x, c("outcome", "parameter"))
  x
}), fill = TRUE)
fwrite(m4_fixed, file.path(outdir, "inla_m4_category_fixed_effects.csv"))
m4_scores <- rbindlist(lapply(names(m4_fits), function(nm) {
  x <- component_score(m4_fits[[nm]])
  x[, outcome := nm]
  setcolorder(x, c("outcome", "log_cpo_sum", "n_cpo", "waic", "dic"))
  x
}))
fwrite(m4_scores, file.path(outdir, "inla_m4_category_scores.csv"))

support <- data.table(
  component = c(
    "bombus_geo", "sc_geo", "sc_bombus_geo",
    "plain_geo", "plain_bombus_geo", "plain_sc_geo", "plain_bombus_sc_geo",
    "form_geo", "form_bombus_geo", "form_sc_geo", "form_bombus_sc_geo",
    "m4_base"
  ),
  n_islands = c(
    nrow(bombus_data), nrow(sc_geo_data), nrow(sc_bombus_data),
    nrow(color_geo_data), nrow(color_bombus_data), nrow(color_sc_data), nrow(color_bombus_sc_data),
    nrow(form_geo_data), nrow(form_bombus_data), nrow(form_sc_data), nrow(form_bombus_sc_data),
    nrow(m4)
  )
)
fwrite(support, file.path(outdir, "inla_component_support.csv"))

meta <- list(
  contract = "v2_inla_m0_m4_all_filled_endpoint_specific_support_v3",
  method = "INLA piecewise latent-Gaussian analysis with grouped binomial-logit-normal responses, endpoint-specific maximum available data, observation-level iid overdispersion, Gaussian logit Bombus deficit, and 10-degree spatial-block random effects",
  analysis_tier = "sensitivity_all",
  distance_scope = "distance_to_continent_km > 0 enforced upstream; zero-distance islands retained only in audit/full-universe input",
  common_support = FALSE,
  support_policy = "each endpoint/model uses all available observations required for that model; model-score comparisons are reported within endpoint rather than summed across endpoints",
  main_interaction = "log distance to continent x log island area",
  m4 = "category-preserving one-vs-rest grouped binomial-logit-normal models using all available trials per category",
  n_base_northern_midlatitude = nrow(base_north),
  n_bombus = nrow(bombus_data),
  n_sc_geo = nrow(sc_geo_data),
  n_color_geo = nrow(color_geo_data),
  n_form_geo = nrow(form_geo_data),
  n_global_m4 = nrow(m4),
  n_spatial_blocks_north = uniqueN(base_north$spatial_id),
  n_spatial_blocks_global = uniqueN(m4$spatial_id)
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)

capture.output(
  lapply(c(components, m4_fits), summary),
  file = file.path(outdir, "inla_model_summaries.txt")
)
