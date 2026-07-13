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

# Standardize once on the full available island universe so regional contrasts share a scale.
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

north_main <- d[
  analysis_regime == "northern_midlatitude" &
    color_trials > 0 & form_trials > 0 & sc_trials > 0 &
    complete.cases(
      z_log_distance_to_continent_km, z_log_island_area_km2,
      z_climate_pc1, z_climate_pc2, z_climate_pc3, z_climate_pc4,
      z_bombus_deficit, z_sc_share, spatial_id
    )
]
if (nrow(north_main) < 30) stop("Too few northern-midlatitude islands on common M0-M3 support")

compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)

fit_bb <- function(lhs_success, lhs_trials, rhs, data) {
  formula <- as.formula(paste(lhs_success, "~", rhs))
  inla(
    formula,
    family = "betabinomial",
    data = data,
    Ntrials = data[[lhs_trials]],
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

# Geographic filtering of the Bombus channel.
components$bombus_geo <- fit_gaussian(
  "bombus_deficit_logit", geo_terms, north_main
)

# Reproductive-assurance alternatives.
components$sc_geo <- fit_bb(
  "sc_successes", "sc_trials", geo_terms, north_main
)
components$sc_bombus_geo <- fit_bb(
  "sc_successes", "sc_trials",
  paste("z_bombus_deficit +", geo_terms), north_main
)

# Final floral responses for M0-M3.
for (response in c("plain", "form")) {
  success <- if (response == "plain") "color_plain" else "form_open_generalized"
  trials <- if (response == "plain") "color_trials" else "form_trials"

  components[[paste0(response, "_geo")]] <- fit_bb(
    success, trials, geo_terms, north_main
  )
  components[[paste0(response, "_bombus_geo")]] <- fit_bb(
    success, trials, paste("z_bombus_deficit +", geo_terms), north_main
  )
  components[[paste0(response, "_sc_geo")]] <- fit_bb(
    success, trials, paste("z_sc_share +", geo_terms), north_main
  )
  components[[paste0(response, "_bombus_sc_geo")]] <- fit_bb(
    success, trials, paste("z_bombus_deficit + z_sc_share +", geo_terms), north_main
  )
}

scenarios <- list(
  M0_geographic_filter = c(
    bombus = "bombus_geo", self = "sc_geo",
    color = "plain_geo", form = "form_geo"
  ),
  M1_bombus_channel = c(
    bombus = "bombus_geo", self = "sc_geo",
    color = "plain_bombus_geo", form = "form_bombus_geo"
  ),
  M2_reproductive_assurance = c(
    bombus = "bombus_geo", self = "sc_geo",
    color = "plain_sc_geo", form = "form_sc_geo"
  ),
  M3_bombus_plus_reproductive_assurance = c(
    bombus = "bombus_geo", self = "sc_bombus_geo",
    color = "plain_bombus_sc_geo", form = "form_bombus_sc_geo"
  )
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

scenario_rows <- rbindlist(lapply(names(scenarios), function(scenario) {
  refs <- scenarios[[scenario]]
  s <- score_rows[match(unname(refs), component)]
  data.table(
    scenario = scenario,
    total_log_cpo = sum(s$log_cpo_sum),
    total_waic = sum(s$waic),
    total_dic = sum(s$dic)
  )
}))
scenario_rows[, delta_log_cpo_from_best := total_log_cpo - max(total_log_cpo)]
setorder(scenario_rows, -total_log_cpo)
fwrite(scenario_rows, file.path(outdir, "inla_m0_m3_scenario_comparison.csv"))

# Approximate posterior products from INLA marginals for the M3 indirect paths.
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
  path = c(
    "Bombus_deficit_to_SC_to_plain",
    "Bombus_deficit_to_SC_to_generalized"
  ),
  mean = c(mean(plain_indirect), mean(form_indirect)),
  q025 = c(quantile(plain_indirect, 0.025), quantile(form_indirect, 0.025)),
  q975 = c(quantile(plain_indirect, 0.975), quantile(form_indirect, 0.975))
)
fwrite(indirect, file.path(outdir, "inla_m3_indirect_effects.csv"))

# M4: retain detailed categories as biologically interpretable one-vs-rest models.
# This avoids collapsing red/pink, blue/purple, tubular, or zygomorphic responses.
m4 <- d[
  analysis_regime %in% c(
    "northern_midlatitude", "tropical", "southern_extratropical"
  ) &
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

meta <- list(
  contract = "v2_inla_m0_m4_main_v1",
  method = "INLA piecewise latent-Gaussian analysis with grouped beta-binomial floral/reproductive responses, Gaussian logit Bombus deficit, and 10-degree spatial-block random effects",
  analysis_tier = "sensitivity_all_all_filled_data_requested_first_run",
  sensitivity_analysis = "deferred by design",
  main_interaction = "log distance to continent x log island area",
  m0_m3_common_support = TRUE,
  m0_m3_scenario_metric = "sum of component log CPO (LPML-like; higher is better), with WAIC and DIC also reported",
  m4 = "category-preserving one-vs-rest beta-binomial models for colour and floral-form categories",
  n_northern_midlatitude_common_support = nrow(north_main),
  n_global_m4 = nrow(m4),
  n_spatial_blocks_north = uniqueN(north_main$spatial_id),
  n_spatial_blocks_global = uniqueN(m4$spatial_id)
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)

capture.output(
  lapply(c(components, m4_fits), summary),
  file = file.path(outdir, "inla_model_summaries.txt")
)
