#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(data.table)
  library(posterior)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: run_bayesian_global_isolation_counterfactual.R <input.csv> <output_dir>")
infile <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
d <- fread(infile)

# This workflow is a noncanonical scalar replication of the canonical INLA analysis.
# Wind pollination is analysed as a whole-flora composition outcome, while floral
# colour, form, and SC remain restricted to animal-mediated species.
n_input <- nrow(d)
n_zero_distance <- sum(is.finite(d$distance_to_continent_km) & d$distance_to_continent_km <= 0, na.rm = TRUE)
d <- d[
  is.finite(distance_to_continent_km) & distance_to_continent_km > 0 &
    is.finite(area_km2) & area_km2 > 0
]
if (nrow(d) < 30) stop("Too few positive-distance islands")

z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) stop("Cannot standardize zero-variance predictor")
  (x - mean(x, na.rm = TRUE)) / s
}

d[, z_log_distance := z(log(distance_to_continent_km))]
d[, z_log_area := z(log(area_km2))]
for (nm in c("climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4")) {
  d[[paste0("z_", nm)]] <- z(d[[nm]])
}
d[, wind_share := fifelse(pollen_vector_trials > 0, n_wind_species / pollen_vector_trials, NA_real_)]
d[, z_bombus_deficit := z(bombus_deficit)]
d[, z_wind_share := z(wind_share)]
d[, z_sc_share := z(sc_share)]
d[, analysis_regime := factor(analysis_regime)]

fit_bb <- function(formula, data, file) {
  if (nrow(data) < 30) stop(paste("Too few rows for", deparse(formula)))
  brm(
    formula = formula,
    data = data,
    family = beta_binomial(link = "logit"),
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1.5), class = "Intercept"),
      prior(exponential(1), class = "phi")
    ),
    chains = 2,
    cores = 2,
    iter = 2000,
    warmup = 1000,
    seed = 20260713,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = save_pars(all = TRUE),
    file = file,
    refresh = 100
  )
}

fit_gaussian <- function(formula, data, file) {
  if (nrow(data) < 30) stop(paste("Too few rows for", deparse(formula)))
  brm(
    formula = formula,
    data = data,
    family = gaussian(),
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1.5), class = "Intercept"),
      prior(exponential(1), class = "sigma")
    ),
    chains = 2,
    cores = 2,
    iter = 2000,
    warmup = 1000,
    seed = 20260713,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = save_pars(all = TRUE),
    file = file,
    refresh = 100
  )
}

# Global falsification equations use maximum available support separately.
base_complete <- complete.cases(d[, .(
  z_log_distance, z_log_area, z_climate_pc1, z_climate_pc2,
  z_climate_pc3, z_climate_pc4, analysis_regime
)])
wind_global <- d[base_complete & pollen_vector_trials > 0]
sc_global <- d[base_complete & sc_trials > 0]
color_global <- d[base_complete & color_trials > 0]
form_global <- d[base_complete & form_trials > 0]

f0_wind <- fit_bb(
  n_wind_species | trials(pollen_vector_trials) ~
    z_log_distance * analysis_regime + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  wind_global,
  file.path(outdir, "f0_global_isolation_wind")
)
f0_sc <- fit_bb(
  sc_successes | trials(sc_trials) ~
    z_log_distance * analysis_regime + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  sc_global,
  file.path(outdir, "f0_global_isolation_sc")
)
f0_color <- fit_bb(
  color_plain | trials(color_trials) ~
    z_log_distance * analysis_regime + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  color_global,
  file.path(outdir, "f0_global_isolation_plain")
)
f0_form <- fit_bb(
  form_open_generalized | trials(form_trials) ~
    z_log_distance * analysis_regime + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  form_global,
  file.path(outdir, "f0_global_isolation_generalized")
)

# Northern-midlatitude path network:
# Isolation -> Bombus deficit -> wind share -> SC share -> floral composition,
# while retaining direct paths from isolation and Bombus deficit to downstream nodes.
# Area and climate PCs are included in every equation as adjustment covariates.
north <- d[analysis_regime == "northern_midlatitude"]
north_base <- complete.cases(north[, .(
  z_log_distance, z_log_area, z_climate_pc1, z_climate_pc2,
  z_climate_pc3, z_climate_pc4, z_bombus_deficit
)])
bombus_north <- north[north_base]
wind_north <- north[north_base & pollen_vector_trials > 0]
sc_north <- north[north_base & sc_trials > 0 & is.finite(z_wind_share)]
color_north <- north[
  north_base & color_trials > 0 & is.finite(z_wind_share) & is.finite(z_sc_share)
]
form_north <- north[
  north_base & form_trials > 0 & is.finite(z_wind_share) & is.finite(z_sc_share)
]

n_bombus <- fit_gaussian(
  z_bombus_deficit ~
    z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  bombus_north,
  file.path(outdir, "northern_isolation_bombus")
)
n_wind <- fit_bb(
  n_wind_species | trials(pollen_vector_trials) ~
    z_bombus_deficit + z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  wind_north,
  file.path(outdir, "northern_bombus_wind")
)
n_sc <- fit_bb(
  sc_successes | trials(sc_trials) ~
    z_bombus_deficit + z_wind_share + z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  sc_north,
  file.path(outdir, "northern_bombus_wind_sc")
)
n_color <- fit_bb(
  color_plain | trials(color_trials) ~
    z_bombus_deficit + z_wind_share + z_sc_share + z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  color_north,
  file.path(outdir, "northern_bombus_wind_sc_plain")
)
n_form <- fit_bb(
  form_open_generalized | trials(form_trials) ~
    z_bombus_deficit + z_wind_share + z_sc_share + z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  form_north,
  file.path(outdir, "northern_bombus_wind_sc_generalized")
)

fits <- list(
  F0_wind = f0_wind,
  F0_SC = f0_sc,
  F0_plain = f0_color,
  F0_generalized = f0_form,
  N_Bombus = n_bombus,
  N_wind = n_wind,
  N_SC = n_sc,
  N_plain = n_color,
  N_generalized = n_form
)
posterior_rows <- rbindlist(lapply(names(fits), function(nm) {
  x <- as.data.table(posterior_summary(fits[[nm]], pars = "^b_"), keep.rownames = "parameter")
  x[, model := nm]
  setcolorder(x, c("model", "parameter"))
  x
}), fill = TRUE)
fwrite(posterior_rows, file.path(outdir, "posterior_fixed_effects.csv"))

support <- data.table(
  equation = c(
    "global_wind", "global_SC", "global_plain", "global_generalized",
    "northern_Bombus", "northern_wind", "northern_SC", "northern_plain", "northern_generalized"
  ),
  n_islands = c(
    nrow(wind_global), nrow(sc_global), nrow(color_global), nrow(form_global),
    nrow(bombus_north), nrow(wind_north), nrow(sc_north), nrow(color_north), nrow(form_north)
  )
)
fwrite(support, file.path(outdir, "equation_support.csv"))

coef_draw <- function(draws, name) {
  if (!name %in% names(draws)) stop(paste("Missing posterior coefficient", name))
  draws[[name]]
}
summarize_draw <- function(path, x) {
  data.table(
    path = path,
    mean = mean(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975)),
    probability_positive = mean(x > 0)
  )
}

# Align posterior draws by row index before multiplying path coefficients.
bombus_draws <- as_draws_df(n_bombus)
wind_draws <- as_draws_df(n_wind)
sc_draws <- as_draws_df(n_sc)
color_draws <- as_draws_df(n_color)
form_draws <- as_draws_df(n_form)
n_draws <- min(
  nrow(bombus_draws), nrow(wind_draws), nrow(sc_draws),
  nrow(color_draws), nrow(form_draws)
)
bombus_draws <- bombus_draws[seq_len(n_draws), ]
wind_draws <- wind_draws[seq_len(n_draws), ]
sc_draws <- sc_draws[seq_len(n_draws), ]
color_draws <- color_draws[seq_len(n_draws), ]
form_draws <- form_draws[seq_len(n_draws), ]

# Core northern path coefficients.
a_i_b <- coef_draw(bombus_draws, "b_z_log_distance")
a_i_w <- coef_draw(wind_draws, "b_z_log_distance")
a_b_w <- coef_draw(wind_draws, "b_z_bombus_deficit")
a_i_s <- coef_draw(sc_draws, "b_z_log_distance")
a_b_s <- coef_draw(sc_draws, "b_z_bombus_deficit")
a_w_s <- coef_draw(sc_draws, "b_z_wind_share")
a_i_p <- coef_draw(color_draws, "b_z_log_distance")
a_b_p <- coef_draw(color_draws, "b_z_bombus_deficit")
a_w_p <- coef_draw(color_draws, "b_z_wind_share")
a_s_p <- coef_draw(color_draws, "b_z_sc_share")

# Direct and mediated effects on plain flowers.
plain_paths <- list(
  summarize_draw("Isolation_to_Plain_direct", a_i_p),
  summarize_draw("Isolation_to_Bombus_to_Plain", a_i_b * a_b_p),
  summarize_draw("Isolation_to_Wind_to_Plain", a_i_w * a_w_p),
  summarize_draw("Isolation_to_SC_to_Plain", a_i_s * a_s_p),
  summarize_draw("Isolation_to_Bombus_to_Wind_to_Plain", a_i_b * a_b_w * a_w_p),
  summarize_draw("Isolation_to_Bombus_to_SC_to_Plain", a_i_b * a_b_s * a_s_p),
  summarize_draw("Isolation_to_Wind_to_SC_to_Plain", a_i_w * a_w_s * a_s_p),
  summarize_draw("Isolation_to_Bombus_to_Wind_to_SC_to_Plain", a_i_b * a_b_w * a_w_s * a_s_p),
  summarize_draw("Bombus_to_Plain_direct", a_b_p),
  summarize_draw("Bombus_to_Wind_to_Plain", a_b_w * a_w_p),
  summarize_draw("Bombus_to_SC_to_Plain", a_b_s * a_s_p),
  summarize_draw("Bombus_to_Wind_to_SC_to_Plain", a_b_w * a_w_s * a_s_p),
  summarize_draw("Wind_to_Plain_direct", a_w_p),
  summarize_draw("Wind_to_SC_to_Plain", a_w_s * a_s_p),
  summarize_draw("SC_to_Plain_direct", a_s_p)
)

isolation_total_plain <-
  a_i_p +
  a_i_b * a_b_p +
  a_i_w * a_w_p +
  a_i_s * a_s_p +
  a_i_b * a_b_w * a_w_p +
  a_i_b * a_b_s * a_s_p +
  a_i_w * a_w_s * a_s_p +
  a_i_b * a_b_w * a_w_s * a_s_p
bombus_total_plain <-
  a_b_p +
  a_b_w * a_w_p +
  a_b_s * a_s_p +
  a_b_w * a_w_s * a_s_p
plain_paths <- c(
  plain_paths,
  list(
    summarize_draw("Isolation_to_Plain_total", isolation_total_plain),
    summarize_draw("Bombus_to_Plain_total", bombus_total_plain)
  )
)
path_effects <- rbindlist(plain_paths)
fwrite(path_effects, file.path(outdir, "northern_path_effects.csv"))

# Keep parallel direct/SC-mediated summaries for generalized floral form.
b_sc_form <- coef_draw(form_draws, "b_z_sc_share")
b_wind_form <- coef_draw(form_draws, "b_z_wind_share")
b_bombus_form <- coef_draw(form_draws, "b_z_bombus_deficit")
generalized_paths <- rbindlist(list(
  summarize_draw("Bombus_to_generalized_direct", b_bombus_form),
  summarize_draw("Bombus_to_SC_to_generalized", a_b_s * b_sc_form),
  summarize_draw("Bombus_to_Wind_to_generalized", a_b_w * b_wind_form),
  summarize_draw("Bombus_to_Wind_to_SC_to_generalized", a_b_w * a_w_s * b_sc_form)
))
fwrite(generalized_paths, file.path(outdir, "northern_generalized_path_effects.csv"))

sampler_diagnostics <- rbindlist(lapply(names(fits), function(nm) {
  np <- nuts_params(fits[[nm]])
  draws <- as_draws_array(fits[[nm]])
  diag <- as.data.table(summarise_draws(draws, rhat, ess_bulk, ess_tail))
  diag <- diag[grepl("^b_", variable)]
  rhat_values <- diag$rhat
  bulk_ess <- diag$ess_bulk
  tail_ess <- diag$ess_tail
  data.table(
    model = nm,
    divergent_transitions = sum(np$Parameter == "divergent__" & np$Value == 1),
    max_treedepth_hits = sum(np$Parameter == "treedepth__" & np$Value >= 15),
    max_rhat = if (length(rhat_values) == 0 || all(is.na(rhat_values))) NA_real_ else max(rhat_values, na.rm = TRUE),
    min_bulk_ess = if (length(bulk_ess) == 0 || all(is.na(bulk_ess))) NA_real_ else min(bulk_ess, na.rm = TRUE),
    min_tail_ess = if (length(tail_ess) == 0 || all(is.na(tail_ess))) NA_real_ else min(tail_ess, na.rm = TRUE)
  )
}), fill = TRUE)
fwrite(sampler_diagnostics, file.path(outdir, "sampler_diagnostics.csv"))
if (any(sampler_diagnostics$divergent_transitions > 0)) warning("Divergent transitions detected; inspect sampler_diagnostics.csv")
if (any(sampler_diagnostics$max_rhat > 1.01, na.rm = TRUE)) warning("R-hat > 1.01 detected; inspect sampler_diagnostics.csv")

capture.output(lapply(fits, summary), file = file.path(outdir, "model_summaries.txt"))

meta <- list(
  contract = "v2_bayesian_scalar_falsification_replication_v5",
  role = "noncanonical scalar path replication of the canonical INLA analysis",
  method = "Global isolation-by-regime falsification plus northern Isolation -> Bombus deficit -> wind share -> SC share -> floral composition path decomposition",
  analysis_tier = "sensitivity_all",
  category_preserving_models = FALSE,
  alternative_pollinator_primary_covariates = FALSE,
  northern_path_network = "Isolation -> Bombus deficit -> wind share -> SC share -> plain/generalized floral composition, with direct downstream paths retained",
  adjustment_covariates = "log island area and climate PCs in every northern equation",
  founder_effect_covariate = "not included: no independently defined founder-effect variable is available in the current locked input",
  wind_pollination_scope = "whole-flora resolved biotic vs abiotic_wind composition; not mixed into animal-mediated floral-trait denominators",
  n_input_islands = n_input,
  n_zero_or_negative_distance_excluded = n_zero_distance,
  n_positive_distance_islands_retained_before_equation_missingness = nrow(d),
  distance_rule = "distance_to_continent_km > 0; log(distance_to_continent_km)",
  support_strategy = "maximum available islands per equation; no forced complete-case intersection across all responses",
  global_falsification = "isolation x analysis_regime for wind share, SC, plain colour, and open/generalized form",
  path_effect_output = "northern_path_effects.csv contains direct, specific indirect, and total effects on plain flowers",
  interpretation_guardrail = "Path products combine coefficients across linked equations and are reported as a scalar replication; the category-preserving INLA analysis remains canonical.",
  sampler = "2 chains x 2000 iterations, 1000 warmup, adapt_delta=0.99, max_treedepth=15"
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
