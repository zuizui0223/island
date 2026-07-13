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

# Noncanonical scalar replication/falsification of the canonical category-preserving INLA analysis.
# Wind pollination is a whole-flora outcome. SC and floral traits remain restricted to animal-mediated species.
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

# Global isolation-by-regime falsification, using maximum support per response.
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

# Fit the same path network independently in the northern focal domain and in
# tropical/southern falsification domains. Area and climate PCs adjust every equation.
fit_regime_network <- function(regime_name, prefix) {
  x <- d[analysis_regime == regime_name]
  ok <- complete.cases(x[, .(
    z_log_distance, z_log_area, z_climate_pc1, z_climate_pc2,
    z_climate_pc3, z_climate_pc4, z_bombus_deficit
  )])
  bombus_dat <- x[ok]
  wind_dat <- x[ok & pollen_vector_trials > 0]
  sc_dat <- x[ok & sc_trials > 0 & is.finite(z_wind_share)]
  color_dat <- x[ok & color_trials > 0 & is.finite(z_wind_share) & is.finite(z_sc_share)]
  form_dat <- x[ok & form_trials > 0 & is.finite(z_wind_share) & is.finite(z_sc_share)]

  fits <- list(
    Bombus = fit_gaussian(
      z_bombus_deficit ~ z_log_distance + z_log_area +
        z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
      bombus_dat, file.path(outdir, paste0(prefix, "_isolation_bombus"))
    ),
    wind = fit_bb(
      n_wind_species | trials(pollen_vector_trials) ~
        z_bombus_deficit + z_log_distance + z_log_area +
        z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
      wind_dat, file.path(outdir, paste0(prefix, "_bombus_wind"))
    ),
    SC = fit_bb(
      sc_successes | trials(sc_trials) ~
        z_bombus_deficit + z_wind_share + z_log_distance + z_log_area +
        z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
      sc_dat, file.path(outdir, paste0(prefix, "_bombus_wind_sc"))
    ),
    plain = fit_bb(
      color_plain | trials(color_trials) ~
        z_bombus_deficit + z_wind_share + z_sc_share + z_log_distance + z_log_area +
        z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
      color_dat, file.path(outdir, paste0(prefix, "_bombus_wind_sc_plain"))
    ),
    generalized = fit_bb(
      form_open_generalized | trials(form_trials) ~
        z_bombus_deficit + z_wind_share + z_sc_share + z_log_distance + z_log_area +
        z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
      form_dat, file.path(outdir, paste0(prefix, "_bombus_wind_sc_generalized"))
    )
  )

  support <- data.table(
    regime = regime_name,
    equation = c("Bombus", "wind", "SC", "plain", "generalized"),
    n_islands = c(nrow(bombus_dat), nrow(wind_dat), nrow(sc_dat), nrow(color_dat), nrow(form_dat))
  )
  list(fits = fits, support = support)
}

networks <- list(
  northern_midlatitude = fit_regime_network("northern_midlatitude", "northern"),
  tropical = fit_regime_network("tropical", "tropical"),
  southern_extratropical = fit_regime_network("southern_extratropical", "southern")
)

fits <- list(
  F0_wind = f0_wind,
  F0_SC = f0_sc,
  F0_plain = f0_color,
  F0_generalized = f0_form
)
for (regime_name in names(networks)) {
  tag <- c(
    northern_midlatitude = "N",
    tropical = "T",
    southern_extratropical = "S"
  )[[regime_name]]
  for (eq in names(networks[[regime_name]]$fits)) {
    fits[[paste0(tag, "_", eq)]] <- networks[[regime_name]]$fits[[eq]]
  }
}

posterior_rows <- rbindlist(lapply(names(fits), function(nm) {
  x <- as.data.table(posterior_summary(fits[[nm]], pars = "^b_"), keep.rownames = "parameter")
  x[, model := nm]
  setcolorder(x, c("model", "parameter"))
  x
}), fill = TRUE)
fwrite(posterior_rows, file.path(outdir, "posterior_fixed_effects.csv"))

support <- rbindlist(c(
  list(data.table(
    regime = "global",
    equation = c("wind", "SC", "plain", "generalized"),
    n_islands = c(nrow(wind_global), nrow(sc_global), nrow(color_global), nrow(form_global))
  )),
  lapply(networks, `[[`, "support")
))
fwrite(support, file.path(outdir, "equation_support.csv"))

coef_draw <- function(draws, name) {
  if (!name %in% names(draws)) stop(paste("Missing posterior coefficient", name))
  draws[[name]]
}

summarize_draw <- function(regime, path, x) {
  data.table(
    regime = regime,
    path = path,
    mean = mean(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975)),
    probability_positive = mean(x > 0)
  )
}

# Compute identical direct, specific indirect, and total effects in each regime.
summarize_network_paths <- function(regime_name, network) {
  b <- as_draws_df(network$fits$Bombus)
  w <- as_draws_df(network$fits$wind)
  s <- as_draws_df(network$fits$SC)
  p <- as_draws_df(network$fits$plain)
  g <- as_draws_df(network$fits$generalized)
  nd <- min(nrow(b), nrow(w), nrow(s), nrow(p), nrow(g))
  b <- b[seq_len(nd), ]
  w <- w[seq_len(nd), ]
  s <- s[seq_len(nd), ]
  p <- p[seq_len(nd), ]
  g <- g[seq_len(nd), ]

  a_i_b <- coef_draw(b, "b_z_log_distance")
  a_i_w <- coef_draw(w, "b_z_log_distance")
  a_b_w <- coef_draw(w, "b_z_bombus_deficit")
  a_i_s <- coef_draw(s, "b_z_log_distance")
  a_b_s <- coef_draw(s, "b_z_bombus_deficit")
  a_w_s <- coef_draw(s, "b_z_wind_share")
  a_i_p <- coef_draw(p, "b_z_log_distance")
  a_b_p <- coef_draw(p, "b_z_bombus_deficit")
  a_w_p <- coef_draw(p, "b_z_wind_share")
  a_s_p <- coef_draw(p, "b_z_sc_share")
  a_i_g <- coef_draw(g, "b_z_log_distance")
  a_b_g <- coef_draw(g, "b_z_bombus_deficit")
  a_w_g <- coef_draw(g, "b_z_wind_share")
  a_s_g <- coef_draw(g, "b_z_sc_share")

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

  isolation_total_generalized <-
    a_i_g +
    a_i_b * a_b_g +
    a_i_w * a_w_g +
    a_i_s * a_s_g +
    a_i_b * a_b_w * a_w_g +
    a_i_b * a_b_s * a_s_g +
    a_i_w * a_w_s * a_s_g +
    a_i_b * a_b_w * a_w_s * a_s_g

  rows <- list(
    summarize_draw(regime_name, "Isolation_to_Bombus_direct", a_i_b),
    summarize_draw(regime_name, "Isolation_to_Wind_direct", a_i_w),
    summarize_draw(regime_name, "Bombus_to_Wind_direct", a_b_w),
    summarize_draw(regime_name, "Isolation_to_SC_direct", a_i_s),
    summarize_draw(regime_name, "Bombus_to_SC_direct", a_b_s),
    summarize_draw(regime_name, "Wind_to_SC_direct", a_w_s),
    summarize_draw(regime_name, "Isolation_to_Plain_direct", a_i_p),
    summarize_draw(regime_name, "Bombus_to_Plain_direct", a_b_p),
    summarize_draw(regime_name, "Wind_to_Plain_direct", a_w_p),
    summarize_draw(regime_name, "SC_to_Plain_direct", a_s_p),
    summarize_draw(regime_name, "Isolation_to_Bombus_to_Plain", a_i_b * a_b_p),
    summarize_draw(regime_name, "Isolation_to_Wind_to_Plain", a_i_w * a_w_p),
    summarize_draw(regime_name, "Isolation_to_SC_to_Plain", a_i_s * a_s_p),
    summarize_draw(regime_name, "Isolation_to_Bombus_to_Wind_to_Plain", a_i_b * a_b_w * a_w_p),
    summarize_draw(regime_name, "Isolation_to_Bombus_to_SC_to_Plain", a_i_b * a_b_s * a_s_p),
    summarize_draw(regime_name, "Isolation_to_Wind_to_SC_to_Plain", a_i_w * a_w_s * a_s_p),
    summarize_draw(regime_name, "Isolation_to_Bombus_to_Wind_to_SC_to_Plain", a_i_b * a_b_w * a_w_s * a_s_p),
    summarize_draw(regime_name, "Bombus_to_Wind_to_Plain", a_b_w * a_w_p),
    summarize_draw(regime_name, "Bombus_to_SC_to_Plain", a_b_s * a_s_p),
    summarize_draw(regime_name, "Bombus_to_Wind_to_SC_to_Plain", a_b_w * a_w_s * a_s_p),
    summarize_draw(regime_name, "Wind_to_SC_to_Plain", a_w_s * a_s_p),
    summarize_draw(regime_name, "Isolation_to_Plain_total", isolation_total_plain),
    summarize_draw(regime_name, "Bombus_to_Plain_total", bombus_total_plain),
    summarize_draw(regime_name, "Isolation_to_Generalized_direct", a_i_g),
    summarize_draw(regime_name, "Bombus_to_Generalized_direct", a_b_g),
    summarize_draw(regime_name, "Wind_to_Generalized_direct", a_w_g),
    summarize_draw(regime_name, "SC_to_Generalized_direct", a_s_g),
    summarize_draw(regime_name, "Isolation_to_Generalized_total", isolation_total_generalized)
  )

  list(
    summary = rbindlist(rows),
    draws = list(
      Isolation_to_Bombus_direct = a_i_b,
      Isolation_to_Wind_direct = a_i_w,
      Bombus_to_Wind_direct = a_b_w,
      Isolation_to_SC_direct = a_i_s,
      Bombus_to_SC_direct = a_b_s,
      Wind_to_SC_direct = a_w_s,
      Isolation_to_Plain_direct = a_i_p,
      Bombus_to_Plain_direct = a_b_p,
      Wind_to_Plain_direct = a_w_p,
      SC_to_Plain_direct = a_s_p,
      Isolation_to_Plain_total = isolation_total_plain,
      Bombus_to_Plain_total = bombus_total_plain,
      Isolation_to_Generalized_total = isolation_total_generalized
    )
  )
}

path_results <- lapply(names(networks), function(regime_name) {
  summarize_network_paths(regime_name, networks[[regime_name]])
})
names(path_results) <- names(networks)
regime_path_effects <- rbindlist(lapply(path_results, `[[`, "summary"))
fwrite(regime_path_effects, file.path(outdir, "regime_path_effects.csv"))

# Explicit falsification contrasts: tropical/southern effects minus the northern focal effect.
contrast_draws <- function(target_regime, path_name) {
  north_draw <- path_results$northern_midlatitude$draws[[path_name]]
  target_draw <- path_results[[target_regime]]$draws[[path_name]]
  n <- min(length(north_draw), length(target_draw))
  target_draw[seq_len(n)] - north_draw[seq_len(n)]
}

contrast_paths <- c(
  "Isolation_to_Bombus_direct",
  "Isolation_to_Wind_direct",
  "Bombus_to_Wind_direct",
  "Isolation_to_SC_direct",
  "Bombus_to_SC_direct",
  "Wind_to_SC_direct",
  "Isolation_to_Plain_direct",
  "Bombus_to_Plain_direct",
  "Wind_to_Plain_direct",
  "SC_to_Plain_direct",
  "Isolation_to_Plain_total",
  "Bombus_to_Plain_total",
  "Isolation_to_Generalized_total"
)

falsification_contrasts <- rbindlist(lapply(
  c("tropical", "southern_extratropical"),
  function(regime_name) rbindlist(lapply(contrast_paths, function(path_name) {
    x <- contrast_draws(regime_name, path_name)
    data.table(
      contrast = paste0(regime_name, "_minus_northern_midlatitude"),
      path = path_name,
      mean_difference = mean(x),
      q025 = unname(quantile(x, 0.025)),
      q975 = unname(quantile(x, 0.975)),
      probability_difference_positive = mean(x > 0)
    )
  }))
))
fwrite(falsification_contrasts, file.path(outdir, "falsification_contrasts_vs_north.csv"))

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
  contract = "v2_bayesian_scalar_falsification_replication_v6",
  role = "noncanonical scalar path replication and cross-regime falsification of the canonical INLA analysis",
  method = "Global isolation-by-regime tests plus matched path networks in northern-midlatitude, tropical, and southern-extratropical regimes",
  analysis_tier = "sensitivity_all",
  category_preserving_models = FALSE,
  alternative_pollinator_primary_covariates = FALSE,
  focal_regime = "northern_midlatitude",
  falsification_regimes = c("tropical", "southern_extratropical"),
  path_network = "Isolation -> Bombus deficit -> wind share -> SC share -> plain/generalized floral composition, with direct downstream paths retained",
  adjustment_covariates = "log island area and climate PCs in every regime-specific equation",
  founder_effect_covariate = "not included: no independently defined founder-effect variable is available in the current locked input",
  wind_pollination_scope = "whole-flora resolved biotic vs abiotic_wind composition; not mixed into animal-mediated floral-trait denominators",
  n_input_islands = n_input,
  n_zero_or_negative_distance_excluded = n_zero_distance,
  n_positive_distance_islands_retained_before_equation_missingness = nrow(d),
  distance_rule = "distance_to_continent_km > 0; log(distance_to_continent_km)",
  support_strategy = "maximum available islands per equation within each regime; no forced complete-case intersection across all responses",
  global_falsification = "isolation x analysis_regime for wind share, SC, plain colour, and open/generalized form",
  regime_path_output = "regime_path_effects.csv contains matched direct, specific indirect, and total effects for northern, tropical, and southern regimes",
  contrast_output = "falsification_contrasts_vs_north.csv contains tropical-minus-northern and southern-minus-northern posterior contrasts for key paths",
  interpretation_guardrail = "Products combine coefficients across linked equations and are scalar path summaries; category-preserving INLA remains canonical.",
  sampler = "2 chains x 2000 iterations, 1000 warmup, adapt_delta=0.99, max_treedepth=15"
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
