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
# It does not fit category-preserving flower models and does not use alternative-
# pollinator guilds as primary covariates.
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
d[, z_bombus_deficit := z(bombus_deficit)]
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

# Each equation uses its own maximum available support. Missing form evidence does
# not remove an otherwise usable island from the colour or SC equations.
base_complete <- complete.cases(d[, .(
  z_log_distance, z_log_area, z_climate_pc1, z_climate_pc2,
  z_climate_pc3, z_climate_pc4, analysis_regime
)])
sc_global <- d[base_complete & sc_trials > 0]
color_global <- d[base_complete & color_trials > 0]
form_global <- d[base_complete & form_trials > 0]

# Global falsification: test whether isolation slopes differ among regimes.
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

# Northern-midlatitude scalar replication of the canonical direct/indirect path.
north <- d[analysis_regime == "northern_midlatitude"]
north_base <- complete.cases(north[, .(
  z_log_distance, z_log_area, z_climate_pc1, z_climate_pc2,
  z_climate_pc3, z_climate_pc4, z_bombus_deficit
)])
sc_north <- north[north_base & sc_trials > 0]
color_north <- north[north_base & color_trials > 0 & is.finite(z_sc_share)]
form_north <- north[north_base & form_trials > 0 & is.finite(z_sc_share)]

n_sc <- fit_bb(
  sc_successes | trials(sc_trials) ~
    z_bombus_deficit + z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  sc_north,
  file.path(outdir, "northern_bombus_sc")
)
n_color <- fit_bb(
  color_plain | trials(color_trials) ~
    z_bombus_deficit + z_sc_share + z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  color_north,
  file.path(outdir, "northern_bombus_sc_plain")
)
n_form <- fit_bb(
  form_open_generalized | trials(form_trials) ~
    z_bombus_deficit + z_sc_share + z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  form_north,
  file.path(outdir, "northern_bombus_sc_generalized")
)

fits <- list(
  F0_SC = f0_sc,
  F0_plain = f0_color,
  F0_generalized = f0_form,
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
    "global_SC", "global_plain", "global_generalized",
    "northern_SC", "northern_plain", "northern_generalized"
  ),
  n_islands = c(
    nrow(sc_global), nrow(color_global), nrow(form_global),
    nrow(sc_north), nrow(color_north), nrow(form_north)
  )
)
fwrite(support, file.path(outdir, "equation_support.csv"))

# Posterior products for the northern SC-mediated paths only.
sc_draws <- as_draws_df(n_sc)
color_draws <- as_draws_df(n_color)
form_draws <- as_draws_df(n_form)
n_draws <- min(nrow(sc_draws), nrow(color_draws), nrow(form_draws))
sc_draws <- sc_draws[seq_len(n_draws), ]
color_draws <- color_draws[seq_len(n_draws), ]
form_draws <- form_draws[seq_len(n_draws), ]
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

b_bombus_sc <- coef_draw(sc_draws, "b_z_bombus_deficit")
b_sc_plain <- coef_draw(color_draws, "b_z_sc_share")
b_sc_form <- coef_draw(form_draws, "b_z_sc_share")
b_bombus_plain <- coef_draw(color_draws, "b_z_bombus_deficit")
b_bombus_form <- coef_draw(form_draws, "b_z_bombus_deficit")

paths <- rbindlist(list(
  summarize_draw("Bombus_deficit_to_SC_to_plain", b_bombus_sc * b_sc_plain),
  summarize_draw("Bombus_deficit_to_SC_to_generalized", b_bombus_sc * b_sc_form),
  summarize_draw("Bombus_deficit_to_plain_direct", b_bombus_plain),
  summarize_draw("Bombus_deficit_to_generalized_direct", b_bombus_form),
  summarize_draw("Bombus_deficit_to_plain_total", b_bombus_plain + b_bombus_sc * b_sc_plain),
  summarize_draw("Bombus_deficit_to_generalized_total", b_bombus_form + b_bombus_sc * b_sc_form)
))
fwrite(paths, file.path(outdir, "northern_direct_indirect_total_effects.csv"))

# Use posterior's stable diagnostic API rather than brms summary column names,
# which differ across brms/posterior versions (for example Rhat vs rhat).
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
  contract = "v2_bayesian_scalar_falsification_replication_v2",
  role = "noncanonical scalar replication of the canonical INLA analysis",
  method = "Global isolation-by-regime falsification plus northern scalar Bombus/SC path replication",
  analysis_tier = "sensitivity_all",
  category_preserving_models = false,
  alternative_pollinator_primary_covariates = false,
  n_input_islands = n_input,
  n_zero_or_negative_distance_excluded = n_zero_distance,
  n_positive_distance_islands_retained_before_equation_missingness = nrow(d),
  distance_rule = "distance_to_continent_km > 0; log(distance_to_continent_km)",
  support_strategy = "maximum available islands per equation; no forced colour+form+SC complete-case intersection",
  global_falsification = "isolation x analysis_regime for SC, plain colour, and open/generalized form",
  northern_replication = "Bombus deficit -> SC and direct/SC-mediated paths to scalar plain/generalized outcomes",
  interpretation_guardrail = "This workflow does not replace the category-preserving INLA analysis and does not estimate alternative-pollinator causal paths.",
  sampler = "2 chains x 2000 iterations, 1000 warmup, adapt_delta=0.99, max_treedepth=15"
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
