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

# Whole-flora pollen-vector outcome.
d[, wind_share := fifelse(pollen_vector_trials > 0, n_wind_species / pollen_vector_trials, NA_real_)]

# Flower-colour contrast: plain = white/green/brown/inconspicuous;
# showy = all resolved non-plain colours.
d[, showy_successes := color_trials - color_plain]
d[, showy_share := fifelse(color_trials > 0, showy_successes / color_trials, NA_real_)]

# Flower-form contrast over the four exhaustive resolved form categories:
# complex = tubular/trumpet + zygomorphic/specialized;
# simple = open/generalized + composite/brush.
d[, complex_successes := form_tubular_trumpet + form_zygomorphic_specialized]
d[, simple_successes := form_open_generalized + form_composite_brush]
stopifnot(all(d$complex_successes + d$simple_successes == d$form_trials, na.rm = TRUE))
d[, complex_share := fifelse(form_trials > 0, complex_successes / form_trials, NA_real_)]

d[, z_bombus_deficit := z(bombus_deficit)]
d[, z_wind_share := z(wind_share)]

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

covars <- c(
  "z_log_distance", "z_log_area", "z_climate_pc1", "z_climate_pc2",
  "z_climate_pc3", "z_climate_pc4"
)

# Northern focal mechanism. SC and floral traits are parallel downstream outcomes;
# SC is deliberately NOT used to explain flower colour or form, avoiding a strong
# causal claim between island-level summaries of traits from the same species pool.
north <- d[analysis_regime == "northern_midlatitude"]
n_base <- complete.cases(north[, ..covars]) & is.finite(north$z_bombus_deficit)
n_bombus_dat <- north[n_base]
n_wind_dat <- north[n_base & pollen_vector_trials > 0]
n_sc_dat <- north[n_base & sc_trials > 0 & is.finite(z_wind_share)]
n_showy_dat <- north[n_base & color_trials > 0 & is.finite(z_wind_share)]
n_complex_dat <- north[n_base & form_trials > 0 & is.finite(z_wind_share)]

n_bombus <- fit_gaussian(
  z_bombus_deficit ~ z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  n_bombus_dat,
  file.path(outdir, "north_isolation_bombus")
)
n_wind <- fit_bb(
  n_wind_species | trials(pollen_vector_trials) ~
    z_bombus_deficit + z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  n_wind_dat,
  file.path(outdir, "north_bombus_wind")
)
n_sc <- fit_bb(
  sc_successes | trials(sc_trials) ~
    z_bombus_deficit + z_wind_share + z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  n_sc_dat,
  file.path(outdir, "north_bombus_wind_sc")
)
n_showy <- fit_bb(
  showy_successes | trials(color_trials) ~
    z_bombus_deficit + z_wind_share + z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  n_showy_dat,
  file.path(outdir, "north_bombus_wind_showy")
)
n_complex <- fit_bb(
  complex_successes | trials(form_trials) ~
    z_bombus_deficit + z_wind_share + z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  n_complex_dat,
  file.path(outdir, "north_bombus_wind_complex")
)

# Tropical + southern counter-domain. No plant-derived pollinator-guild predictor is
# used. The falsification question is simply whether isolation shifts flower colour
# and form in the same direction as in the northern Bombus-loss domain.
south <- d[analysis_regime %in% c("tropical", "southern_extratropical")]
s_base <- complete.cases(south[, ..covars])
s_showy_dat <- south[s_base & color_trials > 0]
s_complex_dat <- south[s_base & form_trials > 0]

s_showy <- fit_bb(
  showy_successes | trials(color_trials) ~
    z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  s_showy_dat,
  file.path(outdir, "south_isolation_showy")
)
s_complex <- fit_bb(
  complex_successes | trials(form_trials) ~
    z_log_distance + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  s_complex_dat,
  file.path(outdir, "south_isolation_complex")
)

fits <- list(
  N_Bombus = n_bombus,
  N_wind = n_wind,
  N_SC = n_sc,
  N_showy = n_showy,
  N_complex = n_complex,
  S_showy = s_showy,
  S_complex = s_complex
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
    "north_Bombus", "north_wind", "north_SC", "north_showy", "north_complex",
    "tropical_southern_showy", "tropical_southern_complex"
  ),
  n_islands = c(
    nrow(n_bombus_dat), nrow(n_wind_dat), nrow(n_sc_dat), nrow(n_showy_dat), nrow(n_complex_dat),
    nrow(s_showy_dat), nrow(s_complex_dat)
  )
)
fwrite(support, file.path(outdir, "equation_support.csv"))

coef_draw <- function(draws, name) {
  if (!name %in% names(draws)) stop(paste("Missing posterior coefficient", name))
  draws[[name]]
}

summarize_draw <- function(domain, outcome, path, x) {
  data.table(
    domain = domain,
    outcome = outcome,
    path = path,
    mean = mean(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975)),
    probability_positive = mean(x > 0)
  )
}

# Northern direct and mediated effects. Because SC is parallel rather than a mediator
# of floral traits, no SC -> flower path products are computed.
b <- as_draws_df(n_bombus)
w <- as_draws_df(n_wind)
s <- as_draws_df(n_sc)
sh <- as_draws_df(n_showy)
cx <- as_draws_df(n_complex)
nd <- min(nrow(b), nrow(w), nrow(s), nrow(sh), nrow(cx))
b <- b[seq_len(nd), ]
w <- w[seq_len(nd), ]
s <- s[seq_len(nd), ]
sh <- sh[seq_len(nd), ]
cx <- cx[seq_len(nd), ]

i_b <- coef_draw(b, "b_z_log_distance")
i_w <- coef_draw(w, "b_z_log_distance")
b_w <- coef_draw(w, "b_z_bombus_deficit")

i_s <- coef_draw(s, "b_z_log_distance")
b_s <- coef_draw(s, "b_z_bombus_deficit")
w_s <- coef_draw(s, "b_z_wind_share")

i_sh <- coef_draw(sh, "b_z_log_distance")
b_sh <- coef_draw(sh, "b_z_bombus_deficit")
w_sh <- coef_draw(sh, "b_z_wind_share")

i_cx <- coef_draw(cx, "b_z_log_distance")
b_cx <- coef_draw(cx, "b_z_bombus_deficit")
w_cx <- coef_draw(cx, "b_z_wind_share")

north_paths <- rbindlist(list(
  summarize_draw("northern_midlatitude", "Bombus_deficit", "Isolation_direct", i_b),

  summarize_draw("northern_midlatitude", "Wind", "Isolation_direct", i_w),
  summarize_draw("northern_midlatitude", "Wind", "Isolation_to_Bombus_to_Wind", i_b * b_w),
  summarize_draw("northern_midlatitude", "Wind", "Isolation_total", i_w + i_b * b_w),
  summarize_draw("northern_midlatitude", "Wind", "Bombus_direct", b_w),

  summarize_draw("northern_midlatitude", "SC", "Isolation_direct", i_s),
  summarize_draw("northern_midlatitude", "SC", "Isolation_to_Bombus_to_SC", i_b * b_s),
  summarize_draw("northern_midlatitude", "SC", "Isolation_to_Wind_to_SC", i_w * w_s),
  summarize_draw("northern_midlatitude", "SC", "Isolation_to_Bombus_to_Wind_to_SC", i_b * b_w * w_s),
  summarize_draw("northern_midlatitude", "SC", "Bombus_direct", b_s),
  summarize_draw("northern_midlatitude", "SC", "Wind_direct", w_s),

  summarize_draw("northern_midlatitude", "Showy", "Isolation_direct", i_sh),
  summarize_draw("northern_midlatitude", "Showy", "Isolation_to_Bombus_to_Showy", i_b * b_sh),
  summarize_draw("northern_midlatitude", "Showy", "Isolation_to_Wind_to_Showy", i_w * w_sh),
  summarize_draw("northern_midlatitude", "Showy", "Isolation_to_Bombus_to_Wind_to_Showy", i_b * b_w * w_sh),
  summarize_draw("northern_midlatitude", "Showy", "Bombus_direct", b_sh),
  summarize_draw("northern_midlatitude", "Showy", "Wind_direct", w_sh),

  summarize_draw("northern_midlatitude", "Complex", "Isolation_direct", i_cx),
  summarize_draw("northern_midlatitude", "Complex", "Isolation_to_Bombus_to_Complex", i_b * b_cx),
  summarize_draw("northern_midlatitude", "Complex", "Isolation_to_Wind_to_Complex", i_w * w_cx),
  summarize_draw("northern_midlatitude", "Complex", "Isolation_to_Bombus_to_Wind_to_Complex", i_b * b_w * w_cx),
  summarize_draw("northern_midlatitude", "Complex", "Bombus_direct", b_cx),
  summarize_draw("northern_midlatitude", "Complex", "Wind_direct", w_cx)
))
fwrite(north_paths, file.path(outdir, "northern_path_effects.csv"))

# Counter-domain isolation effects on the same flower-syndrome axes.
ssh <- as_draws_df(s_showy)
scx <- as_draws_df(s_complex)
ns <- min(nrow(ssh), nrow(scx))
ssh <- ssh[seq_len(ns), ]
scx <- scx[seq_len(ns), ]
s_i_sh <- coef_draw(ssh, "b_z_log_distance")
s_i_cx <- coef_draw(scx, "b_z_log_distance")

counter_paths <- rbindlist(list(
  summarize_draw("tropical_plus_southern", "Showy", "Isolation_direct", s_i_sh),
  summarize_draw("tropical_plus_southern", "Complex", "Isolation_direct", s_i_cx)
))
fwrite(counter_paths, file.path(outdir, "counterdomain_flower_effects.csv"))

# Direct posterior contrasts of the isolation slopes across domains.
contrast_n <- min(length(i_sh), length(s_i_sh), length(i_cx), length(s_i_cx))
contrasts <- rbindlist(list(
  summarize_draw(
    "south_minus_north", "Showy", "Isolation_slope_contrast",
    s_i_sh[seq_len(contrast_n)] - i_sh[seq_len(contrast_n)]
  ),
  summarize_draw(
    "south_minus_north", "Complex", "Isolation_slope_contrast",
    s_i_cx[seq_len(contrast_n)] - i_cx[seq_len(contrast_n)]
  )
))
fwrite(contrasts, file.path(outdir, "counterdomain_contrasts_vs_north.csv"))

sampler_diagnostics <- rbindlist(lapply(names(fits), function(nm) {
  np <- nuts_params(fits[[nm]])
  diag <- as.data.table(summarise_draws(as_draws_array(fits[[nm]]), rhat, ess_bulk, ess_tail))[
    grepl("^b_", variable)
  ]
  data.table(
    model = nm,
    divergent_transitions = sum(np$Parameter == "divergent__" & np$Value == 1),
    max_treedepth_hits = sum(np$Parameter == "treedepth__" & np$Value >= 15),
    max_rhat = if (nrow(diag) == 0 || all(is.na(diag$rhat))) NA_real_ else max(diag$rhat, na.rm = TRUE),
    min_bulk_ess = if (nrow(diag) == 0 || all(is.na(diag$ess_bulk))) NA_real_ else min(diag$ess_bulk, na.rm = TRUE),
    min_tail_ess = if (nrow(diag) == 0 || all(is.na(diag$ess_tail))) NA_real_ else min(diag$ess_tail, na.rm = TRUE)
  )
}), fill = TRUE)
fwrite(sampler_diagnostics, file.path(outdir, "sampler_diagnostics.csv"))
if (any(sampler_diagnostics$divergent_transitions > 0)) warning("Divergent transitions detected; inspect sampler_diagnostics.csv")
if (any(sampler_diagnostics$max_rhat > 1.01, na.rm = TRUE)) warning("R-hat > 1.01 detected; inspect sampler_diagnostics.csv")

capture.output(lapply(fits, summary), file = file.path(outdir, "model_summaries.txt"))

meta <- list(
  contract = "v2_bayesian_noncircular_flower_syndrome_test_v1",
  northern_hypothesis = "Isolation -> Bombus deficit; isolation/Bombus deficit -> wind and SC; isolation/Bombus deficit/wind -> showy colour and complex floral form",
  counterdomain_hypothesis = "Tropical plus southern islands provide a counter-domain for isolation effects on the same showy/plain and complex/simple floral axes",
  colour_axis = "plain = white/green/brown/inconspicuous; showy = all other resolved colours",
  form_axis = "complex = tubular/trumpet + zygomorphic/specialized; simple = open/generalized + composite/brush",
  sc_causal_guardrail = "SC is a parallel response and is not used as a predictor of floral traits",
  pollinator_guild_guardrail = "Plant-derived bird/butterfly/moth/bat guild summaries are excluded from causal flower-trait models",
  adjustment_covariates = "log island area and climate PCs in every equation",
  founder_effect_covariate = "not included: no independent founder-effect variable is available in the locked input",
  analysis_tier = "sensitivity_all",
  n_input_islands = n_input,
  n_zero_or_negative_distance_excluded = n_zero_distance,
  sampler = "2 chains x 2000 iterations, 1000 warmup, adapt_delta=0.99, max_treedepth=15"
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
