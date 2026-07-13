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
d <- d[is.finite(distance_to_continent_km) & distance_to_continent_km > 0 & is.finite(area_km2) & area_km2 > 0]
if (nrow(d) < 30) stop("Too few positive-distance islands")

z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) stop("Cannot standardize zero-variance predictor")
  (x - mean(x, na.rm = TRUE)) / s
}

d[, z_log_distance := z(log(distance_to_continent_km))]
d[, z_log_area := z(log(area_km2))]
for (nm in c("climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4")) d[[paste0("z_", nm)]] <- z(d[[nm]])
d[, wind_share := fifelse(pollen_vector_trials > 0, n_wind_species / pollen_vector_trials, NA_real_)]
d[, vivid_successes := color_trials - color_plain]
d[, tubular_successes := form_tubular_trumpet]
d[, z_bombus_deficit := z(bombus_deficit)]
d[, z_wind_share := z(wind_share)]
d[, z_sc_share := z(sc_share)]
d[, z_mobile_alt := z(showy_alt_guild_share)]

fit_bb <- function(formula, data, file) {
  if (nrow(data) < 30) stop(paste("Too few rows for", deparse(formula)))
  brm(
    formula = formula, data = data, family = beta_binomial(link = "logit"),
    prior = c(prior(normal(0, 1), class = "b"), prior(normal(0, 1.5), class = "Intercept"), prior(exponential(1), class = "phi")),
    chains = 2, cores = 2, iter = 2000, warmup = 1000, seed = 20260713,
    control = list(adapt_delta = 0.99, max_treedepth = 15), save_pars = save_pars(all = TRUE), file = file, refresh = 100
  )
}

fit_gaussian <- function(formula, data, file) {
  if (nrow(data) < 30) stop(paste("Too few rows for", deparse(formula)))
  brm(
    formula = formula, data = data, family = gaussian(),
    prior = c(prior(normal(0, 1), class = "b"), prior(normal(0, 1.5), class = "Intercept"), prior(exponential(1), class = "sigma")),
    chains = 2, cores = 2, iter = 2000, warmup = 1000, seed = 20260713,
    control = list(adapt_delta = 0.99, max_treedepth = 15), save_pars = save_pars(all = TRUE), file = file, refresh = 100
  )
}

covars <- c("z_log_distance", "z_log_area", "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4")

# Northern focal mechanism: isolation -> Bombus deficit -> wind -> SC -> plain flowers.
north <- d[analysis_regime == "northern_midlatitude"]
n_base <- complete.cases(north[, ..covars]) & is.finite(north$z_bombus_deficit)
n_bombus_dat <- north[n_base]
n_wind_dat <- north[n_base & pollen_vector_trials > 0]
n_sc_dat <- north[n_base & sc_trials > 0 & is.finite(z_wind_share)]
n_plain_dat <- north[n_base & color_trials > 0 & is.finite(z_wind_share) & is.finite(z_sc_share)]

n_bombus <- fit_gaussian(
  z_bombus_deficit ~ z_log_distance + z_log_area + z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  n_bombus_dat, file.path(outdir, "north_isolation_bombus")
)
n_wind <- fit_bb(
  n_wind_species | trials(pollen_vector_trials) ~ z_bombus_deficit + z_log_distance + z_log_area + z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  n_wind_dat, file.path(outdir, "north_bombus_wind")
)
n_sc <- fit_bb(
  sc_successes | trials(sc_trials) ~ z_bombus_deficit + z_wind_share + z_log_distance + z_log_area + z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  n_sc_dat, file.path(outdir, "north_bombus_wind_sc")
)
n_plain <- fit_bb(
  color_plain | trials(color_trials) ~ z_bombus_deficit + z_wind_share + z_sc_share + z_log_distance + z_log_area + z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  n_plain_dat, file.path(outdir, "north_bombus_wind_sc_plain")
)

# Tropical + southern alternative-pollinator mechanism:
# isolation -> mobile bird/butterfly/moth guild share -> vivid colour and tubular form.
south <- d[analysis_regime %in% c("tropical", "southern_extratropical")]
s_base <- complete.cases(south[, ..covars]) & is.finite(south$z_mobile_alt)
s_mobile_dat <- south[s_base]
s_vivid_dat <- south[s_base & color_trials > 0]
s_tubular_dat <- south[s_base & form_trials > 0]

s_mobile <- fit_gaussian(
  z_mobile_alt ~ z_log_distance + z_log_area + z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  s_mobile_dat, file.path(outdir, "south_isolation_mobile_alt")
)
s_vivid <- fit_bb(
  vivid_successes | trials(color_trials) ~ z_mobile_alt + z_log_distance + z_log_area + z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  s_vivid_dat, file.path(outdir, "south_mobile_alt_vivid")
)
s_tubular <- fit_bb(
  tubular_successes | trials(form_trials) ~ z_mobile_alt + z_log_distance + z_log_area + z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  s_tubular_dat, file.path(outdir, "south_mobile_alt_tubular")
)

fits <- list(N_Bombus=n_bombus, N_wind=n_wind, N_SC=n_sc, N_plain=n_plain, S_mobile_alt=s_mobile, S_vivid=s_vivid, S_tubular=s_tubular)

posterior_rows <- rbindlist(lapply(names(fits), function(nm) {
  x <- as.data.table(posterior_summary(fits[[nm]], pars = "^b_"), keep.rownames = "parameter")
  x[, model := nm]
  setcolorder(x, c("model", "parameter"))
  x
}), fill = TRUE)
fwrite(posterior_rows, file.path(outdir, "posterior_fixed_effects.csv"))

support <- data.table(
  equation = c("north_Bombus", "north_wind", "north_SC", "north_plain", "south_mobile_alt", "south_vivid", "south_tubular"),
  n_islands = c(nrow(n_bombus_dat), nrow(n_wind_dat), nrow(n_sc_dat), nrow(n_plain_dat), nrow(s_mobile_dat), nrow(s_vivid_dat), nrow(s_tubular_dat))
)
fwrite(support, file.path(outdir, "equation_support.csv"))

coef_draw <- function(draws, name) {
  if (!name %in% names(draws)) stop(paste("Missing posterior coefficient", name))
  draws[[name]]
}
summarize_draw <- function(domain, path, x) data.table(domain=domain, path=path, mean=mean(x), q025=unname(quantile(x, .025)), q975=unname(quantile(x, .975)), probability_positive=mean(x > 0))

# Northern path products.
b <- as_draws_df(n_bombus); w <- as_draws_df(n_wind); s <- as_draws_df(n_sc); p <- as_draws_df(n_plain)
nd <- min(nrow(b), nrow(w), nrow(s), nrow(p)); b <- b[1:nd,]; w <- w[1:nd,]; s <- s[1:nd,]; p <- p[1:nd,]
i_b <- coef_draw(b, "b_z_log_distance"); i_w <- coef_draw(w, "b_z_log_distance"); b_w <- coef_draw(w, "b_z_bombus_deficit")
i_s <- coef_draw(s, "b_z_log_distance"); b_s <- coef_draw(s, "b_z_bombus_deficit"); w_s <- coef_draw(s, "b_z_wind_share")
i_p <- coef_draw(p, "b_z_log_distance"); b_p <- coef_draw(p, "b_z_bombus_deficit"); w_p <- coef_draw(p, "b_z_wind_share"); s_p <- coef_draw(p, "b_z_sc_share")

north_paths <- rbindlist(list(
  summarize_draw("northern_midlatitude", "Isolation_to_Plain_direct", i_p),
  summarize_draw("northern_midlatitude", "Isolation_to_Bombus_to_Plain", i_b*b_p),
  summarize_draw("northern_midlatitude", "Isolation_to_Wind_to_Plain", i_w*w_p),
  summarize_draw("northern_midlatitude", "Isolation_to_SC_to_Plain", i_s*s_p),
  summarize_draw("northern_midlatitude", "Isolation_to_Bombus_to_Wind_to_Plain", i_b*b_w*w_p),
  summarize_draw("northern_midlatitude", "Isolation_to_Bombus_to_SC_to_Plain", i_b*b_s*s_p),
  summarize_draw("northern_midlatitude", "Isolation_to_Wind_to_SC_to_Plain", i_w*w_s*s_p),
  summarize_draw("northern_midlatitude", "Isolation_to_Bombus_to_Wind_to_SC_to_Plain", i_b*b_w*w_s*s_p),
  summarize_draw("northern_midlatitude", "Bombus_to_Plain_direct", b_p),
  summarize_draw("northern_midlatitude", "Bombus_to_Wind_to_Plain", b_w*w_p),
  summarize_draw("northern_midlatitude", "Bombus_to_SC_to_Plain", b_s*s_p),
  summarize_draw("northern_midlatitude", "Bombus_to_Wind_to_SC_to_Plain", b_w*w_s*s_p)
))
fwrite(north_paths, file.path(outdir, "northern_path_effects.csv"))

# Alternative-pollinator direct, indirect, and total effects.
m <- as_draws_df(s_mobile); v <- as_draws_df(s_vivid); t <- as_draws_df(s_tubular)
nd2 <- min(nrow(m), nrow(v), nrow(t)); m <- m[1:nd2,]; v <- v[1:nd2,]; t <- t[1:nd2,]
i_m <- coef_draw(m, "b_z_log_distance")
i_v <- coef_draw(v, "b_z_log_distance"); m_v <- coef_draw(v, "b_z_mobile_alt")
i_t <- coef_draw(t, "b_z_log_distance"); m_t <- coef_draw(t, "b_z_mobile_alt")
alt_paths <- rbindlist(list(
  summarize_draw("tropical_plus_southern", "Isolation_to_Vivid_direct", i_v),
  summarize_draw("tropical_plus_southern", "Isolation_to_MobileAlt_to_Vivid", i_m*m_v),
  summarize_draw("tropical_plus_southern", "Isolation_to_Vivid_total", i_v + i_m*m_v),
  summarize_draw("tropical_plus_southern", "MobileAlt_to_Vivid_direct", m_v),
  summarize_draw("tropical_plus_southern", "Isolation_to_Tubular_direct", i_t),
  summarize_draw("tropical_plus_southern", "Isolation_to_MobileAlt_to_Tubular", i_m*m_t),
  summarize_draw("tropical_plus_southern", "Isolation_to_Tubular_total", i_t + i_m*m_t),
  summarize_draw("tropical_plus_southern", "MobileAlt_to_Tubular_direct", m_t)
))
fwrite(alt_paths, file.path(outdir, "alternative_pollinator_path_effects.csv"))

sampler_diagnostics <- rbindlist(lapply(names(fits), function(nm) {
  np <- nuts_params(fits[[nm]])
  diag <- as.data.table(summarise_draws(as_draws_array(fits[[nm]]), rhat, ess_bulk, ess_tail))[grepl("^b_", variable)]
  data.table(
    model=nm,
    divergent_transitions=sum(np$Parameter == "divergent__" & np$Value == 1),
    max_treedepth_hits=sum(np$Parameter == "treedepth__" & np$Value >= 15),
    max_rhat=if (nrow(diag) == 0 || all(is.na(diag$rhat))) NA_real_ else max(diag$rhat, na.rm=TRUE),
    min_bulk_ess=if (nrow(diag) == 0 || all(is.na(diag$ess_bulk))) NA_real_ else min(diag$ess_bulk, na.rm=TRUE),
    min_tail_ess=if (nrow(diag) == 0 || all(is.na(diag$ess_tail))) NA_real_ else min(diag$ess_tail, na.rm=TRUE)
  )
}), fill=TRUE)
fwrite(sampler_diagnostics, file.path(outdir, "sampler_diagnostics.csv"))
if (any(sampler_diagnostics$divergent_transitions > 0)) warning("Divergent transitions detected; inspect sampler_diagnostics.csv")
if (any(sampler_diagnostics$max_rhat > 1.01, na.rm=TRUE)) warning("R-hat > 1.01 detected; inspect sampler_diagnostics.csv")

capture.output(lapply(fits, summary), file=file.path(outdir, "model_summaries.txt"))
meta <- list(
  contract="v2_bayesian_shortest_two-mechanism_test_v1",
  northern_hypothesis="Island isolation -> Bombus deficit -> wind pollination / self-compatibility -> plain floral syndrome, with direct paths retained",
  tropical_southern_hypothesis="Island isolation -> mobile alternative pollinator guild share (birds, butterflies, moths, bats) -> vivid colours and tubular flowers, with direct isolation effects retained",
  adjustment_covariates="log island area and climate PCs in every equation",
  founder_effect_covariate="not included: no independent founder-effect variable is available in the locked input",
  analysis_tier="sensitivity_all",
  n_input_islands=n_input,
  n_zero_or_negative_distance_excluded=n_zero_distance,
  sampler="2 chains x 2000 iterations, 1000 warmup, adapt_delta=0.99, max_treedepth=15"
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty=TRUE, auto_unbox=TRUE)
