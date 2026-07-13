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

# Isolation is the estimand: zero-distance islands are not treated as zero-dose
# observations. They are excluded before log transformation and standardisation.
n_input <- nrow(d)
n_zero_distance <- sum(is.finite(d$distance_to_continent_km) & d$distance_to_continent_km <= 0, na.rm = TRUE)
d <- d[is.finite(distance_to_continent_km) & distance_to_continent_km > 0 &
         is.finite(area_km2) & area_km2 > 0]
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
d[, z_bombus_deficit := z(bombus_deficit)]
d[, z_sc_share := z(sc_share)]
d[, z_showy_alt_guild_share := z(showy_alt_guild_share)]
d[, z_other_bee_guild_share := z(other_bee_guild_share)]
d[, z_generalist_insect_guild_share := z(generalist_insect_guild_share)]
d[, bombus_primary_regime := as.integer(analysis_regime == "northern_midlatitude")]
d[, analysis_regime := factor(analysis_regime)]

geo <- paste("z_log_distance * z_log_area + z_climate_pc1 + z_climate_pc2 +",
             "z_climate_pc3 + z_climate_pc4 + analysis_regime")
alt <- paste("z_showy_alt_guild_share + z_other_bee_guild_share + z_generalist_insect_guild_share")

fit_bb <- function(formula, data, file) {
  if (nrow(data) < 30) stop(paste("Too few rows for", deparse(formula)))
  brm(
    formula = formula, data = data, family = beta_binomial(link = "logit"),
    prior = c(prior(normal(0, 1), class = "b"),
              prior(normal(0, 1.5), class = "Intercept"),
              prior(exponential(1), class = "phi")),
    chains = 2, cores = 2, iter = 2000, warmup = 1000, seed = 20260713,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = save_pars(all = TRUE), file = file, refresh = 100
  )
}

# Maximal available support per equation. Missing form data no longer removes an
# otherwise usable island from the colour or SC equations, and vice versa.
base_complete <- complete.cases(d[, .(z_log_distance, z_log_area, z_climate_pc1,
  z_climate_pc2, z_climate_pc3, z_climate_pc4, analysis_regime)])
sc_data <- d[base_complete & sc_trials > 0 & is.finite(z_bombus_deficit) &
  is.finite(z_showy_alt_guild_share) & is.finite(z_other_bee_guild_share) &
  is.finite(z_generalist_insect_guild_share)]
color_data <- d[base_complete & color_trials > 0 & is.finite(z_bombus_deficit) &
  is.finite(z_sc_share) & is.finite(z_showy_alt_guild_share) &
  is.finite(z_other_bee_guild_share) & is.finite(z_generalist_insect_guild_share)]
form_data <- d[base_complete & form_trials > 0 & is.finite(z_bombus_deficit) &
  is.finite(z_sc_share) & is.finite(z_showy_alt_guild_share) &
  is.finite(z_other_bee_guild_share) & is.finite(z_generalist_insect_guild_share)]

# M0: global isolation baseline with an explicit isolation x region test.
m0_sc <- fit_bb(as.formula(paste("sc_successes | trials(sc_trials) ~ z_log_distance * analysis_regime + z_log_area + z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4")), sc_data, file.path(outdir, "m0_global_isolation_sc"))
m0_color <- fit_bb(as.formula(paste("color_plain | trials(color_trials) ~ z_log_distance * analysis_regime + z_log_area + z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4")), color_data, file.path(outdir, "m0_global_isolation_color"))
m0_form <- fit_bb(as.formula(paste("form_open_generalized | trials(form_trials) ~ z_log_distance * analysis_regime + z_log_area + z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4")), form_data, file.path(outdir, "m0_global_isolation_form"))

# M1: test whether the Bombus pathway is concentrated in the northern-midlatitude
# focal regime rather than assumed to be a universal global mechanism.
m1_sc <- fit_bb(as.formula(paste("sc_successes | trials(sc_trials) ~ z_bombus_deficit * bombus_primary_regime +", geo)), sc_data, file.path(outdir, "m1_bombus_regime_sc"))
m1_color <- fit_bb(as.formula(paste("color_plain | trials(color_trials) ~ z_bombus_deficit * bombus_primary_regime +", geo)), color_data, file.path(outdir, "m1_bombus_regime_color"))
m1_form <- fit_bb(as.formula(paste("form_open_generalized | trials(form_trials) ~ z_bombus_deficit * bombus_primary_regime +", geo)), form_data, file.path(outdir, "m1_bombus_regime_form"))

# M2: counterfactual alternative channels. showy_alt is the frozen ontology group
# containing birds, butterflies, moths and bats; other bees and generalist insects
# are retained separately.
m2_sc <- fit_bb(as.formula(paste("sc_successes | trials(sc_trials) ~", alt, "+", geo)), sc_data, file.path(outdir, "m2_alternative_channels_sc"))
m2_color <- fit_bb(as.formula(paste("color_plain | trials(color_trials) ~", alt, "+", geo)), color_data, file.path(outdir, "m2_alternative_channels_color"))
m2_form <- fit_bb(as.formula(paste("form_open_generalized | trials(form_trials) ~", alt, "+", geo)), form_data, file.path(outdir, "m2_alternative_channels_form"))

# M3: integrated Bayesian piecewise path model with direct and mediated paths.
m3_sc <- fit_bb(as.formula(paste("sc_successes | trials(sc_trials) ~ z_bombus_deficit * bombus_primary_regime +", alt, "+", geo)), sc_data, file.path(outdir, "m3_integrated_sc"))
m3_color <- fit_bb(as.formula(paste("color_plain | trials(color_trials) ~ z_bombus_deficit * bombus_primary_regime + z_sc_share * bombus_primary_regime +", alt, "+", geo)), color_data, file.path(outdir, "m3_integrated_color"))
m3_form <- fit_bb(as.formula(paste("form_open_generalized | trials(form_trials) ~ z_bombus_deficit * bombus_primary_regime + z_sc_share * bombus_primary_regime +", alt, "+", geo)), form_data, file.path(outdir, "m3_integrated_form"))

fits <- list(M0_SC=m0_sc, M0_color=m0_color, M0_form=m0_form,
  M1_SC=m1_sc, M1_color=m1_color, M1_form=m1_form,
  M2_SC=m2_sc, M2_color=m2_color, M2_form=m2_form,
  M3_SC=m3_sc, M3_color=m3_color, M3_form=m3_form)

posterior_rows <- rbindlist(lapply(names(fits), function(nm) {
  x <- as.data.table(posterior_summary(fits[[nm]], pars = "^b_"), keep.rownames = "parameter")
  x[, model := nm]
  setcolorder(x, c("model", "parameter"))
  x
}), fill = TRUE)
fwrite(posterior_rows, file.path(outdir, "posterior_fixed_effects.csv"))

support <- data.table(
  equation = c("SC", "plain_color", "generalized_form"),
  n_islands = c(nrow(sc_data), nrow(color_data), nrow(form_data)),
  n_northern_midlatitude = c(sum(sc_data$bombus_primary_regime == 1), sum(color_data$bombus_primary_regime == 1), sum(form_data$bombus_primary_regime == 1)),
  n_non_northern_counterfactual = c(sum(sc_data$bombus_primary_regime == 0), sum(color_data$bombus_primary_regime == 0), sum(form_data$bombus_primary_regime == 0))
)
fwrite(support, file.path(outdir, "equation_support.csv"))

# Posterior products for mediated paths. The non-northern effect is the main
# slope; the northern-midlatitude effect adds the binary-regime interaction.
sc_draws <- as_draws_df(m3_sc)
color_draws <- as_draws_df(m3_color)
form_draws <- as_draws_df(m3_form)
n <- min(nrow(sc_draws), nrow(color_draws), nrow(form_draws))
sc_draws <- sc_draws[1:n, ]; color_draws <- color_draws[1:n, ]; form_draws <- form_draws[1:n, ]
coef_draw <- function(draws, name) {
  if (!name %in% names(draws)) stop(paste("Missing posterior coefficient", name))
  draws[[name]]
}
summarize_draw <- function(path, x) data.table(path=path, mean=mean(x),
  q025=unname(quantile(x, 0.025)), q975=unname(quantile(x, 0.975)),
  probability_positive=mean(x > 0))

b_sc_non <- coef_draw(sc_draws, "b_z_bombus_deficit")
b_sc_north <- b_sc_non + coef_draw(sc_draws, "b_z_bombus_deficit:bombus_primary_regime")
b_color_sc_non <- coef_draw(color_draws, "b_z_sc_share")
b_color_sc_north <- b_color_sc_non + coef_draw(color_draws, "b_z_sc_share:bombus_primary_regime")
b_form_sc_non <- coef_draw(form_draws, "b_z_sc_share")
b_form_sc_north <- b_form_sc_non + coef_draw(form_draws, "b_z_sc_share:bombus_primary_regime")

indirect <- rbindlist(list(
  summarize_draw("Bombus_deficit_to_SC_to_plain__non_northern_counterfactual", b_sc_non * b_color_sc_non),
  summarize_draw("Bombus_deficit_to_SC_to_plain__northern_midlatitude", b_sc_north * b_color_sc_north),
  summarize_draw("Bombus_deficit_to_SC_to_generalized__non_northern_counterfactual", b_sc_non * b_form_sc_non),
  summarize_draw("Bombus_deficit_to_SC_to_generalized__northern_midlatitude", b_sc_north * b_form_sc_north),
  summarize_draw("showy_alt_to_SC_to_plain", coef_draw(sc_draws, "b_z_showy_alt_guild_share") * b_color_sc_non),
  summarize_draw("showy_alt_to_SC_to_generalized", coef_draw(sc_draws, "b_z_showy_alt_guild_share") * b_form_sc_non),
  summarize_draw("other_bees_to_SC_to_plain", coef_draw(sc_draws, "b_z_other_bee_guild_share") * b_color_sc_non),
  summarize_draw("generalist_insects_to_SC_to_plain", coef_draw(sc_draws, "b_z_generalist_insect_guild_share") * b_color_sc_non)
))
fwrite(indirect, file.path(outdir, "m3_indirect_effects_by_regime.csv"))
capture.output(lapply(fits, summary), file = file.path(outdir, "model_summaries.txt"))

meta <- list(
  contract = "v2_bayesian_global_isolation_counterfactual_v1",
  method = "Bayesian piecewise path models with regime-moderated Bombus pathway and alternative-pollinator counterfactual",
  analysis_tier = "sensitivity_all_all_filled_data",
  n_input_islands = n_input,
  n_zero_or_negative_distance_excluded = n_zero_distance,
  n_positive_distance_islands_retained_before_equation_missingness = nrow(d),
  distance_rule = "distance_to_continent_km > 0; log(distance_to_continent_km)",
  support_strategy = "maximum available islands per equation; no forced color+form+SC complete-case intersection",
  bombus_scope_test = "z_bombus_deficit x northern_midlatitude indicator",
  alternative_channel_counterfactual = c("showy_alt: birds + butterflies + moths + bats", "other bees", "generalist insects"),
  interpretation_guardrail = "Bombus deficit outside the northern-midlatitude focal regime is a moderated counterfactual; a null or reversed pathway there is evidence against globally universalising the northern-temperate syndrome.",
  sampler = "2 chains x 2000 iterations, 1000 warmup, adapt_delta=0.99, max_treedepth=15"
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty=TRUE, auto_unbox=TRUE)
