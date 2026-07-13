#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(data.table)
  library(posterior)
  library(loo)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: run_bayesian_m0_m4_main.R <input.csv> <output_dir>")
infile <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

d <- fread(infile)

# Main analysis only: no evidence-tier sensitivity analyses here.
# M0-M3 use northern mid-latitude islands, matching the v2 Bombus mechanism scope.
north <- d[analysis_regime == "northern_midlatitude" &
           is.finite(log_distance_to_continent_km) &
           is.finite(log_island_area_km2) &
           is.finite(bombus_deficit)]

z <- function(x) as.numeric(scale(x))
for (nm in c("log_distance_to_continent_km", "log_island_area_km2",
             "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4")) {
  north[[paste0("z_", nm)]] <- z(north[[nm]])
  d[[paste0("z_", nm)]] <- z(d[[nm]])
}
north[, z_bombus_deficit := z(bombus_deficit)]
north[, z_sc_share := z(sc_share)]
d[, z_showy_alt_guild_share := z(showy_alt_guild_share)]
d[, z_other_bee_guild_share := z(other_bee_guild_share)]
d[, z_generalist_insect_guild_share := z(generalist_insect_guild_share)]

g <- "z_log_distance_to_continent_km * z_log_island_area_km2 + climate_pc1 + climate_pc2 + climate_pc3 + climate_pc4"

fit_bb <- function(formula, data, file) {
  brm(
    formula = formula,
    data = data,
    family = beta_binomial(link = "logit"),
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1.5), class = "Intercept"),
      prior(exponential(1), class = "phi")
    ),
    chains = 2, cores = 2, iter = 1200, warmup = 600,
    seed = 20260713, control = list(adapt_delta = 0.95),
    save_pars = save_pars(all = TRUE), file = file,
    refresh = 100
  )
}

fit_beta <- function(formula, data, file) {
  dd <- copy(data)
  lhs <- all.vars(formula)[1]
  dd[[lhs]] <- pmin(pmax(dd[[lhs]], 1e-5), 1 - 1e-5)
  brm(
    formula = formula,
    data = dd,
    family = Beta(link = "logit"),
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1.5), class = "Intercept"),
      prior(exponential(1), class = "phi")
    ),
    chains = 2, cores = 2, iter = 1200, warmup = 600,
    seed = 20260713, control = list(adapt_delta = 0.95),
    save_pars = save_pars(all = TRUE), file = file,
    refresh = 100
  )
}

# M0: geographic filtering, with the key v2 extension distance x island area.
m0_plain <- fit_bb(as.formula(paste("color_plain | trials(color_trials) ~", g)),
                   north[color_trials > 0], file.path(outdir, "m0_plain"))
m0_form <- fit_bb(as.formula(paste("form_open_generalized | trials(form_trials) ~", g)),
                  north[form_trials > 0], file.path(outdir, "m0_form"))

# M1: isolation/area -> Bombus deficit; Bombus deficit adds to floral traits.
m1_bombus <- fit_beta(as.formula(paste("bombus_deficit ~", g)),
                      north, file.path(outdir, "m1_bombus"))
m1_plain <- fit_bb(as.formula(paste("color_plain | trials(color_trials) ~ z_bombus_deficit +", g)),
                   north[color_trials > 0], file.path(outdir, "m1_plain"))
m1_form <- fit_bb(as.formula(paste("form_open_generalized | trials(form_trials) ~ z_bombus_deficit +", g)),
                  north[form_trials > 0], file.path(outdir, "m1_form"))

# M2: Bombus deficit -> reproductive assurance; reproductive assurance -> floral traits.
m2_sc <- fit_bb(as.formula(paste("sc_successes | trials(sc_trials) ~ z_bombus_deficit +", g)),
                north[sc_trials > 0], file.path(outdir, "m2_sc"))
m2_plain <- fit_bb(as.formula(paste("color_plain | trials(color_trials) ~ z_sc_share +", g)),
                   north[color_trials > 0 & is.finite(z_sc_share)], file.path(outdir, "m2_plain"))
m2_form <- fit_bb(as.formula(paste("form_open_generalized | trials(form_trials) ~ z_sc_share +", g)),
                  north[form_trials > 0 & is.finite(z_sc_share)], file.path(outdir, "m2_form"))

# M3: direct Bombus path plus reproductive-assurance path.
m3_plain <- fit_bb(as.formula(paste("color_plain | trials(color_trials) ~ z_bombus_deficit + z_sc_share +", g)),
                   north[color_trials > 0 & is.finite(z_sc_share)], file.path(outdir, "m3_plain"))
m3_form <- fit_bb(as.formula(paste("form_open_generalized | trials(form_trials) ~ z_bombus_deficit + z_sc_share +", g)),
                  north[form_trials > 0 & is.finite(z_sc_share)], file.path(outdir, "m3_form"))

# M4: preserve v2 category information. Multinomial colour and form composition are
# modelled globally as functions of distance x area x pollination regime plus the
# three alternative-pollinator composition axes. This is the replacement/regime test.
m4 <- d[analysis_regime %in% c("northern_midlatitude", "tropical", "southern_extratropical")]
m4[, analysis_regime := factor(analysis_regime)]

m4_color <- brm(
  cbind(color_plain, color_yellow_orange, color_red_pink, color_blue_purple) ~
    z_log_distance_to_continent_km * z_log_island_area_km2 * analysis_regime +
    z_showy_alt_guild_share + z_other_bee_guild_share + z_generalist_insect_guild_share,
  data = m4[color_trials > 0 & complete.cases(z_showy_alt_guild_share, z_other_bee_guild_share, z_generalist_insect_guild_share)],
  family = multinomial(), prior = prior(normal(0, 1), class = "b"),
  chains = 2, cores = 2, iter = 1200, warmup = 600, seed = 20260713,
  control = list(adapt_delta = 0.95), save_pars = save_pars(all = TRUE),
  file = file.path(outdir, "m4_color"), refresh = 100
)

m4_form <- brm(
  cbind(form_open_generalized, form_tubular_trumpet, form_zygomorphic_specialized, form_composite_brush) ~
    z_log_distance_to_continent_km * z_log_island_area_km2 * analysis_regime +
    z_showy_alt_guild_share + z_other_bee_guild_share + z_generalist_insect_guild_share,
  data = m4[form_trials > 0 & complete.cases(z_showy_alt_guild_share, z_other_bee_guild_share, z_generalist_insect_guild_share)],
  family = multinomial(), prior = prior(normal(0, 1), class = "b"),
  chains = 2, cores = 2, iter = 1200, warmup = 600, seed = 20260713,
  control = list(adapt_delta = 0.95), save_pars = save_pars(all = TRUE),
  file = file.path(outdir, "m4_form"), refresh = 100
)

fits <- list(
  M0_plain=m0_plain, M0_form=m0_form,
  M1_bombus=m1_bombus, M1_plain=m1_plain, M1_form=m1_form,
  M2_SC=m2_sc, M2_plain=m2_plain, M2_form=m2_form,
  M3_plain=m3_plain, M3_form=m3_form,
  M4_color=m4_color, M4_form=m4_form
)

summaries <- rbindlist(lapply(names(fits), function(nm) {
  x <- as.data.table(posterior_summary(fits[[nm]], pars = "^b_"), keep.rownames = "parameter")
  x[, model := nm]
  setcolorder(x, c("model", "parameter"))
  x
}), fill = TRUE)
fwrite(summaries, file.path(outdir, "posterior_fixed_effects.csv"))

# Compare only models with identical responses. M4 is a category-preserving extension,
# not forced into the same LOO table as the planned binary contrasts.
loo_rows <- list()
for (response in c("plain", "form")) {
  ms <- list(
    M0 = fits[[paste0("M0_", response)]],
    M1 = fits[[paste0("M1_", response)]],
    M2 = fits[[paste0("M2_", response)]],
    M3 = fits[[paste0("M3_", response)]]
  )
  loos <- lapply(ms, loo)
  cmp <- as.data.table(loo_compare(loos), keep.rownames = "scenario")
  cmp[, response := response]
  loo_rows[[response]] <- cmp
}
fwrite(rbindlist(loo_rows, fill = TRUE), file.path(outdir, "m0_m3_loo_comparison.csv"))

# Posterior products for the main mediated path: Bombus deficit -> SC -> floral trait.
sc_draws <- as_draws_df(m2_sc)
plain_draws <- as_draws_df(m3_plain)
form_draws <- as_draws_df(m3_form)
n <- min(nrow(sc_draws), nrow(plain_draws), nrow(form_draws))
indirect <- data.table(
  path = c("Bombus_deficit_to_SC_to_plain", "Bombus_deficit_to_SC_to_generalized"),
  mean = c(
    mean(sc_draws$b_z_bombus_deficit[1:n] * plain_draws$b_z_sc_share[1:n]),
    mean(sc_draws$b_z_bombus_deficit[1:n] * form_draws$b_z_sc_share[1:n])
  ),
  q025 = c(
    quantile(sc_draws$b_z_bombus_deficit[1:n] * plain_draws$b_z_sc_share[1:n], 0.025),
    quantile(sc_draws$b_z_bombus_deficit[1:n] * form_draws$b_z_sc_share[1:n], 0.025)
  ),
  q975 = c(
    quantile(sc_draws$b_z_bombus_deficit[1:n] * plain_draws$b_z_sc_share[1:n], 0.975),
    quantile(sc_draws$b_z_bombus_deficit[1:n] * form_draws$b_z_sc_share[1:n], 0.975)
  )
)
fwrite(indirect, file.path(outdir, "m2_m3_indirect_effects.csv"))

meta <- list(
  contract = "v2_bayesian_m0_m4_main_v1",
  analysis_tier = "sensitivity_all_all_filled_data_requested_first_run",
  method = "Bayesian piecewise path models; beta-binomial planned contrasts plus multinomial M4 category composition",
  main_interaction = "log distance to continent x log island area",
  n_northern_midlatitude = nrow(north),
  n_global_m4 = nrow(m4),
  sensitivity_analysis = "deferred by design",
  note = "Initial main fit uses 2 chains x 1200 iterations; convergence diagnostics are reported with outputs and should be checked before manuscript claims."
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)

capture.output(lapply(fits, summary), file = file.path(outdir, "model_summaries.txt"))
