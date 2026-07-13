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

# Main analysis only: sensitivity analyses are deliberately deferred.
# M0-M3 use one common northern-midlatitude island set so scenario LOO values
# are comparable. M4 keeps the detailed v2 colour/form categories as a separate
# pollination-regime extension and is not forced into the M0-M3 LOO ranking.
north <- d[
  analysis_regime == "northern_midlatitude" &
    is.finite(log_distance_to_continent_km) &
    is.finite(log_island_area_km2) &
    is.finite(bombus_deficit)
]

z <- function(x) as.numeric(scale(x))
standardize_cols <- c(
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"
)
for (nm in standardize_cols) {
  north[[paste0("z_", nm)]] <- z(north[[nm]])
  d[[paste0("z_", nm)]] <- z(d[[nm]])
}
north[, z_bombus_deficit := z(bombus_deficit)]
north[, z_sc_share := z(sc_share)]
d[, z_showy_alt_guild_share := z(showy_alt_guild_share)]
d[, z_other_bee_guild_share := z(other_bee_guild_share)]
d[, z_generalist_insect_guild_share := z(generalist_insect_guild_share)]

# One common support for every M0-M3 component and response.
north_main <- north[
  color_trials > 0 &
    form_trials > 0 &
    sc_trials > 0 &
    is.finite(z_sc_share) &
    complete.cases(
      z_log_distance_to_continent_km,
      z_log_island_area_km2,
      z_climate_pc1, z_climate_pc2, z_climate_pc3, z_climate_pc4,
      z_bombus_deficit
    )
]
if (nrow(north_main) < 30) stop("Too few northern-midlatitude islands on common M0-M3 support")

geo <- paste(
  "z_log_distance_to_continent_km * z_log_island_area_km2 +",
  "z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4"
)

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

# Fit each unique component once, then map components into pathway scenarios.
# This preserves v1-style summed component LOO while avoiding duplicate MCMC fits.
components <- list()

# Geography -> Bombus is common to all v2 pathway scenarios.
components$bombus_geo <- fit_beta(
  as.formula(paste("bombus_deficit ~", geo)),
  north_main, file.path(outdir, "component_bombus_geo")
)

# Reproductive-assurance alternatives.
components$sc_geo <- fit_bb(
  as.formula(paste("sc_successes | trials(sc_trials) ~", geo)),
  north_main, file.path(outdir, "component_sc_geo")
)
components$sc_bombus_geo <- fit_bb(
  as.formula(paste("sc_successes | trials(sc_trials) ~ z_bombus_deficit +", geo)),
  north_main, file.path(outdir, "component_sc_bombus_geo")
)

# Final floral responses for M0-M3.
for (response in c("plain", "form")) {
  if (response == "plain") {
    lhs <- "color_plain | trials(color_trials)"
  } else {
    lhs <- "form_open_generalized | trials(form_trials)"
  }

  components[[paste0(response, "_geo")]] <- fit_bb(
    as.formula(paste(lhs, "~", geo)),
    north_main, file.path(outdir, paste0("component_", response, "_geo"))
  )
  components[[paste0(response, "_bombus_geo")]] <- fit_bb(
    as.formula(paste(lhs, "~ z_bombus_deficit +", geo)),
    north_main, file.path(outdir, paste0("component_", response, "_bombus_geo"))
  )
  components[[paste0(response, "_sc_geo")]] <- fit_bb(
    as.formula(paste(lhs, "~ z_sc_share +", geo)),
    north_main, file.path(outdir, paste0("component_", response, "_sc_geo"))
  )
  components[[paste0(response, "_bombus_sc_geo")]] <- fit_bb(
    as.formula(paste(lhs, "~ z_bombus_deficit + z_sc_share +", geo)),
    north_main, file.path(outdir, paste0("component_", response, "_bombus_sc_geo"))
  )
}

# M0-M3 pathway definitions. Every scenario has the same four responses and
# exactly the same islands; only the hypothesised predictors differ.
scenarios <- list(
  M0_geographic_filter = list(
    bombus = "bombus_geo",
    self = "sc_geo",
    color = "plain_geo",
    form = "form_geo"
  ),
  M1_bombus_channel = list(
    bombus = "bombus_geo",
    self = "sc_geo",
    color = "plain_bombus_geo",
    form = "form_bombus_geo"
  ),
  M2_reproductive_assurance = list(
    bombus = "bombus_geo",
    self = "sc_geo",
    color = "plain_sc_geo",
    form = "form_sc_geo"
  ),
  M3_bombus_plus_reproductive_assurance = list(
    bombus = "bombus_geo",
    self = "sc_bombus_geo",
    color = "plain_bombus_sc_geo",
    form = "form_bombus_sc_geo"
  )
)

# Component-level posterior summaries.
component_summaries <- rbindlist(lapply(names(components), function(nm) {
  x <- as.data.table(posterior_summary(components[[nm]], pars = "^b_"), keep.rownames = "parameter")
  x[, component := nm]
  setcolorder(x, c("component", "parameter"))
  x
}), fill = TRUE)
fwrite(component_summaries, file.path(outdir, "posterior_fixed_effects.csv"))

# v1-style scenario comparison: sum component ELPD across the same response set.
# Pass newdata explicitly to work around a brms 2.20.x LOO restructuring issue
# where adforms$trials is not recovered correctly for certain formula structures.
component_loos <- lapply(components, function(fit) loo(fit, newdata = fit$data))
loo_component_rows <- list()
scenario_rows <- list()
for (scenario_name in names(scenarios)) {
  refs <- scenarios[[scenario_name]]
  rows <- rbindlist(lapply(names(refs), function(response_name) {
    component_name <- refs[[response_name]]
    x <- component_loos[[component_name]]$estimates
    data.table(
      scenario = scenario_name,
      response = response_name,
      component = component_name,
      elpd_loo = x["elpd_loo", "Estimate"],
      se_elpd_loo = x["elpd_loo", "SE"],
      p_loo = x["p_loo", "Estimate"],
      looic = x["looic", "Estimate"]
    )
  }))
  loo_component_rows[[scenario_name]] <- rows
  scenario_rows[[scenario_name]] <- data.table(
    scenario = scenario_name,
    total_elpd_loo = sum(rows$elpd_loo),
    approx_se_total_elpd_loo = sqrt(sum(rows$se_elpd_loo^2)),
    total_looic = sum(rows$looic)
  )
}
loo_components <- rbindlist(loo_component_rows)
scenario_comparison <- rbindlist(scenario_rows)
scenario_comparison[, delta_elpd_from_best := total_elpd_loo - max(total_elpd_loo)]
setorder(scenario_comparison, -total_elpd_loo)
fwrite(loo_components, file.path(outdir, "m0_m3_component_loo.csv"))
fwrite(scenario_comparison, file.path(outdir, "m0_m3_scenario_comparison_total_elpd.csv"))

# Posterior products for the M3 mediated path: Bombus deficit -> SC -> floral trait.
sc_draws <- as_draws_df(components$sc_bombus_geo)
plain_draws <- as_draws_df(components$plain_bombus_sc_geo)
form_draws <- as_draws_df(components$form_bombus_sc_geo)
n <- min(nrow(sc_draws), nrow(plain_draws), nrow(form_draws))
plain_indirect <- sc_draws$b_z_bombus_deficit[1:n] * plain_draws$b_z_sc_share[1:n]
form_indirect <- sc_draws$b_z_bombus_deficit[1:n] * form_draws$b_z_sc_share[1:n]
indirect <- data.table(
  path = c("Bombus_deficit_to_SC_to_plain", "Bombus_deficit_to_SC_to_generalized"),
  mean = c(mean(plain_indirect), mean(form_indirect)),
  q025 = c(quantile(plain_indirect, 0.025), quantile(form_indirect, 0.025)),
  q975 = c(quantile(plain_indirect, 0.975), quantile(form_indirect, 0.975))
)
fwrite(indirect, file.path(outdir, "m3_indirect_effects.csv"))

# M4: category-preserving pollination-regime extension. This asks whether detailed
# colour/form composition changes with distance x area differently among regions and
# with alternative pollinator-guild composition. It is reported separately because
# its multinomial responses differ from the planned-contrast M0-M3 responses.
m4 <- d[
  analysis_regime %in% c("northern_midlatitude", "tropical", "southern_extratropical") &
    complete.cases(
      z_log_distance_to_continent_km,
      z_log_island_area_km2,
      z_showy_alt_guild_share,
      z_other_bee_guild_share,
      z_generalist_insect_guild_share
    )
]
m4[, analysis_regime := factor(analysis_regime)]

# For multinomial responses, brms creates category-specific distributional
# parameters. A generic class='b' prior does not map to those parameters in all
# brms versions, so M4 deliberately uses brms defaults here. The category model
# itself is unchanged; explicit category-specific priors can be added later using
# get_prior() once the final category contract is frozen.
m4_color <- brm(
  cbind(color_plain, color_yellow_orange, color_red_pink, color_blue_purple) | trials(color_trials) ~
    z_log_distance_to_continent_km * z_log_island_area_km2 * analysis_regime +
    z_showy_alt_guild_share + z_other_bee_guild_share + z_generalist_insect_guild_share,
  data = m4[color_trials > 0],
  family = multinomial(),
  chains = 2, cores = 2, iter = 1200, warmup = 600, seed = 20260713,
  control = list(adapt_delta = 0.95), save_pars = save_pars(all = TRUE),
  file = file.path(outdir, "m4_color"), refresh = 100
)

m4_form <- brm(
  cbind(form_open_generalized, form_tubular_trumpet, form_zygomorphic_specialized, form_composite_brush) | trials(form_trials) ~
    z_log_distance_to_continent_km * z_log_island_area_km2 * analysis_regime +
    z_showy_alt_guild_share + z_other_bee_guild_share + z_generalist_insect_guild_share,
  data = m4[form_trials > 0],
  family = multinomial(),
  chains = 2, cores = 2, iter = 1200, warmup = 600, seed = 20260713,
  control = list(adapt_delta = 0.95), save_pars = save_pars(all = TRUE),
  file = file.path(outdir, "m4_form"), refresh = 100
)

m4_summaries <- rbindlist(list(
  cbind(model = "M4_color", as.data.table(posterior_summary(m4_color, pars = "^b_"), keep.rownames = "parameter")),
  cbind(model = "M4_form", as.data.table(posterior_summary(m4_form, pars = "^b_"), keep.rownames = "parameter"))
), fill = TRUE)
fwrite(m4_summaries, file.path(outdir, "m4_multinomial_fixed_effects.csv"))

all_fits <- c(components, list(M4_color = m4_color, M4_form = m4_form))
capture.output(lapply(all_fits, summary), file = file.path(outdir, "model_summaries.txt"))

meta <- list(
  contract = "v2_bayesian_m0_m4_main_v2",
  analysis_tier = "sensitivity_all_all_filled_data_requested_first_run",
  method = "Bayesian piecewise path models; beta-binomial M0-M3 scenario comparison plus category-preserving multinomial M4",
  main_interaction = "log distance to continent x log island area",
  m0_m3_common_support = TRUE,
  m0_m3_scenario_comparison = "summed component LOO-ELPD across Bombus, self-compatibility, colour, and form responses",
  n_northern_midlatitude_common_support = nrow(north_main),
  n_global_m4 = nrow(m4),
  sensitivity_analysis = "deferred by design",
  note = "Initial main fit uses 2 chains x 1200 iterations. Check R-hat, ESS, divergences, and Pareto-k diagnostics before manuscript claims."
)
write_json(meta, file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
