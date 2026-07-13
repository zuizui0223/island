#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(data.table)
  library(posterior)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: run_bayesian_binary_pathway.R <input.csv> <output_dir>")
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

# Raw geographic distance on a standardized scale with a quadratic term.
# Positive quadratic = accelerating isolation effect at long distance;
# negative quadratic = saturation/turnover.
d[, z_distance := z(distance_to_continent_km)]
d[, z_distance_sq := z_distance^2]
d[, z_log_area := z(log(area_km2))]
for (nm in c("climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4")) {
  d[[paste0("z_", nm)]] <- z(d[[nm]])
}
d[, z_bombus_deficit := z(bombus_deficit)]
d[, z_sc_share := z(sc_share)]

# Functional flower-form regrouping for the binary Bayes layer.
# INLA retains and decomposes all four original form categories.
d[, form_generalist_like := form_open_generalized + form_composite_brush]
if (any(d$form_generalist_like > d$form_trials, na.rm = TRUE)) {
  stop("generalist-like form successes exceed form trials")
}

d[, analysis_regime := factor(
  analysis_regime,
  levels = c(
    "northern_midlatitude",
    "northern_high_latitude",
    "tropical",
    "southern_extratropical"
  )
)]

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

base_complete <- complete.cases(d[, .(
  z_distance, z_distance_sq, z_log_area, z_climate_pc1, z_climate_pc2,
  z_climate_pc3, z_climate_pc4, analysis_regime
)])
sc_global <- d[base_complete & sc_trials > 0]
color_global <- d[base_complete & color_trials > 0]
form_global <- d[base_complete & form_trials > 0]

f_sc <- fit_bb(
  sc_successes | trials(sc_trials) ~
    (z_distance + z_distance_sq) * analysis_regime + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  sc_global,
  file.path(outdir, "falsification_sc")
)
f_plain <- fit_bb(
  color_plain | trials(color_trials) ~
    (z_distance + z_distance_sq) * analysis_regime + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  color_global,
  file.path(outdir, "falsification_plain")
)
f_generalist_like <- fit_bb(
  form_generalist_like | trials(form_trials) ~
    (z_distance + z_distance_sq) * analysis_regime + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  form_global,
  file.path(outdir, "falsification_generalist_like")
)
f_red_pink <- fit_bb(
  color_red_pink | trials(color_trials) ~
    (z_distance + z_distance_sq) * analysis_regime + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  color_global,
  file.path(outdir, "falsification_red_pink")
)
f_yellow_orange <- fit_bb(
  color_yellow_orange | trials(color_trials) ~
    (z_distance + z_distance_sq) * analysis_regime + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  color_global,
  file.path(outdir, "falsification_yellow_orange")
)
f_tubular <- fit_bb(
  form_tubular_trumpet | trials(form_trials) ~
    (z_distance + z_distance_sq) * analysis_regime + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  form_global,
  file.path(outdir, "falsification_tubular_trumpet")
)

north <- d[analysis_regime == "northern_midlatitude"]
north_base <- complete.cases(north[, .(
  z_distance, z_distance_sq, z_log_area, z_climate_pc1, z_climate_pc2,
  z_climate_pc3, z_climate_pc4, z_bombus_deficit
)])
sc_north <- north[north_base & sc_trials > 0]
plain_north <- north[north_base & color_trials > 0 & is.finite(z_sc_share)]
generalist_like_north <- north[north_base & form_trials > 0 & is.finite(z_sc_share)]

n_sc <- fit_bb(
  sc_successes | trials(sc_trials) ~
    z_bombus_deficit + z_distance + z_distance_sq + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  sc_north,
  file.path(outdir, "north_bombus_sc")
)
n_plain <- fit_bb(
  color_plain | trials(color_trials) ~
    z_bombus_deficit + z_sc_share + z_distance + z_distance_sq + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  plain_north,
  file.path(outdir, "north_bombus_sc_plain")
)
n_generalist_like <- fit_bb(
  form_generalist_like | trials(form_trials) ~
    z_bombus_deficit + z_sc_share + z_distance + z_distance_sq + z_log_area +
    z_climate_pc1 + z_climate_pc2 + z_climate_pc3 + z_climate_pc4,
  generalist_like_north,
  file.path(outdir, "north_bombus_sc_generalist_like")
)

fits <- list(
  F_SC = f_sc,
  F_plain = f_plain,
  F_generalist_like = f_generalist_like,
  F_red_pink = f_red_pink,
  F_yellow_orange = f_yellow_orange,
  F_tubular = f_tubular,
  N_SC = n_sc,
  N_plain = n_plain,
  N_generalist_like = n_generalist_like
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
    "global_SC", "global_plain", "global_generalist_like",
    "global_red_pink", "global_yellow_orange", "global_tubular_trumpet",
    "northern_SC", "northern_plain", "northern_generalist_like"
  ),
  n_islands = c(
    nrow(sc_global), nrow(color_global), nrow(form_global),
    nrow(color_global), nrow(color_global), nrow(form_global),
    nrow(sc_north), nrow(plain_north), nrow(generalist_like_north)
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

regime_component <- function(draws, component, regime) {
  base_name <- paste0("b_", component)
  out <- coef_draw(draws, base_name)
  if (regime != "northern_midlatitude") {
    interaction_name <- paste0("b_", component, ":analysis_regime", regime)
    out <- out + coef_draw(draws, interaction_name)
  }
  out
}

regime_distance_table <- function(fit, outcome) {
  draws <- as_draws_df(fit)
  regimes <- c(
    "northern_midlatitude",
    "northern_high_latitude",
    "tropical",
    "southern_extratropical"
  )
  rbindlist(lapply(regimes, function(regime) {
    linear <- regime_component(draws, "z_distance", regime)
    quadratic <- regime_component(draws, "z_distance_sq", regime)
    rbindlist(list(
      cbind(
        outcome = outcome,
        regime = regime,
        component = "linear_distance",
        summarize_draw("distance_linear", linear)
      ),
      cbind(
        outcome = outcome,
        regime = regime,
        component = "quadratic_acceleration",
        summarize_draw("distance_quadratic", quadratic)
      )
    ), fill = TRUE)
  }), fill = TRUE)
}

regime_distance_effects <- rbindlist(list(
  regime_distance_table(f_sc, "self_compatibility"),
  regime_distance_table(f_plain, "plain_flower"),
  regime_distance_table(f_generalist_like, "generalist_like_form"),
  regime_distance_table(f_red_pink, "red_pink_flower"),
  regime_distance_table(f_yellow_orange, "yellow_orange_flower"),
  regime_distance_table(f_tubular, "tubular_trumpet_form")
), fill = TRUE)
fwrite(regime_distance_effects, file.path(outdir, "cross_regime_isolation_slopes.csv"))

sc_draws <- as_draws_df(n_sc)
plain_draws <- as_draws_df(n_plain)
generalist_draws <- as_draws_df(n_generalist_like)
n_draws <- min(nrow(sc_draws), nrow(plain_draws), nrow(generalist_draws))
sc_draws <- sc_draws[seq_len(n_draws), ]
plain_draws <- plain_draws[seq_len(n_draws), ]
generalist_draws <- generalist_draws[seq_len(n_draws), ]

b_bombus_sc <- coef_draw(sc_draws, "b_z_bombus_deficit")
b_sc_plain <- coef_draw(plain_draws, "b_z_sc_share")
b_sc_generalist <- coef_draw(generalist_draws, "b_z_sc_share")
b_bombus_plain <- coef_draw(plain_draws, "b_z_bombus_deficit")
b_bombus_generalist <- coef_draw(generalist_draws, "b_z_bombus_deficit")

paths <- rbindlist(list(
  summarize_draw("Bombus_deficit_to_SC_to_plain", b_bombus_sc * b_sc_plain),
  summarize_draw("Bombus_deficit_to_SC_to_generalist_like", b_bombus_sc * b_sc_generalist),
  summarize_draw("Bombus_deficit_to_plain_direct", b_bombus_plain),
  summarize_draw("Bombus_deficit_to_generalist_like_direct", b_bombus_generalist),
  summarize_draw("Bombus_deficit_to_plain_total", b_bombus_plain + b_bombus_sc * b_sc_plain),
  summarize_draw("Bombus_deficit_to_generalist_like_total", b_bombus_generalist + b_bombus_sc * b_sc_generalist)
))
fwrite(paths, file.path(outdir, "northern_direct_indirect_total_effects.csv"))

sampler_diagnostics <- rbindlist(lapply(names(fits), function(nm) {
  fit <- fits[[nm]]
  np <- nuts_params(fit)
  diag <- as.data.table(posterior::summarise_draws(
    posterior::as_draws_array(fit),
    posterior::rhat,
    posterior::ess_bulk,
    posterior::ess_tail
  ))
  diag <- diag[grepl("^b_", variable)]
  rhat_values <- diag$rhat
  bulk_values <- diag$ess_bulk
  tail_values <- diag$ess_tail
  data.table(
    model = nm,
    divergent_transitions = if (nrow(np)) sum(np$Parameter == "divergent__" & np$Value == 1) else 0L,
    max_treedepth_hits = if (nrow(np)) sum(np$Parameter == "treedepth__" & np$Value >= 15) else 0L,
    max_rhat = if (!length(rhat_values) || all(is.na(rhat_values))) NA_real_ else max(rhat_values, na.rm = TRUE),
    min_bulk_ess = if (!length(bulk_values) || all(is.na(bulk_values))) NA_real_ else min(bulk_values, na.rm = TRUE),
    min_tail_ess = if (!length(tail_values) || all(is.na(tail_values))) NA_real_ else min(tail_values, na.rm = TRUE)
  )
}), fill = TRUE)
fwrite(sampler_diagnostics, file.path(outdir, "sampler_diagnostics.csv"))

capture.output(lapply(fits, summary), file = file.path(outdir, "model_summaries.txt"))
write_json(list(
  contract = "v2_engine_specific_bayesian_binary_pathway_v3",
  role = "binary syndrome inference plus targeted channel-associated trait contrasts; INLA owns full category decomposition",
  analysis_tier = "sensitivity_all",
  binary_outcomes = c(
    "SC_vs_SI",
    "plain_vs_non_plain",
    "generalist_like_open_plus_composite_vs_tubular_plus_zygomorphic"
  ),
  targeted_global_traits = c("red_pink", "yellow_orange", "tubular_trumpet"),
  northern_pathway = "Bombus deficit predicts parallel SC, plain-colour and generalist-like form responses; SC-mediated quantities are descriptive rather than assumed causal",
  falsification = "four-regime raw-distance linear and quadratic effects without global Bombus causal imposition",
  distance_function = "standardized raw distance plus quadratic term; positive quadratic means accelerating isolation effect",
  global_regimes = c(
    "northern_midlatitude",
    "northern_high_latitude",
    "tropical",
    "southern_extratropical"
  ),
  category_models = FALSE,
  alternative_pollinator_primary_covariates = FALSE,
  channel_trait_interpretation = "red/pink, yellow/orange and tubular/trumpet are trait syndromes only; they do not directly measure bird, butterfly or moth abundance",
  n_input_islands = n_input,
  n_zero_or_negative_distance_excluded = n_zero_distance,
  support_strategy = "maximum available support per equation"
), file.path(outdir, "analysis_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
