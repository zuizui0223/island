#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(readr)
  library(posterior)
  library(loo)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript analysis/v2/run_bayesian_m0_m4_primary.R <input.csv> <output_dir>")
}
input_csv <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20260713
options(mc.cores = 2)
set.seed(seed)

dat <- read_csv(input_csv, show_col_types = FALSE)

z <- function(x) as.numeric(scale(x))
clip01 <- function(x, eps = 1e-4) pmin(pmax(x, eps), 1 - eps)

north <- dat %>%
  filter(is_northern_midlatitude) %>%
  mutate(
    z_distance = z(log_distance),
    z_area = z(log_area),
    z_climate1 = z(climate_pc1),
    z_climate2 = z(climate_pc2),
    bombus_logit = qlogis(clip01(bombus_compatibility)),
    z_bombus = z(bombus_logit),
    sc_share = self_compatible_successes / self_compatible_trials,
    z_sc = z(sc_share),
    z_showy_alt = z(showy_alt_share)
  )

required <- c(
  "plain_color_successes", "plain_color_trials",
  "generalized_form_successes", "generalized_form_trials",
  "self_compatible_successes", "self_compatible_trials",
  "z_distance", "z_area", "z_climate1", "z_climate2",
  "z_bombus", "z_sc", "z_showy_alt"
)
north <- north %>% filter(if_all(all_of(required), ~ is.finite(.x)))

if (nrow(north) < 100) stop("Too few complete northern-midlatitude islands for primary analysis")

write_csv(north, file.path(out_dir, "bayesian_m0_m4_model_ready.csv"))

priors_bb <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(gamma(2, 0.1), class = "phi")
)
priors_gaussian <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(exponential(1), class = "sigma")
)

fit_bb <- function(name, formula, data = north) {
  message("Fitting ", name)
  fit <- brm(
    formula = formula,
    data = data,
    family = beta_binomial2(link = "logit"),
    prior = priors_bb,
    chains = 2,
    iter = 1200,
    warmup = 600,
    cores = 2,
    seed = seed,
    refresh = 100,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    save_pars = save_pars(all = TRUE)
  )
  saveRDS(fit, file.path(out_dir, paste0(name, ".rds")))
  fit
}

fit_gaussian <- function(name, formula, data = north) {
  message("Fitting ", name)
  fit <- brm(
    formula = formula,
    data = data,
    family = gaussian(),
    prior = priors_gaussian,
    chains = 2,
    iter = 1200,
    warmup = 600,
    cores = 2,
    seed = seed,
    refresh = 100,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    save_pars = save_pars(all = TRUE)
  )
  saveRDS(fit, file.path(out_dir, paste0(name, ".rds")))
  fit
}

geo <- "z_distance * z_area + z_climate1 + z_climate2"

# Shared exogenous-response model: same in all scenarios, ensuring comparable joint LOO support.
f_bombus <- as.formula(paste("z_bombus ~", geo))
bombus_geo <- fit_gaussian("bombus_geo", f_bombus)

# M0: geographic filtering only.
sc_m0 <- fit_bb("sc_m0", bf(self_compatible_successes | trials(self_compatible_trials) ~ z_distance * z_area + z_climate1 + z_climate2))
plain_m0 <- fit_bb("plain_m0", bf(plain_color_successes | trials(plain_color_trials) ~ z_distance * z_area + z_climate1 + z_climate2))
form_m0 <- fit_bb("form_m0", bf(generalized_form_successes | trials(generalized_form_trials) ~ z_distance * z_area + z_climate1 + z_climate2))

# M1: pollinator-channel pathway to floral composition, without reproductive mediation.
plain_m1 <- fit_bb("plain_m1", bf(plain_color_successes | trials(plain_color_trials) ~ z_bombus + z_distance * z_area + z_climate1 + z_climate2))
form_m1 <- fit_bb("form_m1", bf(generalized_form_successes | trials(generalized_form_trials) ~ z_bombus + z_distance * z_area + z_climate1 + z_climate2))

# M2: reproductive-assurance pathway without Bombus mediation.
plain_m2 <- fit_bb("plain_m2", bf(plain_color_successes | trials(plain_color_trials) ~ z_sc + z_distance * z_area + z_climate1 + z_climate2))
form_m2 <- fit_bb("form_m2", bf(generalized_form_successes | trials(generalized_form_trials) ~ z_sc + z_distance * z_area + z_climate1 + z_climate2))

# M3: full Bombus -> reproductive assurance + direct/indirect floral pathway.
sc_m3 <- fit_bb("sc_m3", bf(self_compatible_successes | trials(self_compatible_trials) ~ z_bombus + z_distance * z_area + z_climate1 + z_climate2))
plain_m3 <- fit_bb("plain_m3", bf(plain_color_successes | trials(plain_color_trials) ~ z_bombus + z_sc + z_distance * z_area + z_climate1 + z_climate2))
form_m3 <- fit_bb("form_m3", bf(generalized_form_successes | trials(generalized_form_trials) ~ z_bombus + z_sc + z_distance * z_area + z_climate1 + z_climate2))

# M4: M3 plus functional replacement by highly mobile/showy alternative guilds.
plain_m4 <- fit_bb("plain_m4", bf(plain_color_successes | trials(plain_color_trials) ~ z_bombus * z_showy_alt + z_sc + z_distance * z_area + z_climate1 + z_climate2))
form_m4 <- fit_bb("form_m4", bf(generalized_form_successes | trials(generalized_form_trials) ~ z_bombus * z_showy_alt + z_sc + z_distance * z_area + z_climate1 + z_climate2))

models <- list(
  bombus_geo = bombus_geo,
  sc_m0 = sc_m0,
  plain_m0 = plain_m0,
  form_m0 = form_m0,
  plain_m1 = plain_m1,
  form_m1 = form_m1,
  plain_m2 = plain_m2,
  form_m2 = form_m2,
  sc_m3 = sc_m3,
  plain_m3 = plain_m3,
  form_m3 = form_m3,
  plain_m4 = plain_m4,
  form_m4 = form_m4
)

# Fixed-effect summaries.
coef_rows <- lapply(names(models), function(nm) {
  fx <- as.data.frame(fixef(models[[nm]], probs = c(0.025, 0.975)))
  fx$predictor <- rownames(fx)
  rownames(fx) <- NULL
  fx$model <- nm
  fx
})
coef_table <- bind_rows(coef_rows) %>%
  rename(estimate = Estimate, estimate_error = Est.Error, q2.5 = Q2.5, q97.5 = Q97.5) %>%
  select(model, predictor, estimate, estimate_error, q2.5, q97.5)
write_csv(coef_table, file.path(out_dir, "bayesian_m0_m4_fixed_effects.csv"))

# Same observed response set for M0-M4 scenario comparison.
loo_cache <- lapply(models, loo)
scenario_components <- list(
  M0 = c("bombus_geo", "sc_m0", "plain_m0", "form_m0"),
  M1 = c("bombus_geo", "sc_m0", "plain_m1", "form_m1"),
  M2 = c("bombus_geo", "sc_m0", "plain_m2", "form_m2"),
  M3 = c("bombus_geo", "sc_m3", "plain_m3", "form_m3"),
  M4 = c("bombus_geo", "sc_m3", "plain_m4", "form_m4")
)
scenario_table <- bind_rows(lapply(names(scenario_components), function(sc) {
  comps <- scenario_components[[sc]]
  elpd <- sum(vapply(comps, function(x) loo_cache[[x]]$estimates["elpd_loo", "Estimate"], numeric(1)))
  se2 <- sum(vapply(comps, function(x) loo_cache[[x]]$estimates["elpd_loo", "SE"]^2, numeric(1)))
  tibble(scenario = sc, total_elpd_loo = elpd, approximate_se = sqrt(se2), components = paste(comps, collapse = ";"))
})) %>% arrange(desc(total_elpd_loo))
scenario_table$delta_elpd_from_best <- scenario_table$total_elpd_loo - max(scenario_table$total_elpd_loo)
write_csv(scenario_table, file.path(out_dir, "bayesian_m0_m4_scenario_comparison.csv"))

# Posterior path products for M3. These are descriptive indirect effects on standardized link scales.
draw_sc <- as_draws_df(sc_m3)
draw_plain <- as_draws_df(plain_m3)
draw_form <- as_draws_df(form_m3)
indirect <- tibble(
  draw = seq_len(min(nrow(draw_sc), nrow(draw_plain), nrow(draw_form))),
  bombus_to_sc = draw_sc$b_z_bombus[seq_len(min(nrow(draw_sc), nrow(draw_plain), nrow(draw_form)))],
  sc_to_plain = draw_plain$b_z_sc[seq_len(min(nrow(draw_sc), nrow(draw_plain), nrow(draw_form)))],
  sc_to_form = draw_form$b_z_sc[seq_len(min(nrow(draw_sc), nrow(draw_plain), nrow(draw_form)))]
) %>% mutate(
  indirect_plain = bombus_to_sc * sc_to_plain,
  indirect_form = bombus_to_sc * sc_to_form
)

summarise_draw <- function(x, path) {
  tibble(
    path = path,
    estimate = mean(x),
    q2.5 = quantile(x, 0.025),
    q50 = quantile(x, 0.5),
    q97.5 = quantile(x, 0.975),
    posterior_positive = mean(x > 0)
  )
}
path_table <- bind_rows(
  summarise_draw(indirect$indirect_plain, "Bombus -> SC -> plain_color"),
  summarise_draw(indirect$indirect_form, "Bombus -> SC -> generalized_form")
)
write_csv(path_table, file.path(out_dir, "bayesian_m3_indirect_effects.csv"))

metadata <- tibble(
  key = c(
    "analysis", "scope", "n_islands", "seed", "chains", "iterations",
    "main_interaction", "floral_scope", "category_policy", "sensitivity_policy"
  ),
  value = c(
    "v2 Bayesian M0-M4 primary analysis",
    "northern_midlatitude, all available sensitivity_all trait fills for first main run",
    as.character(nrow(north)),
    as.character(seed), "2", "1200",
    "standardized log(distance to continent) x standardized log(island area)",
    "animal-pollinated species only (biotic or mixed pollen-vector mode)",
    "original categories preserved separately; planned binary contrasts used for primary path tests",
    "not run in this workflow"
  )
)
write_csv(metadata, file.path(out_dir, "bayesian_m0_m4_metadata.csv"))

message("Bayesian M0-M4 primary analysis completed")
print(scenario_table)
print(path_table)
