#!/usr/bin/env Rscript

# Jointly analyse the support-promoted fine colour x floral-form composition.
# Seven colour groups x seven form groups are available after promotion, but only
# globally observed cells are retained in the fitted composition. Counts still
# exhaust fine_joint_color_form_trials for every island.

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: run_inla_fine_joint_color_form_composition.R <input.csv> <output_dir>")
infile <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260714)

d <- fread(infile)
colors <- c("white", "blue", "yellow", "red", "green", "cream", "rare_other")
forms <- c(
  "tubular", "composite_head", "star", "bell_campanulate",
  "open_radial", "trumpet_funnel", "rare_other"
)
all_categories <- as.vector(outer(colors, forms, paste, sep = "__"))
all_cols <- paste0("fine_joint_", all_categories)
required <- c(
  "island_id", "analysis_regime", "spatial_block", "fine_joint_color_form_trials",
  all_cols,
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
  "bombus_deficit", "sc_share", "wind_mixed_share"
)
missing <- setdiff(required, names(d))
if (length(missing)) stop("missing fine joint-composition columns: ", paste(missing, collapse = ", "))

cell_sum <- rowSums(d[, ..all_cols], na.rm = TRUE)
expected_trials <- fifelse(is.na(d$fine_joint_color_form_trials), 0, d$fine_joint_color_form_trials)
if (any(cell_sum != expected_trials)) stop("fine joint colour-form cells do not sum to fine_joint_color_form_trials")

global_counts <- colSums(d[, ..all_cols], na.rm = TRUE)
observed_cols <- names(global_counts[global_counts > 0])
if (length(observed_cols) < 2) stop("too few observed fine joint categories")
observed_categories <- sub("^fine_joint_", "", observed_cols)

scale_vars <- c(
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
  "bombus_deficit", "sc_share", "wind_mixed_share"
)
north_ref <- d[analysis_regime == "northern_midlatitude"]
scaling <- rbindlist(lapply(scale_vars, function(nm) {
  x <- north_ref[[nm]]
  data.table(variable = nm, center = mean(x, na.rm = TRUE), scale = sd(x, na.rm = TRUE))
}))
if (any(!is.finite(scaling$scale) | scaling$scale <= 0)) stop("invalid fine-joint scaling")
for (i in seq_len(nrow(scaling))) {
  nm <- scaling$variable[i]
  d[[paste0("z_", nm)]] <- (d[[nm]] - scaling$center[i]) / scaling$scale[i]
}
fwrite(scaling, file.path(outdir, "fine_joint_color_form_scaling.csv"))

predictors <- c(
  "z_log_distance_to_continent_km", "z_log_island_area_km2",
  "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
  "z_bombus_deficit", "z_sc_share", "z_wind_mixed_share"
)
support <- d[
  analysis_regime == "northern_midlatitude" &
    fine_joint_color_form_trials > 0 &
    complete.cases(d[, ..predictors])
]
if (nrow(support) < 30) stop("too few northern islands for fine joint colour-form model")

long <- melt(
  support,
  id.vars = c("island_id", "spatial_block", "fine_joint_color_form_trials", predictors),
  measure.vars = observed_cols,
  variable.name = "joint_category",
  value.name = "count"
)
long[, joint_category := factor(sub("^fine_joint_", "", joint_category), levels = observed_categories)]
long[, island_comp_id := as.integer(factor(island_id))]
long[, spatial_id := as.integer(factor(spatial_block))]
long[, obs_id := seq_len(.N)]
long[, log_trials := log(fine_joint_color_form_trials)]

# Structural zero cells are omitted globally, but every observed island composition
# still sums to its complete fine joint trial count.
check <- long[, sum(count), by = island_id]
truth <- support[, .(island_id, fine_joint_color_form_trials)]
check <- merge(check, truth, by = "island_id")
if (any(check$V1 != check$fine_joint_color_form_trials)) stop("observed fine joint cells do not exhaust island trials")

geo_terms <- paste(
  "joint_category:z_log_distance_to_continent_km +",
  "joint_category:z_log_island_area_km2 +",
  "joint_category:z_climate_pc1 + joint_category:z_climate_pc2 +",
  "joint_category:z_climate_pc3 + joint_category:z_climate_pc4"
)
base_terms <- paste(
  "0 + joint_category +", geo_terms,
  "+ f(island_comp_id, model='iid') + f(spatial_id, model='iid') + f(obs_id, model='iid')"
)
models <- list(
  FJ0_geo = base_terms,
  FJ1_bombus = paste(base_terms, "+ joint_category:z_bombus_deficit"),
  FJ2_sc_wind = paste(base_terms, "+ joint_category:z_sc_share + joint_category:z_wind_mixed_share"),
  FJ3_full = paste(
    base_terms,
    "+ joint_category:z_bombus_deficit + joint_category:z_sc_share +",
    "joint_category:z_wind_mixed_share"
  )
)
compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
fits <- list(); score_rows <- list(); fixed_rows <- list()
for (model_name in names(models)) {
  fit <- inla(
    as.formula(paste("count ~", models[[model_name]], "+ offset(log_trials)")),
    family = "poisson", data = long,
    control.compute = compute,
    control.predictor = list(compute = TRUE), verbose = FALSE
  )
  fits[[model_name]] <- fit
  cpo <- fit$cpo$cpo
  cpo <- cpo[is.finite(cpo) & cpo > 0]
  score_rows[[model_name]] <- data.table(
    model = model_name, n_islands = nrow(support),
    n_observed_joint_categories = length(observed_categories),
    n_joint_rows = nrow(long),
    log_cpo_sum = if (length(cpo)) sum(log(cpo)) else NA_real_,
    waic = fit$waic$waic, dic = fit$dic$dic
  )
  fx <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
  fx[, model := model_name]
  fixed_rows[[model_name]] <- fx
}

scores <- rbindlist(score_rows)
scores[, delta_log_cpo_from_best := log_cpo_sum - max(log_cpo_sum, na.rm = TRUE)]
setorder(scores, -log_cpo_sum)
fwrite(scores, file.path(outdir, "fine_joint_color_form_model_comparisons.csv"))

fixed <- rbindlist(fixed_rows, fill = TRUE)
fixed[, excludes_zero := (`0.025quant` > 0 & `0.975quant` > 0) |
  (`0.025quant` < 0 & `0.975quant` < 0)]
fwrite(fixed, file.path(outdir, "fine_joint_color_form_fixed_effects.csv"))

key <- fixed[
  model == "FJ3_full" & grepl("joint_category.*:(z_bombus_deficit|z_sc_share|z_wind_mixed_share)$", parameter)
]
key[, predictor := fifelse(
  grepl("z_bombus_deficit$", parameter), "bombus_deficit",
  fifelse(grepl("z_sc_share$", parameter), "self_compatibility", "wind_mixed")
)]
key[, joint_category := sub("^joint_category([^:]+):.*$", "\\1", parameter)]
fwrite(key, file.path(outdir, "fine_joint_color_form_key_effects.csv"))

category_support <- data.table(
  joint_category = observed_categories,
  total_count = as.numeric(global_counts[observed_cols])
)
fwrite(category_support, file.path(outdir, "fine_joint_color_form_observed_cells.csv"))

meta <- list(
  model = "support_promoted_fine_joint_colour_form_poisson_loglinear",
  hypothesis_role = "primary floral phenotype composition model",
  n_islands = nrow(support),
  n_observed_joint_categories = length(observed_categories),
  categories = observed_categories,
  promotion = list(
    colours = colors,
    forms = forms,
    residual_resolved_states = "rare_other"
  ),
  interpretation = paste(
    "Fine-grained colour-form states are fitted jointly as one composition;",
    "coefficients describe compositional reallocation rather than independent binary responses."
  )
)
write_json(meta, file.path(outdir, "fine_joint_color_form_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(lapply(fits, summary), file = file.path(outdir, "fine_joint_color_form_model_summaries.txt"))
