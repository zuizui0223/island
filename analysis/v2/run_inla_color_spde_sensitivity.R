#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
  library(sf)
  library(jsonlite)
  library(Matrix)
  library(fmesher)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: run_inla_color_spde_sensitivity.R <input.csv> <output_dir>")
infile <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260714)

d <- fread(infile)
colors <- c("white", "blue", "yellow", "red", "green", "cream", "rare_other")
regime_levels <- c(
  "northern_midlatitude", "northern_high_latitude",
  "tropical", "southern_extratropical"
)
d[, analysis_regime := factor(analysis_regime, levels = regime_levels)]

required <- c(
  "island_id", "island_longitude", "island_latitude", "analysis_regime",
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
  "bombus_deficit", "sc_share", "animal_status_observed", "n_wind_species", "n_mixed_species",
  "all_fine_color_trials", "fine_color_trials",
  paste0("all_fine_color_", colors), paste0("fine_color_", colors)
)
missing <- setdiff(required, names(d))
if (length(missing)) stop("missing SPDE sensitivity columns: ", paste(missing, collapse = ", "))

d[, wind_share := fifelse(animal_status_observed > 0, n_wind_species / animal_status_observed, NA_real_)]
d[, mixed_share := fifelse(animal_status_observed > 0, n_mixed_species / animal_status_observed, NA_real_)]

for (nm in c(
  "log_distance_to_continent_km", "log_island_area_km2",
  "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"
)) {
  mu <- mean(d[[nm]], na.rm = TRUE)
  sig <- sd(d[[nm]], na.rm = TRUE)
  if (!is.finite(sig) || sig <= 0) stop("invalid SPDE scale for ", nm)
  d[[paste0("z_", nm)]] <- (d[[nm]] - mu) / sig
}
set(d, j = "z_isolation_sq", value = d$z_log_distance_to_continent_km^2)

north <- d[analysis_regime == "northern_midlatitude"]
mech_candidates <- c("bombus_deficit", "sc_share", "wind_share", "mixed_share")
mech_audit <- rbindlist(lapply(mech_candidates, function(nm) {
  x <- north[[nm]]
  n <- sum(is.finite(x))
  u <- uniqueN(x[is.finite(x)])
  s <- sd(x, na.rm = TRUE)
  ok <- n >= 2 && u >= 2 && is.finite(s) && s > 0
  data.table(variable = nm, n_nonmissing = n, n_unique = u, sd = s, estimable = ok)
}))
for (nm in mech_audit[estimable == TRUE, variable]) {
  mu <- mean(north[[nm]], na.rm = TRUE)
  s <- sd(north[[nm]], na.rm = TRUE)
  d[[paste0("z_", nm)]] <- (d[[nm]] - mu) / s
}
if (!"z_bombus_deficit" %in% names(d)) stop("bombus_deficit is not estimable for SPDE sensitivity")
if (!"z_sc_share" %in% names(d)) stop("sc_share is not estimable for SPDE sensitivity")
poll_controls <- mech_audit[variable %in% c("wind_share", "mixed_share") & estimable == TRUE, variable]

coord_dt <- unique(d[
  is.finite(island_longitude) & is.finite(island_latitude),
  .(island_id, island_longitude, island_latitude)
], by = "island_id")
pts <- st_as_sf(coord_dt, coords = c("island_longitude", "island_latitude"), crs = 4326, remove = FALSE)
pts <- st_transform(pts, 6933)
xy <- st_coordinates(pts)
coord_lookup <- data.table(island_id = coord_dt$island_id, x_spde = xy[, 1], y_spde = xy[, 2])
d <- merge(d, coord_lookup, by = "island_id", all.x = TRUE, sort = FALSE)

mesh <- fmesher::fm_mesh_2d_inla(
  loc = as.matrix(d[is.finite(x_spde) & is.finite(y_spde), .(x_spde, y_spde)]),
  cutoff = 100000,
  max.edge = c(500000, 1500000),
  offset = c(300000, 1500000)
)
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(500000, 0.5),
  prior.sigma = c(1, 0.01)
)

compute <- list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
results_fixed <- list()
results_scores <- list()
mesh_meta <- list(
  n_vertices = mesh$n,
  crs = "EPSG:6933",
  max_edge_m = c(500000, 1500000),
  cutoff_m = 100000,
  prior_range = c(500000, 0.5),
  prior_sigma = c(1, 0.01)
)

make_category_block_A <- function(loc, category_id) {
  A0 <- fmesher::fm_basis(mesh, loc = loc)
  A0 <- as(A0, "dgTMatrix")
  n_spde <- mesh$n
  sparseMatrix(
    i = A0@i + 1L,
    j = A0@j + 1L + (category_id[A0@i + 1L] - 1L) * n_spde,
    x = A0@x,
    dims = c(nrow(loc), n_spde * length(colors))
  )
}

fit_spde_composition <- function(data, endpoint, model_name, prefix, trials, rhs) {
  cols <- paste0(prefix, colors)
  ids <- intersect(c(
    "island_id", "analysis_regime", trials,
    "x_spde", "y_spde",
    "z_log_distance_to_continent_km", "z_isolation_sq", "z_log_island_area_km2",
    "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
    "z_bombus_deficit", "z_sc_share", "z_wind_share", "z_mixed_share"
  ), names(data))
  long <- melt(
    data,
    id.vars = ids,
    measure.vars = cols,
    variable.name = "category",
    value.name = "count"
  )
  long[, category := factor(sub(paste0("^", prefix), "", category), levels = colors)]
  long[, category_id := as.integer(category)]
  long[, obs_id := seq_len(.N)]
  long[, log_trials := log(get(trials))]

  loc <- as.matrix(long[, .(x_spde, y_spde)])
  A <- make_category_block_A(loc, long$category_id)
  spde_index <- inla.spde.make.index(
    "spatial_field",
    n.spde = spde$n.spde * length(colors)
  )

  fixed_data <- data.frame(
    Intercept = 1,
    category = long$category,
    analysis_regime = long$analysis_regime,
    z_log_distance_to_continent_km = long$z_log_distance_to_continent_km,
    z_isolation_sq = long$z_isolation_sq,
    z_log_island_area_km2 = long$z_log_island_area_km2,
    z_climate_pc1 = long$z_climate_pc1,
    z_climate_pc2 = long$z_climate_pc2,
    z_climate_pc3 = long$z_climate_pc3,
    z_climate_pc4 = long$z_climate_pc4,
    z_bombus_deficit = if ("z_bombus_deficit" %in% names(long)) long$z_bombus_deficit else NA_real_,
    z_sc_share = if ("z_sc_share" %in% names(long)) long$z_sc_share else NA_real_,
    z_wind_share = if ("z_wind_share" %in% names(long)) long$z_wind_share else NA_real_,
    z_mixed_share = if ("z_mixed_share" %in% names(long)) long$z_mixed_share else NA_real_,
    obs_id = long$obs_id,
    log_trials = long$log_trials
  )

  stk <- inla.stack(
    data = list(y = long$count),
    A = list(1, A),
    effects = list(fixed_data, spatial_field = spde_index$spatial_field),
    tag = model_name
  )
  stk_data <- inla.stack.data(stk)
  formula <- as.formula(paste(
    "y ~", rhs,
    "+ f(spatial_field, model=generic0, Cmatrix=kronecker(Diagonal(length(colors)), spde$param.inla$M0),",
    "hyper=list(prec=list(prior='loggamma', param=c(1, 0.00005))))",
    "+ f(obs_id, model='iid') + offset(log_trials)"
  ))
  fit <- inla(
    formula,
    family = "poisson",
    data = stk_data,
    control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
    control.compute = compute,
    verbose = FALSE
  )
  fx <- as.data.table(fit$summary.fixed, keep.rownames = "parameter")
  fx[, `:=`(endpoint = endpoint, model = model_name)]
  results_fixed[[length(results_fixed) + 1L]] <<- fx
  cpo <- fit$cpo$cpo
  cpo <- cpo[is.finite(cpo) & cpo > 0]
  results_scores[[length(results_scores) + 1L]] <<- data.table(
    endpoint = endpoint,
    model = model_name,
    n_islands = uniqueN(data$island_id),
    n_rows = nrow(long),
    mesh_vertices = mesh$n,
    log_cpo_sum = if (length(cpo)) sum(log(cpo)) else NA_real_,
    waic = fit$waic$waic,
    dic = fit$dic$dic
  )
}

geo_terms <- paste(
  "category:z_log_island_area_km2 +",
  "category:z_climate_pc1 + category:z_climate_pc2 +",
  "category:z_climate_pc3 + category:z_climate_pc4"
)

global_support <- d[
  !is.na(analysis_regime) & all_fine_color_trials > 0 &
    is.finite(x_spde) & is.finite(y_spde) &
    complete.cases(d[, .(
      z_log_distance_to_continent_km, z_isolation_sq, z_log_island_area_km2,
      z_climate_pc1, z_climate_pc2, z_climate_pc3, z_climate_pc4
    )])
]
global_rhs <- paste(
  "0 + category + category:analysis_regime +",
  "category:analysis_regime:z_log_distance_to_continent_km +",
  "category:analysis_regime:z_isolation_sq +",
  geo_terms
)
fit_spde_composition(
  global_support,
  "all_flora_color",
  "G_global_colour_SPDE",
  "all_fine_color_",
  "all_fine_color_trials",
  global_rhs
)

mech_vars <- c(
  "z_log_distance_to_continent_km", "z_isolation_sq", "z_log_island_area_km2",
  "z_climate_pc1", "z_climate_pc2", "z_climate_pc3", "z_climate_pc4",
  "z_bombus_deficit", "z_sc_share", paste0("z_", poll_controls)
)
north_support <- d[
  analysis_regime == "northern_midlatitude" & fine_color_trials > 0 &
    is.finite(x_spde) & is.finite(y_spde) & complete.cases(d[, ..mech_vars])
]
mechanism_rhs <- paste(
  "0 + category + category:z_log_distance_to_continent_km + category:z_isolation_sq +",
  geo_terms,
  "+ category:z_bombus_deficit + category:z_sc_share",
  if (length(poll_controls)) paste0(" + ", paste0("category:z_", poll_controls, collapse = " + ")) else ""
)
fit_spde_composition(
  north_support,
  "animal_flora_color",
  "N3_colour_SPDE",
  "fine_color_",
  "fine_color_trials",
  mechanism_rhs
)

fixed <- rbindlist(results_fixed, fill = TRUE)
fixed[, excludes_zero := (`0.025quant` > 0 & `0.975quant` > 0) |
  (`0.025quant` < 0 & `0.975quant` < 0)]
scores <- rbindlist(results_scores, fill = TRUE)
key <- fixed[grepl(
  "z_log_distance_to_continent_km|z_isolation_sq|z_bombus_deficit|z_sc_share|z_wind_share|z_mixed_share",
  parameter
)]

fwrite(fixed, file.path(outdir, "color_spde_fixed_effects.csv"))
fwrite(key, file.path(outdir, "color_spde_key_effects.csv"))
fwrite(scores, file.path(outdir, "color_spde_model_scores.csv"))
fwrite(mech_audit, file.path(outdir, "color_spde_mechanism_variable_audit.csv"))
write_json(
  list(
    role = "spatial sensitivity analysis for canonical flower-colour endpoint",
    spatial_model = "category-specific independent spatial fields on a shared EPSG:6933 mesh",
    canonical_model_replaced = FALSE,
    global_population = "all colour-resolved flora",
    northern_mechanism_population = "animal-pollinated colour-resolved northern-midlatitude flora",
    mesh = mesh_meta
  ),
  file.path(outdir, "color_spde_metadata.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
