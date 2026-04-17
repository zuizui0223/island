############################################################
## SEM-ready dataset + Bayesian SEM
## temperate only / insect-pollinated only
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(brms)
  library(posterior)
  library(dplyr)
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
  library(scales)
})

set.seed(1)

############################################################
## 1) Load data
############################################################
infile <- "/Users/rachelzhang/Desktop/list_with_traits_island_distance_bom_with_PDI_reclassified_isolated_insect_trop_temp.csv"
dat <- fread(infile, encoding = "UTF-8")

cat("Rows loaded:", nrow(dat), "\n")

############################################################
## 2) Keep temperate only
############################################################
dat <- dat[climate_zone == "temperate"]

cat("Rows after temperate filter:", nrow(dat), "\n")
cat("Unique islands:", uniqueN(dat$USGS_ISID), "\n")

############################################################
## 3) Safety checks
############################################################
need_cols <- c(
  "pollination_main",
  "flower_color_simple",
  "flower_shape_simple",
  "self_incompatibility",
  "dist_to_continent_km",
  "area_km2",
  "mean_PDI"
)
miss <- setdiff(need_cols, names(dat))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

############################################################
## 4) COLOR (binary)
##    1 = simple (white / green-brown / other)
##    0 = colored (blue / red / yellow)
############################################################
dat[, color_simple_bin := NA_integer_]

dat[flower_color_simple %in% c(
  "white",
  "green_brown_inconspicuous",
  "other"
), color_simple_bin := 1L]

dat[flower_color_simple %in% c(
  "blue_purple",
  "red_pink",
  "yellow_orange"
), color_simple_bin := 0L]

cat("\n[Color check]\n")
print(dat[, .N, by = .(flower_color_simple, color_simple_bin)]
      [order(flower_color_simple)])

if (anyNA(dat$color_simple_bin)) {
  cat("\n[Dropping unknown flower color]\n")
  print(dat[is.na(color_simple_bin),
            .N, by = flower_color_simple][order(-N)][1:min(20, .N)])
  dat <- dat[!is.na(color_simple_bin)]
}

############################################################
## 5) SHAPE (binary)
##    1 = simple/generalist
##    0 = specialized
############################################################
dat[, shape_simple_bin := NA_integer_]

## simple / generalist
dat[flower_shape_simple %in% c(
  "actinomorphic_open",
  "brush_puff",
  "composite_head",
  "reduced_wind"
), shape_simple_bin := 1L]

## specialized
dat[flower_shape_simple %in% c(
  "tubular",
  "zygomorphic"
), shape_simple_bin := 0L]

cat("\n[Shape before exclusion]\n")
print(dat[, .N, by = .(flower_shape_simple, shape_simple_bin)]
      [order(flower_shape_simple)])

dat <- dat[!is.na(shape_simple_bin)]

cat("\n[Shape after exclusion]\n")
print(dat[, .N, by = shape_simple_bin])

############################################################
## 6) SELF-CAPABLE
##    1 = SC / likely_SC
##    0 = SI / likely_SI / obligate_SI
############################################################
dat[, self_incomp_clean := tolower(trimws(self_incompatibility))]
dat[, self_capable := NA_integer_]

dat[self_incomp_clean %in% c("sc", "likely_sc"),
    self_capable := 1L]

dat[self_incomp_clean %in% c("si", "likely_si", "obligate_si"),
    self_capable := 0L]

cat("\n[Self-capable mapping check]\n")
print(dat[, .N, by = .(self_incompatibility, self_capable)]
      [order(self_incompatibility)])

cat("\n[Self-capable summary]\n")
print(dat[, .(
  n_all = .N,
  n_selfcapable_1 = sum(self_capable == 1, na.rm = TRUE),
  n_selfcapable_0 = sum(self_capable == 0, na.rm = TRUE),
  n_na = sum(is.na(self_capable))
)])

############################################################
## 7) Final NA summary before SEM-ready save
############################################################
cat("\n[Final NA summary]\n")
print(dat[, .(
  n_all    = .N,
  na_color = sum(is.na(color_simple_bin)),
  na_shape = sum(is.na(shape_simple_bin)),
  na_self  = sum(is.na(self_capable)),
  na_pdi   = sum(is.na(mean_PDI)),
  na_dist  = sum(is.na(dist_to_continent_km)),
  na_area  = sum(is.na(area_km2))
)])

############################################################
## 8) Save SEM-ready dataset
############################################################
out_semready <- "/Users/rachelzhang/Desktop/temperate_insect_with_PDI_SEMready.csv"
fwrite(dat, out_semready, bom = TRUE)
cat("\nSaved SEM-ready dataset:\n", out_semready, "\n")

############################################################
## 9) Prepare predictors
############################################################
need <- c(
  "dist_to_continent_km",
  "area_km2",
  "mean_PDI",
  "color_simple_bin",
  "shape_simple_bin",
  "self_capable"
)

dat_sem <- copy(dat)
dat_sem <- dat_sem[complete.cases(dat_sem[, ..need])]
dat_sem <- dat_sem[dist_to_continent_km > 0 & area_km2 > 0]

dat_sem[, logdist := log10(dist_to_continent_km)]
dat_sem[, logarea := log10(area_km2)]

dat_sem[, z_logdist := as.numeric(scale(logdist))]
dat_sem[, z_logarea := as.numeric(scale(logarea))]
dat_sem[, z_PDI := as.numeric(scale(mean_PDI))]

dat_sem[, self_capable := as.integer(self_capable)]
dat_sem[, color_simple_bin := as.integer(color_simple_bin)]
dat_sem[, shape_simple_bin := as.integer(shape_simple_bin)]

cat("\nCounts in SEM dataset:\n")
print(table(dat_sem$color_simple_bin))
print(table(dat_sem$shape_simple_bin))
print(table(dat_sem$self_capable))

############################################################
## 10) Bayesian SEM
## Structure:
##   z_PDI ~ z_logdist
##   self_capable ~ z_PDI
##   color_simple_bin ~ z_PDI + z_logarea
##   shape_simple_bin ~ self_capable + z_logarea
############################################################
bf_pdi   <- bf(z_PDI ~ 1 + z_logdist)
bf_sc    <- bf(self_capable ~ 1 + z_PDI, family = bernoulli("logit"))
bf_color <- bf(color_simple_bin ~ 1 + z_PDI + z_logarea, family = bernoulli("logit"))
bf_shape <- bf(shape_simple_bin ~ 1 + self_capable + z_logarea, family = bernoulli("logit"))

pri <- c(
  prior(normal(0, 1), class = "b", resp = "zPDI"),
  prior(normal(0, 1.5), class = "Intercept", resp = "zPDI"),
  prior(exponential(1), class = "sigma", resp = "zPDI"),
  
  prior(normal(0, 1), class = "b", resp = "selfcapable"),
  prior(normal(0, 1.5), class = "Intercept", resp = "selfcapable"),
  
  prior(normal(0, 1), class = "b", resp = "colorsimplebin"),
  prior(normal(0, 1.5), class = "Intercept", resp = "colorsimplebin"),
  
  prior(normal(0, 1), class = "b", resp = "shapesimplebin"),
  prior(normal(0, 1.5), class = "Intercept", resp = "shapesimplebin")
)

ctrl <- list(adapt_delta = 0.95, max_treedepth = 12)

fit_bsem <- brm(
  bf_pdi + bf_sc + bf_color + bf_shape + set_rescor(FALSE),
  data = dat_sem,
  prior = pri,
  chains = 4, cores = 4,
  iter = 3000, warmup = 1500,
  control = ctrl,
  seed = 1
)

print(fit_bsem)

############################################################
## 11) Posterior indirect effects
## ind_dist_color_viaPDI := a1*b1
## ind_dist_shape_viaSC  := a1*c1*s2
############################################################
draws <- as_draws_df(fit_bsem)

a1 <- draws$b_zPDI_z_logdist
c1 <- draws$b_selfcapable_z_PDI
b1 <- draws$b_colorsimplebin_z_PDI
s2 <- draws$b_shapesimplebin_self_capable

ind_dist_color_viaPDI <- a1 * b1
ind_dist_shape_viaSC  <- a1 * c1 * s2

summ <- function(x){
  c(
    mean   = mean(x),
    median = median(x),
    q025   = quantile(x, 0.025),
    q975   = quantile(x, 0.975),
    pd     = mean(x > 0)
  )
}

cat("\n=== Indirect effects (posterior summaries) ===\n")
print(rbind(
  ind_dist_color_viaPDI = summ(ind_dist_color_viaPDI),
  ind_dist_shape_viaSC  = summ(ind_dist_shape_viaSC)
))

############################################################
## 12) Save outputs
############################################################
out_dir <- "/Users/rachelzhang/Desktop/sem_temperate_insect_PDI"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(fit_bsem, file.path(out_dir, "fit_bsem.rds"))

fwrite(as.data.table(summ(ind_dist_color_viaPDI), keep.rownames = "stat"),
       file.path(out_dir, "ind_dist_color_viaPDI_summary.csv"))

fwrite(as.data.table(summ(ind_dist_shape_viaSC), keep.rownames = "stat"),
       file.path(out_dir, "ind_dist_shape_viaSC_summary.csv"))

fwrite(dat_sem, file.path(out_dir, "SEM_analysis_dataset.csv"), bom = TRUE)

cat("\nSaved Bayesian SEM outputs in:\n", out_dir, "\n")

############################################################
## 13) Bayesian SEM Path Diagram
############################################################
coef_names <- grep("^b_", names(draws), value = TRUE)

parse_path <- function(term){
  x <- sub("^b_", "", term)
  sp <- strsplit(x, "_")[[1]]
  to   <- sp[1]
  from <- paste(sp[-1], collapse = "_")
  data.table(from = from, to = to, par = term)
}

edge_map <- rbindlist(lapply(coef_names, parse_path))

summ_fun <- function(x){
  c(
    mean = mean(x),
    q025 = quantile(x, 0.025),
    q975 = quantile(x, 0.975),
    pd   = max(mean(x > 0), mean(x < 0))
  )
}

edge_df <- edge_map[, {
  sx <- summ_fun(draws[[par]])
  .(mean = sx["mean"], q025 = sx["q025"], q975 = sx["q975"], pd = sx["pd"])
}, by = .(from, to, par)]

edge_df[, sig := pd > 0.95]
edge_df <- edge_df[sig == TRUE]

label_map <- c(
  z_logdist        = "Distance to continent",
  zPDI             = "Pollinator deficit (PDI)",
  z_PDI            = "Pollinator deficit (PDI)",
  selfcapable      = "Self-compatible (SC)",
  self_capable     = "Self-compatible (SC)",
  z_logarea        = "Island area",
  colorsimplebin   = "Color simplification",
  color_simple_bin = "Color simplification",
  shapesimplebin   = "Shape simplification",
  shape_simple_bin = "Shape simplification"
)

edge_df[, from_lab := label_map[from]]
edge_df[, to_lab   := label_map[to]]
edge_df <- edge_df[!is.na(from_lab) & !is.na(to_lab)]

print(edge_df)

nodes <- unique(c(edge_df$from_lab, edge_df$to_lab))
make_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

node_txt <- paste(sprintf(
  '%s [label="%s", shape=box, style=rounded, fontname="Helvetica"];',
  make_id(nodes), nodes
), collapse = "\n")

edge_txt <- paste(sprintf(
  '%s -> %s [label="%.2f", color="%s", penwidth=%.2f];',
  make_id(edge_df$from_lab),
  make_id(edge_df$to_lab),
  edge_df$mean,
  ifelse(edge_df$mean >= 0, "#2B8CBE", "#E34A33"),
  scales::rescale(abs(edge_df$mean), to = c(1.5, 4))
), collapse = "\n")

graph_txt <- paste0(
  'digraph G {\n',
  'rankdir=TB;\n',
  'graph [fontsize=16, fontname="Helvetica"];\n',
  'node [fontname="Helvetica"];\n',
  node_txt, '\n',
  edge_txt, '\n',
  '}\n'
)

svg <- export_svg(grViz(graph_txt))
out_png <- file.path(out_dir, "SEM_sigpaths_BAYES_temperate.png")
rsvg_png(charToRaw(svg), out_png, width = 2000, height = 1300)

cat("\nSaved SEM path diagram to:\n", out_png, "\n")
cat("\nDone.\n")