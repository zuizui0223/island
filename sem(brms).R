############################################################
## Bayesian SEM replication of lavaan model (brms)
## Structure:
##   z_PDI ~ z_logdist
##   self_capable ~ z_PDI
##   color_simple_bin ~ z_PDI + z_logarea
##   shape_simple_bin ~ self_capable + z_logarea
## Indirect effects:
##   dist -> color via PDI = a1*b1
##   dist -> shape via PDI->SC = a1*c1*s2
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(brms)
  library(posterior)
})

set.seed(1)


## =========================================================
## 1) Load SEM-ready data
## =========================================================
infile <- file.path(Sys.getenv("HOME"),
                    "Desktop/list_with_traits_island_distance_bom_with_PDI_SEMready.csv")
dat <- fread(infile, encoding = "UTF-8")
cat("Rows loaded:", nrow(dat), "\n")

## =========================================================
## 3) Predictors (same as lavaan)
## =========================================================
need <- c("dist_to_continent_km","area_km2","mean_PDI",
          "color_simple_bin","shape_simple_bin","self_capable")

dat <- dat[complete.cases(dat[, ..need])]
dat <- dat[dist_to_continent_km > 0 & area_km2 > 0]

dat[, logdist := log10(dist_to_continent_km)]
dat[, logarea := log10(area_km2)]

dat[, z_logdist := as.numeric(scale(logdist))]
dat[, z_logarea := as.numeric(scale(logarea))]

## lavaanでは z_PDI は mean_PDI の z
dat[, z_PDI := as.numeric(scale(mean_PDI))]

## brmsのBernoulliは 0/1 numeric でOK
dat[, self_capable := as.integer(self_capable)]
dat[, color_simple_bin := as.integer(color_simple_bin)]
dat[, shape_simple_bin := as.integer(shape_simple_bin)]

cat("\nCounts:\n")
print(table(dat$color_simple_bin))
print(table(dat$shape_simple_bin))
print(table(dat$self_capable))

## =========================================================
## 4) Model specification (multivariate brms)
##    - No PDI -> Shape
##    - Area -> Color & Shape
## =========================================================

bf_pdi   <- bf(z_PDI ~ 1 + z_logdist)
bf_sc    <- bf(self_capable ~ 1 + z_PDI, family = bernoulli("logit"))
bf_color <- bf(color_simple_bin ~ 1 + z_PDI + z_logarea, family = bernoulli("logit"))
bf_shape <- bf(shape_simple_bin ~ 1 + self_capable + z_logarea, family = bernoulli("logit"))

## =========================================================
## 5) Priors (weakly informative, stable)
## =========================================================
pri <- c(
  ## PDI equation (gaussian)
  prior(normal(0, 1), class = "b", resp = "zPDI"),
  prior(normal(0, 1.5), class = "Intercept", resp = "zPDI"),
  prior(exponential(1), class = "sigma", resp = "zPDI"),
  
  ## logistic equations
  prior(normal(0, 1), class = "b", resp = "selfcapable"),
  prior(normal(0, 1.5), class = "Intercept", resp = "selfcapable"),
  
  prior(normal(0, 1), class = "b", resp = "colorsimplebin"),
  prior(normal(0, 1.5), class = "Intercept", resp = "colorsimplebin"),
  
  prior(normal(0, 1), class = "b", resp = "shapesimplebin"),
  prior(normal(0, 1.5), class = "Intercept", resp = "shapesimplebin")
)

ctrl <- list(adapt_delta = 0.95, max_treedepth = 12)

## =========================================================
## 6) Fit
## =========================================================
fit_bsem <- brm(
  bf_pdi + bf_sc + bf_color + bf_shape + set_rescor(FALSE),
  data = dat,
  prior = pri,
  chains = 4, cores = 4,
  iter = 3000, warmup = 1500,
  control = ctrl,
  seed = 1
)

print(fit_bsem)

## =========================================================
## 7) Posterior indirect effects (match lavaan definitions)
##    ind_dist_color_viaPDI := a1*b1
##    ind_dist_shape_viaSC  := a1*c1*s2
## =========================================================
draws <- as_draws_df(fit_bsem)

## 係数名は brms の内部名になるので、まず colnames を確認
## head(grep("b_", names(draws), value=TRUE))

## a1: z_PDI ~ z_logdist
a1 <- draws$b_zPDI_z_logdist

## c1: self_capable ~ z_PDI
c1 <- draws$b_selfcapable_z_PDI

## b1: color ~ z_PDI
b1 <- draws$b_colorsimplebin_z_PDI

## s2: shape ~ self_capable
s2 <- draws$b_shapesimplebin_self_capable

ind_dist_color_viaPDI <- a1 * b1
ind_dist_shape_viaSC  <- a1 * c1 * s2

summ <- function(x){
  c(mean = mean(x),
    median = median(x),
    q025 = quantile(x, 0.025),
    q975 = quantile(x, 0.975),
    pd = mean(x > 0))  ## posterior probability > 0
}

cat("\n=== Indirect effects (posterior summaries) ===\n")
print(rbind(
  ind_dist_color_viaPDI = summ(ind_dist_color_viaPDI),
  ind_dist_shape_viaSC  = summ(ind_dist_shape_viaSC)
))

## （任意）CSV保存
out_dir <- file.path(Sys.getenv("HOME"), "Desktop/sem_PDI_SC_Area_bayes")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fwrite(as.data.table(summ(ind_dist_color_viaPDI), keep.rownames="stat"),
       file.path(out_dir, "ind_dist_color_viaPDI_summary.csv"))
fwrite(as.data.table(summ(ind_dist_shape_viaSC), keep.rownames="stat"),
       file.path(out_dir, "ind_dist_shape_viaSC_summary.csv"))

saveRDS(fit_bsem, file.path(out_dir, "fit_bsem.rds"))
cat("\n✅ Saved bayesian SEM outputs in:", out_dir, "\n")

############################################################
## Bayesian SEM Path Diagram (NO isolated nodes)
############################################################

library(posterior)
library(dplyr)
library(stringr)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(scales)

## =========================================================
## 1) Posterior draws
## =========================================================
draws <- as_draws_df(fit_bsem)

## =========================================================
## 2) 回帰係数を抽出
## =========================================================
coef_names <- grep("^b_", names(draws), value = TRUE)

parse_path <- function(term){
  x <- sub("^b_", "", term)
  sp <- strsplit(x, "_")[[1]]
  to   <- sp[1]
  from <- paste(sp[-1], collapse = "_")
  tibble(from = from, to = to, par = term)
}

edge_map <- bind_rows(lapply(coef_names, parse_path))

## =========================================================
## 3) Posterior summary + 有意パス抽出
## =========================================================
summ_fun <- function(x){
  c(mean = mean(x),
    q025 = quantile(x, 0.025),
    q975 = quantile(x, 0.975),
    pd   = max(mean(x > 0), mean(x < 0)))
}

edge_df <- edge_map %>%
  rowwise() %>%
  mutate(
    mean = summ_fun(draws[[par]])["mean"],
    q025 = summ_fun(draws[[par]])["q025"],
    q975 = summ_fun(draws[[par]])["q975"],
    pd   = summ_fun(draws[[par]])["pd"],
    sig  = pd > 0.95
  ) %>%
  ungroup() %>%
  filter(sig)

## =========================================================
## 4) 表示ラベル（内部名を統合）
## =========================================================
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

edge_df <- edge_df %>%
  mutate(
    from_lab = label_map[from],
    to_lab   = label_map[to]
  ) %>%
  filter(!is.na(from_lab), !is.na(to_lab))

print(edge_df)

## =========================================================
## 5) ★ ノードは「有意エッジに出てきたものだけ」
## =========================================================
nodes <- unique(c(edge_df$from_lab, edge_df$to_lab))

make_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

node_txt <- paste(sprintf(
  '%s [label="%s", shape=box, style=rounded, fontname="Helvetica"];',
  make_id(nodes), nodes
), collapse = "\n")

## =========================================================
## 6) エッジ（青=正, 赤=負）
## =========================================================
edges_txt <- paste(sprintf(
  '%s -> %s [label="%.2f", color="%s", penwidth=%.2f];',
  make_id(edge_df$from_lab),
  make_id(edge_df$to_lab),
  edge_df$mean,
  ifelse(edge_df$mean >= 0, "#2B8CBE", "#E34A33"),
  scales::rescale(abs(edge_df$mean), to = c(1.5, 4))
), collapse = "\n")

## =========================================================
## 7) Graphviz（縦レイアウト）
## =========================================================
graph_txt <- paste0(
  'digraph G {\n',
  'rankdir=TB;\n',
  'graph [fontsize=16, fontname="Helvetica"];\n',
  'node  [fontname="Helvetica"];\n',
  node_txt, '\n',
  edges_txt, '\n',
  '}\n'
)

## =========================================================
## 8) 出力
## =========================================================
svg <- export_svg(grViz(graph_txt))

out_dir <- file.path(Sys.getenv("HOME"), "Desktop/sem_PDI_SC_Area_bayes")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(out_dir, "SEM_sigpaths_BAYES_merged_CLEAN.png")
rsvg_png(charToRaw(svg), out_png, width=2000, height=1300)

cat("\n✅ Left-side isolated nodes REMOVED.\nSaved to:\n", out_png, "\n")