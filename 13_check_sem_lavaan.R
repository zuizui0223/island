############################################################
## SEM: PDI -> Color, Area -> Color/Shape, SC -> Shape
##  - NO: PDI -> Shape
##  - Show significant paths only
##  - Positive = Blue, Negative = Red
############################################################

## =========================================================
## 0) Packages
## =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(lavaan)
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
  library(scales)
})

alpha <- 0.05

## =========================================================
## 1) Load SEM-ready data
## =========================================================
infile <- file.path(Sys.getenv("HOME"),
                    "Desktop/list_with_traits_island_distance_bom_with_PDI_SEMready.csv")
dat <- fread(infile, encoding = "UTF-8")
cat("Rows loaded:", nrow(dat), "\n")

## =========================================================
## 2) Safety checks
## =========================================================
need <- c("dist_to_continent_km","area_km2","mean_PDI",
          "color_simple_bin","shape_simple_bin","self_capable",
          "pollination_main")
miss <- setdiff(need, names(dat))
if (length(miss)) stop("Missing columns in SEMready.csv: ", paste(miss, collapse=", "))


## =========================================================
## 3) Predictors: transform & scale
## =========================================================
dat <- dat[dist_to_continent_km > 0 & area_km2 > 0]

dat[, logdist := log10(dist_to_continent_km)]
dat[, logarea := log10(area_km2)]

dat[, z_logdist := as.numeric(scale(logdist))]
dat[, z_logarea := as.numeric(scale(logarea))]
dat[, z_PDI     := as.numeric(scale(mean_PDI))]
dat[, self_capable := as.numeric(self_capable)]

## =========================================================
## 4) SEM dataset (ordered outcomes only)
## =========================================================
dat_sem <- copy(dat)
dat_sem[, color_simple_bin := ordered(color_simple_bin, levels=c(0,1))]
dat_sem[, shape_simple_bin := ordered(shape_simple_bin, levels=c(0,1))]

cat("\nCounts check:\n")
print(table(dat_sem$color_simple_bin))
print(table(dat_sem$shape_simple_bin))
print(table(dat_sem$self_capable))

## =========================================================
## 5) SEM model
##  - Dist -> PDI
##  - PDI -> SC
##  - PDI -> Color
##  - Area -> Color, Shape
##  - SC -> Shape
##  - NO: PDI -> Shape
## =========================================================
mod_pscs_area <- '
  ## isolation -> pollinator deficit
  z_PDI ~ a1*z_logdist

  ## pollinator deficit -> self-compatibility
  self_capable ~ c1*z_PDI

  ## trait evolution
  color_simple_bin ~ b1*z_PDI + e1*z_logarea
  shape_simple_bin ~ s2*self_capable + e2*z_logarea

  ## indirect effects
  ind_dist_color_viaPDI := a1*b1
  ind_dist_shape_viaSC  := a1*c1*s2
'

## =========================================================
## 6) Fit (WLSMV: robust for categorical)
## =========================================================
fit_sem_wlsmv <- function(model_syntax, data_sem){
  sem(
    model_syntax,
    data = data_sem,
    ordered = c("color_simple_bin","shape_simple_bin"),
    estimator = "WLSMV",
    mimic = "Mplus",
    parameterization = "theta",
    std.lv = TRUE,
    warn = TRUE
  )
}

cat("\n=== Fitting SEM (Area as founder effect) ===\n")
fit_core <- fit_sem_wlsmv(mod_pscs_area, dat_sem)

cat("\n--- Fit measures ---\n")
print(fitMeasures(fit_core, c("cfi","tli","rmsea","srmr")))

## =========================================================
## 7) Extract significant paths
## =========================================================
pe <- parameterEstimates(fit_core, standardized=TRUE)

## direct effects
pe_dir <- pe[pe$op=="~" & pe$pvalue < alpha,
             c("lhs","rhs","label","est","std.all","pvalue")]

## indirect effects
pe_ind <- pe[pe$op==":=" & pe$pvalue < alpha,
             c("lhs","rhs","label","est","std.all","pvalue")]

cat("\n--- Significant direct paths ---\n")
print(pe_dir)

cat("\n--- Significant indirect effects ---\n")
print(pe_ind)

## =========================================================
## 8) Graphviz plot (significant direct paths only)
## =========================================================

make_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

## Node labels
node_labels <- c(
  z_logdist = "Distance to continent",
  z_PDI = "Pollinator deficit (PDI)",
  self_capable = "Self-compatible (SC)",
  z_logarea = "Island area",
  color_simple_bin = "Color simplification",
  shape_simple_bin = "Shape simplification"
)

nodes <- names(node_labels)
node_txt <- paste(sprintf(
  '%s [label="%s", shape=box, style=rounded, fontname="Helvetica"];',
  make_id(nodes),
  node_labels
), collapse = "\n")

## Edges (significant only)
edges_txt <- ""
if(nrow(pe_dir) > 0){
  edges_txt <- paste(sprintf(
    '%s -> %s [label="%.2f", color="%s", penwidth=%.2f];',
    make_id(pe_dir$rhs),
    make_id(pe_dir$lhs),
    pe_dir$std.all,
    ifelse(pe_dir$std.all >= 0, "#2B8CBE", "#E34A33"),
    scales::rescale(abs(pe_dir$std.all), c(1.5, 4))
  ), collapse="\n")
}

graph_txt <- paste0(
  'digraph G {\n',
  'rankdir=LR;\n',
  'graph [fontsize=16, fontname="Helvetica"];\n',
  'node  [fontname="Helvetica"];\n',
  node_txt, '\n',
  edges_txt, '\n',
  '}\n'
)

## =========================================================
## 9) Export
## =========================================================
svg <- export_svg(grViz(graph_txt))

out_dir <- file.path(Sys.getenv("HOME"), "Desktop/sem_PDI_SC_Area")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(out_dir, "SEM_sigpaths_PDI_SC_Area.png")
rsvg_png(charToRaw(svg), out_png, width=2200, height=1300)

## Save tables
write.csv(pe, file.path(out_dir, "lavaan_parameterEstimates_full.csv"), row.names=FALSE)
write.csv(pe_dir, file.path(out_dir, "sig_direct_paths.csv"), row.names=FALSE)
write.csv(pe_ind, file.path(out_dir, "sig_indirect_effects.csv"), row.names=FALSE)

cat("\n✅ Exported:\n", out_png, "\nTables in:\n", out_dir, "\n")