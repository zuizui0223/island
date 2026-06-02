############################################################
## Build SEM-ready dataset
##  - color_simple_bin
##  - shape_simple_bin
##  - self_capable
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

set.seed(1)

## =========================================================
## 1) Load data (isolated + reclassified)
## =========================================================
infile <- file.path(
  Sys.getenv("HOME"),
  "Desktop/list_with_traits_island_distance_bom_with_PDI_reclassified_isolated.csv"
)

dat <- fread(infile, encoding = "UTF-8")
cat("Rows loaded:", nrow(dat), "\n")

## =========================================================
## 2) Safety checks
## =========================================================
need_cols <- c(
  "pollination_main",
  "flower_color_simple",
  "flower_shape_simple",
  "self_incompatibility"
)
miss <- setdiff(need_cols, names(dat))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

## =========================================================
## 3) Pollination: drop unknown
## =========================================================
cat("\n[Pollination before]\n")
print(dat[, .N, by = pollination_main][order(-N)])

dat <- dat[!is.na(pollination_main)]
dat <- dat[pollination_main != "unknown"]

cat("\n[Pollination after removing unknown]\n")
print(dat[, .N, by = pollination_main][order(-N)])

## =========================================================
## 4) COLOR (binary)
##    1 = simple (white / green-brown / other)
##    0 = colored (blue / red / yellow)
## =========================================================
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

## 未分類（空文字・不明）は除外
if (anyNA(dat$color_simple_bin)) {
  cat("\n[Dropping unknown flower color]\n")
  print(dat[is.na(color_simple_bin),
            .N, by = flower_color_simple][order(-N)][1:20])
  dat <- dat[!is.na(color_simple_bin)]
}

## =========================================================
## 5) SHAPE (binary)
##    1 = generalist
##    0 = specialist
##    drop: other / no_flower
## =========================================================
dat[, shape_simple_bin := NA_integer_]

## generalist
dat[flower_shape_simple %in% c(
  "actinomorphic_open",
  "brush_puff",
  "composite_head",
  "reduced_wind"
), shape_simple_bin := 1L]

## specialist
dat[flower_shape_simple %in% c(
  "tubular",
  "zygomorphic"
), shape_simple_bin := 0L]

cat("\n[Shape before exclusion]\n")
print(dat[, .N, by = .(flower_shape_simple, shape_simple_bin)]
      [order(flower_shape_simple)])

## other / no_flower / NA を除外
dat <- dat[!is.na(shape_simple_bin)]

cat("\n[Shape after exclusion]\n")
print(dat[, .N, by = shape_simple_bin])

## =========================================================
## 6) SELF-CAPABLE (strict, SEM-consistent)
##    1 = SC / likely_SC
##    0 = SI / likely_SI / obligate_SI
## =========================================================
dat[, self_incomp_clean := tolower(trimws(self_incompatibility))]
dat[, self_capable := NA_integer_]

## self-compatible
dat[self_incomp_clean %in% c("sc", "likely_sc"),
    self_capable := 1L]

## self-incompatible
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

## =========================================================
## 7) Final NA summary
## =========================================================
cat("\n[Final NA summary]\n")
print(dat[, .(
  n_all    = .N,
  na_color = sum(is.na(color_simple_bin)),
  na_shape = sum(is.na(shape_simple_bin)),
  na_self  = sum(is.na(self_capable))
)])

## =========================================================
## 8) Save SEM-ready dataset
## =========================================================
out_csv <- file.path(
  Sys.getenv("HOME"),
  "Desktop/list_with_traits_island_distance_bom_with_PDI_SEMready.csv"
)
fwrite(dat, out_csv, bom = TRUE)

cat("\n✅ Saved SEM-ready dataset:\n", out_csv, "\n")