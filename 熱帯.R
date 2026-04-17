############################################################
## TROPICAL analysis
## from animal-pollinated master dataset
##
## 方針:
## - 入力は animal_trop_temp ファイル
## - tropical のみ抽出
## - 鳥媒を含む動物媒をそのまま使う
## - PDI は主解析に入れない
## - 熱帯では「温帯型のPDI中心モデルが当てはまりにくい」ことを
##   構成比と島レベル要約からみる
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(brms)　
  library(ggplot2)
})

set.seed(1)

############################################################
## 1) Load data
############################################################
infile <- "/Users/rachelzhang/Desktop/list_with_traits_island_distance_bom_with_PDI_reclassified_isolated_animal_trop_temp.csv"
dat <- fread(infile, encoding = "UTF-8")

cat("Rows loaded:", nrow(dat), "\n")
cat("Columns:\n")
print(names(dat))

############################################################
## 2) Keep tropical only
############################################################
dat <- dat[climate_zone == "tropical"]

cat("Rows after tropical filter:", nrow(dat), "\n")
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
  "species",
  "USGS_ISID"
)
miss <- setdiff(need_cols, names(dat))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

############################################################
## 4) Build helper variables
############################################################
## self-capable
dat[, self_incomp_clean := tolower(trimws(self_incompatibility))]
dat[, self_capable := NA_integer_]

dat[self_incomp_clean %in% c("sc", "likely_sc"),
    self_capable := 1L]

dat[self_incomp_clean %in% c("si", "likely_si", "obligate_si"),
    self_capable := 0L]

## color helpers
dat[, is_white := fifelse(flower_color_simple == "white", 1L, 0L)]
dat[, is_inconspicuous := fifelse(flower_color_simple == "green_brown_inconspicuous", 1L, 0L)]
dat[, is_red_yellow := fifelse(
  flower_color_simple %in% c("red_pink", "yellow_orange"), 1L, 0L
)]

## shape helper
dat[, is_open_shape := fifelse(
  flower_shape_simple == "actinomorphic_open", 1L,
  fifelse(
    flower_shape_simple %in% c(
      "tubular", "zygomorphic", "brush_puff",
      "composite_head", "reduced_wind", "other"
    ),
    0L, NA_integer_
  )
)]

## pollination helper
dat[, is_bird_pollinated := fifelse(pollination_main == "birds", 1L, 0L)]
dat[, is_mixed_pollinated := fifelse(pollination_main == "mixed", 1L, 0L)]

## continuous predictors
dat <- dat[dist_to_continent_km > 0 & area_km2 > 0]
dat[, logdist := log10(dist_to_continent_km)]
dat[, logarea := log10(area_km2)]

cat("\n[Self-capable summary]\n")
print(dat[, .(
  n_all = .N,
  n_SC = sum(self_capable == 1, na.rm = TRUE),
  n_SI = sum(self_capable == 0, na.rm = TRUE),
  n_NA = sum(is.na(self_capable))
)])

############################################################
## 5) Basic composition summaries
############################################################
cat("\n[Flower color summary]\n")
print(dat[, .N, by = flower_color_simple][order(-N)])

cat("\n[Flower shape summary]\n")
print(dat[, .N, by = flower_shape_simple][order(-N)])

cat("\n[Pollination summary]\n")
print(dat[, .N, by = pollination_main][order(-N)])

cat("\n[Bird pollination proportion]\n")
print(dat[, .(
  n_all = .N,
  prop_birds = mean(is_bird_pollinated, na.rm = TRUE),
  prop_mixed = mean(is_mixed_pollinated, na.rm = TRUE)
)])

############################################################
## 6) Compare tropical vs temperate composition
############################################################
all_dat <- fread(infile, encoding = "UTF-8")

comp_color <- all_dat[, .N, by = .(climate_zone, flower_color_simple)]
comp_color[, prop := N / sum(N), by = climate_zone]

comp_shape <- all_dat[, .N, by = .(climate_zone, flower_shape_simple)]
comp_shape[, prop := N / sum(N), by = climate_zone]

comp_poll <- all_dat[, .N, by = .(climate_zone, pollination_main)]
comp_poll[, prop := N / sum(N), by = climate_zone]

############################################################
## 7) Plots
############################################################
p_color <- ggplot(comp_color, aes(x = flower_color_simple, y = prop, fill = climate_zone)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    title = "Flower color composition: tropical vs temperate",
    x = "flower_color_simple",
    y = "proportion"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_shape <- ggplot(comp_shape, aes(x = flower_shape_simple, y = prop, fill = climate_zone)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    title = "Flower shape composition: tropical vs temperate",
    x = "flower_shape_simple",
    y = "proportion"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_poll <- ggplot(comp_poll, aes(x = pollination_main, y = prop, fill = climate_zone)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    title = "Pollination guild composition: tropical vs temperate",
    x = "pollination_main",
    y = "proportion"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_color)
print(p_shape)
print(p_poll)

############################################################
## 8) Island-level summaries
############################################################
island_trop <- dat[, .(
  n_species = uniqueN(species),
  
  prop_white = mean(is_white, na.rm = TRUE),
  prop_inconspicuous = mean(is_inconspicuous, na.rm = TRUE),
  prop_red_yellow = mean(is_red_yellow, na.rm = TRUE),
  
  prop_open_shape = mean(is_open_shape, na.rm = TRUE),
  prop_SC = mean(self_capable, na.rm = TRUE),
  
  prop_birds = mean(is_bird_pollinated, na.rm = TRUE),
  prop_mixed = mean(is_mixed_pollinated, na.rm = TRUE),
  
  logdist = mean(logdist, na.rm = TRUE),
  logarea = mean(logarea, na.rm = TRUE)
), by = USGS_ISID]

cat("\n[Island-level summary]\n")
print(summary(island_trop))

############################################################
## 9) Beta helper
############################################################
## Beta回帰のため 0,1 を少し内側へ
squeeze01 <- function(x, n){
  ((x * (n - 1)) + 0.5) / n
}

n_isl <- nrow(island_trop)

island_trop[, prop_white_b := squeeze01(prop_white, n_isl)]
island_trop[, prop_inconspicuous_b := squeeze01(prop_inconspicuous, n_isl)]
island_trop[, prop_red_yellow_b := squeeze01(prop_red_yellow, n_isl)]
island_trop[, prop_open_shape_b := squeeze01(prop_open_shape, n_isl)]
island_trop[, prop_SC_b := squeeze01(prop_SC, n_isl)]
island_trop[, prop_birds_b := squeeze01(prop_birds, n_isl)]
island_trop[, prop_mixed_b := squeeze01(prop_mixed, n_isl)]

############################################################
## 10) Tropical models WITHOUT PDI
##     ここでは logdist + logarea のみ
############################################################
ctrl <- list(adapt_delta = 0.95, max_treedepth = 12)

fit_white_trop <- brm(
  prop_white_b ~ logdist + logarea,
  data = island_trop,
  family = Beta(),
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000,
  seed = 1,
  control = ctrl
)

fit_redyellow_trop <- brm(
  prop_red_yellow_b ~ logdist + logarea,
  data = island_trop,
  family = Beta(),
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000,
  seed = 1,
  control = ctrl
)

fit_SC_trop <- brm(
  prop_SC_b ~ logdist + logarea,
  data = island_trop[is.finite(prop_SC_b) & !is.na(prop_SC_b)],
  family = Beta(),
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000,
  seed = 1,
  control = ctrl
)

fit_birds_trop <- brm(
  prop_birds_b ~ logdist + logarea,
  data = island_trop,
  family = Beta(),
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000,
  seed = 1,
  control = ctrl
)

fit_mixed_trop <- brm(
  prop_mixed_b ~ logdist + logarea,
  data = island_trop,
  family = Beta(),
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000,
  seed = 1,
  control = ctrl
)

cat("\n=== Tropical model summaries ===\n")
print(summary(fit_white_trop))
print(summary(fit_redyellow_trop))
print(summary(fit_SC_trop))
print(summary(fit_birds_trop))
print(summary(fit_mixed_trop))

############################################################
## 11) Save outputs
############################################################
out_dir <- "/Users/rachelzhang/Desktop/tropical_animalpoll_noPDI"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(fit_white_trop, file.path(out_dir, "fit_white_trop.rds"))
saveRDS(fit_redyellow_trop, file.path(out_dir, "fit_redyellow_trop.rds"))
saveRDS(fit_SC_trop, file.path(out_dir, "fit_SC_trop.rds"))
saveRDS(fit_birds_trop, file.path(out_dir, "fit_birds_trop.rds"))
saveRDS(fit_mixed_trop, file.path(out_dir, "fit_mixed_trop.rds"))

fwrite(island_trop, file.path(out_dir, "tropical_island_level_summary.csv"), bom = TRUE)
fwrite(comp_color, file.path(out_dir, "compare_color_trop_vs_temp.csv"), bom = TRUE)
fwrite(comp_shape, file.path(out_dir, "compare_shape_trop_vs_temp.csv"), bom = TRUE)
fwrite(comp_poll, file.path(out_dir, "compare_pollination_trop_vs_temp.csv"), bom = TRUE)

cat("\nSaved tropical outputs in:\n", out_dir, "\n")
cat("\nDone.\n")