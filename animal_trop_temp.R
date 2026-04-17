############################################################
## FULL SCRIPT
## 目的:
## 1) 再分類済み・isolated済みデータを読む
## 2) 風媒・自殖・不明を除き、動物媒を残す
## 3) 島を tropical / temperate に分類
## 4) 植物データへ climate_zone を付与
## 5) 保存
##
## 入力:
## - /Users/rachelzhang/Desktop/island/list_with_traits_island_distance_bom_with_PDI_reclassified_isolated.csv
## - /Users/rachelzhang/Desktop/island/data_global_islands/global_islands_over20km2.gpkg
##
## 出力:
## - /Users/rachelzhang/Desktop/island_zone_table.csv
## - /Users/rachelzhang/Desktop/list_with_traits_island_distance_bom_with_PDI_reclassified_isolated_animal_trop_temp.csv
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(dplyr)
})

############################################################
## 1) Paths
############################################################
in_file <- "/Users/rachelzhang/Desktop/island/list_with_traits_island_distance_bom_with_PDI_reclassified_isolated.csv"
island_gpkg <- "/Users/rachelzhang/Desktop/island/data_global_islands/global_islands_over20km2.gpkg"

out_zone_file <- "/Users/rachelzhang/Desktop/island_zone_table.csv"
out_file <- "/Users/rachelzhang/Desktop/list_with_traits_island_distance_bom_with_PDI_reclassified_isolated_animal_trop_temp.csv"

############################################################
## 2) Read reclassified isolated dataset
############################################################
cat("Reading data...\n")
dat <- fread(in_file, encoding = "UTF-8")

cat("Rows loaded:", nrow(dat), "\n")
cat("Columns:\n")
print(names(dat))

stopifnot("USGS_ISID" %in% names(dat))
stopifnot("pollination_main" %in% names(dat))

############################################################
## 3) Keep animal-pollinated only
##
## pollination_main の定義:
## self, wind, birds, bumblebees, bees, moths,
## butterflies, flies, mixed, unknown
##
## ここでは除外:
##   - wind
##   - self
##   - unknown
##
## 残す:
##   - birds
##   - bumblebees
##   - bees
##   - moths
##   - butterflies
##   - flies
##   - mixed
############################################################
drop_poll <- c(
  "wind",
  "self",
  "unknown"
)

dat_animal <- dat[!pollination_main %in% drop_poll]

cat("Rows after animal-pollinated filter:", nrow(dat_animal), "\n")
cat("Unique islands after animal-pollinated filter:", uniqueN(dat_animal$USGS_ISID), "\n")
cat("Pollination summary after filter:\n")
print(dat_animal[, .N, by = pollination_main][order(-N)])

############################################################
## 4) Read island polygons
############################################################
cat("Reading island polygons...\n")
islands_sf <- st_read(island_gpkg, quiet = TRUE)
stopifnot("USGS_ISID" %in% names(islands_sf))

## 動物媒植物が存在する島だけ残す
target_ids <- unique(dat_animal$USGS_ISID)
islands_sel <- islands_sf %>%
  dplyr::filter(USGS_ISID %in% target_ids)

cat("Selected islands in polygon file:", nrow(islands_sel), "\n")

############################################################
## 5) Assign tropical / temperate by representative point
############################################################
cat("Fixing geometry and assigning climate zones...\n")

if (is.na(st_crs(islands_sel))) {
  st_crs(islands_sel) <- 4326
}

## geometry 修正
islands_sel <- st_make_valid(islands_sel)

## 緯度経度へ
islands_ll <- st_transform(islands_sel, 4326)

## s2 を一時的にオフ
old_s2 <- sf::sf_use_s2(FALSE)

## centroid より安全な representative point
pts <- st_point_on_surface(islands_ll)
coords <- st_coordinates(pts)

## s2 を元に戻す
sf::sf_use_s2(old_s2)

island_zone_dt <- data.table(
  USGS_ISID = islands_ll$USGS_ISID,
  centroid_lon = coords[, 1],
  centroid_lat = coords[, 2]
)

## 北回帰線・南回帰線
tropic_limit <- 23.4366

island_zone_dt[, climate_zone := fifelse(
  abs(centroid_lat) < tropic_limit,
  "tropical",
  "temperate"
)]

## 念のため重複除去
island_zone_dt <- unique(island_zone_dt, by = "USGS_ISID")

cat("Climate zone counts:\n")
print(island_zone_dt[, .N, by = climate_zone])

############################################################
## 6) Join climate zone back to plant data
############################################################
cat("Joining climate zones...\n")

dat_animal_zone <- merge(
  dat_animal,
  island_zone_dt,
  by = "USGS_ISID",
  all.x = TRUE
)

cat("Rows after join:", nrow(dat_animal_zone), "\n")
cat("Rows with climate zone:", sum(!is.na(dat_animal_zone$climate_zone)), "\n")
cat("Rows without climate zone:", sum(is.na(dat_animal_zone$climate_zone)), "\n")

cat("Climate zone summary in plant data:\n")
print(dat_animal_zone[, .N, by = climate_zone][order(climate_zone)])

############################################################
## 7) Save
############################################################
fwrite(island_zone_dt, out_zone_file, bom = TRUE)
cat("Saved island zone table to:\n", out_zone_file, "\n")

fwrite(dat_animal_zone, out_file, bom = TRUE)
cat("Saved final filtered dataset to:\n", out_file, "\n")

############################################################
## 8) Optional quick subsets
############################################################
dat_all  <- dat_animal_zone[!is.na(climate_zone)]
dat_temp <- dat_animal_zone[climate_zone == "temperate"]
dat_trop <- dat_animal_zone[climate_zone == "tropical"]

cat("\nQuick summary:\n")
cat("All animal-pollinated rows with zone:", nrow(dat_all), "\n")
cat("Temperate rows:", nrow(dat_temp), "\n")
cat("Tropical rows:", nrow(dat_trop), "\n")

cat("\nDone.\n")