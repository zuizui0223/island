## =========================================================
## 0. Packages
## =========================================================
library(data.table)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

## =========================================================
## 1. Force UTF-8 locale（文字化け防止）
## =========================================================
Sys.setlocale("LC_ALL", "en_US.UTF-8")

## =========================================================
## 2. File paths
## =========================================================
traits_file     <- file.path(Sys.getenv("HOME"), "Desktop/traits.csv")
list_file       <- file.path(Sys.getenv("HOME"), "Desktop/list.csv")
islandinfo_file <- file.path(Sys.getenv("HOME"), "Desktop/islandinfo.csv")
gpkg_file       <- "/Users/zui/data_global_islands/global_islands_over20km2.gpkg"

## 出力
out_file <- file.path(Sys.getenv("HOME"),
                      "Desktop/list_with_traits_island_distance_bom.csv")

## =========================================================
## 3. Read data (UTF-8 / BOM)
## =========================================================
cat("Reading input tables...\n")
traits_df <- fread(traits_file,     encoding = "UTF-8")
list_df   <- fread(list_file,       encoding = "UTF-8")
island_df <- fread(islandinfo_file, encoding = "UTF-8")

setDT(traits_df); setDT(list_df); setDT(island_df)

## =========================================================
## 4. Text normalization (species)
## =========================================================
clean_text <- function(x){
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- tolower(x)
  x <- trimws(x)
  x
}

traits_df[, species_clean := clean_text(species)]
list_df[,   species_clean := clean_text(species)]

## =========================================================
## 5. Join traits to list (exact species match)
## =========================================================
setkey(traits_df, species_clean)
setkey(list_df,   species_clean)

list_df <- traits_df[list_df, on = "species_clean"]

cat("After species match:\n")
cat("Rows with traits:", sum(!is.na(list_df$flower_color)), "\n")

## =========================================================
## 6. Genus-level fill for remaining NAs
## =========================================================
list_df[is.na(flower_color), genus := tstrsplit(species_clean, " ")[[1]]]
traits_df[, genus := tstrsplit(species_clean, " ")[[1]]]

setkey(traits_df, genus)
setkey(list_df,   genus)

list_df[is.na(flower_color),
        c("flower_color","flower_shape","pollination_guild",
          "pollination_notes","mating_system","self_incompatibility",
          "evidence_type","confidence") :=
          traits_df[.SD, on = "genus", mult = "first",
                    .(flower_color, flower_shape, pollination_guild,
                      pollination_notes, mating_system, self_incompatibility,
                      evidence_type, confidence)],
        .SDcols = "genus"]

cat("After genus match:\n")
cat("Rows with traits:", sum(!is.na(list_df$flower_color)), "\n")

## 後始末
list_df[, c("species_clean","genus") := NULL]

## =========================================================
## 7. Join islandinfo by island_id
## =========================================================
setkey(list_df, island_id)
setkey(island_df, island_id)

list_df <- island_df[list_df, on = "island_id"]

cat("Rows with island info:",
    sum(!is.na(list_df[[colnames(island_df)[2]]])), "\n")

## =========================================================
## 8. Compute distance to continent (BOUNDARY–BOUNDARY)
## =========================================================
cat("Reading island polygons...\n")
islands_sf <- st_read(gpkg_file, quiet = TRUE)

cat("Downloading land polygons...\n")
land_sf <- ne_download(scale = 110, type = "land",
                       category = "physical", returnclass = "sf")

## 投影（メートル）
islands_sf <- st_transform(islands_sf, 3857)
land_sf    <- st_transform(land_sf,    3857)

cat("Computing boundary-to-boundary distances...\n")

# 境界↔境界の距離
dist_mat <- st_distance(islands_sf, land_sf)

# 各島の最小距離（m）
min_dist_m <- apply(dist_mat, 1, min)

# km へ変換
islands_sf$dist_to_continent_km <- as.numeric(min_dist_m) / 1000

## =========================================================
## 9. Build distance table by USGS_ISID
## =========================================================
island_dist_df <- as.data.table(islands_sf)
if ("geometry" %in% names(island_dist_df)) island_dist_df[, geometry := NULL]

stopifnot("USGS_ISID" %in% names(island_dist_df))

# 同一 USGS_ISID が複数ある場合は最小距離
island_dist_unique <- island_dist_df[
  , .(dist_to_continent_km = min(dist_to_continent_km, na.rm = TRUE)),
  by = USGS_ISID
]

cat("Distance summary (new):\n")
print(summary(island_dist_unique$dist_to_continent_km))

## =========================================================
## 10. Join distance to main table by USGS_ISID
## =========================================================
stopifnot("USGS_ISID" %in% names(list_df))

setkey(island_dist_unique, USGS_ISID)
setkey(list_df,           USGS_ISID)

list_df <- island_dist_unique[list_df, on = "USGS_ISID"]

## =========================================================
## 11. Check
## =========================================================
cat("Rows:", nrow(list_df), "\n")
cat("Distance non-NA:", sum(!is.na(list_df$dist_to_continent_km)), "\n")
cat("Range (km):\n")
print(range(list_df$dist_to_continent_km, na.rm = TRUE))

## =========================================================
## 12. Export (UTF-8 BOM for Excel)
## =========================================================
fwrite(list_df, out_file, bom = TRUE)
cat("Saved:", out_file, "\n")