library(rgbif)

# Bombus の taxonKey（GBIF backbone）
key <- name_backbone(name = "Bombus", rank = "genus")$usageKey
key


# ---- 1. ダウンロードクエリを発行 ----
# 学名 Bombus*、自然個体（HUMAN_OBSERVATION, PRESERVED_SPECIMEN）を対象
req <- occ_download(
  pred("taxonKey", key),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_in("basisOfRecord", 
          c("HUMAN_OBSERVATION", "OBSERVATION", "PRESERVED_SPECIMEN")),
  format = "SIMPLE_CSV"
)

req

# ステータス確認
occ_download_meta(req)

# 完了まで待つ
occ_download_wait(req)

# ダウンロードして読み込む
d_file <- occ_download_get(req)
bombus_df <- occ_download_import(d_file)

dim(bombus_df)
head(bombus_df)

library(rgbif)

# Bombus の taxonKey（GBIF backbone）
key <- name_backbone(name = "Apis", rank = "genus")$usageKey
key


# ---- 1. ダウンロードクエリを発行 ----
# 学名 Bombus*、自然個体（HUMAN_OBSERVATION, PRESERVED_SPECIMEN）を対象
req <- occ_download(
  pred("taxonKey", key),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_in("basisOfRecord", 
          c("HUMAN_OBSERVATION", "OBSERVATION", "PRESERVED_SPECIMEN")),
  format = "SIMPLE_CSV"
)

req

# ステータス確認
occ_download_meta(req)

# 完了まで待つ
occ_download_wait(req)

# ダウンロードして読み込む
a_file <- occ_download_get(req)
apis_df <- occ_download_import(a_file)

dim(apis_df)
head(apis_df)

library(dplyr)
library(sf)

# 座標が必要
bombus_clean <- bombus_df %>%
  filter(!is.na(decimalLongitude),
         !is.na(decimalLatitude)) %>%
  distinct(decimalLongitude, decimalLatitude, species, .keep_all = TRUE)

thin_sf <- function(df, km = 5){
  
  pts <- st_as_sf(df,
                  coords = c("decimalLongitude", "decimalLatitude"),
                  crs = 4326) |> 
    st_transform(3857)
  
  coords <- st_coordinates(pts)
  
  df$grid_x <- floor(coords[,1] / (km*1000))
  df$grid_y <- floor(coords[,2] / (km*1000))
  
  df_thin <- df %>%
    group_by(species, grid_x, grid_y) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(-grid_x, -grid_y)
  
  return(df_thin)
}
bombus_thinned_list <- bombus_clean %>%
  group_split(species) %>%
  lapply(function(df){
    thin_sf(df, km = 5)
  })

bombus_thin <- bind_rows(bombus_thinned_list)

cat("Bombus レコード数（thin後）:", nrow(bombus_thin), "\n")
apis_clean <- apis_df %>%
  filter(!is.na(decimalLongitude),
         !is.na(decimalLatitude)) %>%
  distinct(decimalLongitude, decimalLatitude, species, .keep_all = TRUE)

apis_thin <- thin_sf(apis_clean, km = 5)

cat("Apis レコード数（thin後）:", nrow(apis_thin), "\n")
write.csv(bombus_thin, "~/Desktop/bombus_occ_thinned.csv", row.names = FALSE)
write.csv(apis_thin,   "~/Desktop/apis_occ_thinned.csv",   row.names = FALSE)



library(sf)
library(dplyr)

sf_use_s2(FALSE)

# ---- 1. 島データを読み込み & 修復 ----
islands <- st_read("/Users/zui/data_global_islands/global_islands_over20km2.gpkg")
islands <- st_make_valid(islands)
bombus <- read.csv("~/Desktop/bombus_occ_thinned.csv")
apis   <- read.csv("~/Desktop/apis_occ_thinned.csv")
# ---- 2. pollinator データを sf に変換 ----
bombus_sf <- st_as_sf(bombus,
                      coords = c("decimalLongitude", "decimalLatitude"),
                      crs = 4326)
bombus_sf <- st_make_valid(bombus_sf)

apis_sf <- st_as_sf(apis,
                    coords = c("decimalLongitude", "decimalLatitude"),
                    crs = 4326)
apis_sf <- st_make_valid(apis_sf)

# ---- 3. intersect ----
int_bombus <- st_intersects(islands, bombus_sf)
int_apis   <- st_intersects(islands, apis_sf)

# ---- 4. richness 計算 ----
Bombus_richness <- sapply(int_bombus, function(idx)
  if(length(idx)==0) 0 else length(unique(bombus_sf$species[idx])))

Apis_richness <- sapply(int_apis, function(idx)
  if(length(idx)==0) 0 else length(unique(apis_sf$species[idx])))

# ---- 5. 島データに追加 ----
islands$Bombus_richness <- Bombus_richness
islands$Apis_richness   <- Apis_richness





