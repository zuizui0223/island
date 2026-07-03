library(sf)
library(dplyr)

# --- ① 等面積投影で面積計算 ---
big_proj <- st_transform(big, 6933)
big_proj$area_km2 <- as.numeric(st_area(big_proj)) / 1e6

# --- ② 20 km² 以上の島を抽出 ---
islands <- big_proj %>% filter(area_km2 >= 20)

# --- ③ GBIF 用の座標系（WGS84）へ戻す ---
islands <- st_transform(islands, 4326)

# --- ④ 保存するディレクトリを作成 ---
dir.create("data_global_islands", showWarnings = FALSE)

# --- ⑤ GeoPackage（推奨）で保存 ---
st_write(
  islands,
  "data_global_islands/global_islands_over20km2.gpkg",
  delete_dsn = TRUE
)

# --- ⑥ CSV（属性のみ）も保存 ---
islands_df <- islands %>% 
  st_drop_geometry()

write.csv(
  islands_df,
  "data_global_islands/global_islands_over20km2.csv",
  row.names = FALSE
)

# --- ⑦ 確認メッセージ ---
cat("保存完了：島数 =", nrow(islands), "\n")



