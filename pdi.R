
library(data.table)
library(sf)
library(dplyr)

## =========================================================
## 1. ファイルパス
## =========================================================
list_file  <- file.path(Sys.getenv("HOME"), "Desktop/list_with_traits_island_distance_bom.csv")
island_gpkg <- "/Users/zui/data_global_islands/global_islands_over20km2.gpkg"

## =========================================================
## 2. 植物リスト（島ID付き）を読む
## =========================================================
cat("Reading plant list...\n")
island_list <- fread(list_file, encoding = "UTF-8")

# 列名チェック
print(names(island_list))

# USGS_ISID があることを確認
stopifnot("USGS_ISID" %in% names(island_list))

cat("Number of rows in list:", nrow(island_list), "\n")
cat("Number of unique islands in list:", length(unique(island_list$USGS_ISID)), "\n")

## =========================================================
## 3. 島ポリゴンを読む
## =========================================================
cat("Reading island polygons...\n")
islands_sf <- st_read(island_gpkg, quiet = TRUE)

# 列名チェック
print(names(islands_sf))

stopifnot("USGS_ISID" %in% names(islands_sf))

cat("Total islands in shapefile:", nrow(islands_sf), "\n")

## =========================================================
## 4. ★ USGS_ISID で島をフィルタ ★
## =========================================================
cat("Filtering islands by USGS_ISID in plant list...\n")

island_ids <- unique(island_list$USGS_ISID)

islands_sel_sf <- islands_sf %>%
  dplyr::filter(USGS_ISID %in% island_ids)

cat("Selected islands:", nrow(islands_sel_sf), "\n")

library(terra)

## sf → terra vector
islands_vect <- vect(islands_sel_sf)

## 島だけの環境
env_island <- mask(env_global, islands_vect)

## 保存
out_env_island <- file.path(Sys.getenv("HOME"), "Desktop/islandbee/env_island_vif.tif")
writeRaster(env_island, out_env_island, overwrite = TRUE)

cat("Saved env_island:", out_env_island, "\n")

## =========================================================
## 0) Packages
## =========================================================
library(data.table)
library(terra)
library(sf)

## =========================================================
## 1) Paths
## =========================================================
bombus_file <- file.path(Sys.getenv("HOME"), "Desktop/bombus_occ_thinned.csv")

# すでに作成済みのラスタ
env_global_file <- file.path(Sys.getenv("HOME"), "Desktop/islandbee/env_stack_vif.tif")
env_island_file <- file.path(Sys.getenv("HOME"), "Desktop/islandbee/env_island_vif.tif")

# 出力
out_csv <- file.path(Sys.getenv("HOME"), "Desktop/islandbee/island_PDI_bombus_FAST.csv")

## =========================================================
## 2) Read environment
## =========================================================
cat("Reading environment...\n")
env_global <- rast(env_global_file)
env_island <- rast(env_island_file)

## =========================================================
## 3) Read Bombus occurrences
## =========================================================
cat("Reading Bombus occurrences...\n")
bombus_raw <- fread(bombus_file, encoding="UTF-8")

## =========================================================
## 4) Occurrence preprocessing
## =========================================================
prep_occurrences <- function(df, lon_col, lat_col, sp_col="species",
                             min_n=30, max_per_species=2000){
  dt <- as.data.table(df)
  stopifnot(all(c(lon_col, lat_col, sp_col) %in% names(dt)))
  setnames(dt, c(lon_col, lat_col, sp_col), c("x","y","species"))
  
  dt <- dt[is.finite(x) & is.finite(y)]
  dt[, species := trimws(as.character(species))]
  dt <- dt[!is.na(species) & species != "" & species != "Bombus"]
  dt <- dt[!grepl("\\bsp\\.?\\b", species, ignore.case=TRUE)]
  dt <- dt[!grepl("cf\\.|aff\\.", species, ignore.case=TRUE)]
  
  sp_n <- dt[, .N, by=species]
  keep <- sp_n[N >= min_n, species]
  dt <- dt[species %in% keep]
  
  cat("After n >=", min_n, ": records =", nrow(dt),
      ", species =", uniqueN(dt$species), "\n")
  
  set.seed(1)
  dt[, rnd := runif(.N), by=species]
  setorder(dt, species, rnd)
  dt <- dt[, head(.SD, max_per_species), by=species]
  dt[, rnd := NULL]
  
  cat("After cap:", max_per_species, ": records =", nrow(dt),
      ", species =", uniqueN(dt$species), "\n")
  
  dt[, .(species, x, y)]
}

cat("Preprocessing Bombus...\n")
bombus_occ <- prep_occurrences(bombus_raw,
                               lon_col="decimalLongitude",
                               lat_col="decimalLatitude",
                               min_n=30,
                               max_per_species=2000)

## =========================================================
## 5) Build points vector
## =========================================================
occ_vect <- vect(bombus_occ, geom=c("x","y"), crs=crs(env_global))

## =========================================================
## 6) FAST extract with progress (chunked)
## =========================================================
extract_with_progress <- function(env, pts_vect, chunk_size = 5000, cores = 2){
  n <- nrow(pts_vect)
  idx <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  res_list <- vector("list", length(idx))
  
  cat("Total points:", n, "\n")
  cat("Chunk size:", chunk_size, "\n")
  cat("Chunks:", length(idx), "\n\n")
  
  t0 <- Sys.time()
  
  for(i in seq_along(idx)){
    id <- idx[[i]]
    cat(sprintf("Chunk %d / %d (rows %d–%d)... ",
                i, length(idx), min(id), max(id)))
    
    t1 <- Sys.time()
    res_list[[i]] <- extract(env, pts_vect[id, ], cores = cores)
    cat("done in",
        round(difftime(Sys.time(), t1, units="secs"), 2),
        "sec\n")
  }
  
  cat("\nTotal time:",
      round(difftime(Sys.time(), t0, units="mins"), 2),
      "minutes\n")
  
  rbindlist(res_list)
}

cat("Extracting env at occurrences (GLOBAL, fast)...\n")
occ_env <- extract_with_progress(
  env         = env_global,
  pts_vect   = occ_vect,
  chunk_size = 5000,   # 3000–10000で調整可
  cores      = 2
)

## species を戻す
occ_env[, species := bombus_occ$species]

## extract() が返す ID 列を削除
if("ID" %in% names(occ_env)) occ_env[, ID := NULL]

## =========================================================
## 7) Build niches (5–95%) = "world training"
## =========================================================
cat("Building species niches (5–95%)...\n")

occ_long <- melt(occ_env,
                 id.vars="species",
                 variable.name="var",
                 value.name="val")

niche_dt <- occ_long[
  , .(q05 = quantile(val, 0.05, na.rm=TRUE),
      q95 = quantile(val, 0.95, na.rm=TRUE)),
  by = .(species, var)
]

species_list <- unique(niche_dt$species)
cat("Species used:", length(species_list), "\n")

out_niche_rds <- file.path(Sys.getenv("HOME"), "Desktop/islandbee/bombus_niche_5_95.rds")

saveRDS(niche_dt, out_niche_rds)

cat("Saved niche_dt (RDS) to:", out_niche_rds, "\n")


library(data.table)
library(terra)
library(sf)

## =========================================================
## 1) Paths
## =========================================================
env_island_file <- file.path(Sys.getenv("HOME"), "Desktop/islandbee/env_island_vif.tif")
niche_rds_file  <- file.path(Sys.getenv("HOME"), "Desktop/islandbee/bombus_niche_5_95.rds")
island_gpkg     <- "/Users/zui/data_global_islands/global_islands_over20km2.gpkg"

# ★ 別CSV（新方式のPDI）
out_csv <- file.path(Sys.getenv("HOME"), "Desktop/islandbee/island_PDI_bombus_CONTINUOUS.csv")

# 途中経過ログ（再開用・安全）
log_file <- file.path(Sys.getenv("HOME"), "Desktop/islandbee/_pdi_progress.log")

## =========================================================
## 2) Read env_island
## =========================================================
cat("Reading env_island...\n")
env_island <- rast(env_island_file)

## =========================================================
## 3) Read niche (RDS)
## =========================================================
cat("Reading niche_dt from RDS...\n")
niche_dt <- readRDS(niche_rds_file)
setDT(niche_dt)

species_list <- unique(niche_dt$species)
cat("Species loaded:", length(species_list), "\n")

# 変数名
env_vars <- unique(niche_dt$var)

## speciesごとにニッチ情報をリスト化（高速化）
niche_list <- lapply(species_list, function(sp){
  rr <- niche_dt[species == sp]
  list(
    lo = rr$q05,
    hi = rr$q95
  )
})
names(niche_list) <- species_list

## =========================================================
## 4) Read island polygons
## =========================================================
cat("Reading island polygons...\n")
islands_sf  <- st_read(island_gpkg, quiet=TRUE)
islands_vect <- vect(islands_sf)
crs(islands_vect) <- crs(env_island)

## =========================================================
## 5) 植物が存在する島IDを高速取得（rasterize）
## =========================================================
cat("Rasterizing island IDs...\n")
tmpl <- env_island[[1]]
island_id_ras <- rasterize(islands_vect, tmpl, field="USGS_ISID", touches=TRUE)

cat("Extracting plant-present island IDs (FAST)...\n")
plant_island_ids <- unique(values(island_id_ras, na.rm=TRUE))
plant_island_ids <- plant_island_ids[!is.na(plant_island_ids)]
cat("Plant-present islands:", length(plant_island_ids), "\n")

# 対象島のみ
islands_sel <- islands_vect[islands_vect$USGS_ISID %in% plant_island_ids, ]
cat("Islands selected:", nrow(islands_sel), "\n")

## =========================================================
## 6) 連続型ニッチ：ESI → PDI = 1 - ESI
## =========================================================

## ---- 単一島のESI/PDIを計算する関数
calc_ESI_one_island <- function(isl_vals, niche_list, env_vars){
  
  # 行列化（数値演算を高速に）
  isl_mat <- as.matrix(isl_vals[, env_vars, drop=FALSE])
  ncell   <- nrow(isl_mat)
  nvar    <- ncol(isl_mat)
  
  out <- vector("list", length(niche_list))
  names(out) <- names(niche_list)
  
  for(sp in names(niche_list)){
    
    lo <- niche_list[[sp]]$lo
    hi <- niche_list[[sp]]$hi
    mid <- (lo + hi) / 2
    half_range <- (hi - lo) / 2
    
    # 連続型適合度（中心からの距離）
    s <- 1 - abs(isl_mat - matrix(mid, ncell, nvar, byrow=TRUE)) /
      matrix(half_range, ncell, nvar, byrow=TRUE)
    s[s < 0] <- 0  # 範囲外は0に切る
    
    S_cell <- rowMeans(s)
    ESI <- mean(S_cell, na.rm=TRUE)
    PDI <- 1 - ESI
    
    out[[sp]] <- data.table(
      species = sp,
      island_ESI = ESI,
      island_PDI = PDI
    )
  }
  
  rbindlist(out)
}

## =========================================================
## 7) 島×種：逐次保存（絶対落ちない設計）
## =========================================================
cat("Computing island-level PDI (continuous)...\n")

N_CELL_LIMIT <- 5000   # 島ごとの最大セル数（速度×精度のトレードオフ）

# 出力初期化
if (file.exists(out_csv)) file.remove(out_csv)
if (file.exists(log_file)) file.remove(log_file)

# 途中再開用：既に処理済みの島
done_ids <- character(0)

for(ii in seq_len(nrow(islands_sel))){
  
  isl_id <- as.character(islands_sel$USGS_ISID[ii])
  if(isl_id %in% done_ids) next
  
  poly_i <- islands_sel[ii, ]
  
  tryCatch({
    
    ## --- 島ごとに切り出し（空間演算はここ1回）
    isl_env <- mask(crop(env_island, poly_i), poly_i)
    isl_vals <- values(isl_env, na.rm=TRUE, dataframe=TRUE)
    if(nrow(isl_vals) == 0){
      cat("Skipping island (no cells):", isl_id, "\n")
      next
    }
    
    ## --- 巨大島はサンプリング
    if(nrow(isl_vals) > N_CELL_LIMIT){
      set.seed(1)
      isl_vals <- isl_vals[sample.int(nrow(isl_vals), N_CELL_LIMIT), ]
    }
    
    ## --- 種ごとのESI/PDI（連続）
    res_i <- calc_ESI_one_island(isl_vals, niche_list, env_vars)
    res_i[, USGS_ISID := isl_id]
    
    ## --- 逐次保存
    fwrite(res_i, out_csv, append=TRUE)
    fwrite(data.table(USGS_ISID=isl_id), log_file, append=TRUE)
    
    if(ii %% 10 == 0)
      cat("Processed islands:", ii, "/", nrow(islands_sel), "\n")
    
  }, finally = {
    # メモリ掃除
    rm(list=ls(pattern="isl_"))
    gc()
    terra::tmpFiles(remove=TRUE)
  })
}

cat("Finished.\nSaved to:\n", out_csv, "\n")

library(data.table)

# CSVの読み込み
res <- fread("~/Desktop/islandbee/island_PDI_bombus_CONTINUOUS.csv")

# 数値型に変換（念のため）
res[, island_ESI := as.numeric(island_ESI)]
res[, island_PDI := as.numeric(island_PDI)]

# RDSとして保存
out_rds <- "~/Desktop/islandbee/island_PDI_bombus_CONTINUOUS.rds"
saveRDS(res, out_rds)

cat("Saved RDS to:", out_rds, "\n")


library(data.table)

## =========================================================
## 1) ファイルパス
## =========================================================
list_file <- file.path(Sys.getenv("HOME"), "Desktop/list_with_traits_island_distance_bom.csv")
pdi_file  <- file.path(Sys.getenv("HOME"), "Desktop/islandbee/island_PDI_bombus_CONTINUOUS.csv")

# 出力
out_file  <- file.path(Sys.getenv("HOME"), "Desktop/list_with_traits_island_distance_bom_with_PDI.csv")

## =========================================================
## 2) 植物リストを読む
## =========================================================
cat("Reading plant list...\n")
plant_dt <- fread(list_file, encoding="UTF-8")
stopifnot("USGS_ISID" %in% names(plant_dt))

cat("Plants:", nrow(plant_dt), "\n")
cat("Unique islands (plants):", uniqueN(plant_dt$USGS_ISID), "\n")

## =========================================================
## 3) PDIデータを読む
## =========================================================
cat("Reading PDI results...\n")
pdi_dt <- fread(pdi_file)

# 念のため数値化
pdi_dt[, island_ESI := as.numeric(island_ESI)]
pdi_dt[, island_PDI := as.numeric(island_PDI)]

cat("PDI rows:", nrow(pdi_dt), "\n")
cat("Unique islands (PDI):", uniqueN(pdi_dt$USGS_ISID), "\n")

## =========================================================
## 4) 島ごとにPDIを要約（平均）
## =========================================================
cat("Aggregating PDI by island...\n")

island_pdi_mean <- pdi_dt[
  , .(
    mean_PDI = mean(island_PDI, na.rm=TRUE),
    mean_ESI = mean(island_ESI, na.rm=TRUE),
    n_species = sum(!is.na(island_PDI))
  ),
  by = USGS_ISID
]

cat("Islands with PDI:", nrow(island_pdi_mean), "\n")

## =========================================================
## 5) 植物リストに結合
## =========================================================
cat("Joining PDI to plant list...\n")

plant_with_pdi <- merge(
  plant_dt,
  island_pdi_mean,
  by = "USGS_ISID",
  all.x = TRUE
)

## 確認
cat("Rows after join:", nrow(plant_with_pdi), "\n")
cat("Plants with PDI:", sum(!is.na(plant_with_pdi$mean_PDI)), "\n")

## =========================================================
## 6) 保存
## =========================================================
fwrite(plant_with_pdi, out_file)

cat("Saved merged file to:\n", out_file, "\n")

