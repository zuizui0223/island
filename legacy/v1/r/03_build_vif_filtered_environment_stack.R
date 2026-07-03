## =========================================================
## 0. Packages
## =========================================================
library(terra)
library(usdm)
library(data.table)

## =========================================================
## 1. Path
## =========================================================
env_dir <- file.path(Sys.getenv("HOME"), "Desktop/islandbee")

out_stack_vif <- file.path(env_dir, "env_stack_vif.tif")
out_vif_csv   <- file.path(env_dir, "vif_results.csv")

## =========================================================
## 2. Read all TIFs (no resample, no project)
## =========================================================
tif_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(tif_files) > 1)

cat("Number of raster layers:", length(tif_files), "\n")

ras_list <- lapply(tif_files, rast)

## 確認
lapply(ras_list, function(r){
  c(res = paste(res(r), collapse=","),
    ext = paste(ext(r), collapse=","),
    crs = as.character(crs(r)))
})

## =========================================================
## 3. Find common overlapping extent
##    （解像度・グリッドは変えない。重なる範囲だけ使う）
## =========================================================
common_ext <- Reduce(intersect, lapply(ras_list, ext))
cat("Common extent:\n")
print(common_ext)

## 各ラスタを「共通部分だけ」に crop
ras_crop <- lapply(ras_list, function(r){
  crop(r, common_ext)
})

## =========================================================
## 4. Create virtual stack (still different resolution OK)
## =========================================================
env_stack <- rast(ras_crop)
print(env_stack)

## =========================================================
## 5. Prepare data for VIF
##    - 重なった領域からランダムサンプル
## =========================================================
set.seed(1)

sample_n <- 50000   # 必要に応じて調整
vals <- spatSample(env_stack, size = sample_n, na.rm = TRUE, as.df = TRUE)

cat("Sampled rows for VIF:", nrow(vals), "\n")

## =========================================================
## 6. VIF stepwise selection
## =========================================================
vif_res <- vifstep(vals, th = 10)

print(vif_res)

selected_vars <- vif_res@results$Variables
cat("Selected variables after VIF:\n")
print(selected_vars)

## 保存
fwrite(as.data.table(vif_res@results), out_vif_csv)
cat("Saved VIF table:", out_vif_csv, "\n")

## =========================================================
## 7. Create final stack (VIF後)
##    ※ 解像度は各レイヤ元のまま
## =========================================================
env_stack_vif <- env_stack[[selected_vars]]
print(env_stack_vif)

writeRaster(env_stack_vif, out_stack_vif, overwrite = TRUE)
cat("Saved VIF-filtered stack:", out_stack_vif, "\n")

## =========================================================
## 8. Final check
## =========================================================
cat("Final number of layers:", nlyr(env_stack_vif), "\n")
cat("Layer names:\n")
print(names(env_stack_vif))

