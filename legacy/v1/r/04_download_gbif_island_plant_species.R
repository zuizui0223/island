## ============================================================
## 0. パッケージ読み込み
## ============================================================
library(sf)
library(dplyr)
library(httr)
library(jsonlite)

## ============================================================
## 1. 島ポリゴンの準備（GlbIslands.gdb → 面積 ≥ 20 km²）
## ============================================================

# 1-1) GDB から読み込み（geometry が入った sf として読める前提）
big <- st_read("/Users/zui/Desktop/GlbIslands.gdb", layer = "BigIslands")

# 1-2) 面積計算のため等面積投影に変換
#     EPSG:6933 = World Cylindrical Equal Area
big_proj <- st_transform(big, 6933)

# 1-3) 面積（km²）を計算
big_proj$area_km2 <- as.numeric(st_area(big_proj)) / 1e6

# 1-4) 20 km² 以上の島だけ抽出
islands <- big_proj %>% filter(area_km2 >= 20)

# 1-5) GBIF 用に WGS84 に戻す
islands <- st_transform(islands, 4326)

# 確認
n_islands <- nrow(islands)
cat("対象島数:", n_islands, "\n")
print(islands[1, c("area_km2")])

## ============================================================
## 2. 1島分の種リストを facet API で取得する関数
## ============================================================
#  - geometry: 各島の bbox（POLYGON WKT）
#  - taxonKey: 7707728 = angiosperms
#  - facet=scientificName で種リストだけ取得
#  - limit=0 + facetLimit=2000 で、多くても 2000 種まで

get_species_facet_safe <- function(poly,
                                   sleep_sec    = 1.0,
                                   facet_limit  = 2000) {
  
  # レート制限を少し回避
  Sys.sleep(sleep_sec)
  
  # --- bbox → WKT POLYGON ---
  bb <- sf::st_bbox(poly)
  wkt <- sprintf(
    "POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))",
    bb["xmin"], bb["ymin"],
    bb["xmax"], bb["ymin"],
    bb["xmax"], bb["ymax"],
    bb["xmin"], bb["ymax"],
    bb["xmin"], bb["ymin"]
  )
  
  # --- URL を安全に組み立てる ---
  base_url <- "https://api.gbif.org/v1/occurrence/search"
  
  # httr::modify_url でクエリを指定（自動 URL エンコード）
  url <- httr::modify_url(
    base_url,
    query = list(
      geometry         = wkt,
      taxonKey         = 7707728,          # angiosperms
      hasCoordinate    = "true",
      year             = "2010,2025",
      limit            = 0,                # レコード本体は返さない
      facet            = "scientificName",
      facetLimit       = facet_limit
    )
  )
  
  # --- HTTP GET ---
  resp <- try(httr::GET(url), silent = TRUE)
  if (inherits(resp, "try-error")) {
    message("  ❌ HTTP error (GET に失敗)")
    return(NULL)
  }
  
  status <- httr::status_code(resp)
  
  # レート制限（429）の場合は NULL 返して、上位ループで空扱いにする
  if (status == 429L) {
    message("  ⚠️ HTTP 429 Too Many Requests → NULL 返す")
    return(NULL)
  }
  
  if (status >= 400L) {
    message("  ⚠️ HTTP status ", status, " → NULL 返す")
    return(NULL)
  }
  
  # --- JSON をパース（簡約を切って構造をそのまま保持） ---
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  res <- try(jsonlite::fromJSON(txt, simplifyVector = FALSE), silent = TRUE)
  if (inherits(res, "try-error")) {
    message("  ❌ JSON parse error")
    return(NULL)
  }
  
  # --- facets 構造の安全チェック ---
  if (is.null(res$facets)) {
    # この bbox にはオカレンス自体がない
    return(NULL)
  }
  
  f <- res$facets
  
  # facets が list でない場合（まれ）も防御
  if (!is.list(f) || length(f) == 0) {
    return(NULL)
  }
  
  # 最初の facet（scientificName）を取り出す
  # 期待される構造: f[[1]]$counts は list( list(name=..., count=...), ... )
  if (is.null(f[[1]]$counts)) {
    return(NULL)
  }
  
  counts <- f[[1]]$counts
  
  # counts が data.frame の場合と list の場合を両方ケア
  species_vec <- NULL
  
  if (is.data.frame(counts)) {
    # col "name" を期待
    if ("name" %in% names(counts)) {
      species_vec <- counts$name
    }
  } else if (is.list(counts)) {
    # list の場合は各要素の $name を取り出す
    species_vec <- vapply(counts, function(z) {
      if (is.list(z) && !is.null(z$name)) {
        as.character(z$name)
      } else {
        NA_character_
      }
    }, FUN.VALUE = character(1))
  }
  
  if (is.null(species_vec)) {
    return(NULL)
  }
  
  # NA を除いてユニーク種名
  species_vec <- unique(species_vec[!is.na(species_vec) & species_vec != ""])
  
  return(species_vec)
}

## ============================================================
## 3. 全島ループ + オートセーブ
## ============================================================

save_path <- "/Users/zui/Desktop/island_flowering_species_facet_progress.rds"
n_islands <- nrow(islands)

if (file.exists(save_path)) {
  results <- readRDS(save_path)
} else {
  results <- vector("list", n_islands)
}

for (i in seq_len(n_islands)) {
  cat("\n🌴 Island", i, "/", n_islands, "\n")
  
  # 既に埋まっていればスキップ
  if (!is.null(results[[i]])) {
    cat("   → already done, skip\n")
    next
  }
  
  sp <- get_species_facet_safe(islands[i, ], sleep_sec = 0.8, facet_limit = 2000)
  
  if (is.null(sp)) {
    cat("   → NULL / no data → empty\n")
    results[[i]] <- character(0)
  } else {
    cat("   → OK:", length(sp), "species\n")
    results[[i]] <- sp
  }
  
  if (i %% 5 == 0) {
    cat("💾 Autosaving at island", i, "\n")
    saveRDS(results, save_path)
  }
}

saveRDS(results, save_path)