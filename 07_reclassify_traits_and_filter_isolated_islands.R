############################################################
## 0) Packages
############################################################
library(data.table)
library(stringr)
library(brms)
library(ggplot2)

############################################################
## 1) Load raw data
############################################################
in_raw <- file.path(Sys.getenv("HOME"),
                    "Desktop/list_with_traits_island_distance_bom_with_PDI.csv")
dat <- fread(in_raw, encoding = "UTF-8")

cat("Rows loaded:", nrow(dat), "\n")

## Remove suspicious non-island records
if ("USGS_ISID" %in% names(dat)) {
  dat <- dat[USGS_ISID != 0]
}

############################################################
## 2) Flower color reclassification (strict, review-proof)
############################################################
dat[, flower_color_clean := tolower(trimws(flower_color))]
dat[, flower_color_simple := NA_character_]

dat[, has_flower := !str_detect(
  flower_color_clean,
  regex("none|flowerless|no_flower|catkin|cone|spadix|unknown",
        ignore_case = TRUE)
)]

table(dat$has_flower)


dat[, flower_color_simple := NA_character_]

## green / brown / inconspicuous
dat[has_flower == TRUE &
      str_detect(flower_color_clean,
                 regex("green|brown|inconspicuous|dull",
                       ignore_case=TRUE)),
    flower_color_simple := "green_brown_inconspicuous"]

## white
dat[has_flower == TRUE &
      str_detect(flower_color_clean,
                 regex("white|cream|pale|whitish",
                       ignore_case=TRUE)),
    flower_color_simple := "white"]

## yellow / orange
dat[has_flower == TRUE &
      str_detect(flower_color_clean,
                 regex("yellow|orange|gold",
                       ignore_case=TRUE)),
    flower_color_simple := "yellow_orange"]

## red / pink
dat[has_flower == TRUE &
      str_detect(flower_color_clean,
                 regex("red|pink|magenta|scarlet|crimson",
                       ignore_case=TRUE)),
    flower_color_simple := "red_pink"]

## blue / purple
dat[has_flower == TRUE &
      str_detect(flower_color_clean,
                 regex("blue|purple|violet|lilac|lavender",
                       ignore_case=TRUE)),
    flower_color_simple := "blue_purple"]

## other
dat[has_flower == TRUE & is.na(flower_color_simple),
    flower_color_simple := "other"]

dat_color <- dat[has_flower == TRUE]

cat("Removed non-flowering taxa:",
    nrow(dat) - nrow(dat_color), "\n")

print(table(dat_color$flower_color_simple))

############################################################
## Pollination guild -> pollination_main
## ルール：
## 1) "butterflies" を含むものは全て蝶媒にする（混合より優先）
## 2) その後に混合や他ギルドを判定
############################################################

dat[, pollination_guild_clean := tolower(str_squish(pollination_guild))]
pg <- dat$pollination_guild_clean
dat[, pollination_main := NA_character_]

## ---------------------------------------------------------
## 1) BUTTERFLIES 最優先
##    → 混合表記でもすべて蝶媒に分類
## ---------------------------------------------------------
dat[str_detect(pg, "butterfly|butterflies|lepidoptera"),
    pollination_main := "butterflies"]

## ---------------------------------------------------------
## 2) 明示的な混合（ただし蝶媒はすでに除外済み）
## ---------------------------------------------------------
dat[is.na(pollination_main) &
      str_detect(pg, "\\+|/| and | & |;|mixed"),
    pollination_main := "mixed"]

## ---------------------------------------------------------
## 3) 単一ギルド
## ---------------------------------------------------------

## self / 非送粉
dat[is.na(pollination_main) &
      str_detect(pg, "self|autogam|cleistogam|apomictic|flowerless|no_flower|spore"),
    pollination_main := "self"]

## wind / water
dat[is.na(pollination_main) &
      str_detect(pg, "wind|anemoph|water|hydroph"),
    pollination_main := "wind"]

## birds
dat[is.na(pollination_main) &
      str_detect(pg, "bird|birds|hummingbird|sunbird"),
    pollination_main := "birds"]

## bumblebees（bees より先）
dat[is.na(pollination_main) &
      str_detect(pg, "bumblebee|bumble"),
    pollination_main := "bumblebees"]

## bees
dat[is.na(pollination_main) &
      str_detect(pg, "\\bbee\\b|bees|euglossine|andrena|solitary"),
    pollination_main := "bees"]

## moths
dat[is.na(pollination_main) &
      str_detect(pg, "moth|moths|hawkmoth|nocturnal"),
    pollination_main := "moths"]

## flies
dat[is.na(pollination_main) &
      str_detect(pg, "fly|flies|gnat|midge|diptera"),
    pollination_main := "flies"]

## insects（広義）→ mixed
dat[is.na(pollination_main) &
      str_detect(pg, "\\binsect\\b|insects"),
    pollination_main := "mixed"]

## ---------------------------------------------------------
## 4) unknown
## ---------------------------------------------------------
dat[is.na(pollination_main) | pg == "",
    pollination_main := "unknown"]
############################################################
## Force into main pollination guild (9 categories)
############################################################

## typo
dat[pollination_main == "unknown" &
      pollination_guild_clean == "files",
    pollination_main := "flies"]

## beetles / wasps / ants / bats / fig wasps -> mixed
dat[pollination_main == "unknown" &
      pollination_guild_clean %in% c("beetles","wasps","ants","bats"),
    pollination_main := "mixed"]

dat[pollination_main == "unknown" &
      str_detect(pollination_guild_clean, "^fig_wasp"),
    pollination_main := "mixed"]

cat("\nPollination guild summary:\n")
print(table(dat$pollination_main))

############################################################
## 4) Mating system -> mating_simple
############################################################
dat[, mating_simple := NA_character_]

## obligate outcrossing
dat[str_detect(mating_system,
               regex("obligate_outcross|outcrossing|self_incompatible", ignore_case=TRUE)),
    mating_simple := "obligate_outcrossing"]

## mixed mating
dat[str_detect(mating_system, regex("mixed", ignore_case=TRUE)),
    mating_simple := "mixed_mating"]

## mainly selfing
dat[str_detect(mating_system,
               regex("mainly_self|mostly_self", ignore_case=TRUE)),
    mating_simple := "mainly_selfing"]

## obligate selfing
dat[str_detect(mating_system, regex("obligate_self|selfing", ignore_case=TRUE)) &
      is.na(mating_simple),
    mating_simple := "obligate_selfing"]

## unknown
dat[is.na(mating_simple) | mating_system == "" |
      str_detect(mating_system, regex("flowerless|no_flower|none", ignore_case=TRUE)),
    mating_simple := "unknown"]

cat("\nMating system summary:\n")
print(table(dat$mating_simple))

############################################################
## 5) Self-incompatibility -> self_incomp_simple
############################################################
dat[, self_incomp_simple := NA_character_]

dat[str_detect(self_incompatibility,
               regex("^SI$|obligate_SI", ignore_case=TRUE)),
    self_incomp_simple := "SI"]
dat[str_detect(self_incompatibility,
               regex("^SC$", ignore_case=TRUE)),
    self_incomp_simple := "SC"]
dat[str_detect(self_incompatibility,
               regex("likely_SI", ignore_case=TRUE)),
    self_incomp_simple := "likely_SI"]
dat[str_detect(self_incompatibility,
               regex("likely_SC", ignore_case=TRUE)),
    self_incomp_simple := "likely_SC"]
dat[is.na(self_incomp_simple) | self_incompatibility == "",
    self_incomp_simple := "unknown"]

cat("\nSelf-incompatibility summary:\n")
print(table(dat$self_incomp_simple))

############################################################
## 6) Flower SHAPE（カテゴリー分類・完成版）
############################################################

dat[, flower_shape_clean := tolower(trimws(flower_shape))]
dat[, flower_shape_simple := NA_character_]

## ---------------------------------------------------------
## 6-1) no_flower（最優先・非被子植物や無花）
## ---------------------------------------------------------
dat[str_detect(flower_shape_clean,
               "flowerless|no_flower|none|fern|gymnosperm|cone|spore|catkin|spadix|syconium"),
    flower_shape_simple := "no_flower"]

cat("\n[STEP 6-1] no_flower:\n")
print(table(dat$flower_shape_simple, useNA = "ifany"))

## ---------------------------------------------------------
## 6-2) reduced_wind（風媒・縮退花）
## ---------------------------------------------------------
dat[str_detect(flower_shape_clean,
               "wind|anemoph|apetalous|minute|tiny|inconspicuous|grass|spikelet|sedg|rush|cyperaceae") &
      is.na(flower_shape_simple),
    flower_shape_simple := "reduced_wind"]

cat("\n[STEP 6-2] reduced_wind:\n")
print(table(dat$flower_shape_simple, useNA = "ifany"))

## ---------------------------------------------------------
## 6-3) composite_head（頭状花）
## ---------------------------------------------------------
dat[str_detect(flower_shape_clean,
               "capitul|composite|head|aster|daisy|radiate") &
      is.na(flower_shape_simple),
    flower_shape_simple := "composite_head"]

cat("\n[STEP 6-3] composite_head:\n")
print(table(dat$flower_shape_simple, useNA = "ifany"))

## ---------------------------------------------------------
## 6-4) brush_puff（ブラシ・パフ型）
## ---------------------------------------------------------
dat[str_detect(flower_shape_clean,
               "brush|puff|powder|pom|stamen|globose") &
      is.na(flower_shape_simple),
    flower_shape_simple := "brush_puff"]

cat("\n[STEP 6-4] brush_puff:\n")
print(table(dat$flower_shape_simple, useNA = "ifany"))

## ---------------------------------------------------------
## 6-5) zygomorphic（左右相称・特殊化）
## ---------------------------------------------------------
dat[str_detect(flower_shape_clean,
               "zygomorphic|bilabiate|bilateral|papilion|orchid|labiate|pea") &
      is.na(flower_shape_simple),
    flower_shape_simple := "zygomorphic"]

cat("\n[STEP 6-5] zygomorphic:\n")
print(table(dat$flower_shape_simple, useNA = "ifany"))

## ---------------------------------------------------------
## 6-6) tubular（管状）
## ---------------------------------------------------------
dat[str_detect(flower_shape_clean,
               "tubular|trumpet|funnel|salverform|pipe|urceolate|urn") &
      is.na(flower_shape_simple),
    flower_shape_simple := "tubular"]

cat("\n[STEP 6-6] tubular:\n")
print(table(dat$flower_shape_simple, useNA = "ifany"))

## ---------------------------------------------------------
## 6-7) actinomorphic_open（放射相称・開放）
## ---------------------------------------------------------
dat[str_detect(flower_shape_clean,
               "actinomorphic|radial|rotate|open|bowl|cup|star|campanulate|umbel|panicle|raceme|cluster|cruciform|crucifer|petaled|rosaceous|geranium|bell") &
      is.na(flower_shape_simple),
    flower_shape_simple := "actinomorphic_open"]

cat("\n[STEP 6-7] actinomorphic_open:\n")
print(table(dat$flower_shape_simple, useNA = "ifany"))

## ---------------------------------------------------------
## 6-8) other（形態として定義不能な残り）
## ---------------------------------------------------------
dat[is.na(flower_shape_simple),
    flower_shape_simple := "other"]

cat("\n[FINAL] Flower shape categories:\n")
print(table(dat$flower_shape_simple))
############################################################
## 7) Exclude non-flowering records
############################################################
n_before <- nrow(dat)
dat <- dat[flower_shape_simple != "no_flower"]
cat("\nExcluded no_flower:", n_before - nrow(dat), "\n")
cat("Remaining rows:", nrow(dat), "\n")

############################################################
## 8) Save reclassified file
############################################################
out_recl <- file.path(Sys.getenv("HOME"),
                      "Desktop/list_with_traits_island_distance_bom_with_PDI_reclassified.csv")
fwrite(dat, out_recl, bom = TRUE)
cat("\nSaved reclassified file to:\n", out_recl, "\n")

library(data.table)

# 読み込み（すでに dat があるなら不要）
dat <- fread("~/Desktop/list_with_traits_island_distance_bom_with_PDI_reclassified.csv")

## =========================================================
## 1) 不正・非島データの除外
## =========================================================

# USGS_ISID = 0 を除外
dat_clean <- dat[USGS_ISID != 0]

# 距離 = 0 km（大陸接続・ポリゴン接触）を除外
dat_clean <- dat_clean[dist_to_continent_km > 0]

## =========================================================
## 2) 確認
## =========================================================
cat("Before filtering:", nrow(dat), "\n")
cat("After filtering :", nrow(dat_clean), "\n")

summary(dat_clean$dist_to_continent_km)

## =========================================================
## 3) 保存（解析用ファイル）
## =========================================================
out_file2 <- file.path(Sys.getenv("HOME"),
                       "Desktop/list_with_traits_island_distance_bom_with_PDI_reclassified_isolated.csv")

fwrite(dat_clean, out_file2, bom = TRUE)
cat("Saved isolated-island dataset to:\n", out_file2, "\n")