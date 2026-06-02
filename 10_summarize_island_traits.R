############################################################
## Save island × species data to CSV
############################################################

## =========================================================
## 0. Package
## =========================================================
suppressPackageStartupMessages({
  library(data.table)
})

## =========================================================
## 1. Load RDS
## =========================================================
save_path <- "/Users/zui/Desktop/island_flowering_species_facet_progress.rds"
results <- readRDS(save_path)

## =========================================================
## 2. Clean species lists
## =========================================================
results_clean <- lapply(results, function(x) {
  x <- unique(x)
  x <- x[!is.na(x)]
  x
})

## =========================================================
## 3. Island × Species long table
## =========================================================
dt_island_species <- rbindlist(
  lapply(seq_along(results_clean), function(i) {
    data.table(
      island_id = i,
      species   = results_clean[[i]]
    )
  })
)

dt_island_species <- unique(dt_island_species)

## =========================================================
## 4. Save CSV
## =========================================================
out_csv <- "/Users/zui/Desktop/island_species_longtable.csv"

fwrite(
  dt_island_species,
  file = out_csv
)

## =========================================================
## 5. Confirmation
## =========================================================
cat("CSV saved to:\n", out_csv, "\n")
cat("Rows:", nrow(dt_island_species), "\n")
cat("Islands:", uniqueN(dt_island_species$island_id), "\n")
cat("Species:", uniqueN(dt_island_species$species), "\n")