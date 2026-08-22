#!/usr/bin/env Rscript

# Export exact-name accepted-species native range sizes from the WCVP data
# package plus the WGSRPD level-3 geometry used by WCVP/POWO. This is a
# corroboration layer for endemism, not a replacement taxonomy for the island
# master.

suppressPackageStartupMessages({
  library(rWCVP)
  library(rWCVPdata)
  library(dplyr)
  library(readr)
  library(sf)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: export_wcvp_native_ranges.R ISLAND_TAXA_CSV OUTPUT_DIR")
}
input_path <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target <- readr::read_csv(input_path, show_col_types = FALSE) %>%
  transmute(accepted_species = trimws(as.character(accepted_species))) %>%
  filter(accepted_species != "") %>%
  distinct()

names_tbl <- rWCVPdata::wcvp_names
dist_tbl <- rWCVPdata::wcvp_distributions

required_names <- c("plant_name_id", "taxon_name", "taxon_rank", "taxon_status")
missing_names <- setdiff(required_names, names(names_tbl))
if (length(missing_names)) stop(paste("WCVP names missing:", paste(missing_names, collapse=", ")))
required_dist <- c("plant_name_id", "area_code_l3", "introduced", "extinct", "location_doubtful")
missing_dist <- setdiff(required_dist, names(dist_tbl))
if (length(missing_dist)) stop(paste("WCVP distributions missing:", paste(missing_dist, collapse=", ")))

accepted <- names_tbl %>%
  mutate(
    taxon_name = trimws(as.character(taxon_name)),
    taxon_rank = tolower(as.character(taxon_rank)),
    taxon_status = tolower(as.character(taxon_status))
  ) %>%
  filter(taxon_rank == "species", taxon_status == "accepted") %>%
  inner_join(target, by = c("taxon_name" = "accepted_species")) %>%
  select(accepted_species = taxon_name, plant_name_id) %>%
  distinct()

native <- dist_tbl %>%
  inner_join(accepted, by = "plant_name_id") %>%
  mutate(
    introduced = as.numeric(introduced),
    extinct = as.numeric(extinct),
    location_doubtful = as.numeric(location_doubtful),
    area_code_l3 = trimws(as.character(area_code_l3))
  ) %>%
  filter(introduced == 0, extinct == 0, location_doubtful == 0, area_code_l3 != "") %>%
  distinct(accepted_species, plant_name_id, area_code_l3)

summary <- native %>%
  arrange(accepted_species, area_code_l3) %>%
  group_by(accepted_species, plant_name_id) %>%
  summarise(
    n_native_l3 = n_distinct(area_code_l3),
    native_l3_codes = paste(sort(unique(area_code_l3)), collapse = "|"),
    .groups = "drop"
  )

all_matched <- accepted %>%
  left_join(summary, by = c("accepted_species", "plant_name_id")) %>%
  mutate(
    n_native_l3 = ifelse(is.na(n_native_l3), 0L, as.integer(n_native_l3)),
    native_l3_codes = ifelse(is.na(native_l3_codes), "", native_l3_codes)
  )

readr::write_csv(all_matched, file.path(out_dir, "wcvp_native_range_summary.csv.gz"))
st_write(rWCVP::wgsrpd3, file.path(out_dir, "wgsrpd_level3.gpkg"), delete_dsn = TRUE, quiet = TRUE)

manifest <- list(
  contract = "wcvp_native_range_exact_name_v1",
  n_target_species = nrow(target),
  n_exact_accepted_species = nrow(accepted),
  n_species_with_native_range = sum(all_matched$n_native_l3 > 0),
  n_single_native_l3_species = sum(all_matched$n_native_l3 == 1),
  endemism_policy = paste(
    "A single WCVP native TDWG-L3 region is corroboration only; final endemic",
    "classification also requires an independent GIFT source-reported endemic claim",
    "and focal-island membership in the same TDWG-L3 region."
  )
)
jsonlite::write_json(manifest, file.path(out_dir, "wcvp_native_range_manifest.json"), auto_unbox=TRUE, pretty=TRUE)
print(jsonlite::toJSON(manifest, auto_unbox=TRUE, pretty=TRUE))
