#!/usr/bin/env Rscript

# Export source-backed island floristic status and geometries from GIFT.
#
# This script intentionally exports raw checklist/status fields plus metadata.
# It does not decide the final native/endemic/introduced class; that fail-closed
# decision is made downstream after spatial matching to the frozen GSHHG islands.

suppressPackageStartupMessages({
  library(GIFT)
  library(sf)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: fetch_gift_floristic_status.R OUTPUT_DIR [GIFT_VERSION]")
}

output_dir <- args[[1]]
gift_version <- if (length(args) >= 2) args[[2]] else "latest"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Fetching GIFT island angiosperm checklists, version=", gift_version)

# We request all suitable island-like regions first so the metadata and original
# checklist status columns are retained. The primary spatial matcher later keeps
# only entity_class == 'Island'; Island Group / Island Part are preserved here for
# audit rather than silently discarded.
checklists <- GIFT_checklists(
  taxon_name = "Angiospermae",
  complete_taxon = TRUE,
  floristic_group = "all",
  complete_floristic = FALSE,
  geo_type = "Island",
  suit_geo = TRUE,
  namesmatched = TRUE,
  list_set_only = FALSE,
  GIFT_version = gift_version
)

list_meta <- checklists$lists
species_rows <- checklists$checklists
if (is.null(list_meta) || nrow(list_meta) == 0) {
  stop("GIFT returned no island checklist metadata")
}
if (is.null(species_rows) || nrow(species_rows) == 0) {
  stop("GIFT returned no island checklist species rows")
}

write.csv(
  list_meta,
  file.path(output_dir, "gift_island_list_metadata.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  species_rows,
  file.path(output_dir, "gift_island_checklist_rows.csv"),
  row.names = FALSE,
  na = ""
)

# Only exact Island polygons are materialized for primary matching. Group and part
# metadata remain in gift_island_list_metadata.csv and therefore remain auditable.
if (!("entity_class" %in% colnames(list_meta))) {
  stop("GIFT checklist metadata lacks entity_class")
}
exact_islands <- list_meta[list_meta$entity_class == "Island", , drop = FALSE]
entity_ids <- unique(exact_islands$entity_ID)
entity_ids <- entity_ids[!is.na(entity_ids)]
if (length(entity_ids) == 0) {
  stop("GIFT returned no exact Island entity_ID values")
}

shapes <- GIFT_shapes(entity_ID = entity_ids, GIFT_version = gift_version)
if (is.null(shapes) || nrow(shapes) == 0) {
  stop("GIFT returned no exact-island geometries")
}

st_write(
  shapes,
  dsn = file.path(output_dir, "gift_exact_island_shapes.gpkg"),
  layer = "gift_islands",
  delete_dsn = TRUE,
  quiet = TRUE
)

versions <- tryCatch(GIFT_versions(), error = function(e) NULL)
if (!is.null(versions)) {
  write.csv(
    versions,
    file.path(output_dir, "gift_versions.csv"),
    row.names = FALSE,
    na = ""
  )
}

manifest <- list(
  requested_version = gift_version,
  n_list_metadata_rows = nrow(list_meta),
  n_checklist_rows = nrow(species_rows),
  n_exact_island_entities = length(unique(entity_ids)),
  n_exact_island_shapes = nrow(shapes),
  policy = paste(
    "Raw GIFT status fields and references are exported without status inference;",
    "only exact Island geometries are materialized for primary GSHHG matching."
  )
)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required to write the manifest")
}
jsonlite::write_json(
  manifest,
  file.path(output_dir, "gift_floristic_status_export_manifest.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

message(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE))
