#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
out_dir <- "data/v2/staging/traits/bulk/bien"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

shard_index <- as.integer(Sys.getenv("BIEN_SHARD_INDEX", "0"))
shard_count <- as.integer(Sys.getenv("BIEN_SHARD_COUNT", "1"))
if (is.na(shard_index) || is.na(shard_count) || shard_count < 1 || shard_index < 0 || shard_index >= shard_count) {
  stop("Invalid BIEN shard configuration")
}
suffix <- if (shard_count > 1) sprintf("_shard_%02d_of_%02d", shard_index + 1, shard_count) else ""

if (!requireNamespace("BIEN", quietly = TRUE)) {
  install.packages("BIEN", dependencies = NA)
}

exports <- getNamespaceExports("BIEN")
required <- c("BIEN_trait_list", "BIEN_trait_trait")
missing <- setdiff(required, exports)
if (length(missing)) stop("Required BIEN API functions missing: ", paste(missing, collapse = ", "))

trait_list <- BIEN::BIEN_trait_list()
write.csv(trait_list, file.path(out_dir, paste0("trait_list", suffix, ".csv")), row.names = FALSE, na = "")

# BIEN 1.2.8 documents BIEN_trait_list() and BIEN_trait_trait(); identify the
# returned trait-name column from the actual response rather than assuming a schema.
trait_col <- intersect(c("trait_name", "trait", "Trait", "traitName"), names(trait_list))[1]
if (is.na(trait_col)) {
  char_cols <- names(trait_list)[vapply(trait_list, is.character, logical(1))]
  if (!length(char_cols)) stop("Could not identify trait-name column in BIEN_trait_list()")
  trait_col <- char_cols[1]
}
traits_all <- unique(trimws(as.character(trait_list[[trait_col]])))
traits_all <- traits_all[nzchar(traits_all) & !is.na(traits_all)]
traits <- traits_all[(seq_along(traits_all) - 1L) %% shard_count == shard_index]
if (!length(traits)) stop("BIEN shard contains no traits")
message(sprintf("BIEN shard %d/%d: %d of %d traits", shard_index + 1, shard_count, length(traits), length(traits_all)))

master <- read.csv("data/v2/staging/gbif/collected/island_taxa.csv", stringsAsFactors = FALSE, check.names = FALSE)
master_col <- "accepted_species"
if (!master_col %in% names(master)) stop("accepted_species missing from island master")
master_names <- unique(trimws(as.character(master[[master_col]])))
master_names <- master_names[nzchar(master_names) & !is.na(master_names)]

pick_name_col <- function(x) {
  preferred <- c("scrubbed_species_binomial", "accepted_species", "species", "species_binomial", "taxon_name", "scientific_name")
  hit <- preferred[preferred %in% names(x)]
  if (length(hit)) return(hit[1])
  candidates <- names(x)[grepl("species|taxon|scientific", names(x), ignore.case = TRUE)]
  if (length(candidates)) return(candidates[1])
  NA_character_
}

summary_rows <- list()
matched_chunks <- list()
for (i in seq_along(traits)) {
  tr <- traits[[i]]
  message(sprintf("BIEN shard %d/%d trait %d/%d: %s", shard_index + 1, shard_count, i, length(traits), tr))
  dat <- tryCatch(BIEN::BIEN_trait_trait(trait = tr), error = function(e) e)
  if (inherits(dat, "error")) {
    summary_rows[[length(summary_rows)+1]] <- data.frame(trait = tr, records = 0, candidate_taxa = 0, exact_master_matches = 0, error = conditionMessage(dat))
    next
  }
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  name_col <- pick_name_col(dat)
  if (is.na(name_col)) {
    summary_rows[[length(summary_rows)+1]] <- data.frame(trait = tr, records = nrow(dat), candidate_taxa = 0, exact_master_matches = 0, error = paste("No taxon column; columns:", paste(names(dat), collapse = "|")))
    next
  }
  taxa <- trimws(as.character(dat[[name_col]]))
  keep <- !is.na(taxa) & taxa %in% master_names
  matched <- dat[keep, , drop = FALSE]
  if (nrow(matched)) {
    matched$bien_query_trait <- tr
    matched_chunks[[length(matched_chunks)+1]] <- matched
  }
  summary_rows[[length(summary_rows)+1]] <- data.frame(
    trait = tr,
    records = nrow(dat),
    candidate_taxa = length(unique(taxa[nzchar(taxa) & !is.na(taxa)])),
    exact_master_matches = length(unique(taxa[keep])),
    error = ""
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(out_dir, paste0("trait_coverage", suffix, ".csv")), row.names = FALSE, na = "")

if (length(matched_chunks)) {
  all_names <- unique(unlist(lapply(matched_chunks, names)))
  matched_chunks <- lapply(matched_chunks, function(x) {
    miss <- setdiff(all_names, names(x)); for (m in miss) x[[m]] <- NA
    x[all_names]
  })
  matched_all <- do.call(rbind, matched_chunks)
  gz <- gzfile(file.path(out_dir, paste0("master_matched_trait_records", suffix, ".csv.gz")), "wt")
  write.csv(matched_all, gz, row.names = FALSE, na = "")
  close(gz)
  matched_unique_taxa <- length(unique(trimws(as.character(matched_all[[pick_name_col(matched_all)]]))))
  matched_records <- nrow(matched_all)
} else {
  matched_unique_taxa <- 0
  matched_records <- 0
}

report <- list(
  source = "BIEN",
  access = "official CRAN BIEN package querying BIEN database",
  package_version = as.character(utils::packageVersion("BIEN")),
  shard_index = shard_index,
  shard_count = shard_count,
  traits_available = length(traits_all),
  traits_queried = length(traits),
  successful_trait_queries = sum(summary_df$error == ""),
  failed_trait_queries = sum(summary_df$error != ""),
  total_records_returned = sum(summary_df$records),
  unique_exact_master_matches = matched_unique_taxa,
  matched_records_saved = matched_records
)
json_value <- function(v) {
  if (is.character(v)) return(paste0('"', gsub('"', '\\\\"', v), '"'))
  if (is.logical(v)) return(tolower(as.character(v)))
  as.character(v)
}
json <- paste0("{\n", paste(sprintf('  "%s": %s', names(report), vapply(report, json_value, character(1))), collapse = ",\n"), "\n}\n")
writeLines(json, file.path(out_dir, paste0("coverage_report", suffix, ".json")))
print(report)
