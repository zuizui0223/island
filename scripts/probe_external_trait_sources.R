#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

out_dir <- "data/v2/staging/traits/bulk/source_probe"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

packages <- c("BIEN", "BROT", "GIFT")
summary_rows <- list()

safe_json <- function(x, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite", quiet = TRUE)
  }
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}

for (pkg in packages) {
  installed <- requireNamespace(pkg, quietly = TRUE)
  install_error <- NULL
  if (!installed) {
    tryCatch(
      install.packages(pkg, quiet = TRUE),
      error = function(e) install_error <<- conditionMessage(e)
    )
    installed <- requireNamespace(pkg, quietly = TRUE)
  }

  info <- list(
    package = pkg,
    installed = installed,
    install_error = install_error,
    version = if (installed) as.character(utils::packageVersion(pkg)) else NULL,
    exports = character(),
    trait_exports = character(),
    datasets = character(),
    extracted_tables = character()
  )

  if (installed) {
    exports <- getNamespaceExports(pkg)
    info$exports <- exports
    info$trait_exports <- exports[grepl("trait|growth|flower|pollin|reproduct|phenolog", exports, ignore.case = TRUE)]

    datasets <- tryCatch(utils::data(package = pkg)$results, error = function(e) NULL)
    if (!is.null(datasets) && nrow(datasets) > 0) {
      dataset_names <- unique(as.character(datasets[, "Item"]))
      dataset_names <- sub("\\s+\\(.*$", "", dataset_names)
      info$datasets <- dataset_names

      for (dataset_name in dataset_names) {
        env <- new.env(parent = emptyenv())
        loaded <- tryCatch({
          utils::data(list = dataset_name, package = pkg, envir = env)
          TRUE
        }, error = function(e) FALSE)
        if (!loaded || !exists(dataset_name, envir = env, inherits = FALSE)) next
        obj <- get(dataset_name, envir = env, inherits = FALSE)
        if (!is.data.frame(obj)) next
        path <- file.path(out_dir, paste0(tolower(pkg), "__", dataset_name, ".csv"))
        utils::write.csv(obj, path, row.names = FALSE, na = "")
        info$extracted_tables <- c(info$extracted_tables, path)
      }
    }

    trait_export_rows <- lapply(info$trait_exports, function(fn) {
      formals_text <- tryCatch(
        paste(names(formals(getExportedValue(pkg, fn))), collapse = ","),
        error = function(e) ""
      )
      data.frame(package = pkg, function_name = fn, arguments = formals_text, stringsAsFactors = FALSE)
    })
    if (length(trait_export_rows) > 0) {
      utils::write.csv(
        do.call(rbind, trait_export_rows),
        file.path(out_dir, paste0(tolower(pkg), "_trait_exports.csv")),
        row.names = FALSE
      )
    }
  }

  safe_json(info, file.path(out_dir, paste0(tolower(pkg), "_probe.json")))
  summary_rows[[pkg]] <- data.frame(
    source = pkg,
    installed = installed,
    version = if (installed) info$version else "",
    trait_exports = length(info$trait_exports),
    datasets = length(info$datasets),
    extracted_tables = length(info$extracted_tables),
    install_error = if (is.null(install_error)) "" else install_error,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
utils::write.csv(summary_df, file.path(out_dir, "source_probe_summary.csv"), row.names = FALSE)
print(summary_df)
