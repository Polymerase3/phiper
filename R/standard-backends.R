#' @title Read and register "long" phiper data into a DuckDB-backed database
#'
#' @description This internal function ingests one or more data files (Parquet or CSV)
#' specified by `cfg$data_long_path` into a single DuckDB view named
#' `data_long`, applying user-provided column mappings (`colmap`) to
#' rename each source column to the standard PHIPER names. The resulting
#' `phip_data` object contains a lazy DuckDB table that can be queried
#' with dplyr without loading the full dataset into R until explicitly
#' collected.
#'
#' @param cfg Named list, must contain element `data_long_path` pointing
#'   to either a single file or a directory of files. Supported file
#'   extensions are `.parquet`, `.parq`, `.pq`, and `.csv`.
#' @param colmap Named character list mapping **standard** PHIPER column
#'   names (e.g. `"sample_id"`, `"peptide_id"`, ...) to the **actual**
#'   column names found in the source files.
#'
#' @return A `phip_data` S3/S4 object (depending on your package
#'   implementation) whose `data_long` slot is a `dplyr::tbl_dbi`
#'   representing the union of all source files. Calculations against
#'   `data_long` remain lazy until `collect()` is called.
#'
#' @details
#' - If `cfg$data_long_path` is a **directory**, all matching files
#'   within it are UNION ALL'ed.
#' - Parquet files are read via `parquet_scan()`, CSV via
#'   `read_csv_auto()`.
#' - Column renaming is performed in SQL with `AS`, so no R-level
#'   `rename()` calls are needed.
#' - A DuckDB **VIEW** called `data_long` is created (dropped if it
#'   existed previously) for downstream queries.
#'
#' @keywords internal
.standard_read_duckdb_backend <- function(cfg, colmap) {
  rlang::check_installed(c("duckdb", "DBI", "dbplyr"),
    reason = "duckdb backend"
  )

  ## 0. open a DuckDB connection -------------------------------------------
  ## keep it persistent only for this R session, but ON DISK (so it can spill)
  cache_dir <- withr::local_tempdir("phiper_cache", .local_envir = globalenv())
  duckdb_file <- file.path(cache_dir, "phip_cache.duckdb")
  tmp_dir <- file.path(cache_dir, "tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  ## set threads + enable spilling via memory_limit & temp_directory
  con <- DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = duckdb_file,
    config = list(
      threads        = as.character(cfg$n_cores),
      memory_limit   = cfg$duckdb_mem %||% "32GB",
      temp_directory = tmp_dir
    )
  )

  ## -- 1. Determine files to load ---------------------------------------------
  info <- file.info(cfg$data_long_path)

  files <- if (isTRUE(info$isdir)) {
    list.files(cfg$data_long_path,
      pattern = "\\.(parquet|parq|pq|csv)$", full.names = TRUE
    )
  } else {
    cfg$data_long_path
  }

  stopifnot(length(files) > 0)

  ## -- 2. Helper: load a single file into a (disk-backed) TABLE ---------------
  load_file <- function(path, tbl_name) {
    path_q <- gsub("'", "''", path)
    if (grepl("\\.csv$", path, ignore.case = TRUE)) {
      DBI::dbExecute(
        con,
        sprintf(
          "CREATE TABLE %s AS
             SELECT * FROM read_csv_auto('%s', HEADER=TRUE);",
          tbl_name, path_q
        )
      )
    } else {
      DBI::dbExecute(
        con,
        sprintf(
          "CREATE TABLE %s AS
             SELECT * FROM parquet_scan('%s');",
          tbl_name, path_q
        )
      )
    }
  }

  tbl_names <- sprintf("f_%d", seq_along(files))
  Map(load_file, files, tbl_names)

  ## -- 3. UNION ALL --> raw_combined ------------------------------------------
  union_sql <- paste(
    sprintf("SELECT * FROM %s", tbl_names),
    collapse = " UNION ALL "
  )

  # 2. Execute either a materialized table or a view based on cfg$table_type
  if (cfg$materialise_table) {
    # Create or replace a real table (snapshot persisted to disk)
    DBI::dbExecute(
      con,
      sprintf(
        "CREATE OR REPLACE TABLE raw_combined AS %s;",
        union_sql
      )
    )
  } else {
    # Create or replace a view (virtual table; query re-runs on each SELECT)
    DBI::dbExecute(
      con,
      sprintf(
        "CREATE OR REPLACE VIEW raw_combined AS %s;",
        union_sql
      )
    )
  }

  ## -- 4. Standardise column names inside DuckDB ------------------------------
  .rename_to_standard_inplace(
    con = con,
    tbl = "raw_combined",
    colname_map = colmap
  )

  # return "BASE TABLE", "VIEW", or NA if it does not exist
  obj_type <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT table_type
       FROM information_schema.tables
      WHERE table_schema = current_schema
        AND table_name   = %s",
      DBI::dbQuoteString(con, "raw_combined")
    )
  )$table_type[1]

  if (!is.na(obj_type) && toupper(obj_type) != "VIEW") {
    DBI::dbExecute(con, "ANALYZE raw_combined;")
  } else {
    message("Skipping ANALYZE - raw_combined is a view.")
  }

  ## -- 5. Clean up intermediate tables ----------------------------------------
  ## we actually want to clean up the tables only, when the main table is
  ## materialised. When it's not materialised, the view will take a look on the
  ## original files/tables, so we can not delete them --> we have to have
  ## something to reference to
  if (!is.na(obj_type) && toupper(obj_type) != "VIEW") {
    invisible(lapply(
      tbl_names,
      function(tn) {
        DBI::dbExecute(
          con,
          sprintf("DROP TABLE %s;", tn)
        )
      }
    ))
  }

  invisible(con) # return the open connection
}
