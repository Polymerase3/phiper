# ---------------------------------------------------------------------------
#  Internal helper: build counts_final in an in-memory DuckDB connection
# ---------------------------------------------------------------------------
#  * No rows are collected to R.
#  * The caller receives a ready-to-use DBI connection with
#      counts_final   – long table (peptide × sample) + metadata
#      samples2/3     – cleaned sample metadata
#  * Caller is responsible for DBI::dbDisconnect().
# ---------------------------------------------------------------------------
.read_duckdb_backend <- function(cfg,
                                 meta) {
  rlang::check_installed(c("duckdb", "DBI", "dbplyr"),
    reason = "duckdb backend"
  )

  cache_dir <- withr::local_tempdir("phiper_cache") # optional name-prefix
  duckdb_file <- file.path(cache_dir, "phip_cache.duckdb")

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = duckdb_file)

  q <- function(x) DBI::dbQuoteString(con, x)
  qi <- function(x) DBI::dbQuoteIdentifier(con, x)

  ## file reader based on the file extension (.csv/.parquet)
  ## returns a character string (SQL query)
  duckdb_load_sql <- function(path, table_name, header = TRUE) {
    if (is.null(path)) {
      return(NULL)
    } # fallback to NULL

    ext <- paste0(
      tolower(strsplit(basename(path), "\\.", fixed = FALSE)[[1]][-1]),
      collapse = "."
    )

    is_parquet <- grepl("^parq(uet)?(\\.|$)|^pq(\\.|$)", ext)
    reader_fun <- if (is_parquet) "parquet_scan" else "read_csv_auto"
    qpath <- sprintf("'%s'", gsub("'", "''", path))

    if (reader_fun == "parquet_scan") {
      sprintf(
        "CREATE TEMP TABLE %s AS SELECT * FROM parquet_scan(%s);",
        table_name, qpath
      )
    } else {
      hdr_flag <- if (isTRUE(header)) "HEADER=TRUE" else "HEADER=FALSE"
      sprintf(
        "CREATE TEMP TABLE %s AS SELECT * FROM read_csv_auto(%s,%s);",
        table_name, qpath, hdr_flag
      )
    }
  }

  # -----------------------------------------------------------------------
  # 1. OPTIONAL MATRICES ---------------------------------------------------
  load_and_unpivot <- function(file, wide_name, long_name, value_col) {
    if (is.null(file)) {
      return(NULL)
    }
    DBI::dbExecute(con, duckdb_load_sql(file, wide_name))

    first_col <- DBI::dbGetQuery(
      con,
      sprintf(
        "SELECT column_name
                 FROM duckdb_columns()
                WHERE table_name = '%s'
             ORDER BY column_index LIMIT 1;",
        wide_name
      )
    )$column_name

    samp_cols <- DBI::dbGetQuery(
      con,
      sprintf(
        "SELECT column_name
                 FROM duckdb_columns()
                WHERE table_name = '%s' AND column_name <> '%s'
             ORDER BY column_index;",
        wide_name, first_col
      )
    )$column_name

    DBI::dbExecute(
      con,
      sprintf(
        "CREATE TEMP TABLE %s AS
           SELECT %s AS peptide_id, sample_id, %s
             FROM %s
             UNPIVOT (%s FOR sample_id IN (%s));",
        long_name, qi(first_col), value_col, wide_name, value_col,
        paste(qi(samp_cols), collapse = ", ")
      )
    )
    long_name
  }

  ## create a list of matrices (or actually names, as the tables are registered
  ## in the temporary duckdb dir)
  long_tbls <- list(
    tbl_counts = load_and_unpivot(
      cfg$exist_file, "counts_wide",
      "counts_long", "exist"
    ),
    tbl_fc = load_and_unpivot(
      cfg$fold_change_file, "fold_change_wide",
      "fold_change_long", "fold_change"
    ),
    tbl_inp = load_and_unpivot(
      cfg$input_file, "input_wide",
      "input_long", "input_count"
    ),
    tbl_hit = load_and_unpivot(
      cfg$hit_file, "hit_wide",
      "hit_long", "hit_count"
    )
  )

  ## remove the NULLs
  long_tbls <- Filter(Negate(is.null), long_tbls)

  stopifnot(length(long_tbls) > 0) # at least one matrix must be present

  base_tbl <- long_tbls[[1]]

  ## ensure `counts` exists and is the base -------------------------------
  DBI::dbExecute(
    con,
    sprintf("CREATE TEMP TABLE final_long AS SELECT * FROM %s;", base_tbl)
  )

  ## join additional matrices --------------------------------------------
  for (tbl in long_tbls[-1]) {
    value_col <- DBI::dbGetQuery(
      con,
      sprintf(
        "SELECT column_name
                 FROM duckdb_columns()
                WHERE table_name = '%s'
                  AND column_name NOT IN ('peptide_id','sample_id')
             LIMIT 1;",
        tbl
      )
    )$column_name

    DBI::dbExecute(
      con,
      sprintf(
        "CREATE OR REPLACE TEMP TABLE final_long AS
           SELECT f.*, t.%s
             FROM final_long f
             LEFT JOIN %s t USING (peptide_id, sample_id);",
        qi(value_col), tbl
      )
    )

    ## clean after merging --> removing unnecessary large tables, everything is
    ## in final_long now
    DBI::dbExecute(con, sprintf("DROP TABLE %s;", tbl))
  }

  # -----------------------------------------------------------------------
  # 2. samples metadata ----------------------------------------------------
  duckdb::duckdb_register(con, "samples_raw", meta$samples)

  first_col_samples <- DBI::dbGetQuery(
    con,
    "SELECT column_name
       FROM duckdb_columns()
      WHERE table_name = 'samples_raw'
   ORDER BY column_index LIMIT 1;"
  )$column_name

  DBI::dbExecute(
    con,
    sprintf(
      "CREATE TEMP TABLE samples2 AS
         SELECT %s AS sample_id, * EXCLUDE (%s)
           FROM samples_raw;",
      qi(first_col_samples), qi(first_col_samples)
    )
  )

  DBI::dbExecute(
    con,
    "CREATE OR REPLACE TEMP TABLE final_long AS
       SELECT f.*, s.* EXCLUDE sample_id
         FROM final_long f
         LEFT JOIN samples2 s USING (sample_id);"
  )

  # -----------------------------------------------------------------------

  if (!is.null(meta$comparisons)) {
    duckdb::duckdb_register(con, "comparisons", meta$comparisons)
  }

  invisible(con)
}

#' @title Build a **memory-backed** `phip_data` object (no DuckDB / Arrow)
#'
#' @description Reads up to four wide matrices (`exist`, `fold_change`, `input`, `hit`)
#' directly into R, pivots each to long format with **data.table::melt()**
#' (faster than `stats::reshape()`), left-joins the available matrices on
#' `(sample_id, peptide_id)`, attaches the cleaned sample metadata produced by
#' `legacy_prepare_metadata()`, and returns an in-memory `phip_data`.
#' Heavy validation of matrix content is deliberately **omitted**—that should be
#' performed by the `phip_data` class validator.
#'
#' @param cfg       Named list from `legacy_resolve_paths()`; must contain the
#'                  four `*_file` elements.
#' @param meta      List returned by `legacy_prepare_metadata()` containing
#'                  `samples` (data-frame) and `comparisons` (data-frame or
#'                  `NULL`).
#'
#' @return A `phip_data` object whose `data_long` slot is a tibble in RAM.
#' @keywords internal
.read_memory_backend <- function(cfg, meta) {
  # ---------- tiny helper ----------------------------------------------------
  ## safe fallback when file_name is NULL
  .read_or_null <- function(path) if (is.null(path)) NULL else .auto_read(path)

  ## safe fallback when NULL
  .process_matrix <- function(mat, value_col) {
    if (is.null(mat)) {
      return(NULL)
    }

    data.table::setDT(mat) # in-place, no copy
    id_col <- names(mat)[1]
    valcols <- names(mat)[-1]

    long <- data.table::melt(
      mat,
      id.vars = id_col,
      measure.vars = valcols,
      variable.name = "sample_id",
      value.name = value_col,
      variable.factor = FALSE
    )

    data.table::setnames(long, id_col, "peptide_id")
    long[]
  }

  # ---------- 1.  read & pivot matrices --------------------------------------
  mats <- list(
    counts      = .process_matrix(.read_or_null(cfg$exist_file), "exist"),
    fold_change = .process_matrix(.read_or_null(cfg$fold_change_file), "fold_change"),
    input       = .process_matrix(.read_or_null(cfg$input_file), "input_count"),
    hit         = .process_matrix(.read_or_null(cfg$hit_file), "hit_count")
  )

  ## filter out NULLs
  mats <- Filter(Negate(is.null), mats)

  # ---------- 2.  progressively join value columns ---------------------------
  out <- mats[[1L]]
  for (tbl in mats[-1]) {
    vcol <- setdiff(names(tbl), c("sample_id", "peptide_id"))
    out <- dplyr::left_join(out,
      tbl[, c("sample_id", "peptide_id", vcol)],
      by = c("sample_id", "peptide_id")
    )
  }

  ## clear memory
  rm(mats)
  invisible(gc())

  # ---------- 3.  attach sample metadata -------------------------------------
  out <- merge(out, meta$samples, by = "sample_id")

  # ---------- 4.  finish ------------------------------------------------------
  new_phip_data(
    data_long = tibble::as_tibble(out),
    comparisons = if (is.null(meta$comparisons)) NULL else tibble::as_tibble(meta$comparisons),
    backend = "memory",
    peptide_library = cfg$peptide_library
  )
}
