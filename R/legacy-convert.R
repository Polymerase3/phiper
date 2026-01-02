#' @title Convert legacy Carlos-style input to a modern **phip_data** object
#'
#' @description
#' `phip_convert_legacy()` ingests the original three-file PhIP-Seq input
#' (binary *exist* matrix, *samples* metadata, optional *timepoints* map) plus
#' an optional *comparisons* file.
#' Paths can be supplied directly or via a single YAML config; explicit
#' arguments always override the YAML.  The function normalises the chosen
#' DuckDB storage, validates every file, and returns a ready-to-use
#' `phip_data` object.
#'
#' @details
#' **Validation rules**
#' *1 – exist CSV*
#' * First column **must** be `peptide_id` and unique.
#' * Remaining columns are `sample_id`s found in the samples file.
#' * Values allowed: `0`, `1`, or `NA` – anything else aborts.
#'
#' *2 – samples CSV*
#' * First column **must** be `sample_id`, unique.
#' * Extra columns kept only if listed in `extra_cols`.
#' * If dummy group columns are referenced by `comparisons_file`, each row’s
#'   dummy sum must equal **1**.
#'
#' *3 – timepoints CSV* (optional, longitudinal)
#' * First column **must** be `ind_id` (subject).
#' * Other columns are time-point names; cells are `sample_id` or `NA`.
#' * Column names must match `timepoint` values in the data; every `sample_id`
#'   appears at most once.
#'
#' *4 – comparisons CSV* (optional)
#' * Columns required: `comparison`, `group1`, `group2`, `variable`.
#' * Labels in `group1`/`group2` must exist in the derived `group` column or the
#'   `timepoint` column (for longitudinal data).
#'
#' Files failing any rule trigger an informative `.chk_cond()` error.
#'
#' @param exist_file       Path to the **exist** CSV (peptide x sample binary
#'   matrix). *Required unless given in `config_yaml`.*
#' @param fold_change_file       Path to the **fold_change** CSV (peptide x
#'   sample numeric matrix). *Required unless given in `config_yaml`.*
#' @param input_file,hit_file    Paths to the **raw_counts** CSV (peptide x
#'   sample integer matrix). *Required unless given in `config_yaml`.*
#' @param samples_file     Path to the **samples** CSV (sample metadata).
#'   *Required unless given in `config_yaml`.*
#' @param timepoints_file  Path to the **timepoints** CSV (subject <-> sample
#'   mapping). Optional for cross-sectional data.
#' @param extra_cols       Character vector of extra metadata columns to retain.
#' @param comparisons_file Path to a **comparisons** CSV. Optional.
#' @param output_dir       *Deprecated.* Ignored with a warning.
#' @param peptide_library logical, defining if the `peptide_library` is to be
#'    downloaded from the official `phiper` GitHub
#' @param config_yaml      Optional YAML file containing any of the above
#'   parameters (see example).
#' @param n_cores Integer >= 1. Number of CPU threads DuckDB may use while
#'   reading and writing files.
#'
#' @param materialise_table Logical. If `FALSE` the result is registered as a
#'   **view**; if `TRUE` the table is fully **materialised** and stored on disk,
#'   trading higher load time and storage for faster repeated queries.
#' @return A validated `phip_data` object whose `data_long` slot is backed by a
#'   DuckDB connection.
#'
#' @examples
#' \dontrun{
#' ## 1. Direct-path usage
#' pd <- phip_convert_legacy(
#'   exist_file = "legacy/exist.csv",
#'   samples_file = "legacy/samples.csv",
#'   timepoints_file = "legacy/timepoints.csv",
#'   comparisons_file = "legacy/comparisons.csv"
#' )
#'
#' ## 2. YAML-driven usage (explicit args override YAML)
#' # --- config/legacy_config.yaml ---
#' # exist_file:       data/exist.csv
#' # samples_file:     meta/samples.csv
#' # timepoints_file:  meta/timepoints.csv
#' # comparisons_file: meta/comparisons.csv
#' # extra_cols: [sex, age]
#' # -------------------------------
#'
#' pd <- phip_convert_legacy(
#'   config_yaml = "config/legacy_config.yaml"
#' )
#' }
#'
#' @export

phip_convert_legacy <- function(
  exist_file = NULL,
  fold_change_file = NULL,
  samples_file = NULL,
  input_file = NULL,
  hit_file = NULL,
  timepoints_file = NULL,
  extra_cols = NULL,
  comparisons_file = NULL,
  output_dir = NULL, # hard deprecation
  peptide_library = TRUE,
  n_cores = 8,
  materialise_table = TRUE,
  config_yaml = NULL
) {
  #' @importFrom rlang .data

  # ------------------------------------------------------------------
  # 1. arg checks
  # ------------------------------------------------------------------
  chk::chk_numeric(n_cores)
  chk::chk_flag(materialise_table)

  # ------------------------------------------------------------------
  # 2. resolving the paths to absolute
  # ------------------------------------------------------------------
  cfg <- .resolve_paths(
    exist_file = exist_file,
    fold_change_file = fold_change_file,
    samples_file = samples_file,
    input_file = input_file,
    hit_file = hit_file,
    timepoints_file = timepoints_file,
    extra_cols = extra_cols,
    comparisons_file = comparisons_file,
    output_dir = output_dir,
    backend = "duckdb",
    peptide_library = peptide_library,
    config_yaml = config_yaml,
    n_cores = n_cores,
    materialise_table = materialise_table
  )

  # ------------------------------------------------------------------
  # 3. prepare the metadata
  # ------------------------------------------------------------------
  meta_list <- .legacy_prepare_metadata(
    cfg$samples_file,
    cfg$comparisons_file,
    cfg$timepoints_file,
    cfg$extra_cols
  )

  # ------------------------------------------------------------------
  # 4. create the phip_data object (DuckDB)
  # ------------------------------------------------------------------
  con <- .read_duckdb_backend(cfg, meta_list)

  ## duckdb-specific code
  long <- dplyr::tbl(con, "final_long")

  comps <- if (DBI::dbExistsTable(con, "comparisons")) {
    dplyr::tbl(con, "comparisons") |> dplyr::collect()
  }

  # returning the phip_data object
  new_phip_data(
    data_long = long,
    comparisons = comps,
    peptide_library = cfg$peptide_library,
    backend = "duckdb",
    meta = list(con = con)
  )
}

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
