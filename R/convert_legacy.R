#' @title Convert legacy Carlos-style input to a modern **phip_data** object
#'
#' @description
#' `phip_convert_legacy()` ingests the original three-file PhIP-Seq input
#' (binary *exist* matrix, *samples* metadata, optional *timepoints* map) plus
#' an optional *comparisons* file.
#' Paths can be supplied directly or via a single YAML config; explicit
#' arguments always override the YAML.  The function normalises the chosen
#' storage `backend`, validates every file, and returns a ready-to-use
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
#' @param exist_file       Path to the **exist** CSV (peptide × sample binary
#'   matrix). *Required unless given in `config_yaml`.*
#' @param samples_file     Path to the **samples** CSV (sample metadata).
#'   *Required unless given in `config_yaml`.*
#' @param timepoints_file  Path to the **timepoints** CSV (subject <-> sample
#'   mapping). Optional for cross-sectional data.
#' @param extra_cols       Character vector of extra metadata columns to retain.
#' @param comparisons_file Path to a **comparisons** CSV. Optional.
#' @param output_dir       *Deprecated.* Ignored with a warning.
#' @param backend          Storage backend: `"arrow"`, `"duckdb"`, or
#'   `"memory"`. Defaults to `"duckdb"`.
#' @param config_yaml      Optional YAML file containing any of the
#'   above parameters (see example).
#'
#' @return A validated `phip_data` object whose `data_long` slot is backed by
#'   a tibble (memory), a DuckDB connection, or an Arrow dataset, depending on
#'   `backend`.
#'
#' @examples
#' \dontrun{
#' ## 1. Direct-path usage
#' pd <- phip_convert_legacy(
#'   exist_file = "legacy/exist.csv",
#'   samples_file = "legacy/samples.csv",
#'   timepoints_file = "legacy/timepoints.csv",
#'   comparisons_file = "legacy/comparisons.csv",
#'   backend = "duckdb"
#' )
#'
#' ## 2. YAML-driven usage (explicit args override YAML)
#' # --- config/legacy_config.yaml ---
#' # exist_file:       data/exist.csv
#' # samples_file:     meta/samples.csv
#' # timepoints_file:  meta/timepoints.csv
#' # comparisons_file: meta/comparisons.csv
#' # extra_cols: [sex, age]
#' # backend: duckdb
#' # -------------------------------
#'
#' pd <- phip_convert_legacy(
#'   config_yaml = "config/legacy_config.yaml",
#'   backend     = "arrow" # overrides YAML backend
#' )
#' }
#'
#' @export



phip_convert_legacy <- function(
    exist_file = NULL,
    samples_file = NULL,
    timepoints_file = NULL,
    extra_cols = NULL,
    comparisons_file = NULL,
    output_dir = NULL, # hard deprecation
    backend = NULL,
    config_yaml = NULL) {
  #' @importFrom rlang .data

  # ------------------------------------------------------------------
  # db-backend: default to "duckdb" if user supplies nothing
  # ------------------------------------------------------------------
  backend_choices <- c("arrow", "duckdb", "memory")

  backend <- if (is.null(backend)) {
    "duckdb" # implicit default
  } else {
    match.arg(backend, choices = backend_choices)
  }

  # ------------------------------------------------------------------
  # figure out the base directory
  # ------------------------------------------------------------------
  base_dir <- if (!is.null(config_yaml)) {
    fs::path_dir(fs::path_abs(config_yaml)) # dir that holds config.yaml
  } else {
    ## check if exist_file and samples file provided (bare minimum)
    .chk_cond(
      is.null(exist_file) || is.null(samples_file),
      "The `exist_file` and `samples_file` arguments can not be null,
              when no `config_yaml` provided."
    )

    fs::path_dir(fs::path_abs(exist_file)) # dir that holds exist_file
  }

  # ---------------------------------------------------------------------------
  # pick up values from YAML, if provided
  # ---------------------------------------------------------------------------
  if (!is.null(config_yaml)) {
    # validate the path to config.yaml with custom function
    .chk_path(
      config_yaml,
      "config_yaml",
      c("yml", "yaml")
    )

    # check dependencies
    rlang::check_installed("yaml")

    # read yaml
    cfg <- yaml::read_yaml(config_yaml)
  } else {
    # safe placeholder
    cfg <- NULL
  }

  # helper to fetch the paths and abort when path not found
  fetch <- function(arg,
                    key,
                    validator = NULL,
                    optional = FALSE,
                    ...) {
    # if both NULL, then %||% returns NULL
    val <- cfg[[key]] %||% arg # YAML first, then explicit arg

    # abort if val is NULL (required arg not found)
    # actually a fallback, works only for comparisons file, the other two
    # (samples and exist) are checked earlier
    if (is.null(val) && !optional) {
      chk::abort_chk(sprintf("Missing required argument `%s`: not defined in
                               the YAML config or as an explicit path. Aborting
                               data load.", key))
    }

    # if the validator if .chk_path, expand the path to absolute
    if (!is.null(val) && identical(validator, .chk_path)) {
      if (!fs::is_absolute_path(val)) {
        val <- fs::path_abs(val, start = base_dir)
      }
    }

    # run the validator only if supplied and value is not NULL
    if (!is.null(val) && is.function(validator)) {
      validator(val, key, ...)
    }

    return(val)
  }

  # ---------------------------------------------------------------------------
  # load the arguments from .yaml or use direct path definitions and
  # quick validate required inputs
  # ---------------------------------------------------------------------------
  exist_file <- fetch(exist_file, "exist_file", .chk_path, extension = c("csv", "parquet"))
  samples_file <- fetch(samples_file, "samples_file",
    .chk_path,
    extension = "csv"
  )
  timepoints_file <- fetch(timepoints_file, "timepoints_file", .chk_path,
    optional = TRUE, extension = "csv"
  )
  extra_cols <- fetch(extra_cols, "extra_cols", optional = TRUE)
  comparisons_file <- fetch(comparisons_file, "comparisons_file", .chk_path,
    extension = "csv", optional = TRUE
  )
  output_dir <- fetch(output_dir, "output_dir", optional = TRUE)

  backend <- fetch(backend, "backend", chk::chk_string)

  # validate the tables itself --> the helper functions are defined at the
  # bottom of the script

  # the legacy workflow has an important limitation: the comparisons specified
  # in the comparisons_file are suitable only for the longitudinal data -->
  # this means that both timepoints_file and comparisons_file have to be
  # provided at the same time --> log error if only one provided

  ## ==> this is actually wrong, solution: see long comment around line 180 with
  ## XXX at the beginning

  # Warn about deprecation
  .chk_cond(!is.null(output_dir),
    error = FALSE,
    "'output_dir' is deprecated."
  )

  ##############################################################################
  #  L E G A C Y   B R I D G E
  ##############################################################################
  ## ---- 1.  samples & contrasts (both small) ---------------------------------
  samples <- .auto_read_csv(samples_file)

  if (!is.null(comparisons_file)) {
    # XXX
    # soooo, if i get it right, the groups in Carlos's script have to be defined
    # in the comparisons file and in the metadata as columns; so my idea to
    # solve this and allow comparisons for both long and cross-sectional data is
    # to first pull all unique values from the group1 and group2 cols in the
    # comparisons file, then select respective column from the metadata table.
    # They are dummy-coded so we can merge them into one single column: group

    comparisons <- .auto_read_csv(comparisons_file)

    ## unique levels
    unique_comp_levels <- unique(c(comparisons$group1, comparisons$group2))

    ## select the columns from samples, which match the unique_comp_levels
    extra_cols <- c(extra_cols, unique_comp_levels)

    ## default name
    comparisons$variable <- "group"
  } else {
    comparisons <- NULL
  }

  # rename and keep only requested extra cols (if present)
  names(samples)[1] <- "sample_id"
  keep_cols <- c("sample_id", intersect(extra_cols, names(samples)))
  samples <- samples[keep_cols]

  # handle the grouping variables
  if (!is.null(comparisons_file)) {
    # run the check ------------------------------------------------------------
    .chk_cond(
      any(rowSums(samples[, unique_comp_levels]) != 1),
      sprintf(
        "The grouping columns in the samples_file have to be mutually
        exclusive; %d row%s violate this (%s).",
        sum(rowSums(samples[, unique_comp_levels]) != 1),
        ifelse(sum(rowSums(samples[, unique_comp_levels]) != 1) == 1, "", "s"),
        paste(which(rowSums(samples[, unique_comp_levels]) != 1),
          collapse = ", "
        )
      ),
      error = TRUE
    )

    ## recode the dummy vars
    which_max <- max.col(samples[, unique_comp_levels], ties.method = "first")
    samples$group <- names(samples[, unique_comp_levels])[which_max]

    ## delete the dummy cols
    samples <- samples[, !(names(samples) %in% unique_comp_levels)]
  }

  ## ---- 2.  counts  ----------------------------------------------------------
  # The steps needed to curate the old-style data into the new, long format are:
  #   1.) Pivot the exist_file to long format
  #   2.) Attach the metadata
  #   3.) Convert the sample_id to subject_id + timepoint (cause sample_id in
  #      the mock files right now is a mixture of subject x timepoint - because
  #      of that we currently need the samples2ind table to encode th IDs)
  #   4.) Construct the phip_data S3 object with comparisons

  if (backend == "memory") {
    # Step 1.) =================================================================
    counts <- .auto_read_csv(exist_file)
    .validate_exist(counts)

    # the result of this is a full crossover of peptide x sample_ID --> it
    # contains all possible combinations of the both variables and additionally
    # a binary variable "present" defining if the combination was detected or
    # not
    counts <- stats::reshape(
      counts,
      direction = "long",
      varying   = names(counts)[-1], # all sample columns
      v.names   = "present", # new column with 0/1
      timevar   = "sample_id", # new column with sample name
      times     = names(counts)[-1] # keep the original names
    )

    counts <- counts[, -ncol(counts)] # remove the last col
    names(counts)[1] <- "peptide_id" # rename
    row.names(counts) <- NULL # tidy row-names

    # Step 2.) =================================================================
    counts <- merge(counts, samples, by = "sample_id")

    # Step 3.) - is actually optional ==========================================
    if (!is.null(timepoints_file)) {
      # load the timepoints
      timepoints <- .auto_read_csv(timepoints_file)

      # reshape the timepoints to long
      timepoints <- stats::reshape(
        timepoints,
        direction = "long",
        varying   = names(timepoints)[-1], # all sample columns
        v.names   = "sample_id",
        times     = names(timepoints)[-1], # keep original names
        timevar   = "timepoint",
        idvar     = "subject_id"
      )

      # remove the last variable, rename the first and reset the rownames
      timepoints <- timepoints[, -ncol(timepoints)]
      names(timepoints)[names(timepoints) == "ind_id"] <- "subject_id"
      row.names(timepoints) <- NULL

      # filter the NAs out
      timepoints <- timepoints[!is.na(timepoints$sample_id), ]

      # ---------- 1. reference sets -------------------------------------------
      valid_tp <- unique(timepoints$timepoint) # all time-points
      has_group <- "group" %in% names(counts) # did we already make a group col?
      valid_grp <- if (has_group) unique(counts$group) else character()

      valid_vals <- union(valid_tp, valid_grp)
      # what comparisons may legally refer to?

      # ---------- 2. sanity-check comparisons ---------------------------------
      bad <- setdiff(
        unique(c(comparisons$group1, comparisons$group2)),
        valid_vals
      )

      .chk_cond(
        length(bad) > 0L,
        sprintf(
          "Found comparison(s) referring to unknown group / timepoint: %s.\n
          Valid values are: %s",
          paste(bad, collapse = ", "),
          paste(valid_vals, collapse = ", ")
        )
      )

      # ---------- 3. merge or drop redundant column ---------------------------
      # add the time-point info
      counts <- merge(timepoints,
        counts,
        by = "sample_id"
      )

      # if 'group' duplicates 'timepoint' row-by-row, drop it to keep the table
      # tidy
      if ("group" %in% names(counts) &&
        identical(counts$group, counts$timepoint)) {
        counts$group <- NULL
        if (!is.null(comparisons_file)) {
          comparisons$variable <- "timepoint"
        }
      }
    } else {
      # rename the sample_id to subject_id if data not longitudinal
      colnames(counts)[1] <- "subject_id"
    }

    new_phip_data(
      data_long = counts,
      comparisons = comparisons,
      backend = "memory"
    )
  } else if (backend == "duckdb") {
    # --------------------------------------------------------------------------
    ## -------------------------------------------------------------------------
    ## DuckDB path -------------------------------------------------------------
    ## -------------------------------------------------------------------------

    # using the helper (a lot of code is repeated in duckdb and arrow, so i
    # decided to export it into a separate internale helper to reuse it)
    con <- .build_counts_duckdb(
      exist_file, samples,
      timepoints_file, comparisons, extra_cols
    )
    counts_tbl <- dplyr::tbl(con, "counts_final")
# print(head(counts_tbl, 5))
    if (DBI::dbExistsTable(con, "comparisons")) {
      comparisons <- dplyr::tbl(con, "comparisons") |> dplyr::collect()
    } else {
      comparisons <- NULL
    }

    # returning the phip_data object
    new_phip_data(
      data_long   = counts_tbl,
      comparisons = comparisons,
      backend     = "duckdb",
      meta        = list(con = con)
    )
  } else if (backend == "arrow") {
    # --------------------------------------------------------------------------
    ## -------------------------------------------------------------------------
    ## Arrow path --------------------------------------------------------------
    ## -------------------------------------------------------------------------
    ## check dependency
    rlang::check_installed(c("arrow"), reason = "arrow backend")

    # same as up - use helper to create the data
    con <- .build_counts_duckdb(
      exist_file, samples,
      timepoints_file, comparisons, extra_cols
    )

    # arrow-specific code, create tempdir to store the data
    arrow_dir <- fs::path(
      tempdir(),
      sprintf(
        "phip_arrow_%s",
        format(
          Sys.time(),
          "%Y%m%d%H%M%OS6"
        )
      )
    )
    dir.create(arrow_dir)

    # store the data as .parquet (more efficient than plain .csv)
    DBI::dbExecute(
      con,
      sprintf(
        "COPY counts_final TO %s (FORMAT 'parquet', PER_THREAD_OUTPUT TRUE);",
        DBI::dbQuoteString(con, arrow_dir)
      )
    )

    counts_ds <- arrow::open_dataset(arrow_dir)

    if (DBI::dbExistsTable(con, "comparisons")) {
      comparisons <- dplyr::tbl(con, "comparisons") |> dplyr::collect()
    } else {
      comparisons <- NULL
    }

    # returning the phip_data object
    new_phip_data(
      data_long   = counts_ds,
      comparisons = comparisons,
      backend     = "arrow",
      meta        = list(parquet_dir = arrow_dir)
    )
  }
}

# ------------------------------------------------------------------------------
#  helper: fastest possible CSV reader with delimiter sniffing -
#   variation on the Carlos's function
# ------------------------------------------------------------------------------
.auto_read_csv <- function(path,
                           ...) {
  # read header line only
  hdr <- readLines(path, n = 1L, warn = FALSE)

  # count delimiters without stringr (base R only)
  n_comma <- lengths(regmatches(hdr, gregexpr(",", hdr, fixed = TRUE)))
  n_semi <- lengths(regmatches(hdr, gregexpr(";", hdr, fixed = TRUE)))

  sep <- if (n_semi > n_comma) ";" else ","

  # prefer data.table::fread() if available
  if (requireNamespace("data.table", quietly = TRUE)) {
    rlang::check_installed("data.table")
    data.table::fread(path,
      sep = sep, data.table = FALSE,
      check.names = FALSE, showProgress = FALSE, ...
    )
  } else {
    utils::read.csv(path,
      header = TRUE, check.names = FALSE, sep = sep,
      stringsAsFactors = FALSE, ...
    )
  }
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
.build_counts_duckdb <- function(exist_file,
                                 samples,
                                 timepoints_file = NULL,
                                 comparisons = NULL,
                                 extra_cols = character()) {
  rlang::check_installed(c("duckdb", "DBI", "dbplyr"),
    reason = "duckdb backend"
  )
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")

  q <- function(x) DBI::dbQuoteString(con, x) # safe path quoting

  # ---- 1. counts_wide  ----
  DBI::dbExecute(
    con,
    sprintf(
      "CREATE TEMP TABLE counts_wide AS
      SELECT * FROM parquet_scan(%s);",
    # SELECT * FROM read_csv_auto(%s, HEADER=TRUE);",
    #
      q(exist_file)
    )
  )

  # ---- validate the original wide table ----
  counts_wide_tbl <- dplyr::tbl(con, "counts_wide")
  .validate_exist(counts_wide_tbl)

  first_col_counts <- DBI::dbGetQuery(
    con,
    "SELECT column_name
       FROM duckdb_columns()
      WHERE table_name = 'counts_wide'
      ORDER BY column_index
      LIMIT 1;"
  )$column_name

  samp_cols <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT column_name
         FROM duckdb_columns()
        WHERE table_name = 'counts_wide'
          AND column_name <> '%s'
        ORDER BY column_index;",
      first_col_counts
    )
  )$column_name

  ## cast the exists binary to integer

  # 1. Get all the column names
  all_cols <- DBI::dbGetQuery(con, "PRAGMA table_info('counts_wide');")$name

  # 2. Identify the first (character) column
  char_col  <- all_cols[1]

  # 3. The rest should become integer
  int_cols  <- setdiff(all_cols, char_col)

  # 4. Loop and alter each one in place
  for (col in int_cols) {
    sql <- sprintf(
      "ALTER TABLE counts_wide
       ALTER COLUMN %s TYPE INTEGER
       USING CAST(%s AS INTEGER);",
      DBI::dbQuoteIdentifier(con, col),
      DBI::dbQuoteIdentifier(con, col)
    )
    DBI::dbExecute(con, sql)
  }

  # pivot
  DBI::dbExecute(
    con,
    sprintf(
      "CREATE TEMP TABLE counts AS
         SELECT %s AS peptide_id, sample_id, present
           FROM counts_wide
           UNPIVOT (present FOR sample_id IN (%s));",
      DBI::dbQuoteIdentifier(con, first_col_counts),
      paste(DBI::dbQuoteIdentifier(con, samp_cols), collapse = ", ")
    )
  )

  # ---- 2. samples metadata  ----
  duckdb::duckdb_register(con, "samples_raw", samples)

  first_col_samples <- DBI::dbGetQuery(
    con,
    "SELECT column_name
       FROM duckdb_columns()
      WHERE table_name = 'samples_raw'
      ORDER BY column_index
      LIMIT 1;"
  )$column_name

  DBI::dbExecute(
    con,
    sprintf(
      "CREATE TEMP TABLE samples2 AS
         SELECT %s AS sample_id,
                * EXCLUDE (%s)
           FROM samples_raw;",
      DBI::dbQuoteIdentifier(con, first_col_samples),
      DBI::dbQuoteIdentifier(con, first_col_samples)
    )
  )

  DBI::dbExecute(
    con,
    sprintf(
      "CREATE TEMP TABLE counts2 AS
         SELECT c.*, s.* EXCLUDE sample_id
           FROM counts c
           LEFT JOIN %s s USING (sample_id);",
      "samples2"
    )
  )

  # ---- 3. optional timepoints  ----
  if (!is.null(timepoints_file)) {
    DBI::dbExecute(
      con,
      sprintf(
        "CREATE TEMP TABLE tp_wide AS
           SELECT * FROM read_csv_auto(%s, HEADER=TRUE);",
        q(timepoints_file)
      )
    )

    tp_cols <- DBI::dbGetQuery(
      con,
      "SELECT column_name
         FROM duckdb_columns()
        WHERE table_name = 'tp_wide' AND column_name <> 'ind_id'
        ORDER BY column_index;"
    )$column_name

    DBI::dbExecute(
      con,
      sprintf(
        "CREATE TEMP TABLE timepoints AS
           SELECT ind_id AS subject_id, timepoint, sample_id
             FROM tp_wide
             UNPIVOT (sample_id FOR timepoint IN (%s))
            WHERE sample_id IS NOT NULL;",
        paste(DBI::dbQuoteIdentifier(con, tp_cols), collapse = ", ")
      )
    )

    ## sanity check on comparisons
    # make a dplyr tbl so we can use distinct()/pull()
    timepoints_tbl <- dplyr::tbl(con, "timepoints")

    ##     a) unique time-points
    .data <- rlang::.data ## to silence the lintr note: no visible binding
    ## for the variable timepoint or .data

    valid_tp <- timepoints_tbl |>
      dplyr::distinct(.data$timepoint) |>
      dplyr::pull(.data$timepoint)

    ##     b) unique groups (only if the column exists in counts2)
    has_group <- DBI::dbGetQuery(
      con,
      "SELECT EXISTS (
         SELECT 1
           FROM duckdb_columns()
          WHERE table_name = 'counts2' AND column_name = 'group'
       ) AS has_group;"
    )$has_group

    valid_grp <- if (has_group) {
      dplyr::tbl(con, "counts2") |>
        dplyr::distinct(.data$group) |>
        dplyr::pull()
    } else {
      character()
    }

    valid_vals <- union(valid_tp, valid_grp)

    ## -----------------------------------------------------------------------
    ## 5.  merge time-points into counts and keep tidy (= R’s merge())
    ## -----------------------------------------------------------------------
    DBI::dbExecute(
      con,
      "CREATE OR REPLACE TEMP TABLE counts_final AS
       SELECT tp.*, c2.* EXCLUDE sample_id
         FROM timepoints tp
         JOIN counts2   c2 USING (sample_id);"
    )

    ## -----------------------------------------------------------------------
    ## 6.  drop `group` when it is *exactly* the same as `timepoint`
    ## -----------------------------------------------------------------------
    if (has_group) {
      identical_cols <- DBI::dbGetQuery(
        con,
        "SELECT COUNT(*) AS n_diff
         FROM counts_final
        WHERE \"group\" IS DISTINCT FROM timepoint;"
      )$n_diff == 0L

      if (identical_cols) {
        DBI::dbExecute(con, "ALTER TABLE counts_final DROP COLUMN \"group\";")
        comparisons$variable <- "timepoint"
      }
    }
  } else {
    DBI::dbExecute(
      con,
      "CREATE TEMP TABLE counts_final AS
         SELECT sample_id AS subject_id,
                * EXCLUDE sample_id
           FROM counts2;"
    )
  }
  if(!is.null(comparisons)) {
    duckdb::duckdb_register(con, "comparisons", comparisons)
  }


  # ---- return the live connection ------------------------------------------
  invisible(con)
}

## ---- valdiate  exist_file ---------------------------------------------------
### actually the only two requirements for the exist_file are, that the first
### column consists ONLY of unique IDs, and the values in the rest of the
### columns are only 0's and 1's
#' @keywords internal
.validate_exist <- function(table) {
  cols <- colnames(table)
  id_col <- cols[1]
  other_cols <- cols[-1]
  id_sym <- rlang::sym(id_col)

  # quick validate if at least two rows in the data
  .data <- rlang::.data ## to silence the linting and check errors
  n_rows <- table |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::pull(.data$n)

  .chk_cond(
    n_rows < 1,
    "The `exist_file` has no rows! No peptides are specified."
  )

  # quick validate if at least two cols in the data
  #print(head(table, 5))
  .chk_cond(
    length(cols) < 2,
    "The `exist_file` has only one column! No subjects are specified."
  )
  # 1) missing IDs
  na_ids <- table |>
    dplyr::summarise(missing = sum(is.na(!!id_sym))) |>
    dplyr::pull(missing)
  .chk_cond(
    na_ids > 0L,
    sprintf("Column '%s' must not contain NA values.", id_col)
  )

  # 2) duplicate IDs
  dup_cnt <- table |>
    dplyr::count(!!id_sym, name = "n") |>
    dplyr::filter(.data$n > 1L) |>
    dplyr::summarise(total = sum(.data$n, na.rm = TRUE)) |>
    dplyr::pull(.data$total)

  if (is.na(dup_cnt)) dup_cnt <- 0

  .chk_cond(
    dup_cnt > 0L,
    sprintf("Found %d duplicate IDs in column '%s'.", dup_cnt, id_col)
  )

  # 3) rest only 0/1 or NA
  if (length(other_cols) > 0) {
    bad_row <- table |>
      dplyr::filter(
        dplyr::if_any(
          dplyr::all_of(other_cols),
          ~ !(.x %in% c(0, 1) | is.na(.x))
        )
      ) |>
      dplyr::collect()
    if (nrow(bad_row) > 0) {
      bad_vals <- bad_row[other_cols]

      bad_info <- mapply(
        FUN = function(value, col_name) {
          if (!is.na(value) && !(value %in% c(0, 1))) {
            sprintf("'%s'=%s", col_name, value)
          } else {
            NA_character_
          }
        },
        value = bad_vals,
        col_name = names(bad_vals),
        USE.NAMES = FALSE
      )

      .chk_cond(
        TRUE,
        sprintf("Invalid value in exist table: %s", stats::na.omit(bad_info)[1])
      )
    }
  }

  invisible(TRUE)
}
