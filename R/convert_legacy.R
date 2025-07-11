#' @title Convert the Carlos's legacy 3-file PhIP-seq input into a new phiper
#' phip_dataset S3 object
#'
#' @description This version of `phip_convert_legacy()` is backward‑compatible
#'   with the old interface but can also take a single **YAML** configuration
#'   file instead of a long list of arguments.  Explicit function arguments
#'   always win over values supplied through the YAML file.
#'
#' @param exist_csv       Path to the existing‑peptides CSV created by the old
#'                        pipeline.
#' @param samples_csv     Path to the sample‑annotation CSV.
#' @param timepoint_csv   Path to the sample‑timepoint CSV.
#' @param extra_cols      A character vector of additional column names to keep
#'                        when converting.
#' @param comparisons_csv Path to the comparisons CSV.
#' @param output_dir      (Deprecated) directory in which to place converted
#'                        output files. Use the working directory instead.
#' @param backend         One of "arrow", "duckdb", or "memory".
#' @param config_yaml     **Optional.** Path to a YAML file that specifies any or
#'                        all of the above parameters (see example below).
#'
#' ```yaml
#' # Example config.yaml
#' # it follows the naming convention of previous Carlos's workflow
#' exist_file: data/existing.csv
#' samples_file: meta/samples.csv
#' timepoints_file: meta/timepoints.csv
#' extra_cols: [sex, age]
#' comparisons_file: meta/comparisons.csv
#' output_dir: out          # optional / deprecated
#' ```
#'
#' @return fff
#' @export

phip_convert_legacy <- function(
    exist_file       = NULL,
    samples_file     = NULL,
    timepoints_file  = NULL,
    extra_cols       = NULL,
    comparisons_file = NULL,
    output_dir       = NULL,   # hard deprecation
    backend          = NULL,
    config_yaml      = NULL
) {
  # ------------------------------------------------------------------
  # db-backend: default to "duckdb" if user supplies nothing
  # ------------------------------------------------------------------
  backend_choices <- c("arrow", "duckdb", "memory")

  backend <- if (is.null(backend)) {
    "duckdb"                                 # implicit default
  } else {
    match.arg(backend, choices = backend_choices)
  }

  # ------------------------------------------------------------------
  # figure out the base directory
  # ------------------------------------------------------------------
  base_dir <- if (!is.null(config_yaml)) {
    fs::path_dir(fs::path_abs(config_yaml))     # dir that holds config.yaml
  } else {
    getwd()                                 # wherever caller is
  }

  # ---------------------------------------------------------------------------
  # pick up values from YAML, if provided
  # ---------------------------------------------------------------------------
  if (!is.null(config_yaml)) {
    # validate the path to config.yaml with custom function
    .chk_path(config_yaml,
              "config_yaml",
              c("yml", 'yaml'))

    # read yaml
    cfg <- yaml::read_yaml(config_yaml)

  }

  # helper to fetch the paths and abort when path not found

  fetch <- function(arg,
                    key,
                    validator = NULL,
                    optional = FALSE,
                    ...) {
    # if both NULL, then %||% returns NULL
    val <- cfg[[key]] %||% arg          # YAML first, then explicit arg

    # abort if val is NULL (required arg not found)
    if (is.null(val) && !optional)
      chk::abort_chk(sprintf("Missing required argument `%s`: not defined in
                               the YAML config or as an explicit path. Aborting
                               data load.", key))

    # if the validator if .chk_path, expand the path to absolute
    if (!is.null(val) && identical(validator, .chk_path)) {
      if (!fs::is_absolute_path(val)) {
        val <- fs::path_abs(val, start = base_dir)
      }
    }

    # run the validator only if supplied and value is not NULL
    if (!is.null(val) && is.function(validator))
      validator(val, key, ...)

    return(val)
  }

  # ---------------------------------------------------------------------------
  # load the arguments from .yaml or use direct path definitions and
  # validate required inputs
  # ---------------------------------------------------------------------------
  exist_file       <- fetch(exist_file, "exist_file", .chk_path,extension="csv")
  samples_file     <- fetch(samples_file, "samples_file",
                            .chk_path, extension="csv")
  timepoints_file  <- fetch(timepoints_file,"timepoints_file",.chk_path,
                            optional=TRUE, extension="csv")
  extra_cols       <- fetch(extra_cols,"extra_cols",optional=TRUE)
  comparisons_file <- fetch(comparisons_file,"comparisons_file",.chk_path,
                            extension="csv")
  output_dir       <- fetch(output_dir,"output_dir",optional=TRUE)
  backend          <- fetch(backend,"backend",chk::chk_string)

  # Warn about deprecation
  .chk_cond(!is.null(output_dir), error = FALSE,
            "'output_dir' is deprecated.")

  ##############################################################################
  #  L E G A C Y   B R I D G E
  ##############################################################################
  ## ---- 1.  samples & contrasts (both small) ---------------------------------
  samples     <- .auto_read_csv(samples_file)
  comparisons <- .auto_read_csv(comparisons_file)

  # rename and keep only requested extra cols (if present)
  names(samples)[1] <- "sample_id"
  keep_cols <- c("sample_id", intersect(extra_cols, names(samples)))
  samples   <- samples[keep_cols]

  ## ---- 2.  counts  ----------------------------------------------------------
  # The steps needed to curate the old-style data into the new, long format are:
  #   1.) Pivot the exist_file to long format
  #   2.) Attach the metadata
  #   3.) Convert the sample_id to subject_id + timepoint (cause sample_id in
  #      the mock files right now is a mixture of subject x timepoint - because
  #      of that we currently need the samples2ind table to encode th IDs)
  #   4.) Construct the phip_data S3 object with comparisons

  # choose fast reader if it exists, else base
  header  <- .auto_read_csv(exist_file, nrows = 0)
  samp_id <- names(header)[-1]           # all sample columns

  if (backend == "memory") {
    # Step 1.) =================================================================
    counts <- .auto_read_csv(exist_file)

    # the result of this is a full crossover of peptide x sample_ID --> it
    # contains all possible combinations of the both variables and additionally
    # a binary variable "present" defining if the combination was detected or
    # not
    counts <- stats::reshape(
      counts,
      direction = "long",
      varying   = names(counts)[-1],          # all sample columns
      v.names   = "present",                  # new column with 0/1
      timevar   = "sample_id",                # new column with sample name
      times     = names(counts)[-1]          # keep the original names
    )

    counts <- counts[, -ncol(counts)]                  # remove the last col
    names(counts)[names(counts) == "V1"] <- "peptide"  # rename
    row.names(counts) <- NULL                          # tidy row-names

    # Step 2.) =================================================================
    counts <- merge(counts, samples, by = "sample_id")

    # Step 3.) - is actually optional ==========================================
    if(!is.null(timepoints_file)) {
      # load the timepoints
      timepoints <- .auto_read_csv(timepoints_file)

      # reshape the timepoints to long
      timepoints <- reshape(
        timepoints,
        direction = "long",
        varying   = names(timepoints)[-1],          # all sample columns
        v.names   = "sample_id",
        times     = names(timepoints)[-1],          # keep original names
        timevar   = "timepoint",
        idvar     = "subject_id"
      )

      # remove the last variable, rename the first and reset the rownames
      timepoints <- timepoints[, -ncol(timepoints)]
      names(timepoints)[names(timepoints) == "ind_id"] <- "subject_id"
      row.names(timepoints) <- NULL

      # filter the NAs out
      timepoints <- timepoints[!is.na(timepoints$sample_id), ]

      # now merge the timepoints with the counts table
      counts <- merge(timepoints,
                          counts,
                          by = "sample_id")

      print(counts)

    } else {
      # rename the sample_id to subject_id if data not longitudinal
      colnames(counts)[1] <- "subject_id"
    }
   # print(counts)

    new_phip_data(counts    = counts,
                  samples   = samples,
                  contrasts = contrasts,
                  backend   = "memory")

  } else if (backend == "duckdb") {
    # --------------------------------------------------------------------------
    ## -------------------------------------------------------------------------
    ## DuckDB path -------------------------------------------------------------
    ## -------------------------------------------------------------------------
    rlang::check_installed(c("duckdb", "DBI"), reason = "data processing")
    con <- DBI::dbConnect(duckdb::duckdb(),
                          dbdir = ":memory:")

    on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

    # helper – quote paths safely ----------------------------------------------
    q  <- function(x) DBI::dbQuoteString(con, x)

    # ---- 3‑a  counts wide → temp table ---------------------------------------
    DBI::dbExecute(
      con,
      sprintf(
        "CREATE TEMP TABLE counts_wide AS
        SELECT * FROM read_csv_auto(%s, HEADER=TRUE);",
        q(exist_file)
      )
    )

    samp_cols <- DBI::dbGetQuery(
      con,
      "SELECT column_name
      FROM duckdb_columns()
      WHERE table_name = 'counts_wide' AND column_name <> 'V1'
      ORDER BY column_index;"
    )$column_name

    DBI::dbExecute(
      con,
      sprintf(
        "CREATE TEMP TABLE counts AS\n         SELECT V1 AS peptide, sample_id, present FROM counts_wide\n         UNPIVOT (present FOR sample_id IN (%s));",
        paste(DBI::dbQuoteIdentifier(con, samp_cols), collapse = ", ")
      )
    )

    # ---- 3‑b  samples metadata into DuckDB -------------------------------
    DBI::dbExecute(
      con,
      sprintf(
        "CREATE TEMP TABLE samples AS\n         SELECT * FROM read_csv_auto(%s, HEADER=TRUE);",
        q(samples_file)
      )
    )

    # rename 1st col to sample_id and keep only requested columns -----------
    first_col <- DBI::dbGetQuery(
      con,
      "SELECT column_name FROM duckdb_columns()             \n       WHERE table_name='samples' ORDER BY column_index     \n       LIMIT 1;"
    )$column_name

    DBI::dbExecute(
      con,
      sprintf(
        "CREATE TEMP TABLE samples2 AS\n         SELECT %s AS sample_id, * EXCLUDE (%s) FROM samples;",
        DBI::dbQuoteIdentifier(con, first_col),
        DBI::dbQuoteIdentifier(con, first_col)
      )
    )

    if (length(extra_cols) > 0) {
      keep_sql <- paste(DBI::dbQuoteIdentifier(con, c("sample_id", extra_cols)),
                        collapse = ", ")
      DBI::dbExecute(con, sprintf("CREATE TEMP TABLE samples3 AS\n                                  SELECT %s FROM samples2;", keep_sql))
      samples_tbl <- "samples3"
    } else {
      samples_tbl <- "samples2"
    }

    # ---- 3‑c  join counts + samples --------------------------------------
    DBI::dbExecute(
      con,
      sprintf(
        "CREATE TEMP TABLE counts2 AS\n         SELECT c.*, s.* EXCLUDE sample_id\n         FROM counts c LEFT JOIN %s s USING (sample_id);",
        samples_tbl
      )
    )

    # ---- 3‑d  optional timepoints ----------------------------------------
    if (!is.null(timepoints_file)) {

      DBI::dbExecute(
        con,
        sprintf(
          "CREATE TEMP TABLE tp_wide AS\n           SELECT * FROM read_csv_auto(%s, HEADER=TRUE);",
          q(timepoints_file)
        )
      )

      tp_cols <- DBI::dbGetQuery(
        con,
        "SELECT column_name FROM duckdb_columns()                \n         WHERE table_name='tp_wide' AND column_name <> 'ind_id'  \n         ORDER BY column_index;"
      )$column_name

      DBI::dbExecute(
        con,
        sprintf(
          "CREATE TEMP TABLE timepoints AS\n           SELECT ind_id AS subject_id, timepoint, sample_id\n           FROM tp_wide\n           UNPIVOT (sample_id FOR timepoint IN (%s))\n           WHERE sample_id IS NOT NULL;",
          paste(DBI::dbQuoteIdentifier(con, tp_cols), collapse = ", ")
        )
      )

      DBI::dbExecute(
        con,
        "CREATE TEMP TABLE counts_final AS\n         SELECT tp.*, c2.* EXCLUDE sample_id\n         FROM timepoints tp JOIN counts2 c2 USING (sample_id);"
      )
    } else {
      DBI::dbExecute(
        con,
        "CREATE TEMP TABLE counts_final AS\n         SELECT sample_id AS subject_id, * EXCLUDE sample_id FROM counts2;"
      )
    }

    # ---- 3‑e  collect back to R (optional) -------------------------------
    counts <- DBI::dbGetQuery(con, "SELECT * FROM counts_final")
    return(new_phip_data(counts, samples, comparisons, backend = "duckdb"))
  } else {                       # backend == "arrow"
    if (!requireNamespace("arrow",  quietly = TRUE))
      stop("backend = 'arrow' but package 'arrow' is not installed.")

    dir_out <- persist_to %||% file.path(tempdir(), "arrow_phip")
    if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

    for (s in samp_id) {
      tmp <- fread(exist_csv,
                   select    = c(1, which(names(header) == s)),
                   col.names = c("peptide_id", "count"))
      tmp <- tmp[tmp$count > 0, , drop = FALSE]
      if (nrow(tmp)) {
        tmp$sample_id <- s
        tmp <- merge(tmp, samples, by = "sample_id", all.x = TRUE)
        arrow::write_parquet(tmp, file.path(dir_out, paste0(s, ".parquet")))
      }
    }

    dataset <- arrow::open_dataset(dir_out)

    new_phip_data(counts    = dataset,
                  samples   = samples,
                  contrasts = contrasts,
                  backend   = "arrow",
                  meta      = list(parquet_dir = dir_out))
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
  n_comma <- lengths(regmatches(hdr, gregexpr(",",  hdr, fixed = TRUE)))
  n_semi  <- lengths(regmatches(hdr, gregexpr(";", hdr, fixed = TRUE)))

  sep <- if (n_semi > n_comma) ";" else ","

  # prefer data.table::fread() if available
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, sep = sep, data.table = FALSE,
                      check.names = FALSE, showProgress = FALSE, ...)
  } else {
    read.csv(path, header = TRUE, check.names = FALSE, sep = sep,
             stringsAsFactors = FALSE, ...)
  }
}
