# ==============================================================================
# phiperio logging utilities (ASCII only; based on the chk and cli packages)
# ==============================================================================
# ---- user-tweakable globals (set via options() in .onLoad or zzz.R) ----------
# options(
#   phiperio.log.verbose   = TRUE,
#   phiperio.log.time_fmt  = "%Y-%m-%d %H:%M:%S",
#   phiperio.log.width     = getOption("width", 80)
# )

# ==============================================================================
# LOGGING: OPTIONS + FORMATTING
# ==============================================================================
#' @title Internal helper: .ph_opt
#' @description Read a phiperio logging option with a fallback default.
#' @keywords internal
.ph_opt <- function(key,
                    default = NULL) {
  getOption(paste0("phiperio.log.", key), default)
}

#' @title Internal helper: .ph_now
#' @description Return the current time formatted for log prefixes.
#' @keywords internal
.ph_now <- function() {
  format(Sys.time(), .ph_opt("time_fmt", "%H:%M:%S"))
}

#' @title Internal helper: .ph_base_prefix
#' @description Build the base log prefix for a given level label.
#' @keywords internal
.ph_base_prefix <- function(level = "INFO") {
  sprintf("[%s] %-5s ", .ph_now(), toupper(level)[1])
}

# wraps the text nicely, regardless of the console width
#' @title Internal helper: .ph_wrap
#' @description Wrap text to the configured width, preserving the prefix.
#' @keywords internal
.ph_wrap <- function(text,
                     prefix) {
  w <- .ph_opt("width", getOption("width", 80))
  # strwrap: 'initial' for first line, 'prefix' for continuations
  strwrap(text, width = w, initial = prefix, prefix = strrep(
    " ",
    nchar(prefix)
  ))
}

# Compose multi-depth message lines
# currently the maximal supported log depth is 3:
# headline (required), step (optional), bullets (optional chr vec)
#' @title Internal helper: .ph_compose_lines
#' @description Compose multi-line log output with headline, step, and bullets.
#' @keywords internal
.ph_compose_lines <- function(level,
                              headline,
                              step = NULL,
                              bullets = NULL) {
  base <- .ph_base_prefix(level)
  stepP <- paste0(strrep(" ", nchar(base)), "-> ")
  bullP <- paste0(strrep(" ", nchar(base)), "  - ")

  out <- character(0)
  if (!is.null(headline) && nzchar(headline)) {
    out <- c(out, .ph_wrap(headline, base))
  }
  if (!is.null(step) && nzchar(step)) {
    out <- c(out, .ph_wrap(step, stepP))
  }
  if (length(bullets)) {
    for (b in bullets) {
      if (isTRUE(is.na(b)) || !nzchar(b)) next
      out <- c(out, .ph_wrap(b, bullP))
    }
  }
  out
}

# ==============================================================================
# LOGGING: EMITTERS
# ==============================================================================
## monitor task progress
#' @title Internal helper: .ph_log_info
#' @description Emit an INFO log block if verbose logging is enabled.
#' @keywords internal
.ph_log_info <- function(headline,
                         step = NULL,
                         bullets = NULL,
                         verbose = .ph_opt("verbose", TRUE)) {
  if (!isTRUE(verbose)) {
    return(invisible(character()))
  }
  lines <- .ph_compose_lines("INFO", headline, step, bullets)
  cat(paste0(lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}

## monitor task progression
#' @title Internal helper: .ph_log_ok
#' @description Emit an OK log block if verbose logging is enabled.
#' @keywords internal
.ph_log_ok <- function(headline,
                       step = NULL,
                       bullets = NULL,
                       verbose = .ph_opt("verbose", TRUE)) {
  if (!isTRUE(verbose)) {
    return(invisible(character()))
  }
  lines <- .ph_compose_lines("OK", headline, step, bullets)
  cat(paste0(lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}

# Warnings/errors via chk, but formatted to match the style of the logger
#' @title Internal helper: .ph_warn
#' @description Emit a WARN log block using chk when available.
#' @keywords internal
.ph_warn <- function(headline, step = NULL, bullets = NULL, ...) {
  lines <- .ph_compose_lines("WARN", headline, step, bullets)
  msg <- paste(lines, collapse = "\n")
  if (requireNamespace("chk", quietly = TRUE)) {
    chk::wrn(msg, ...) # single string -> respects \n
  } else {
    warning(msg, call. = FALSE, ...) # fallback if chk not installed
  }
  invisible(lines)
}

#' @title Internal helper: .ph_abort
#' @description Emit an ERROR log block and abort execution.
#' @keywords internal
.ph_abort <- function(headline, step = NULL, bullets = NULL, ...) {
  lines <- .ph_compose_lines("ERROR", headline, step, bullets)
  msg <- paste(lines, collapse = "\n")
  if (requireNamespace("chk", quietly = TRUE)) {
    chk::abort_chk(msg, ...) # single string -> respects \n
  } else {
    stop(msg, call. = FALSE, ...) # fallback if chk not installed
  }
}

# ==============================================================================
# LOGGING: TIMING
# ==============================================================================
## many tasks in phiperio can be long/take a while; it was important to have the
## infos on timing - this func wraps a task to get a start/end pair in the same
## style
#' @title Internal helper: .ph_with_timing
#' @description Conditionally raise a formatted warning or error.
#' @keywords internal
.ph_with_timing <- function(headline,
                            step = NULL,
                            bullets = NULL,
                            expr,
                            verbose = .ph_opt("verbose", TRUE)) {
  t0 <- Sys.time()
  .ph_log_info(headline = headline, step = step, bullets = bullets,
               verbose = verbose)

  res <- tryCatch(
    {
      force(expr)
    }, # evaluate user's code
    finally = {
      dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)
      .ph_log_ok(
        headline = paste0(headline, " - done"),
        step     = sprintf("elapsed: %ss", dt),
        verbose  = verbose
      )
    }
  )

  res
}

# ==============================================================================
# VALIDATION / CHECKS (unified with phiperio logger)
# ==============================================================================

# original conditional helper to not break down older code
# upgraded to the unified phiperio style
#' @title Internal helper: .ph_check_cond
#' @description Execute an expression with start/stop timing logs.
#' @keywords internal
.ph_check_cond <- function(condition,
                           error_message,
                           error = TRUE,
                           step = NULL,
                           bullets = NULL,
                           ...) {
  # log nopthing
  if (!isTRUE(condition)) {
    return(invisible(FALSE))
  }

  # print error and abort exec
  if (isTRUE(error)) {
    .ph_abort(headline = error_message, step = step, bullets = bullets, ...)
  } else {
    # print warning and go on
    .ph_warn(headline = error_message, step = step, bullets = bullets, ...)
  }
  invisible(TRUE)
}

#' @title Internal helper: .ph_check_pd
#' @description Validate that an object is a phip_data instance.
#' @keywords internal
.ph_check_pd <- function(obj) {
  .ph_check_cond(
    !inherits(obj, "phip_data"),
    "`x` must be a <phip_data> object."
  )
}

# -- check if filename has given extension -------------------------------------
# comes in handy when loading .csv or .parquet; provide filename and vector of
# extensions to check (eg c(".csv", ".parquet"))
#' @title Internal helper: .ph_check_extension
#' @description Validate a filename extension against an allowed set.
#' @keywords internal
.ph_check_extension <- function(name,
                                x_name,
                                ext_vec) {
  if (is.null(ext_vec) || !length(ext_vec)) {
    return(invisible(TRUE))
  }

  base <- basename(name %||% "") # extracting filename from paths
  parts <- strsplit(base, "\\.", fixed = FALSE)[[1]] # names + complex ext

  # taking last ext after . (eg .tar.gz --> .gz)
  ext <- if (length(parts) > 1L) {
    tolower(paste(parts[-1L], collapse = "."))
  } else {
    ""
  }

  norm <- function(x) sub("^\\.+", "", tolower(x)) # normalize ext
  got <- if (nzchar(ext)) norm(ext) else "<none>"
  ok <- nzchar(ext) && got %in% norm(ext_vec) # final ext check

  if (!ok) {
    .ph_abort(
      headline = sprintf("Invalid file extension for `%s`.", x_name),
      step = sprintf("validating path: %s", name),
      bullets = c(
        sprintf("got: %s", .ph_add_quotes(got, 2L)),
        sprintf(
          "allowed: %s",
          .ph_word_list(.ph_add_quotes(norm(ext_vec), 2L), and_or = "or")
        )
      )
    )
  }
  invisible(TRUE)
}

# -- check if NULL and replace with default when TRUE (warn in unified style) --
#' @title Internal helper: .ph_check_null_default
#' @description Replace NULL with a default, logging a warning.
#' @keywords internal
.ph_check_null_default <- function(x,
                                   x_name,
                                   method,
                                   default) {
  if (is.null(x)) {
    # format the default for print
    fmt <- function(v) {
      if (is.character(v) && length(v) == 1L) {
        return(.ph_add_quotes(v, 2L))
      }
      if (is.atomic(v) && length(v) == 1L) {
        return(as.character(v))
      }
      sprintf("<%s>", paste(class(v), collapse = "/"))
    }

    # generate warning and the replace
    .ph_warn(
      headline = sprintf("Argument `%s` missing; using default.", x_name),
      step     = sprintf("method: %s", .ph_add_quotes(method, 2L)),
      bullets  = sprintf("default: %s", fmt(default))
    )
    x <- default
  }
  x
}

# -- validate path to file -----------------------------------------------------
#' @title Internal helper: .ph_check_path
#' @description Validate a file or directory path (and optional extension).
#' @keywords internal
.ph_check_path <- function(path,
                           arg_name,
                           extension,
                           is_dir = FALSE) {
  ## error when path not a string
  .ph_check_cond(
    !chk::vld_string(path),
    sprintf("`%s` must be a character scalar.", arg_name),
    step    = "path validation",
    bullets = sprintf("got class: %s", paste(class(path), collapse = "/"))
  )

  ## error when path does not exist
  if (is_dir) {
    .ph_check_cond(
      !chk::vld_dir(path),
      sprintf("Folder for `%s` does not exist.", arg_name),
      step    = "path validation",
      bullets = sprintf("path: %s", path)
    )
  } else {
    .ph_check_cond(
      !chk::vld_file(path),
      sprintf("File for `%s` does not exist.", arg_name),
      step    = "path validation",
      bullets = sprintf("path: %s", path)
    )
  }

  # optionally extension check if provided
  if (!missing(extension) && length(extension)) {
    ## error when both is_dir and extension are given
    .ph_check_cond(
      is_dir,
      sprintf("Can't check if `%s` is both a valid direcotry and has a
                certain extension", arg_name),
      step    = "path validation",
      bullets = sprintf("path: %s", path)
    )

    .ph_check_extension(
      path,
      arg_name,
      extension
    )
  }

  invisible(TRUE)
}

# ==============================================================================
# STRING FORMATTING
# ==============================================================================
# -- clean wordlists for message generation ------------------------------------
# for multiple arguments/values
#' @title Internal helper: .ph_word_list
#' @description Build a human-readable list from a character vector, with
#'   optional quoting and conjunction handling.
#' @keywords internal
.ph_word_list <- function(word_list = NULL,
                          and_or = "and",
                          is_are = FALSE,
                          quotes = FALSE) {
  # Make "a and b" / "a, b, and c"; optionally append "is/are".
  word_list <- setdiff(word_list, c(NA_character_, ""))

  if (is.null(word_list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }

  word_list <- .ph_add_quotes(word_list, quotes)

  len_wl <- length(word_list)

  if (len_wl == 1L) {
    out <- word_list
    if (is_are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }

  if (is.null(and_or) || isFALSE(and_or)) {
    out <- paste(word_list, collapse = ", ")
  } else {
    and_or <- match.arg(and_or, c("and", "or"))
    if (len_wl == 2L) {
      out <- sprintf("%s %s %s", word_list[1L], and_or, word_list[2L])
    } else {
      out <- sprintf(
        "%s, %s %s",
        paste(word_list[-len_wl], collapse = ", "),
        and_or, word_list[len_wl]
      )
    }
  }

  if (is_are) out <- sprintf("%s are", out)
  attr(out, "plural") <- TRUE
  out
}

# -- quoting helper (unified error style) --------------------------------------
# define number of quotes you want --> for printing logs/messages/warnings
# or define the quotes itself as a string
#' @title Internal helper: .ph_add_quotes
#' @description Wrap character values in quotes for log output. Supports
#'   FALSE/TRUE, 0/1/2, or a single-character quote string.
#' @keywords internal
.ph_add_quotes <- function(x,
                           quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }
  if (isTRUE(quotes)) quotes <- '"'

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, quotes))
  }

  if (!chk::vld_count(quotes) || quotes > 2) {
    .ph_abort(
      headline = "Invalid `quotes` argument.",
      step = "formatting .ph_add_quotes()",
      bullets = c(
        "allowed: FALSE, TRUE, 0, 1, 2, or a single-character string",
        sprintf("got class: %s", paste(class(quotes), collapse = "/"))
      )
    )
  }

  if (quotes == 0L) {
    return(x)
  }
  if (quotes == 1L) {
    return(sprintf("'%s'", x))
  }
  sprintf('"%s"', x)
}

# ==============================================================================
# OPERATORS
# ==============================================================================
# -- not-in operator -----------------------------------------------------------
`%nin%` <- function(x, inx) {
  !(x %in% inx)
}

# -- NULL-coalescing helper ----------------------------------------------------
`%||%` <- function(x, y) if (!is.null(x)) x else y


# ==============================================================================
# PATH RESOLUTION (FAST-FAIL)
# ==============================================================================
#' @title Resolve legacy-import paths and perform fast-fail argument checks
#'
#' @description Combines explicit arguments with a YAML config (if given),
#'   expands every relative path to an absolute path (relative paths are
#'   evaluated against `dirname(config_yaml)` (!!!) when YAML is used, otherwise
#'   against the directory that contains the first supplied data matrix (!!!)),
#'   and returns a fully populated list of file locations and options ready for
#'   downstream conversion. Only cheap, load-blocking checks are done here:
#'
#' * `input_file` and `hit_file` must be supplied together or both omitted.
#'
#' * At least one matrix source (`exist_file`, `fold_change_file`, or the
#'   `input_file` + `hit_file` pair) must be present.
#'
#' * Deprecated `output_dir` triggers a soft warning.
#'
#'  All deeper table-content validation is deferred to `phip_data` class
#'  validation.
#'
#' @param exist_file,fold_change_file,input_file,hit_file,samples_file,timepoints_file
#'  Character paths (relative or absolute) to the respective CSV/Parquet inputs.
#'  `NULL` means "not supplied".
#' @param extra_cols Character vector of extra metadata columns to keep; may be
#'   `NULL`.
#' @param output_dir Ignored (soft-deprecated).
#' @param config_yaml Optional path to a YAML file whose keys mirror the
#'   function arguments; relative paths inside the YAML are resolved against the
#'   YAML’s own directory.
#'
#' @return A named list with absolute paths, `extra_cols`, and `base_dir`;
#'   suitable for downstream helper functions.
#'
#' @keywords internal
.ph_resolve_paths <- function(
    exist_file = NULL,
    fold_change_file = NULL,
    samples_file = NULL,
    input_file = NULL,
    hit_file = NULL,
    timepoints_file = NULL,
    extra_cols = NULL,
    output_dir = NULL, # deprecated
    data_long_path = NULL,
    peptide_library = TRUE,
    n_cores = NULL,
    materialise_table = NULL,
    auto_expand = NULL,
    sample_id_from_filenames = NULL,
    config_yaml = NULL
) {
  ## ------------------------------------------------------------------------ ##
  ## 1.  locate base directory & read yaml (if any provided)                  ##
  ## ------------------------------------------------------------------------ ##
  is_abs_path <- function(p) {
    grepl("^(/|[A-Za-z]:[/\\\\])", p)
  }
  abs_path <- function(p, start = NULL) {
    if (is.null(p)) return(p)
    if (!is.null(start) && !is_abs_path(p)) p <- file.path(start, p)
    normalizePath(p, winslash = "/", mustWork = FALSE)
  }

  # Determine base_dir depending on which input source was provided
  base_dir <- if (!is.null(config_yaml)) {
    # 1) If a YAML config was given, take its folder
    dirname(abs_path(config_yaml))
  } else if (!is.null(data_long_path)) {
    # 2) If a data_long_path directory was given - use it
    dirname(abs_path(data_long_path))
  } else {
    # 3) Otherwise require at least a samples_file (or exist_file)
    .ph_check_cond(
      is.null(samples_file) && is.null(exist_file),
      "When neither 'config_yaml' nor 'data_long_path' is provided,
      you must supply 'samples_file' or 'exist_file'."
    )
    # pick samples_file if present, else exist_file, then take its parent folder
    dirname(abs_path(samples_file %||% exist_file))
  }

  yaml_cfg <- if (!is.null(config_yaml)) {
    ## validate the file extension
    .ph_check_path(config_yaml, "config_yaml", c("yml", "yaml"))

    ## read the yamlW
    rlang::check_installed("yaml")
    yaml::read_yaml(config_yaml)
  } else {
    ## safe fallback
    NULL
  }

  ## ------------------------------------------------------------------------ ##
  ## 2.  helper to merge yaml + explicit args, validate & absolutise          ##
  ## ------------------------------------------------------------------------ ##
  fetch <- function(arg,
                    key,
                    validate = NULL,
                    optional = FALSE,
                    absolutize = FALSE,
                    ...) {
    # safe fallback to NULL --> if both NULL, then %||% returns NULL
    val <- yaml_cfg[[key]] %||% arg # yaml first, then explicit

    # required argument have to be provided!
    .ph_check_cond(
      is.null(val) && !optional,
      sprintf("Missing required argument '%s' in YAML or call.", key)
    )

    # if the validator is .ph_check_path or absolute == TRUE, expand the path
    # to absolute
    if ((!is.null(val) && identical(validate, .ph_check_path)) ||
        (!is.null(val) && absolutize)) {
      if (!is_abs_path(val)) {
        val <- abs_path(basename(val), start = base_dir)
      } else {
        val <- abs_path(val)
      }
    }

    # perform the custom validation if specified
    if (!is.null(val) && is.function(validate)) {
      validate(val, key, ...)
    }

    # return
    val
  }

  ## ------------------------------------------------------------------------ ##
  ## 3.  resolve every supported argument                                     ##
  ## ------------------------------------------------------------------------ ##
  samples_required <- !is.null(data_long_path)

  cfg <- list(
    exist_file = fetch(exist_file,
                       "exist_file",
                       .ph_check_path,
                       optional = TRUE,
                       extension = c("csv", "parquet", "parq", "pq")
    ),
    fold_change_file = fetch(fold_change_file,
                             "fold_change_file",
                             .ph_check_path,
                             optional = TRUE,
                             extension = c("csv", "parquet", "parq", "pq")
    ),
    input_file = fetch(input_file,
                       "input_file",
                       .ph_check_path,
                       optional = TRUE,
                       extension = c("csv", "parquet", "parq", "pq")
    ),
    hit_file = fetch(hit_file,
                     "hit_file",
                     .ph_check_path,
                     optional = TRUE,
                     extension = c("csv", "parquet", "parq", "pq")
    ),
    samples_file = fetch(samples_file,
                         "samples_file",
                         .ph_check_path,
                         optional = samples_required,
                         extension = c("csv", "parquet", "parq", "pq")
    ),
    timepoints_file = fetch(timepoints_file,
                            "timepoints_file",
                            .ph_check_path,
                            optional = TRUE,
                            extension = c("csv", "parquet", "parq", "pq")
    ),
    extra_cols = fetch(extra_cols,
                       "extra_cols",
                       optional = TRUE
    ),
    output_dir = fetch(output_dir,
                       "output_dir",
                       optional = TRUE
    ),
    data_long_path = fetch(data_long_path,
                           "data_long_path",
                           optional = !samples_required,
                           absolutize = TRUE
    ),
    peptide_library = peptide_library,
    n_cores = n_cores,
    materialise_table = materialise_table,
    auto_expand = auto_expand,
    sample_id_from_filenames = sample_id_from_filenames,
    base_dir = base_dir # for downstream helpers
  )

  ## ------------------------------------------------------------------------ ##
  ## 4.  fast-fail rules that really must hold before heavy work              ##
  ## ------------------------------------------------------------------------ ##
  #  rule 1: input_file and hit_file must be provided together ----------
  .ph_check_cond(
    xor(is.null(cfg$input_file), is.null(cfg$hit_file)),
    "Arguments 'input_file' and 'hit_file' must be provided together."
  )

  # Rule 2a: if data_long_path is provided, it must be the ONLY file argument
  if (!is.null(cfg$data_long_path)) {
    others_supplied <- any(
      !is.null(cfg$exist_file),
      !is.null(cfg$fold_change_file),
      !is.null(cfg$input_file),
      !is.null(cfg$hit_file)
    )
    .ph_check_cond(
      others_supplied,
      "When 'data_long_path' is supplied, do not supply 'exist_file',
      'fold_change_file', 'input_file', or 'hit_file'."
    )
  } else {
    # Rule 2b: if data_long_path is NOT provided,
    #          require at least one of the other file arguments
    all_null <- with(
      cfg,
      is.null(exist_file) &&
        is.null(fold_change_file) &&
        is.null(input_file) &&
        is.null(hit_file)
    )
    .ph_check_cond(
      all_null,
      paste0(
        "Supply at least one of:\n",
        "  * 'exist_file'\n",
        "  * 'fold_change_file'\n",
        "  * both 'input_file' and 'hit_file'"
      )
    )
  }

  #  deprecation notice -------------------------------------------------
  .ph_check_cond(!is.null(cfg$output_dir),
                 error = FALSE,
                 "'output_dir' is deprecated and will be ignored."
  )

  # validate the tables itself --> the logic has been moved entirely to the
  # phip_data class validator

  cfg
}

# ==============================================================================
# DATABASE HELPERS
# ==============================================================================

# Ensure peptide_library is queryable from the SAME connection as data_long.
# First tries a zero-copy DuckDB ATTACH; falls back to copy_to() as a temp
# table when the two connections are incompatible.
#' @keywords internal
.ph_peplib_on_main <- function(x, schema_alias = "peplib") {
  main_con <- dbplyr::remote_con(x$data_long)
  pep_con <- if (!is.null(x$meta$peptide_con)) {
    x$meta$peptide_con
  } else {
    dbplyr::remote_con(x$peptide_library)
  }

  # -- fast path: both DuckDB -> ATTACH the peptide db file ------------------
  if (inherits(main_con, "duckdb_connection") &&
      inherits(pep_con, "duckdb_connection")) {
    pep_db_path <- try(pep_con@driver@dbdir, silent = TRUE)
    if (!inherits(pep_db_path, "try-error") &&
        is.character(pep_db_path) && nzchar(pep_db_path)) {
      try(
        DBI::dbExecute(
          main_con,
          sprintf("ATTACH '%s' AS %s;", pep_db_path, schema_alias)
        ),
        silent = TRUE
      )

      base_name <- tryCatch(
        {
          nm <- dbplyr::remote_name(x$peptide_library)
          if (is.null(nm) || !nzchar(nm)) "peptide_meta"
          else sub("^.*\\.", "", nm)
        },
        error = function(e) "peptide_meta"
      )

      try_tbl <- function(sql_expr) {
        tryCatch(dplyr::tbl(main_con, dbplyr::sql(sql_expr)),
                 error = function(e) NULL)
      }

      for (q in c(
        sprintf("SELECT * FROM %s.%s",      schema_alias, base_name),
        sprintf("SELECT * FROM %s.main.%s", schema_alias, base_name)
      )) {
        out <- try_tbl(q)
        if (!is.null(out)) return(out)
      }
    }
  }

  # -- fallback: collect peptide library and copy into main connection --------
  peplib_local <- dplyr::collect(x$peptide_library)
  tmp_name <- paste0("peptide_meta_tmp_", as.integer(Sys.time()))
  dplyr::copy_to(main_con, peplib_local, tmp_name,
                 temporary = TRUE, overwrite = TRUE)
}

# ==============================================================================
# Internal helper: order-agnostic pair filter for ph_prev_result / data.frame
# ==============================================================================

#' @title Internal helper: .ph_filter_pairs
#' @description Filter prevalence result rows by group pair(s), rank, features,
#'   universe, and optional significance thresholds. Preserves `ph_prev_result`
#'   class and `prev_meta` attribute when the input carries them.
#' @keywords internal
.ph_filter_pairs <- function(
  df,
  gA, gB,
  ranks = NULL,
  features = NULL,
  features_regex = FALSE,
  group_universe = NULL,
  universe_regex = FALSE,
  col_rank = "rank",
  col_feature = "feature",
  col_g1 = "group1",
  col_g2 = "group2",
  col_groupcol = "group_col",
  p_raw_max = NULL,
  q_bh_max = NULL,
  q_wbh_max = NULL,
  passed_only = FALSE,
  drop_na = TRUE,
  keep_cols = NULL
) {
  normalize_pairs <- function(gA, gB) {
    if (is.character(gA) && is.character(gB) && length(gA) == 1L && length(gB) == 1L) {
      return(list(c(gA, gB)))
    }
    if (is.list(gA) && missing(gB)) {
      ok <- vapply(gA, function(x) is.character(x) && length(x) == 2L, logical(1))
      if (!all(ok)) stop(".ph_filter_pairs: all list elements must be length-2 character vectors.")
      return(gA)
    }
    if (is.data.frame(gA) && missing(gB)) {
      if (!all(c("gA", "gB") %in% names(gA))) stop(".ph_filter_pairs: data.frame must have columns 'gA' and 'gB'.")
      return(split(as.matrix(gA[, c("gA", "gB"), drop = FALSE]), seq_len(nrow(gA))))
    }
    stop(".ph_filter_pairs: provide either scalars gA,gB; a list of pairs; or a data.frame with columns gA,gB.")
  }
  pairs_list <- normalize_pairs(gA, gB)

  needed <- c(col_g1, col_g2)
  if (!all(needed %in% names(df))) stop(".ph_filter_pairs: df is missing required columns: ", paste(setdiff(needed, names(df)), collapse = ", "))

  dt <- as.data.frame(df)

  if (!is.null(ranks)) {
    if (!col_rank %in% names(dt)) stop(".ph_filter_pairs: column '", col_rank, "' not found.")
    dt <- dt[as.character(dt[[col_rank]]) %in% as.character(ranks), , drop = FALSE]
  }

  if (!is.null(features)) {
    if (!col_feature %in% names(dt)) stop(".ph_filter_pairs: column '", col_feature, "' not found.")
    fvec <- as.character(dt[[col_feature]]); fvec[is.na(fvec)] <- ""
    if (isTRUE(features_regex)) {
      keep <- rep(FALSE, nrow(dt))
      for (p in as.character(features)) keep <- keep | grepl(p, fvec, ignore.case = TRUE, perl = TRUE)
      dt <- dt[keep, , drop = FALSE]
    } else {
      dt <- dt[tolower(fvec) %in% tolower(as.character(features)), , drop = FALSE]
    }
  }

  if (!is.null(group_universe) && col_groupcol %in% names(dt)) {
    gvec <- as.character(dt[[col_groupcol]]); gvec[is.na(gvec)] <- ""
    if (isTRUE(universe_regex)) {
      keep <- rep(FALSE, nrow(dt))
      for (p in as.character(group_universe)) keep <- keep | grepl(p, gvec, ignore.case = TRUE, perl = TRUE)
      dt <- dt[keep, , drop = FALSE]
    } else {
      dt <- dt[tolower(gvec) %in% tolower(as.character(group_universe)), , drop = FALSE]
    }
  }

  c1 <- dt[[col_g1]]; c2 <- dt[[col_g2]]
  ok <- if (drop_na) (!is.na(c1) & !is.na(c2)) else rep(TRUE, nrow(dt))
  lo <- ifelse(c1 <= c2, c1, c2); hi <- ifelse(c1 <= c2, c2, c1)
  keep_any <- rep(FALSE, nrow(dt))
  for (pp in pairs_list) {
    tgt <- sort(as.character(pp))
    keep_any <- keep_any | (ok & lo == tgt[1] & hi == tgt[2])
  }
  dt <- dt[keep_any, , drop = FALSE]

  if (!is.null(p_raw_max) && "p_raw" %in% names(dt)) {
    dt <- dt[is.na(dt[["p_raw"]]) | dt[["p_raw"]] <= p_raw_max, , drop = FALSE]
  }
  if (!is.null(q_bh_max) && "p_adj_rank" %in% names(dt)) {
    dt <- dt[is.na(dt[["p_adj_rank"]]) | dt[["p_adj_rank"]] <= q_bh_max, , drop = FALSE]
  }
  if (!is.null(q_wbh_max) && "p_adj_rank_wbh" %in% names(dt)) {
    dt <- dt[is.na(dt[["p_adj_rank_wbh"]]) | dt[["p_adj_rank_wbh"]] <= q_wbh_max, , drop = FALSE]
  }
  if (isTRUE(passed_only)) {
    pass_cols <- intersect(c("passed_rank_bh", "passed_rank_wbh"), names(dt))
    if (length(pass_cols)) {
      pass_any <- Reduce(`|`, lapply(pass_cols, function(cc) ifelse(is.na(dt[[cc]]), FALSE, dt[[cc]])))
      dt <- dt[pass_any, , drop = FALSE]
    }
  }

  if (!is.null(keep_cols)) {
    dt <- dt[, intersect(keep_cols, names(dt)), drop = FALSE]
  }

  res <- as.data.frame(dt)

  if (inherits(df, "ph_prev_result")) {
    meta <- attr(df, "prev_meta") %||% list()
    meta$subsetted <- TRUE
    meta$subset_n <- nrow(res)
    attr(res, "prev_meta") <- meta
    class(res) <- class(df)
  }

  res
}
