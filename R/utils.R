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

# ==============================================================================
# STRING FORMATTING
# ==============================================================================
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
# -- NULL-coalescing helper ----------------------------------------------------
`%||%` <- function(x, y) if (!is.null(x)) x else y


# ==============================================================================
# DATABASE HELPERS
# ==============================================================================

# Ensure peptide_library is queryable from the SAME connection as data_long.
# First tries a zero-copy DuckDB ATTACH; falls back to copy_to() as a temp
# table when the two connections are incompatible. The ATTACH is READ_ONLY: a
# writable attach of phiperio's cache file corrupts tables that
# get_peptide_library() later writes to it in the same R session.
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
          sprintf("ATTACH '%s' AS %s (READ_ONLY);", pep_db_path, schema_alias)
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

# Fetch the phiperio peptide libraries the given peptide IDs belong to, detected
# from their prefixes (agilent_, humanProteome_, icam_, ...). Returns NULL when
# no ID matches a known library.
#' @keywords internal
.ph_library_for_peptides <- function(peptide_ids) {
  libraries <- detect_peptide_libraries(peptide_ids)
  if (length(libraries) == 0L) return(NULL)
  get_peptide_library(libraries)
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
