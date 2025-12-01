#' @importFrom utils head
#' @exportS3Method head phip_data
head.phip_data <- function(x, ...) {
  tryCatch(
    {
      tbl <- x$data_long
      if (inherits(tbl, c("tbl_dbi", "arrow_dplyr_query"))) {
        tbl |>
          utils::head(...) |>
          dplyr::collect()
      } else {
        utils::head(tbl, ...)
      }
    },
    error = function(e) tibble::tibble(.error = "<head access failed>")
  )
}

#' @exportS3Method dim phip_data
dim.phip_data <- function(x) {
  n_cols  = ncol(x$data_long)

  .data <- rlang::.data
  n_rows = tryCatch(
    if (inherits(x$data_long, c("tbl_dbi", "arrow_dplyr_query"))) {
      x$data_long |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::pull(.data$n)
    } else {
      nrow(x$data_long)
    },
    error = function(e) NA_integer_
  )

  c(n_rows, n_cols)
}

#' @exportS3Method print phip_data
print.phip_data <- function(x, ...) {
  cat(cli::rule(
    left  = "<phip_data>",
    right = paste0("backend: ", x$backend)
  ), "\n\n")

  # ---- counts preview -------------------------------------------------------
  cat(cli::col_cyan("counts (first 5 rows):"), "\n")
  print(head(x, 5))

  ## --- size summary ----------------------------------------------------------
  x_dim <- dim(x)

  cat("\n")
  cat(sprintf(
    "table size: %s rows x %s columns\n\n",
    format(x_dim[1], big.mark = ","), x_dim[2]
  ))

  # ---- contrasts ------------------------------------------------------------
  if (!is.null(x$comparisons)) {
    cat(cli::col_cyan("contrasts:"), "\n")
    cat(paste0(knitr::kable(x$comparisons, format = "simple"), collapse = "\n"))
    cat("\n\n")
  }

  # ---- peptide library preview ---------------------------------------------
  cat(cli::col_cyan("peptide library preview (first 5 rows):"), "\n")

  lib <- x$peptide_library

  # If the DuckDB connection was closed, reopen the cached table
  if (inherits(lib, "tbl_dbi")) {
    if (!DBI::dbIsValid(lib$src$con)) {
      lib <- get_peptide_meta()
      x$peptide_library <- lib # keep the fresh handle for later prints
    }
  }

  if (!is.null(lib)) {
    show_cols <- intersect(c("peptide_id", "pos", "len_seq"), colnames(lib))

    lib_preview <- tryCatch(
      {
        lib |>
          utils::head(5) |> # LIMIT 5
          dplyr::collect() |> # into R
          dplyr::select(dplyr::all_of(show_cols))
      },
      error = function(e) tibble::tibble(.error = "<preview failed>")
    )

    print(lib_preview)

    # summary: extra-column count + overall dim ---------------------------------
    n_cols <- length(colnames(lib))
    extra_nc <- n_cols - length(show_cols)

    .data <- rlang::.data ## to silence the R CMD CHECK and lintr
    n_rows <- tryCatch(
      lib |> dplyr::summarise(n = dplyr::n()) |> dplyr::pull(.data$n),
      error = function(e) NA_integer_
    )

    cat(sprintf("... plus %d more columns\n\n", extra_nc))
    cat(sprintf(
      "library size: %s rows x %s columns\n\n",
      format(n_rows, big.mark = ","), n_cols
    ))
  }

  # ---- meta flags -----------------------------------------------------------
  cat(cli::col_cyan("meta flags:"), "\n")
  for (nm in names(x$meta)) {
    val <- x$meta[[nm]]
    val_txt <- if (is.atomic(val)) {
      paste(val, collapse = ", ")
    } else {
      paste0("<", class(val)[1], ">")
    }
    cat(sprintf("  %-15s %s\n", paste0(nm, ":"), val_txt))
  }
  cat("\n")

  invisible(x)
}


################################################################################
## plain accessors (cause no S3 generics) --------------------------------------
################################################################################
.check_pd <- function(obj) {
  .chk_cond(
    !inherits(obj, "phip_data"),
    "`x` must be a <phip_data> object."
  )
}

#' @title Retrieve the main PhIP-Seq counts table
#'
#' @description Quick accessor for the `data_long` slot of a **phip_data**
#' object.
#'
#' @param x A valid `phip_data` object.
#'
#' @return A tibble or lazy table with one row per peptide × sample pair.
#' @export
get_counts <- function(x) {
  .check_pd(x)
  x$data_long
}

#' @title Retrieve the comparisons definition table
#'
#' @description Returns the two-way comparison specifications stored inside a
#' **phip_data** object.
#'
#' @inheritParams get_counts
#' @return A tibble with columns: `comparison`, `group1`, `group2`,
#'   and `variable`.
#' @export
get_comparisons <- function(x) {
  .check_pd(x)
  x$comparisons
}

#' @title Report the storage backend in use
#'
#' @description Tells you whether the counts are held in `"duckdb"`, `"arrow"`,
#' or `"memory"`.
#'
#' @inheritParams get_counts
#' @return A single string (`"duckdb"`, `"arrow"`, or `"memory"`).
#' @export
get_backend <- function(x) {
  .check_pd(x)
  x$backend
}

#' @title Retrieve the metadata list
#'
#' @description Accesses the `meta` slot, which holds flags such as whether the
#'  table is a full peptide × sample grid, the available outcome columns, etc.
#'
#' @inheritParams get_counts
#' @return A named list.
#' @export
get_meta <- function(x) {
  .check_pd(x)
  x$meta
}

#' @title Retrieve the peptide-library annotation table
#'
#' @description Quick accessor for the `peptide_library` slot of a **phip_data**
#' object.
#'
#' @inheritParams get_counts
#' @return A tibble with one row per `peptide_id` and associated annotation.
#' @export
get_peptide_library <- function(x) {
  .check_pd(x)
  x$peptide_library
}

################################################################################
## dplyr verb wrappers  --------------------------------------------------------
## (dplyr already provides the generics, we only add methods)
################################################################################

# helper to copy x, modify counts, and return
.modify_pd <- function(.data, new_counts) {
  .data$data_long <- new_counts
  .data
}

#' @importFrom dplyr filter
#' @exportS3Method filter phip_data
filter.phip_data <- function(.data, ..., .preserve = FALSE) {
  .modify_pd(
    .data,
    dplyr::filter(.data$data_long, ..., .preserve = .preserve)
  )
}

#' @importFrom dplyr select
#' @exportS3Method select phip_data
select.phip_data <- function(.data, ...) {
  .modify_pd(.data, dplyr::select(.data$data_long, ...))
}

#' @importFrom dplyr mutate
#' @exportS3Method mutate phip_data
mutate.phip_data <- function(.data, ...) {
  .modify_pd(.data, dplyr::mutate(.data$data_long, ...))
}

#' @importFrom dplyr arrange
#' @exportS3Method arrange phip_data
arrange.phip_data <- function(.data, ...) {
  .modify_pd(.data, dplyr::arrange(.data$data_long, ...))
}

#' @importFrom dplyr summarise
#' @exportS3Method summarise phip_data
summarise.phip_data <- function(.data, ..., .groups = NULL) {
  dplyr::summarise(.data$data_long, ..., .groups = .groups)
}

#' @importFrom dplyr collect
#' @exportS3Method collect phip_data
collect.phip_data <- function(x, ...) {
  dplyr::collect(x$data_long, ...)
}

#' @importFrom dplyr group_by
#' @exportS3Method group_by phip_data
group_by.phip_data <- function(.data, ..., .add = FALSE,
                               .drop = dplyr::group_by_drop_default(.data$data_long)) {
  .modify_pd(
    .data,
    dplyr::group_by(.data$data_long, ..., .add = .add, .drop = .drop)
  )
}

#' @importFrom dplyr distinct
#' @exportS3Method distinct phip_data
distinct.phip_data <- function(.data, ..., .keep_all = FALSE) {
  dplyr::distinct(.data$data_long, ..., .keep_all=.keep_all)
}

#' @importFrom dplyr ungroup
#' @exportS3Method ungroup phip_data
ungroup.phip_data <- function(x, ...) {
  .modify_pd(x, dplyr::ungroup(x$data_long, ...))
}

################################################################################
## merge() + join wrappers for phip_data  --------------------------------------
################################################################################

# Utilities ---------------------------------------------------------------

.extract_data_long <- function(obj) {
  if (inherits(obj, "phip_data")) obj$data_long else obj
}

###############################################################################
##  Pretty, CLI-styled yes/no prompt  --------------------------------------
###############################################################################
.cli_yesno <- function(question,
                       yes = c("y", "yes"),
                       no = c("n", "no")) {
  yes <- tolower(yes)
  no <- tolower(no)

  repeat {
    cli::cli_text(
      "{.question {question}} {.dim [y/n]} ",
      .envir = list(question = question)
    )
    answer <- tolower(trimws(readline(cli::style_dim("→ "))))
    if (answer %in% yes) {
      return(TRUE)
    }
    if (answer %in% no) {
      return(FALSE)
    }
    cli::cli_alert_danger("Please answer {.strong y} or {.strong n}.")
  }
}

###############################################################################
##  CLI-aware warning helper  ----------------------------------------------
###############################################################################
.cli_warn <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_alert_warning(msg)
  } else {
    warning(msg, call. = FALSE)
  }
}

# Main merge() method -----------------------------------------------------

#' Merge or join a `phip_data` object
#'
#' @param x      A `phip_data` object.
#' @param y      A data-frame–like object *or* another `phip_data`.
#' @param ...    Arguments forwarded to either [base::merge()] or the chosen
#'               **dplyr** join (e.g. `by =`, `suffix =`, etc.).
#' @param confirm Logical.  When `TRUE` *and* `type = "base"` *and* the session
#'               is interactive, the user is asked to confirm.  Set to `FALSE`
#'               to skip the prompt (use sparingly—OOM risk remains).
#'
#' @return A new `phip_data` whose `data_long` contains the merged / joined
#'         tibble.
#' @exportS3Method merge phip_data
merge.phip_data <- function(x, y,
                            ...,
                            confirm = interactive()) {
  y <- .extract_data_long(y)

  # -----------------------------------------------------------------------
  #  Base merge (potentially memory-hungry) ------------------------------
  # -----------------------------------------------------------------------
  if (confirm && interactive()) {
    .cli_warn("`merge()` copies both tables in full; this may exhaust RAM.")
    ok <- .cli_yesno("Proceed with base::merge()?")
    if (!isTRUE(ok)) {
      chk::abort_chk("Merge aborted.  Use `dplyr` (or another join) for
             a memory-efficient alternative.",
        call. = FALSE
      )
    }
  }

  merged_tbl <- base::merge(x$data_long, y, ...)

  .modify_pd(x, merged_tbl)
}

#' dplyr joins for `phip_data`
#'
#' Cienkie wrappery: stosują `dplyr::<join>` do `x$data_long`
#' po uprzednim rozwinięciu `y` (jeśli to `phip_data`), a wynik
#' zawijają z powrotem przez `.modify_pd()`.
#'
#' @param x A `phip_data` object.
#' @param y A `phip_data` or a data frame / tbl.
#' @param ... Passed to the corresponding `dplyr::<join>` function.
#' @return A `phip_data` object with updated `data_long`.
#'
#' @importFrom dplyr left_join right_join inner_join full_join semi_join anti_join
# helpers (internal)
.extract_data_long <- function(y) if (inherits(y, "phip_data")) y$data_long else y

#' @export
#' @method left_join phip_data
left_join.phip_data <- function(x, y, ...) {
  y <- .extract_data_long(y)
  .modify_pd(x, dplyr::left_join(x$data_long, y, ...))
}

#' @export
#' @method right_join phip_data
right_join.phip_data <- function(x, y, ...) {
  y <- .extract_data_long(y)
  .modify_pd(x, dplyr::right_join(x$data_long, y, ...))
}

#' @export
#' @method inner_join phip_data
inner_join.phip_data <- function(x, y, ...) {
  y <- .extract_data_long(y)
  .modify_pd(x, dplyr::inner_join(x$data_long, y, ...))
}

#' @export
#' @method full_join phip_data
full_join.phip_data <- function(x, y, ...) {
  y <- .extract_data_long(y)
  .modify_pd(x, dplyr::full_join(x$data_long, y, ...))
}

#' @export
#' @method semi_join phip_data
semi_join.phip_data <- function(x, y, ...) {
  y <- .extract_data_long(y)
  .modify_pd(x, dplyr::semi_join(x$data_long, y, ...))
}

#' @export
#' @method anti_join phip_data
anti_join.phip_data <- function(x, y, ...) {
  y <- .extract_data_long(y)
  .modify_pd(x, dplyr::anti_join(x$data_long, y, ...))
}

###############################################################################
##  Pretty, CLI-styled yes/no prompt  --------------------------------------
###############################################################################

#' @title Ensure an existence flag (all ones) on `data_long`
#'
#' @description Appends/overwrites a column (default: "exist") filled with 1L on
#'   the lazy `data_long` backend. Preserves laziness; no collection is forced.
#'
#' @param phip_data A <phip_data> object.
#' @param exist_col Name of the existence column to append/overwrite.
#' @param overwrite If FALSE and the column exists, abort with a phiper-style error.
#' @return Modified <phip_data> with updated `data_long`.
#' @examples
#' pd <- phip_load_example_data()
#' pd <- add_exist(pd, overwrite = TRUE) # overwrites if present
#' @export add_exist
add_exist <- function(phip_data,
                      exist_col = "exist",
                      overwrite = FALSE) {
  x <- phip_data
  stopifnot(inherits(x, "phip_data"))
  .data <- rlang::.data

  .ph_with_timing(
    headline = "Ensuring existence flag on data_long",
    step = sprintf(
      "column: %s; overwrite: %s",
      add_quotes(exist_col, 1L), as.character(overwrite)
    ),
    expr = {
      tbl <- x$data_long

      if (exist_col %in% colnames(tbl)) {
        if (!isTRUE(overwrite)) {
          .ph_abort(
            headline = "Existence column already present.",
            step = "input validation",
            bullets = c(
              sprintf("column: %s", add_quotes(exist_col, 2L)),
              "set overwrite = TRUE to replace existing values"
            )
          )
        } else {
          .ph_warn(
            headline = "Overwriting existing existence flag.",
            step     = "adding existence indicator",
            bullets  = sprintf("column: %s", add_quotes(exist_col, 2L))
          )
        }
      } else {
        .ph_log_info(
          "Adding existence flag column",
          bullets = sprintf("column: %s", add_quotes(exist_col, 2L))
        )
      }

      # lazy mutate; stays in DuckDB/Arrow without materialising
      tbl_new <- dplyr::mutate(tbl, !!rlang::sym(exist_col) := 1L)

      x_new <- .modify_pd(x, tbl_new)
      # mark availability of an existence flag; don't touch full_cross/prop here
      x_new$meta$exist <- TRUE
      x_new
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

################################################################################
## helper to close duckdb connection -------------------------------
################################################################################

#' @title Disconnect backend database connections
#'
#' @description Closes any open database connections held inside a
#' `phip_data` object-both the `data_long` backend (e.g. DuckDB,
#' Arrow) and the optional `peptide_library` connection.
#'
#' Call this when you want to release resources explicitly; it is
#' also safe to rely on garbage collection to close the connections
#' automatically at the end of the R session.
#'
#' @param x A valid `phip_data` object.
#'
#' @return The input `phip_data` object, invisibly, with its
#'   connections closed.
#'
#' @export

disconnect <- function(x) {
  # close data_long backend, if any
  if (!is.null(x$meta$con) && DBI::dbIsValid(x$meta$con)) {
    DBI::dbDisconnect(x$meta$con, shutdown = TRUE)
  }
  # close peptide library connection
  if (!is.null(x$meta$peptide_con) && DBI::dbIsValid(x$meta$peptide_con)) {
    DBI::dbDisconnect(x$meta$peptide_con, shutdown = TRUE)
  }

  invisible(x)
}
