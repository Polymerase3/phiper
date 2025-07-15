################################################################################
## phip_data S3 class  ---------------------------------------------------------
################################################################################

#' phip_data constructor function
#' @export
new_phip_data <- function(data_long,
                          comparisons,
                          backend = c("memory", "duckdb", "arrow"),
                          peptide_library = NULL,
                          meta = list()
                          ) {
  backend <- match.arg(backend)

  # quick sanity check
  if (!is.null(comparisons)) {
    .chk_cond(
      !inherits(comparisons, "data.frame"),
      "`comparisons` must be a data.frame (or NULL)."
    )
  } else {
    comparisons <- tibble::tibble(         # empty placeholder
      comparison = character(),
      group1     = character(),
      group2     = character(),
      variable   = character()
    )
  }

  # --------------------------------------------------------------------------
  # Download the peptide metadata library
  # --------------------------------------------------------------------------
  peptide_library <- get_peptide_meta()

  # --------------------------------------------------------------------------
  # Scan column names for automatic meta flags
  # --------------------------------------------------------------------------
  # colnames() works for tibble, tbl_dbi, and arrow_dplyr_query
  cols <- colnames(data_long)
  standard_cols <- c("subject_id", "sample_id", "timepoint", "peptide_id",
                     "present", "fold_change", "counts_control", "counts_hit")

  meta$longitudinal <- all(c("timepoint", "sample_id") %in% cols)
  meta$exist        <- "present" %in% cols
  meta$fold_change  <- "fold_change" %in% cols
  meta$raw_counts   <- all(c("counts_control", "counts_hit") %in% cols)
  meta$extra_cols   <- cols[cols %nin% standard_cols]
  meta$peptide_con  <- attr(peptide_library, "duckdb_con")

  # check if the data is a full_cross --> that means all peptide x unique_id
  # combinations are present in the data
  is_full_cross <- function(tbl,
                            peptide = "peptide_id",
                            sample  = "sample_id") {

    pep_sym <- rlang::sym(peptide)
    smp_sym <- rlang::sym(sample)

    dims <- tbl %>%
      dplyr::summarise(
        n_pep = dplyr::n_distinct(!!pep_sym),
        n_smp = dplyr::n_distinct(!!smp_sym),
        n_obs = dplyr::n()
      ) %>%
      dplyr::collect()

    dims$n_obs == dims$n_pep * dims$n_smp
  }

  # pick the correct column ----------------------------------------------
  sample_dim <- if (meta$longitudinal) "sample_id" else "subject_id"

  meta$full_cross <- is_full_cross(data_long,
                                   peptide = "peptide_id",
                                   sample  = sample_dim)

  # --------------------------------------------------------------------------

  # define the object
  obj <- structure(
    list(
      data_long       = data_long, # lazy tbl or tibble
      comparisons     = tibble::as_tibble(comparisons),
      backend         = backend,
      peptide_library = peptide_library,
      meta            = meta
    ),
    class = "phip_data"
  )

  # --------------------------------------------------------------------------
  # Validate the objects used to construct the phip_data
  # --------------------------------------------------------------------------
  validate_phip_data(obj)   # will stop on error, warn on non-fatal issues

  # return clean validated phip_data
  return(obj)
}

#' @export
print.phip_data <- function(x, ...) {
  cat(cli::rule(
    left  = "<phip_data>",
    right = paste0("backend: ", x$backend)
  ), "\n\n")

  # ---- counts preview -------------------------------------------------------
  prev <- tryCatch({
    tbl <- x$data_long
    if (inherits(tbl, c("tbl_dbi", "arrow_dplyr_query"))) {
      tbl %>% head(5) %>% dplyr::collect()
    } else {
      utils::head(tbl, 5)
    }
  }, error = function(e) tibble::tibble(.error = "<preview failed>"))

  cat(cli::col_cyan("counts (first 5 rows):"), "\n")
  print(prev)

  ## --- size summary ----------------------------------------------------------
  n_cols <- ncol(x$data_long)

  n_rows <- tryCatch(
    if (inherits(x$data_long, c("tbl_dbi", "arrow_dplyr_query"))) {
      x$data_long %>% dplyr::summarise(n = dplyr::n()) %>% dplyr::pull(n)
    } else {
      nrow(x$data_long)
    },
    error = function(e) NA_integer_
  )
  cat("\n")
  cat(sprintf("table size: %s rows × %s columns\n\n",
              format(n_rows, big.mark = ","), n_cols))

  # ---- contrasts ------------------------------------------------------------
  cat(cli::col_cyan("contrasts:"), "\n")
  cat(paste0(knitr::kable(x$comparisons, format = "simple"), collapse = "\n"))
  cat("\n\n")

  # ---- peptide library preview ---------------------------------------------
  cat(cli::col_cyan("peptide library preview (first 5 rows):"), "\n")

  lib <- x$peptide_library

  # If the DuckDB connection was closed, reopen the cached table
  if (inherits(lib, "tbl_dbi")) {
    if (!DBI::dbIsValid(lib$src$con)) {
      lib <- get_peptide_meta()
      x$peptide_library <- lib    # keep the fresh handle for later prints
    }
  }

  show_cols <- intersect(c("peptide_id", "pos", "len_seq"), colnames(lib))

  lib_preview <- tryCatch({
    lib %>%
      head(5) %>%                # LIMIT 5
      dplyr::collect() %>%        # into R
      dplyr::select(dplyr::all_of(show_cols))
  }, error = function(e) tibble::tibble(.error = "<preview failed>"))

  print(lib_preview)

  # summary: extra-column count + overall dim ---------------------------------
  n_cols   <- length(colnames(lib))
  extra_nc <- n_cols - length(show_cols)

  n_rows <- tryCatch(
    lib %>% dplyr::summarise(n = dplyr::n()) %>% dplyr::pull(n),
    error = function(e) NA_integer_
  )

  cat(sprintf("... plus %d more columns\n\n", extra_nc))
  cat(sprintf("library size: %s rows × %s columns\n\n",
              format(n_rows, big.mark = ","), n_cols))


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

#' @export
get_counts <- function(x) {
  .check_pd(x)
  x$data_long
}
#' @export
get_comparisons <- function(x) {
  .check_pd(x)
  x$comparisons
}
#' @export
get_backend <- function(x) {
  .check_pd(x)
  x$backend
}
#' @export
get_meta <- function(x) {
  .check_pd(x)
  x$meta
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

################################################################################
## helper to close duckdb connection -------------------------------
################################################################################

#' Disconnect a phip_data backend (if any)
#' @export
disconnect.phip_data <- function(x, ...) {
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

################################################################################
## validator function for phip_data -------------------------------
################################################################################
# helper: verify we got a phip_data object
.check_pd <- function(obj) {
  .chk_cond(
    !inherits(obj, "phip_data"),
    "`x` must be a <phip_data> object."
  )
}

#' Internal validator for <phip_data> objects
#' @keywords internal
validate_phip_data <- function(x,
                               na_warn_thresh = 0.50,   # warn if >50 % NA / zero
                               auto_expand   = TRUE) {  # fill missing grid?
  .check_pd(x)                    # existing helper

  tbl      <- x$data_long
  cols     <- colnames(tbl)
  reserved <- c("subject_id","sample_id","timepoint",
                "peptide_id","present","fold_change",
                "counts_control","counts_hit")

  ## ---------------------------------------------- 1  STRUCTURE
  is_long  <- all(c("subject_id","timepoint") %in% cols)
  need     <- if (is_long) {
    c("subject_id","timepoint","sample_id","peptide_id")
  } else {
    c("subject_id","peptide_id")
  }

  miss <- setdiff(need, cols)
  .chk_cond(length(miss) > 0,
            sprintf("Missing mandatory column(s): %s",
                    paste(miss, collapse =", ")))

  ## ---------------------------------------------- 2  OUTCOME FAMILY
  have <- list(
    present     = "present" %in% cols,
    fold_change = "fold_change" %in% cols,
    raw_counts  = all(c("counts_control","counts_hit") %in% cols)
  )
  k <- sum(unlist(have))
  .chk_cond(k == 0,
            "No outcome column found (need present /
            fold_change / raw_counts).")

  ## --------------------------------------------- 3  RESERVED-NAME COLLISIONS
  overlap <- intersect(x$meta$extra_cols, reserved)
  .chk_cond(length(overlap) > 0,
            sprintf("extra_cols overlap with reserved names: %s",
                    paste(overlap, collapse=", ")))

  ## --------------------------------------------- 4  ATOMIC COLUMNS
  sample0 <- if (inherits(tbl, "tbl_lazy")) {
    # LIMIT 0 → fetch only the schema, then drop Arrow classes
    tbl %>% head(0) %>% dplyr::collect() %>% as.data.frame()
  } else {
    as.data.frame(tbl)
  }

  bad_atomic <- names(sample0)[vapply(sample0, is.list, logical(1))]

  .chk_cond(
    length(bad_atomic) > 0,
    sprintf("Non-atomic (list) columns found: %s",
            paste(bad_atomic, collapse = ", "))
  )

  ## --------------------------------------------- 5  UNIQUENESS
  if (is_long) {
    dup <- tbl %>%
      dplyr::count(subject_id, timepoint, peptide_id) %>%
      dplyr::filter(n > 1) %>% head(1) %>% dplyr::collect()
    .chk_cond(nrow(dup) > 0,
              "Each (subject_id, peptide_id, timepoint) must
              map to exactly one sample_id.")
  } else {
    dup <- tbl %>%
      dplyr::count(subject_id, peptide_id) %>%
      dplyr::filter(n > 1) %>% head(1) %>% dplyr::collect()
    .chk_cond(nrow(dup) > 0, "Duplicate sample_id values found.")
  }

  ## ---------------------------------------------- 6 VALUE RANGES
  if (have$present) {
    bad <- tbl %>%
      dplyr::filter(!is.na(present) & !(present %in% c(0, 1))) %>%
      head(1) %>% dplyr::collect()
    .chk_cond(nrow(bad) > 0, "`present` must be 0/1/NA.")
  }

  if (have$fold_change) {
    ok <- tbl %>%
      dplyr::summarise(all_finite = all(is.finite(fold_change))) %>%
      dplyr::pull(all_finite)
    .chk_cond(!ok, "`fold_change` contains Inf/-Inf or NA.")
  }

  if (have$raw_counts) {
    neg <- tbl %>%
      dplyr::filter(counts_control < 0 | counts_hit < 0) %>%
      head(1) %>% dplyr::collect()
    .chk_cond(nrow(neg) > 0, "Raw counts must be non-negative.")
  }

  ## ---------------------------------------------- 7  SPARSITY WARNING
  if (have$present) {
    prop_na <- tbl %>%
      dplyr::summarise(
        p = sum(dplyr::if_else(is.na(present), 1, 0)) / dplyr::n()
      ) %>%
      dplyr::pull(p)

    if (prop_na > na_warn_thresh)
      cli::cli_warn(sprintf("present column is %.0f %% NA.", prop_na * 100))
  }

  if (have$raw_counts) {
    prop_zero <- tbl %>%
      dplyr::summarise(
        p = sum(dplyr::if_else(counts_control == 0 & counts_hit == 0, 1, 0)) /
          dplyr::n()
      ) %>%
      dplyr::pull(p)

    if (prop_zero > na_warn_thresh)
      cli::cli_warn(sprintf("Both count columns are 0 in %.0f %% of rows.",
                            prop_zero * 100))
  }

  ## ---------------------------------------------- 8  PEPTIDE-ID COVERAGE
  missing_in_lib <- setdiff(
    tbl %>% dplyr::distinct(peptide_id) %>% dplyr::pull(),
    x$peptide_library %>% dplyr::distinct(peptide_id) %>% dplyr::pull()
  )
  missing_in_lib <- missing_in_lib[order(missing_in_lib)]

  .chk_cond(
    length(missing_in_lib) > 0,
    sprintf("peptide_id not found in peptide_library (e.g. %s)",
            missing_in_lib[1]),
    error = FALSE                    # emit warning instead of abort
  )

  ## ---------------------------------------------- 9  COMPARISONS TABLE
  cmp <- x$comparisons
  if (nrow(cmp) > 0) {
    allowed_cmp_cols <- c("comparison","group1","group2","variable")
    extra_cmp <- setdiff(colnames(cmp), allowed_cmp_cols)
    .chk_cond(length(extra_cmp) > 0,
              sprintf("Unexpected columns in comparisons: %s",
                      paste(extra_cmp, collapse = ", ")))

    valid_labels <- if (is_long) {
      # which of the two exist?
      valid_cols <- dplyr::intersect(c("timepoint", "group"),
                              colnames(tbl))

      tbl %>% # keep only the existing ones
        dplyr::select(dplyr::all_of(valid_cols)) %>%
        dplyr::distinct() %>%           # unique row-wise combinations
        dplyr::collect() %>%            # pull the tiny result into R
        unlist(use.names = FALSE) %>%   # flatten to one vector
        unique()                        # final de-dup

    } else {
      tbl %>% dplyr::distinct(group) %>% dplyr::pull(group)
    }

    bad_cmp <- unique(c(setdiff(cmp$group1, valid_labels),
                        setdiff(cmp$group2, valid_labels)))
    .chk_cond(length(bad_cmp) > 0,
              sprintf("comparisons refer to unknown group labels (e.g. %s)",
                      bad_cmp[1]))
  }

  ## ---------------------------------------------- 10 FULL GRID
  # decide which column indexes samples / subjects
  # decide which column indexes samples / subjects
  sample_dim <- if (is_long) "sample_id" else "subject_id"
  sample_sym <- rlang::sym(sample_dim)   # <- turn into a symbol

  dims <- tbl %>%
    dplyr::summarise(
      n_pep = dplyr::n_distinct(peptide_id),
      n_smp = dplyr::n_distinct(!!sample_sym),   # tidy-eval here
      n_obs = dplyr::n()
    ) %>%
    dplyr::collect()

  expect <- dplyr::pull(dims, n_pep) * dplyr::pull(dims, n_smp)
  if (dplyr::pull(dims, n_obs) != expect) {
    cli::cli_warn(sprintf(
      "Counts table is not a full peptide × sample grid (%s / %s rows).",
      dplyr::pull(dims, n_obs), expect))
    if (isTRUE(auto_expand)) {
      full_grid <- tidyr::crossing(
        sample_id  = tbl %>% dplyr::distinct(sample_id) %>% dplyr::collect() %>% dplyr::pull(),
        peptide_id = tbl %>% dplyr::distinct(peptide_id) %>% dplyr::collect() %>% dplyr::pull()
      )
      tbl <- full_grid %>% dplyr::left_join(tbl,
                                            by = c("sample_id","peptide_id"))
      x$data_long       <- tbl
      x$meta$full_cross <- FALSE
    }
  } else {
    x$meta$full_cross <- TRUE
  }

  invisible(x)
}
