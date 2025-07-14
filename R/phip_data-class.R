################################################################################
## phip_data S3 class  ---------------------------------------------------------
################################################################################

#' phip_data constructor function
#' @export
new_phip_data <- function(data_long,
                          comparisons,
                          backend = c("memory", "duckdb", "arrow"),
                          meta = list()) {

  backend <- match.arg(backend)

  # quick sanity check
  .chk_cond(!inherits(comparisons, "data.frame"),
            "The `comparisons` table has to be a data.frame!")

  # --------------------------------------------------------------------------
  # Scan column names for automatic meta flags
  # --------------------------------------------------------------------------
  # colnames() works for tibble, tbl_dbi, and arrow_dplyr_query
  cols <- colnames(data_long)

  meta$longitudinal <- "timepoint"  %in% cols
  meta$fold_change  <- "fold_change" %in% cols

  # --------------------------------------------------------------------------

  # define the object
  structure(
    list(
      data_long   = data_long,              # lazy tbl or tibble
      comparisons = tibble::as_tibble(comparisons),
      backend     = backend,
      meta        = meta
    ),
    class = "phip_data"
  )
}


################################################################################
## pretty printing for phip_data  ----------------------------------------------
################################################################################

#' @export
print.phip_data <- function(x, ...) {
  cat(cli::rule(left = "<phip_data>", right = paste0("backend: ", x$backend)),
      "\n", sep = "")

  # ---- counts preview -------------------------------------------------------
  prev <- tryCatch(
    {
      if (inherits(x$data_long, "tbl_dbi")) {
        dplyr::collect(head(x$data_long, 5L))
      } else {
        head(x$data_long, 5L)
      }
    },
    error = function(e) tibble::tibble(.error = "<preview failed>")
  )

  cat(cli::col_cyan("counts (first 5 rows):"), "\n")
  print(prev, n = 5, width = Inf)
  cat("\n")

  # ---- contrasts ------------------------------------------------------------
  cat(cli::col_cyan("contrasts:"), "\n")
  print(knitr::kable(x$comparisons, format = "simple"), row.names = FALSE)
  cat("\n")

  invisible(x)
}


################################################################################
## plain accessors (cause no S3 generics) --------------------------------------
################################################################################

.check_pd <- function(obj) {
  .chk_cond(!inherits(obj, "phip_data"),
                 "`x` must be a <phip_data> object.")
}

#' @export
get_counts <- function(x) { .check_pd(x); x$data_long }
#' @export
get_comparisons <- function(x) { .check_pd(x); x$comparisons }
#' @export
get_backend <- function(x) { .check_pd(x); x$backend }
#' @export
get_meta <- function(x) { .check_pd(x); x$meta }


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
  .modify_pd(.data,
             dplyr::filter(.data$data_long, ..., .preserve = .preserve))
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
  if (!is.null(x$meta$con) && DBI::dbIsValid(x$meta$con)) {
    DBI::dbDisconnect(x$meta$con, shutdown = TRUE)
  }
  invisible(x)
}
