#' @title Convert raw PhIP-Seq output into a `phip_data` object
#'
#' @description `phip_convert()` ingests a "long" table of PhIPsSeq read counts /
#' enrichment statistics, optionally expands it to the full
#' `sample_id x peptide_id` grid, and registers the result in DuckDB.
#' The function returns a fully initialised **`phip_data`** object that can be
#' queried with the tidy API used throughout the package.
#'
#' @param data_long_path Character scalar. File or directory containing the
#'   *long-format* PhIP-Seq data. Allowed extensions are **`.csv`** and
#'   **`.parquet`**. Directories are treated as partitions of a parquet set.
#' @param sample_id,peptide_id,subject_id,timepoint,exist,fold_change,counts_input,counts_hit
#' Optional character strings. Supply these only if your column names differ
#' from the defaults (`"sample_id"`, `"peptide_id"`, `"subject_id"`,
#' `"timepoint"`, `"exist"`, `"fold_change"`, `"counts_input"`,
#' `"counts_hit"`). Each argument should contain the *name* of the column in the
#' incoming data; `NULL` lets the default stand.
#'
#' @param n_cores Integer >= 1. Number of CPU threads DuckDB may use while
#'   reading and writing files.
#'
#' @param materialise_table Logical. If `FALSE` the result is registered as a
#'   **view**; if `TRUE` the table is fully **materialised** and stored on disk,
#'   trading higher load time and storage for faster repeated queries.
#'
#' @param auto_expand Logical. If `TRUE` and the incoming data are **not** a
#'   complete Cartesian product of `sample_id x peptide_id`, missing
#'   combinations are generated:
#'   * Columns that are constant within each `sample_id` (metadata) are copied
#'     to the new rows.
#'   * Non-recyclable measurement columns (`fold_change`, `exist`,
#'     `counts_input`, `counts_hit`, etc.) are initialised to 0.
#'   The expanded table replaces the original *in place*.
#'
#' @param peptide_library Logical. If `TRUE` (default) `phip_convert()` will
#'   attempt to locate and attach the matching peptide-library metadata for
#'   downstream annotation. Set to `FALSE` to skip this step.
#'
#' @return An S3 object of class **`phip_data`** containing:
#' \describe{
#'   \item{`data_long`}{The (possibly expanded) long-format table.}
#'   \item{`comparisons`}{A tibble of pre-computed group comparisons or
#'     `NULL` if none were supplied.}
#'   \item{`peptide_library`}{Loaded peptide-library metadata (if
#'     `peptide_library = TRUE`).}
#'   \item{`backend`}{Character vector indicating the storage engine.}
#'   \item{`meta`}{List with DuckDB connection handles.}
#' }
#'
#' @details
#' *Paths are resolved to absolute form* before any work begins, and explicit
#' checks confirm existence as well as extension validity.
#' @examples
#' # Basic import, auto-detecting default column names
#' phip_obj <- phip_convert(
#'   data_long_path = phip_example_path("phip_mixture"),
#'   n_cores = 4,
#'   materialise_table = TRUE
#' )
#'
#' \dontrun{
#' # Import a CSV and rename columns
#' phip_mem <- phip_convert(
#'   data_long_path = "data/phip_long.csv",
#'   sample_id      = "sample",
#'   peptide_id     = "pep"
#' )
#' }
#'
#' @seealso
#' * `new_phip_data()` for the object constructor.
#' * `dplyr::tbl()` to query DuckDB tables lazily.
#'
#' @export

phip_convert <- function(
  data_long_path,
  sample_id = NULL,
  peptide_id = NULL,
  subject_id = NULL,
  timepoint = NULL,
  exist = NULL,
  fold_change = NULL,
  counts_input = NULL,
  counts_hit = NULL,
  n_cores = 8,
  materialise_table = TRUE,
  auto_expand = FALSE,
  peptide_library = TRUE
) {
  # ------------------------------------------------------------------
  # 1. arg checks
  # ------------------------------------------------------------------
  chk::chk_numeric(n_cores)
  chk::chk_flag(materialise_table)

  # ------------------------------------------------------------------
  # 2. resolving the data_long_file path to absolute
  # ------------------------------------------------------------------
  ## check if the data_long_path provided
  .chk_cond(
    missing(data_long_path) || !nzchar(data_long_path),
    "'data_long_path' must be provided and non-empty,
            no default is set."
  )

  ## check if this is a directory or a file
  info <- file.info(data_long_path)

  ## pre-check if exists
  .chk_cond(
    all(is.na(info)),
    sprintf("Path '%s' does not exist", data_long_path)
  )

  if (!info$isdir || is.na(info$isdir)) {
    .chk_path(data_long_path,
      "data_long_path",
      extension = c("csv", "parquet", "parq", "pq")
    )
  }

  ## resolve the path to data_long_path
  cfg <- .resolve_paths(
    data_long_path = data_long_path,
    backend = "duckdb",
    peptide_library = peptide_library,
    n_cores = n_cores,
    materialise_table = materialise_table,
    auto_expand = auto_expand
  )

  ## filter the NULLs
  cfg <- Filter(Negate(is.null), cfg)

  # ------------------------------------------------------------------
  # 3. mapping the column names
  # ------------------------------------------------------------------
  colname_map <- list(
    sample_id    = sample_id %||% "sample_id",
    peptide_id   = peptide_id %||% "peptide_id",
    subject_id   = subject_id %||% "subject_id",
    timepoint    = timepoint %||% "timepoint",
    exist        = exist %||% "exist",
    fold_change  = fold_change %||% "fold_change",
    counts_input = counts_input %||% "counts_input",
    counts_hit   = counts_hit %||% "counts_hit"
  )

  # ------------------------------------------------------------------
  # 4. create the phip_data object (DuckDB)
  # ------------------------------------------------------------------
  # using the helper to keep the backend setup consistent
  con <- .standard_read_duckdb_backend(cfg, colname_map)

  ## duckdb-specific code
  long <- dplyr::tbl(con, "raw_combined")

  # returning the phip_data object
  new_phip_data(
    data_long = long,
    comparisons = NULL,
    peptide_library = cfg$peptide_library,
    backend = "duckdb",
    auto_expand = cfg$auto_expand,
    materialise_table = cfg$materialise_table,
    meta = list(con = con)
  )
}
