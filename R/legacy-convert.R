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
