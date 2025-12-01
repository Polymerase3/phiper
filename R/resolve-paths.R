#' @title Resolve legacy-import paths and perform fast-fail argument checks
#'
#' @description Combines explicit arguments with a YAML config (if given),
#'   expands every relative path to an absolute path (relative paths are
#'   evaluated against `dirname(config_yaml)` (!!!) when YAML is used, otherwise
#'   against the directory that contains the first supplied data matrix (!!!)),
#'   and returns a fully populated list of file locations and options ready for
#'   backend builders. Only cheap, load-blocking checks are done here:
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
#' @param exist_file,fold_change_file,input_file,hit_file,samples_file,timepoints_file,comparisons_file
#'  Character paths (relative or absolute) to the respective CSV/Parquet inputs.
#'  `NULL` means "not supplied".
#' @param extra_cols Character vector of extra metadata columns to keep; may be
#'   `NULL`.
#' @param output_dir Ignored (soft-deprecated).
#' @param backend Desired backend string (`"duckdb"`, `"arrow"`, `"memory"`);
#'   may come from YAML or explicit arg.
#' @param config_yaml Optional path to a YAML file whose keys mirror the
#'   function arguments; relative paths inside the YAML are resolved against the
#'   YAML’s own directory.
#'
#' @return A named list with absolute paths, `backend`, `extra_cols`, and
#'   `base_dir`; suitable for downstream helper functions.
#'
#' @keywords internal

.resolve_paths <- function(
  exist_file = NULL,
  fold_change_file = NULL,
  samples_file = NULL,
  input_file = NULL,
  hit_file = NULL,
  timepoints_file = NULL,
  extra_cols = NULL,
  comparisons_file = NULL,
  output_dir = NULL, # deprecated
  data_long_path = NULL,
  backend = NULL,
  peptide_library = TRUE,
  n_cores = NULL,
  materialise_table = NULL,
  auto_expand = NULL,
  config_yaml = NULL
) {
  ## ------------------------------------------------------------------------ ##
  ## 1.  locate base directory & read yaml (if any provided)                  ##
  ## ------------------------------------------------------------------------ ##
  # Determine base_dir depending on which input source was provided
  base_dir <- if (!is.null(config_yaml)) {
    # 1) If a YAML config was given, take its folder
    fs::path_dir(fs::path_abs(config_yaml))
  } else if (!is.null(data_long_path)) {
    # 2) If a data_long_path directory was given - use it
    fs::path_dir(fs::path_abs(data_long_path))
  } else {
    # 3) Otherwise require at least a samples_file (or exist_file)
    .chk_cond(
      is.null(samples_file) && is.null(exist_file),
      "When neither 'config_yaml' nor 'data_long_path' is provided,
      you must supply 'samples_file' or 'exist_file'."
    )
    # pick samples_file if present, else exist_file, then take its parent folder
    fs::path_dir(fs::path_abs(samples_file %||% exist_file))
  }

  yaml_cfg <- if (!is.null(config_yaml)) {
    ## validate the file extension
    .chk_path(config_yaml, "config_yaml", c("yml", "yaml"))

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
    .chk_cond(
      is.null(val) && !optional,
      sprintf("Missing required argument '%s' in YAML or call.", key)
    )

    # if the validator is .chk_path or absolute == TRUE, expand the path
    # to absolute
    if ((!is.null(val) && identical(validate, .chk_path)) ||
      (!is.null(val) && absolutize)) {
      if (!fs::is_absolute_path(val)) {
        val <- basename(val)
        val <- fs::path_abs(val, start = base_dir)
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
      .chk_path,
      optional = TRUE,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    fold_change_file = fetch(fold_change_file,
      "fold_change_file",
      .chk_path,
      optional = TRUE,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    input_file = fetch(input_file,
      "input_file",
      .chk_path,
      optional = TRUE,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    hit_file = fetch(hit_file,
      "hit_file",
      .chk_path,
      optional = TRUE,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    samples_file = fetch(samples_file,
      "samples_file",
      .chk_path,
      optional = samples_required,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    timepoints_file = fetch(timepoints_file,
      "timepoints_file",
      .chk_path,
      optional = TRUE,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    comparisons_file = fetch(comparisons_file,
      "comparisons_file",
      .chk_path,
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
    backend = fetch(
      backend,
      "backend",
      chk::chk_string
    ),
    peptide_library = peptide_library,
    n_cores = n_cores,
    materialise_table = materialise_table,
    auto_expand = auto_expand,
    base_dir = base_dir # for downstream helpers
  )

  ## ------------------------------------------------------------------------ ##
  ## 4.  fast-fail rules that really must hold before heavy work              ##
  ## ------------------------------------------------------------------------ ##
  #  rule 1: input_file and hit_file must be provided together ----------
  .chk_cond(
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
    .chk_cond(
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
    .chk_cond(
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
  .chk_cond(!is.null(cfg$output_dir),
    error = FALSE,
    "'output_dir' is deprecated and will be ignored."
  )

  # validate the tables itself --> the logic has been moved entirely to the
  # phip_data class validator

  cfg
}

#' @title Path to example PhIP-Seq datasets shipped with phiper
#'
#' @param name Character scalar. Name of the example dataset.
#'   Currently supported: `"phip_mixture"`.
#'
#' @return A character scalar with an absolute path to the file.
#'
#' @examples
#' sim_path <- phip_example_path("phip_mixture")
#' # phip_obj <- phip_convert(sim_path)
#'
#' @export
phip_example_path <- function(name = c("phip_mixture")) {
  name <- match.arg(name)
  fname <- switch(
    name,
    phip_mixture = "phip_mixture.parquet"
  )

  path <- system.file("extdata", fname, package = "phiper")
  if (path == "") {
    stop("File ", fname, " not found in extdata/", call. = FALSE)
  }
  path
}

#' @title Load Example PhIP-Seq Dataset as <phip_data>
#'
#' @description
#' Convenience helper to quickly load a shipped example dataset ("phip_mixture") into a `<phip_data>` object,
#' suitable for downstream analysis and visualization. This function wraps \code{\link{phip_convert}},
#' automatically supplying the correct parameters for the included example data.
#'
#' @param name Character scalar. Name of the shipped example dataset. Currently supported: \code{"phip_mixture"}.
#'
#' @return A `<phip_data>` object created from the chosen example dataset.
#'
#' @examples
#' # Load the example data shipped with the package:
#' ex <- phip_load_example_data()
#' # ex is now a <phip_data> object ready for analysis
#'
#' # Specify the dataset name explicitly (currently only "phip_mixture" is available)
#' ex2 <- phip_load_example_data("phip_mixture")
#'
#' # Use with plotting functions
#' p = plot_enrichment_counts(ex, group_cols = "timepoint")
#'
#' @export
phip_load_example_data <- function(name = c("phip_mixture")) {
  name <- match.arg(name)
  phip_convert(
    data_long_path = phip_example_path(name),
    backend = "duckdb",
    peptide_library = TRUE,
    subject_id = "subject_id",
    peptide_id = "peptide_id",
    sample_id  = "sample_id",
    exist      = "exist",
    timepoint  = "time",
    fold_change= "fold_change",
    materialise_table = TRUE,
    auto_expand = FALSE,
    n_cores = 5
  )
}
