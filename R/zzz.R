#' phiper
#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  # Register bundled Montserrat and activate showtext rendering
  ttf <- system.file("fonts", "Montserrat-Regular.ttf", package = pkgname)
  if (nzchar(ttf) && file.exists(ttf)) {
    bold   <- system.file("fonts", "Montserrat-Bold.ttf",   package = pkgname)
    italic <- system.file("fonts", "Montserrat-Italic.ttf", package = pkgname)
    tryCatch(
      sysfonts::font_add(
        family  = "Montserrat",
        regular = ttf,
        bold    = if (nzchar(bold)   && file.exists(bold))   bold   else NULL,
        italic  = if (nzchar(italic) && file.exists(italic)) italic else NULL
      ),
      error = function(e) NULL
    )
  }
  tryCatch(showtext::showtext_auto(enable = TRUE), error = function(e) NULL)

  # Set package-wide ggplot2 defaults
  ggplot2::theme_set(theme_phip())
  options(
    ggplot2.discrete.colour = phip_palette,
    ggplot2.discrete.fill   = phip_palette
  )
}

## usethis namespace: start
#' @useDynLib phiper, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel setThreadOptions
#' @importFrom dplyr %>% select
#' @importFrom rlang .data
#' @importFrom utils head
#' @importFrom showtext showtext_auto
#' @importFrom sysfonts font_add
#' @import phiperio
## usethis namespace: end
NULL

# Re-export phiperio functions used in phiper's examples / user-facing API,
# so they are available after library(phiper) without a separate
# library(phiperio) call.

#' @title Load Example PhIP-Seq Dataset as <phip_data>
#'
#' @description Convenience helper to quickly load a shipped example dataset
#'   (`"phip_mixture"`) into a `<phip_data>` object, suitable for downstream
#'   analysis and visualization. This function wraps
#'   \code{\link[phiperio]{convert_standard}}, automatically supplying the
#'   correct parameters for the included example data.
#'
#' @param name Character scalar. Name of the shipped example dataset.
#'   Currently supported: \code{"phip_mixture"}, \code{"small_mixture"}.
#'
#' @return A `<phip_data>` object created from the chosen example dataset.
#'
#' @examples
#' # Load the example data shipped with the package:
#' ex <- load_example_data()
#' # ex is now a <phip_data> object ready for analysis
#'
#' # Specify the dataset name explicitly
#' ex2 <- load_example_data("small_mixture")
#'
#' @importFrom phiperio load_example_data
#' @export
load_example_data <- phiperio::load_example_data

#' @title Path to Example PhIP-Seq Datasets
#'
#' @description Return the path to an example dataset shipped with phiperio,
#'   suitable for use with \code{\link{load_example_data}} or
#'   \code{\link[phiperio]{convert_standard}}.
#'
#' @param name Character scalar. Name of the example dataset.
#'   Currently supported: \code{"phip_mixture"}.
#'
#' @return A character scalar with an absolute path to the file.
#'
#' @examples
#' sim_path <- get_example_path("phip_mixture")
#' # phip_obj <- convert_standard(sim_path)
#'
#' @importFrom phiperio get_example_path
#' @export
get_example_path <- phiperio::get_example_path

# There is a note for not-declared global variables in R CMD CHECK. This is a
# little bit tricky when using dplyr to manipulate the data, as you use the
# unquoted names of the variables to access the columns in the data. R CMD CHECK
# sees it as a global variable and cries, that you use it, but haven't defined
# it. It is a pretty well known bug, i bypassed it a little bit differently in
# vecmatch, but actually defining all the variables as global will also do the
# trick + you have all the vars in one place.
utils::globalVariables(c(
  "%>%",
  ".",
  "..keep_cols",
  ".abundance",
  ".env",
  ".exist",
  ".x",
  ":=",
  "H_ln",
  "N",
  "N1",
  "N2",
  "N_1",
  "N_2",
  "N_paired",
  "PAIRS",
  "POOL",
  "T_obs",
  "T_obs_stand",
  "Z_from_p",
  "b",
  "category",
  "category_rank_bh",
  "category_rank_wbh",
  "category_use",
  "cohort",
  "d1",
  "d2",
  "delta",
  "delta_jit",
  "design",
  "distance",
  "dummy",
  "example",
  "feature",
  "frac_delta_pos",
  "frac_delta_pos_w",
  "group",
  "group1",
  "group2",
  "group_col",
  "group_value",
  "group_value_1",
  "group_value_2",
  "idx_hits_g1",
  "idx_hits_g2",
  "k",
  "k_levels",
  "log2ratio",
  "m_eff",
  "m_r",
  "max_delta",
  "n",
  "n1",
  "n1_paired",
  "n1_tot",
  "n1_val",
  "n2",
  "n2_paired",
  "n2_tot",
  "n2_val",
  "n_dups",
  "n_peptides",
  "n_peptides_used",
  "n_pos",
  "n_present",
  "n_present_1",
  "n_present_2",
  "n_subjects_paired",
  "nlog10p",
  "p1_val",
  "p2_val",
  "padj_bh",
  "p_adj_rank",
  "p_adj_rank_wbh",
  "p_perm",
  "p_raw",
  "passed_rank_bh",
  "passed_rank_wbh",
  "pct1val",
  "pct2val",
  "pep_key",
  "peptide_id",
  "percent1",
  "percent2",
  "percent_1",
  "percent_2",
  "pooled",
  "pooled_jit",
  "present",
  "present_1",
  "present_2",
  "prop",
  "prop1",
  "prop1_eps",
  "prop1_paired",
  "prop2",
  "prop2_eps",
  "prop2_paired",
  "prop_1",
  "prop_2",
  "rank_val",
  "ratio",
  "richness",
  "sample_id",
  "score",
  "shannon",
  "simpson",
  "subject_id",
  "time_cont",
  "v",
  "view",
  "x",
  "y",
  "T_null_mean",
  "T_null_sd",
  "cohens_d",
  "diagonal",
  "delta_ratio",
  "fill_d",
  "label",
  "n01",
  "n10",
  "p_adj",
  "p_value",
  "sig",
  "stars",
  "test"
))
