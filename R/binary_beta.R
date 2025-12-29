#' @title Compute Pairwise Sample Distances
#'
#' @description This function builds a sample-by-feature abundance matrix from a
#' \code{phip_data} object (using \code{ps$data_long}), optionally normalizes
#' the matrix, and then computes pairwise distances between samples.
#'
#' The normalized abundance matrix used for distance calculation is attached
#' to the returned \code{dist} object as attribute \code{"abundances"}.
#'
#' Note: this function collects \code{ps$data_long} into memory and pivots to a
#' wide matrix in r. This can be large for big cohorts and/or large peptide
#' sets.
#'
#' @param ps input data. either:
#'   \itemize{
#'     \item a \code{phip_data} object, in which case \code{ps$data_long} is
#'     used, or
#'     \item a \code{data_long}-like table (a data.frame or dplyr \code{tbl})
#'       containing at least \code{sample_id}, \code{peptide_id}, and the column
#'       given by \code{value_col}.
#'   }
#'
#' @param value_col character scalar. Name of the abundance column in
#' \code{ps$data_long}. If \code{NULL}, the function tries (in order)
#' \code{exist},  \code{counts_hit}, \code{counts_input}, \code{fold_change}.
#'
#' @param method_normalization character scalar. normalization applied to the
#' abundance matrix before distance computation. one of:
#' \itemize{
#'   \item \code{"auto"}: uses \code{"none"} for binary (0/1) data, otherwise
#'     uses \code{"relative"}.
#'   \item \code{"relative"}: divide each row by its row sum.
#'   \item \code{"hellinger"}: \code{sqrt(relative)}.
#'   \item \code{"log"}: \code{log1p} transform.
#'   \item \code{"none"}: no normalization.
#' }
#'
#' @param distance character scalar. distance method name. the string is
#'   lowercased internally.
#'
#'   supported methods depend on which packages are installed:
#'
#'   - fast path (if package 'parallelDist' is installed):
#'     * "bray" (bray-curtis). Computed via threaded l1 distances and
#'       normalization (equivalent to bray-curtis on the normalized matrix).
#'     * "euclidean"
#'     * "minkowski"
#'     * "manhattan"
#'     * "canberra"
#'     * "binary"
#'     * "maximum" (maximum/supremum/chebyshev distance). Note: 'parallelDist'
#'       documents this as method "maximum"; passing "chebyshev" may fail unless
#'       you map it to "maximum" before calling parDist().
#'     * "cosine"
#'
#'   - fallback path (requires package 'vegan'). Any method supported by
#'   vegan::vegdist(), partial match allowed:
#'     "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
#'     "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#'     "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger",
#'     "aitchison", "robust.aitchison".
#'
#'   if 'parallelDist' is installed but the requested method is not in the fast
#'   list above, the function falls back to vegan::vegdist().
#'
#' @param n_threads integer scalar. Number of cpu threads passed to
#' \code{parallelDist::parDist(threads = ...)}.
#'
#' @return a \code{dist} object of pairwise sample distances. The attribute
#' \code{"abundances"} contains the normalized abundance matrix used for the
#' calculation (rows are samples, columns are features).
#'
#' @examples
#' \donttest{
#' # build an example <phip_data> object from the package example dataset
#' phip_path <- phip_example_path()
#'
#' ps <- phip_convert(
#'   data_long_path    = phip_path,
#'   backend           = "duckdb",
#'   peptide_library   = TRUE,
#'   subject_id        = "subject_id",
#'   peptide_id        = "peptide_id",
#'   sample_id         = "sample_id",
#'   exist             = "exist",
#'   timepoint         = "timepoint_factor",
#'   fold_change       = "fold_change",
#'   materialise_table = TRUE,
#'   auto_expand       = TRUE,
#'   n_cores           = 2
#' )
#'
#' # small subset for speed: 5 peptides at time t1
#' keep_pep <- c("16627", "5243", "24799", "16196", "18003")
#' dat_cols <- dplyr::tbl_vars(ps$data_long)
#' tp_col <- "time"
#'
#' ps_small <- ps %>%
#'   dplyr::filter(
#'     peptide_id %in% keep_pep,
#'     !!rlang::sym(tp_col) == "T1"
#'   ) %>%
#'   dplyr::collect()
#'
#' # compute distances (needs either 'parallelDist' or 'vegan')
#'   val_col <- if ("counts_hits" %in% dplyr::tbl_vars(ps_small)) {
#'     "counts_hits"
#'   } else {
#'     "exist"
#'   }
#'
#'   d <- compute_distance(
#'     ps_small,
#'     value_col = val_col,
#'     method_normalization = "hellinger",
#'     distance = "bray",
#'     n_threads = 2L
#'   )
#'
#'   a <- attr(d, "abundances")
#'   a[1:min(5, nrow(a)), 1:min(5, ncol(a)), drop = FALSE]
#' }
#' @export
compute_distance <- function(ps,
                             value_col = NULL,
                             method_normalization = c("auto", "relative",
                                                      "hellinger", "log",
                                                      "none"),
                             distance = "bray",
                             n_threads = 1L) {
  # ----------------------------------------------------------------------------
  # input validation (chk)
  # ----------------------------------------------------------------------------
  if (!is.null(value_col)) chk::chk_string(value_col)
  chk::chk_character(method_normalization)
  chk::chk_string(distance)
  chk::chk_count(n_threads)
  chk::chk_gt(n_threads, 0)

  method_normalization <- match.arg(method_normalization)
  dist_method <- tolower(distance)

  # if ps is <phip_data>, overwrite ps with ps$data_long;
  # otherwise treat ps as data_long
  if (inherits(ps, "phip_data")) {
    chk::chk_not_null(ps$data_long)
    ps <- ps$data_long
  }

  # ----------------------------------------------------------------------------
  # 1) column selection and basic structure checks
  # ----------------------------------------------------------------------------
  dat <- ps
  dat_cols <- dplyr::tbl_vars(dat)

  # decide which abundance column to use
  if (is.null(value_col)) {
    candidates <- c("exist", "counts_hits", "counts_control", "fold_change")
    hit <- candidates[candidates %in% dat_cols]
    if (length(hit) == 0L) {
      .ph_abort(
        paste0(
          "could not infer an abundance column in `ps$data_long`. ",
          "tried: ", paste(candidates, collapse = ", "),
          ". please specify `value_col` explicitly."
        ),
        step = "compute_distance"
      )
    }
    value_col <- hit[1L]
    .ph_log_info(
      paste0("auto-detected `value_col = \"", value_col,
             "\"` from `ps$data_long`."),
      step = "compute_distance"
    )
  }

  if (!value_col %in% dat_cols) {
    .ph_abort(
      paste0("column `", value_col, "` not found in `ps`."),
      step = "compute_distance"
    )
  }

  required_cols <- c("sample_id", "peptide_id", value_col)
  missing_cols  <- setdiff(required_cols, dat_cols)
  if (length(missing_cols) > 0L) {
    .ph_abort(
      paste0("missing required column(s) in `ps`: ",
             paste(missing_cols, collapse = ", ")),
      step = "compute_distance"
    )
  }

  .ph_log_info(
    paste0("building abundance matrix from `ps` using `",
           value_col, "`."),
    step = "compute_distance"
  )

  value_sym <- rlang::sym(value_col)

  # ----------------------------------------------------------------------------
  # 2) collect only needed columns and pivot in r
  # ----------------------------------------------------------------------------
  .ph_log_info("collecting long table (sample_id, peptide_id, value).",
               step = "compute_distance")

  dat_small <- dat |>
    dplyr::select(sample_id, peptide_id, !!value_sym) |>
    dplyr::collect()

  # basic sanity checks on collected data
  if (anyNA(dat_small$sample_id) || anyNA(dat_small$peptide_id)) {
    .ph_abort(
      "`ps` contains missing values in `sample_id` and/or
      `peptide_id`.",
      step = "compute_distance"
    )
  }

  # duplicates will break pivot_wider (or create list-cols); require uniqueness
  if (anyDuplicated(dat_small[, c("sample_id", "peptide_id")]) > 0L) {
    .ph_abort(
      "found duplicated (sample_id, peptide_id) pairs in `ps`.",
      step = "compute_distance"
    )
  }

  # replace nas with 0 in abundance column
  dat_small[[value_col]][is.na(dat_small[[value_col]])] <- 0

  # ensure abundance is numeric
  if (!is.numeric(dat_small[[value_col]])) {
    .ph_abort(
      paste0("`", value_col, "` must be numeric after collect()."),
      step = "compute_distance"
    )
  }

  .ph_log_info("pivoting to wide abundance matrix in r.",
               step = "compute_distance")

  wide_df <- dat_small |>
    tidyr::pivot_wider(
      id_cols     = sample_id,
      names_from  = peptide_id,
      values_from = !!value_sym,
      values_fill = 0
    )

  if (!"sample_id" %in% names(wide_df)) {
    .ph_abort("failed to construct wide abundance table (no `sample_id`).",
              step = "compute_distance")
  }

  mat <- wide_df |>
    tibble::column_to_rownames("sample_id") |>
    as.matrix()

  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    .ph_abort(
      "abundance matrix is empty after reshaping. check filters and
      `value_col`.",
      step = "compute_distance"
    )
  }

  .ph_log_info(
    paste0("abundance matrix has ", nrow(mat), " samples and ",
           ncol(mat), " features after preprocessing."),
    step = "compute_distance"
  )

  # ---------------------------------------------------------------------------
  # 3) normalization
  # ---------------------------------------------------------------------------
  if (identical(method_normalization, "auto")) {
    vals <- mat[!is.na(mat)]
    is_binary_data <- length(vals) > 0L && all(vals == 0 | vals == 1)
    method_normalization <- if (is_binary_data) "none" else "relative"
    .ph_log_info(
      paste0("auto normalization selected -> using ", method_normalization),
      step = "compute_distance"
    )
  }

  norm_mat <- switch(
    method_normalization,
    "none" = mat,
    "relative" = {
      rs <- rowSums(mat, na.rm = TRUE)
      rs[rs == 0] <- 1
      mat / rs
    },
    "hellinger" = {
      rs <- rowSums(mat, na.rm = TRUE)
      rs[rs == 0] <- 1
      sqrt(mat / rs)
    },
    "log" = log1p(mat)
  )

  # ---------------------------------------------------------------------------
  # 4) distance computation
  # ---------------------------------------------------------------------------
  if (dist_method == "chebyshev") dist_method <- "maximum"

  .ph_log_info(
    paste0("computing distance: ", dist_method),
    step = "compute_distance"
  )

  dist_obj <- NULL

  pd_methods <- c(
    "bray", "euclidean", "minkowski", "manhattan", "canberra",
    "binary", "maximum", "cosine"
  )

  vegan_methods <- c(
    "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
    "jaccard", "gower", "altgower", "morisita", "horn", "mountford", "raup",
    "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger",
    "aitchison", "robust.aitchison"
  )

  use_pd <- dist_method %in% pd_methods

  if (use_pd) {
    if (!rlang::is_installed("parallelDist")) {
      .ph_abort(
        paste0(
          "`distance = \"", distance, "\"` requires the suggested package
          'parallelDist'. ",
          "Please install 'parallelDist' or choose a vegan::vegdist() method."
        ),
        step = "compute_distance"
      )
    }

    dist_obj <- parallelDist::parDist(
      norm_mat,
      method  = dist_method,
      threads = n_threads
    )

  } else {
    dist_obj <- vegan::vegdist(norm_mat, method = dist_method)
  }

  .ph_log_info("distance matrix computation complete.",
               step = "compute_distance")

  # attach normalized abundance matrix as attribute
  attr(dist_obj, "abundances") <- norm_mat

  dist_obj
}

#' @title Principal Components Analysis (PCoA) on a Distance Matrix
#' @description Performs PCoA on a distance matrix (typically from
#' \code{compute_distance()}), optionally correcting for negative eigenvalues,
#' and returns coordinates, eigenvalues, variance explained, and feature
#' loadings.
#'
#' @param dist_obj a \code{dist} object (for example returned by
#'   \code{compute_distance()}). The normalized abundance matrix used to compute
#'   the distances is attached as attribute \code{"abundances"} (numeric
#'   matrix with samples in rows and features in columns). If missing, feature
#'   loadings are skipped.
#' @param neg_correction character scalar. Method for adjusting negative
#'   eigenvalues (if any). One of \code{"none"}, \code{"lingoes"}, or
#'   \code{"cailliez"}. Default is \code{"none"}.
#' @param n_axes integer scalar. Number of pcoa axes to return in the sample
#'   scores. Must be > 0. Internally, \code{k = min(n_axes, n_samples - 1)}.
#' @param top_features integer scalar. Number of features to keep per axis when
#'   reporting loadings. Features are selected by taking the union of the top
#'   \code{top_features} features (by absolute loading) for each returned axis.
#'   Must be > 0.
#'
#' @return a list of class \code{"beta_pcoa"} with elements:
#' \itemize{
#'   \item \code{sample_coords}: tibble with \code{sample_id} and columns
#'     \code{PCoA1}, \code{PCoA2}, ... up to \code{n_axes} (or fewer if
#'     \code{n_samples - 1} is smaller).
#'   \item \code{eigenvalues}: numeric vector of eigenvalues from the pcoa.
#'   \item \code{var_explained}: one-row tibble with percent variance explained
#'     by the returned axes and \code{%Other}. percentages are computed from the
#'     sum of positive eigenvalues.
#'   \item \code{feature_loadings}: tibble of feature loadings for the returned
#'     axes (empty if \code{"abundances"} is missing or cannot be aligned).
#' }
#'
#' @details
#' Negative eigenvalues indicate that the distances are not perfectly euclidean.
#' If \code{neg_correction} is \code{"lingoes"} or \code{"cailliez"}, a
#' correction is applied via \code{vegan::wcmdscale(add = ...)}.
#'
#' Feature loadings are computed as abundance-weighted averages of sample
#' scores:
#' \code{t(X) %*% U / colSums(X)}, where \code{X} is the abundance matrix and
#' \code{U} are the sample coordinates.
#'
#' @examples
#' \donttest{
#' # compute a distance matrix with an attached abundance matrix
#' # build an example <phip_data> object from the package example dataset
#' phip_path <- phip_example_path()
#'
#' ps <- phip_convert(
#'   data_long_path    = phip_path,
#'   backend           = "duckdb",
#'   peptide_library   = TRUE,
#'   subject_id        = "subject_id",
#'   peptide_id        = "peptide_id",
#'   sample_id         = "sample_id",
#'   exist             = "exist",
#'   timepoint         = "timepoint_factor",
#'   fold_change       = "fold_change",
#'   materialise_table = TRUE,
#'   auto_expand       = TRUE,
#'   n_cores           = 2
#' )
#'
#' # small subset for speed: 5 peptides at time t1
#' keep_pep <- c("16627", "5243", "24799", "16196", "18003")
#' dat_cols <- dplyr::tbl_vars(ps$data_long)
#' tp_col <- "time"
#'
#' ps_small <- ps %>%
#'   dplyr::filter(
#'     peptide_id %in% keep_pep,
#'     !!rlang::sym(tp_col) == "T1"
#'   ) %>%
#'   dplyr::collect()
#'
#' # compute distances (needs either 'parallelDist' or 'vegan')
#'   val_col <- if ("counts_hits" %in% dplyr::tbl_vars(ps_small)) {
#'     "counts_hits"
#'   } else {
#'     "exist"
#'   }
#'
#'   d <- compute_distance(
#'     ps_small,
#'     value_col = val_col,
#'     method_normalization = "hellinger",
#'     distance = "bray",
#'     n_threads = 2L
#'   )
#'
#'   pcoa_res <- compute_pcoa(dist_bc, neg_correction = "none", n_axes = 3L)
#'   pcoa_res$sample_coords
#'   pcoa_res$var_explained
#'   pcoa_res$feature_loadings
#' }
#' @export
compute_pcoa <- function(dist_obj,
                         neg_correction = c("none", "lingoes", "cailliez"),
                         n_axes = 5L,
                         top_features = 30L) {
  # ----------------------------------------------------------------------------
  # 1) input validation
  # ----------------------------------------------------------------------------
  chk::chk_s3_class(dist_obj, "dist")
  neg_correction <- match.arg(neg_correction)
  chk::chk_count(n_axes)
  chk::chk_gt(n_axes, 0)
  chk::chk_count(top_features)
  chk::chk_gt(top_features, 0)

  n <- attr(dist_obj, "Size")
  if (is.null(n) || n < 2L) {
    .ph_abort("`dist_obj` must contain at least 2 samples.",
              step = "compute_pcoa")
  }

  labels <- attr(dist_obj, "Labels")
  if (is.null(labels) || length(labels) != n) {
    labels <- as.character(seq_len(n))
  }

  # ----------------------------------------------------------------------------
  # 2) pcoa computation
  # ----------------------------------------------------------------------------
  .ph_log_info(
    "performing principal coordinates analysis",
    step = "compute_pcoa",
    bullets = if (identical(neg_correction, "none")) NULL else
      paste("using", neg_correction, "correction")
  )

  # cmdscale requires k in [1, n - 1]
  k_cmd <- min(n_axes, n - 1L)

  pcoa_fit <- if (identical(neg_correction, "none")) {
    stats::cmdscale(dist_obj, eig = TRUE, k = k_cmd)
  } else {
    vegan::wcmdscale(dist_obj, eig = TRUE, k = k_cmd, add = neg_correction)
  }

  eig_vals <- pcoa_fit$eig
  if (is.null(eig_vals)) eig_vals <- numeric(0L)
  eig_vals <- as.numeric(eig_vals)

  coords <- as.matrix(pcoa_fit$points)
  if (is.null(coords)) {
    coords <- matrix(0, nrow = n, ncol = 0L)
  }

  # enforce stable axis and sample naming
  if (ncol(coords) > 0L) {
    colnames(coords) <- paste0("PCoA", seq_len(ncol(coords)))
  }
  if (nrow(coords) == n) {
    rownames(coords) <- labels
  }

  # ----------------------------------------------------------------------------
  # 3) sample coordinates (first n_axes, or fewer)
  # ----------------------------------------------------------------------------
  k_use <- min(n_axes, ncol(coords))
  coords_k <- if (k_use > 0L) {
    coords[, seq_len(k_use), drop = FALSE]
  } else {
    matrix(
      numeric(0L),
      nrow = n,
      ncol = 0L,
      dimnames = list(labels, character(0L))
    )
  }

  sample_coords <- tibble::as_tibble(coords_k, rownames = "sample_id")

  # ----------------------------------------------------------------------------
  # 4) variance explained (based on positive eigenvalues)
  # ----------------------------------------------------------------------------
  pos_eig <- pmax(eig_vals, 0)
  sum_pos <- sum(pos_eig, na.rm = TRUE)

  if (sum_pos > 0 && k_use > 0L) {
    pct_axes <- 100 * pos_eig[seq_len(k_use)] / sum_pos
    pct_other <- if (length(pos_eig) > k_use) {
      100 * sum(pos_eig[(k_use + 1L):length(pos_eig)], na.rm = TRUE) / sum_pos
    } else {
      0
    }
  } else {
    pct_axes <- rep(NA_real_, k_use)
    pct_other <- NA_real_
  }

  names(pct_axes) <- paste0("%PCoA", seq_len(k_use))
  var_explained <- tibble::as_tibble_row(
    c(as.list(round(pct_axes, 3)), `%Other` = round(pct_other, 3))
  )

  # ----------------------------------------------------------------------------
  # 5) feature loadings (requires abundances attribute)
  # ----------------------------------------------------------------------------
  feature_loadings <- tibble::tibble()

  X <- attr(dist_obj, "abundances")
  if (is.null(X)) {
    .ph_warn(
      "no 'abundances' attribute found on `dist_obj`;
      skipping feature loadings.",
      step = "compute_pcoa"
    )
  } else if (k_use < 1L) {
    # no axes, nothing to load
    feature_loadings <- tibble::tibble()
  } else {
    X <- as.matrix(X)

    coords_ids <- rownames(coords)
    X_ids <- rownames(X)

    if (is.null(coords_ids) || is.null(X_ids)) {
      .ph_warn(
        "row names missing in coordinates or 'abundances'; cannot align samples
        for feature loadings.",
        step = "compute_pcoa"
      )
    } else {
      common_ids <- intersect(coords_ids, X_ids)

      if (length(common_ids) < 2L) {
        .ph_warn(
          "insufficient overlap between distance labels and abundance rows;
          skipping feature loadings.",
          step = "compute_pcoa"
        )
      } else {
        # compute loadings for all returned axes
        ax_idx <- seq_len(k_use)
        U <- coords[common_ids, ax_idx, drop = FALSE]
        Xsub <- X[common_ids, , drop = FALSE]

        # keep features with nonzero total weight
        w <- colSums(Xsub, na.rm = TRUE)
        keep_feats <- which(w > 0)

        if (length(keep_feats) > 0L) {
          # feature x axis loadings: weighted mean of sample scores
          S <- t(Xsub[, keep_feats, drop = FALSE]) %*% U
          S <- sweep(S, 1, w[keep_feats], "/")

          rownames(S) <- colnames(Xsub)[keep_feats]
          colnames(S) <- colnames(U)

          load_tbl <- tibble::as_tibble(S, rownames = "feature")

          # select up to top_features per axis, then take the union
          ax_names <- colnames(U)
          top_list <- unique(unlist(lapply(seq_along(ax_names), function(j) {
            ord <- order(abs(S[, j]), decreasing = TRUE)
            head(rownames(S)[ord], top_features)
          })))

          feature_loadings <- dplyr::filter(load_tbl,
                                            .data$feature %in% top_list)
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 6) return
  # ---------------------------------------------------------------------------
  result <- list(
    sample_coords    = sample_coords,
    eigenvalues      = eig_vals,
    var_explained    = var_explained,
    feature_loadings = feature_loadings
  )
  class(result) <- "beta_pcoa"

  .ph_log_info("pcoa analysis complete.", step = "compute_pcoa")
  result
}

#' @title Constrained Ordination (db-rda / cap) on Distance Matrix
#'
#' @description Performs distance-based redundancy analysis (constrained pcoa,
#'   a.k.a. cap) on a distance matrix using \pkg{vegan}::\code{capscale}, with
#'   optional negative eigenvalue correction. Returns constrained sample scores,
#'   eigenvalues, variance partitioning, and feature loadings.
#'
#' @param dist_obj A \code{dist} object returned by \code{compute_distance()}.
#'   The normalized abundance matrix used to compute the distances is expected
#'   to be attached as an attribute \code{"abundances"} (samples in rows,
#'   features in columns).
#' @param ps A \code{phip_data} object or a table providing sample-level
#'   metadata. This table must contain \code{sample_id} and all variables
#'   referenced on the right-hand side of \code{formula}. Note that variable
#'   detection uses \code{all.vars(terms(formula))}, so transformed terms like
#'   \code{log(age)} are supported as long as \code{age} exists.
#' @param formula An R formula specifying the constraints (independent
#'   variables) for the ordination, e.g. \code{~ sex + age}. Do not include a
#'   response on the left-hand side; the distance matrix is provided via
#'   \code{dist_obj}.
#' @param neg_correction One of \code{"none"}, \code{"lingoes"},
#'   \code{"cailliez"}. Method for negative eigenvalue correction. Default is
#'   \code{"none"}. Passed to the \code{add} argument of
#'   \code{vegan::capscale()}.
#' @param top_features Integer scalar. Number of top features to return in
#'   loadings (selected per constrained axis by absolute loading, then unioned).
#'   Default is 30.
#'
#' @return A list of class \code{"beta_capscale"} with elements:
#' \item{sample_coords}{Tibble of sample scores on constrained axes
#'   (\code{CAP1}, \code{CAP2}, ...). Contains \code{sample_id} and coordinates.}
#' \item{eigenvalues}{Numeric vector of eigenvalues of the constrained axes.}
#' \item{variance_partition}{Tibble with total inertia and inertia partitioned
#'   into constrained and unconstrained components, with their proportion of total.}
#' \item{feature_loadings}{Tibble of top feature loadings for constrained axes
#'   (possibly empty if the \code{"abundances"} attribute is missing or cannot be aligned).
#'   To limit runtime/memory, loadings are computed for at most 10 constrained axes.}
#' \item{cap_model}{The full \code{vegan::capscale} model object.}
#'
#' @examples
#' \donttest{
#' phip_path <- phip_example_path()
#'
#' ps <- phip_convert(
#'   data_long_path    = phip_path,
#'   backend           = "duckdb",
#'   peptide_library   = TRUE,
#'   subject_id        = "subject_id",
#'   peptide_id        = "peptide_id",
#'   sample_id         = "sample_id",
#'   exist             = "exist",
#'   timepoint         = "timepoint_factor",
#'   fold_change       = "fold_change",
#'   materialise_table = TRUE,
#'   auto_expand       = TRUE,
#'   n_cores           = 2
#' )
#'
#' # small subset for speed: 5 peptides at time t1
#' keep_pep <- c("16627", "5243", "24799", "16196", "18003")
#' dat_cols <- dplyr::tbl_vars(ps$data_long)
#' tp_col <- if ("time" %in% dat_cols) "time" else "timepoint_factor"
#'
#' ps_small <- ps |>
#'   dplyr::filter(
#'     peptide_id %in% keep_pep,
#'     !!rlang::sym(tp_col) == "T1"
#'   ) |>
#'   dplyr::collect()
#'
#' val_col <- if ("fold_change" %in% dplyr::tbl_vars(ps_small)) {
#'   "fold_change"
#' } else {
#'   "exist"
#' }
#'
#' dist_bc <- compute_distance(
#'   ps_small,
#'   value_col = val_col,
#'   method_normalization = "hellinger",
#'   distance = "bray",
#'   n_threads = 2L
#' )
#'
#' # pick a simple constraint that exists in the example data (fallback order)
#' dat <- ps_small
#' cand <- c("group", "big_group", "type_person", "sex", "age")
#' rhs_var <- cand[cand %in% dplyr::tbl_vars(dat)][1]
#'
#' cap_res <- compute_capscale(
#'   dist_bc,
#'   ps      = ps_small,
#'   formula = stats::as.formula(paste0("~ ", rhs_var)),
#'   neg_correction = "none",
#'   top_features   = 30L
#' )
#'
#' cap_res$variance_partition
#' cap_res$sample_coords
#' cap_res$feature_loadings
#' }
#' @export
compute_capscale <- function(dist_obj,
                             ps,
                             formula,
                             neg_correction = c("none", "lingoes", "cailliez"),
                             top_features = 30L) {
  # ----------------------------------------------------------------------------
  # 0) input validation
  # ----------------------------------------------------------------------------
  chk::chk_s3_class(dist_obj, "dist")
  chk::chk_true(inherits(formula, "formula"))
  neg_correction <- match.arg(neg_correction)
  chk::chk_count(top_features)
  chk::chk_gt(top_features, 0)

  # ---------------------------------------------------------------------------
  # 1) distance + labels
  # ---------------------------------------------------------------------------
  d <- dist_obj
  labels <- attr(d, "Labels")
  n <- attr(d, "Size")

  # abundance matrix from dist attributes (may be null)
  X_all <- attr(dist_obj, "abundances")

  # ---------------------------------------------------------------------------
  # 2) metadata from ps + variable checks + alignment + na handling
  # ---------------------------------------------------------------------------
  dat <- if("phip_data" %in% class(ps)) ps$data_long else ps
  if (is.null(dat)) {
    .ph_abort("`ps` is missing. cannot construct metadata.",
              step = "compute_capscale")
  }

  dat_cols <- dplyr::tbl_vars(dat)
  if (!"sample_id" %in% dat_cols) {
    .ph_abort("`ps$data_long` must contain a `sample_id` column.",
              step = "compute_capscale")
  }

  # required variables from formula rhs (supports transformations like log(age))
  vars_needed <- all.vars(stats::terms(formula))
  vars_needed <- setdiff(vars_needed, ".") # defensive

  if (length(vars_needed) == 0L) {
    .ph_abort(
      "No constraints provided in formula (rhs is empty).
      Use compute_pcoa() for unconstrained ordination.",
      step = "compute_capscale"
    )
  }

  missing_vars <- setdiff(vars_needed, dat_cols)
  if (length(missing_vars) > 0L) {
    .ph_abort(
      paste0(
        "the following variables from the formula are missing in `ps`: ",
        paste(missing_vars, collapse = ", ")
      ),
      step = "compute_capscale"
    )
  }

  .ph_log_info(
    "building metadata from `ps$data_long`.",
    step = "compute_capscale"
  )

  # sample-level metadata (one row per sample_id)
  meta_all <- dat |>
    dplyr::select(sample_id, dplyr::all_of(vars_needed)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect() |>
    as.data.frame()

  if (nrow(meta_all) == 0L) {
    .ph_abort(
      "constructed metadata has zero rows. check that `ps` is not empty.",
      step = "compute_capscale"
    )
  }

  rownames(meta_all) <- meta_all$sample_id

  # align metadata to distance order
  if (!is.null(labels)) {
    idx <- match(labels, meta_all$sample_id)
    missing_samples <- labels[is.na(idx)]
    if (length(missing_samples) > 0L) {
      .ph_abort(
        paste0(
          "the following samples from `dist_obj` are missing in `ps`: ",
          paste(missing_samples, collapse = ", ")
        ),
        step = "compute_capscale"
      )
    }
    meta_sub <- meta_all[idx, , drop = FALSE]
    rownames(meta_sub) <- labels
  } else {
    meta_sub <- meta_all
    labels <- rownames(meta_sub)
    n <- length(labels)

    .ph_warn(
      "no labels found in `dist_obj`; assuming metadata row order matches
      the distance order.",
      step = "compute_capscale"
    )

    # critical: propagate assumed labels to the dist object so downstream
    # indexing works
    if (!is.null(attr(d, "Size")) && attr(d, "Size") != n) {
      .ph_abort(
        paste0(
          "cannot align distance and metadata without labels: dist size is ",
          attr(d, "Size"),
          " but metadata has ", n, " samples."
        ),
        step = "compute_capscale"
      )
    }
    attr(d, "Labels") <- labels
  }

  # drop samples with any missing values in required variables
  rhs_df <- meta_sub[, vars_needed, drop = FALSE]
  keep <- stats::complete.cases(rhs_df)

  if (!all(keep)) {
    dropped <- sum(!keep)
    .ph_log_info(
      paste0("dropping ", dropped,
             " samples with missing values in constrained variables."),
      step = "compute_capscale"
    )
  }

  meta_df <- meta_sub[keep, , drop = FALSE]
  if (nrow(meta_df) == 0L) {
    .ph_abort(
      "all samples have missing values in constrained variables; cannot fit
      cap.",
      step = "compute_capscale"
    )
  }

  keep_labels <- rownames(meta_df)

  # ----------------------------------------------------------------------------
  # 3) subset distance and abundances to complete-case samples
  # ----------------------------------------------------------------------------
  mat_d <- as.matrix(d)
  mat_d_sub <- mat_d[keep_labels, keep_labels, drop = FALSE]
  d <- stats::as.dist(mat_d_sub)
  n <- attr(d, "Size")
  labels <- attr(d, "Labels")

  X_sub <- NULL
  if (!is.null(X_all)) {
    X_all <- as.matrix(X_all)
    if (!is.null(rownames(X_all))) {
      X_sub <- X_all[keep_labels, , drop = FALSE]
    } else {
      .ph_warn(
        "abundance matrix has no row names; cannot align precisely with
        samples.",
        step = "compute_capscale"
      )
    }
  }

  # ----------------------------------------------------------------------------
  # 4) build capscale formula: d_resp ~ rhs
  # ----------------------------------------------------------------------------
  d_resp <- d
  cap_formula <- stats::update.formula(formula, d_resp ~ .)
  environment(cap_formula) <- environment()

  add_arg <- if (identical(neg_correction, "none")) FALSE else neg_correction

  .ph_log_info(
    "fitting constrained ordination (cap/db-rda)",
    step = "compute_capscale",
    bullets = c(
      paste0("formula: ", paste(deparse(formula), collapse = " ")),
      if (!identical(add_arg, FALSE)) paste0("neg_correction: ", add_arg)
    )
  )

  # ----------------------------------------------------------------------------
  # 5) fit capscale model
  # ----------------------------------------------------------------------------
  cap_fit <- vegan::capscale(cap_formula, data = meta_df, add = add_arg)

  # ----------------------------------------------------------------------------
  # 6) sample scores on constrained axes
  # ----------------------------------------------------------------------------
  rank_constrained <- cap_fit$CCA$rank
  if (is.null(rank_constrained)) rank_constrained <- 0L

  if (rank_constrained > 0L) {
    site_scores <- vegan::scores(
      cap_fit,
      display = "sites",
      choices = seq_len(rank_constrained)
    )
    pts <- as.matrix(site_scores)
  } else {
    pts <- matrix(0, nrow = n, ncol = 0L,
                  dimnames = list(labels, character(0L)))
  }

  if (ncol(pts) > 0L && is.null(colnames(pts))) {
    colnames(pts) <- paste0("CAP", seq_len(ncol(pts)))
  }
  if (!is.null(labels) && nrow(pts) == length(labels)) {
    rownames(pts) <- labels
  }

  sample_coords <- tibble::as_tibble(pts, rownames = "sample_id")

  # ----------------------------------------------------------------------------
  # 7) eigenvalues and variance partition
  # ----------------------------------------------------------------------------
  eig_constrained <- cap_fit$CCA$eig
  if (is.null(eig_constrained)) eig_constrained <- numeric()

  tot_inertia <- cap_fit$tot.chi
  cons_inertia <- cap_fit$CCA$tot.chi
  if (is.null(cons_inertia)) cons_inertia <- sum(cap_fit$CCA$eig %||% 0)

  uncon_inertia <- cap_fit$CA$tot.chi
  if (is.null(uncon_inertia)) uncon_inertia <- sum(cap_fit$CA$eig %||% 0)

  variance_partition <- tibble::tibble(
    component  = c("Total", "Constrained", "Unconstrained"),
    inertia    = c(tot_inertia, cons_inertia, uncon_inertia),
    proportion = c(
      1,
      if (is.finite(tot_inertia) && tot_inertia > 0)
        cons_inertia / tot_inertia else NA_real_,
      if (is.finite(tot_inertia) && tot_inertia > 0)
        uncon_inertia / tot_inertia else NA_real_
    )
  )

  # ---------------------------------------------------------------------------
  # 8) feature loadings on constrained axes (using X_sub)
  # ---------------------------------------------------------------------------
  feature_loadings <- tibble::tibble()

  if (!is.null(X_sub) && rank_constrained > 0L) {
    coords_ids <- rownames(pts)
    X_ids <- rownames(X_sub)

    if (is.null(coords_ids) || is.null(X_ids)) {
      .ph_warn(
        "row names missing in sample scores or abundance matrix; cannot align
        for feature loadings.",
        step = "compute_capscale"
      )
    } else {
      common_ids <- intersect(coords_ids, X_ids)

      if (length(common_ids) < 2L) {
        .ph_warn(
          "insufficient overlap between distance labels and abundance rows;
          skipping feature loadings.",
          step = "compute_capscale"
        )
      } else {
        # to limit runtime/memory, compute loadings for at most 10 constrained
        # axes
        ax_idx <- seq_len(min(rank_constrained, 10L))
        U <- pts[common_ids, ax_idx, drop = FALSE]
        Xsub <- X_sub[common_ids, , drop = FALSE]

        w <- colSums(Xsub, na.rm = TRUE)
        keep_feats <- which(w > 0)

        if (length(keep_feats) > 0L) {
          S <- t(Xsub[, keep_feats, drop = FALSE]) %*% U
          S <- sweep(S, 1, w[keep_feats], "/")

          rownames(S) <- colnames(Xsub)[keep_feats]
          colnames(S) <- colnames(U)

          load_tbl <- tibble::as_tibble(S, rownames = "feature")

          ax_names <- colnames(U)
          top_list <- unique(unlist(lapply(seq_along(ax_names), function(j) {
            ord <- order(abs(S[, j]), decreasing = TRUE)
            head(rownames(S)[ord], top_features)
          })))

          feature_loadings <- dplyr::filter(load_tbl,
                                            .data$feature %in% top_list)
        }
      }
    }
  }

  # ----------------------------------------------------------------------------
  # 9) return
  # ----------------------------------------------------------------------------
  result <- list(
    sample_coords      = sample_coords,
    eigenvalues        = as.numeric(eig_constrained),
    variance_partition = variance_partition,
    feature_loadings   = feature_loadings,
    cap_model          = cap_fit
  )
  class(result) <- "beta_capscale"

  .ph_log_info("cap analysis complete.", step = "compute_capscale")

  result
}

#' @title PERMANOVA with Global and Post-hoc Tests on Beta Diversity
#' @description Performs PERMANOVA (adonis2) on a distance matrix for overall group/time effects,
#' and optionally conducts post-hoc pairwise or contrast tests (e.g., between each pair of groups, etc.).
#' Supports stratified permutations for repeated measures.
#'
#' @param dist_obj A \code{dist} object of distances between samples
#'   (e.g., output of \code{compute_distance()}).
#' @param ps A \code{phip_data} object or a table providing sample-level
#'   metadata. This table must contain \code{sample_id} and the
#'   columns specified in \code{group_col}, \code{time_col}, and optionally
#'   \code{subject_col}.
#' @param group_col Name of the grouping column in \code{ps}
#'   (between-subject factor). Use \code{NULL} if no group factor.
#' @param time_col Name of the time factor column in \code{ps}
#'   (within-subject factor for longitudinal data). Use \code{NULL} if not
#'   applicable. This should be a \emph{categorical} factor for this function
#'   (continuous time not supported).
#' @param subject_col Name of the subject identifier column in \code{ps}
#'   (for repeated measures). Default \code{"subject_id"}. If this column is
#'   present and \code{time_col} is provided, permutations will be stratified by
#'   subject.
#' @param permutations Number of permutations for significance testing
#'   (default 999).
#' @param contrasts Character vector specifying which post-hoc contrasts to perform.
#'   Options: \code{"none"} (default) for no post-hoc tests,
#'   \code{"pairwise"} for all pairwise comparisons,
#'   \code{"each_vs_rest"} for each level vs. the rest,
#'   and \code{"baseline"} for a specific baseline comparison.
#' @param baseline_level If \code{contrasts} includes \code{"baseline"},
#'   specify the factor level to use as baseline. For a time factor, if
#'   \code{baseline_level} is \code{NULL}, the first level of \code{time_col}
#'   will be used as baseline.
#'
#' @return A tibble with columns:
#' \item{scope}{The scope of the test (e.g., \code{"global"}, \code{"group_pairwise"},
#'   \code{"time_pairwise"}, \code{"each_vs_rest"}, \code{"baseline"}).}
#' \item{contrast}{Description of the contrast (e.g., \code{"<global>"} for overall test,
#'   \code{"A vs B"} for pairwise group comparisons, \code{"X vs other"} for each-vs-rest, etc.).}
#' \item{term}{The term being tested. For global tests, this will be the factor name (or interaction).
#'   For post-hoc tests, it may be \code{"group"} or \code{"time"} indicating which factor is being contrasted.}
#' \item{F_stat}{F statistic of the PERMANOVA test (for global tests and some contrasts where applicable).}
#' \item{R2}{R-squared (variance explained) for the term (global tests).}
#' \item{p_value}{Permutation p-value for the test.}
#' \item{n_perm}{Number of permutations used.}
#'
#' @details
#' Global PERMANOVA uses a model with main effects of \code{group_col} and
#' \code{time_col} and their interaction (if both are provided and have >1 level).
#' Samples with missing values in the constrained variables (group/time and,
#' when used for stratification, subject) are dropped \emph{before} fitting the
#' model, and the distance matrix is subset to the remaining samples so that
#' distances and metadata are always aligned.
#'
#' Post-hoc contrasts (\code{"pairwise"}, \code{"each_vs_rest"}, \code{"baseline"})
#' follow the same logic as described in the detailed comments in the source,
#' using \code{adonis2} with appropriate subsetting and, where applicable,
#' subject stratification.
#'
#' @examples
#' \dontrun{
#'   permanova_res <- compute_permanova(
#'     dist_bc,
#'     ps        = ps,
#'     group_col = "type_person",
#'     time_col  = "timepoint"
#'   )
#'
#'   permanova_res2 <- compute_permanova(
#'     dist_bc,
#'     ps        = ps,
#'     group_col = "type_person",
#'     time_col  = "timepoint",
#'     contrasts = c("pairwise", "baseline"),
#'     baseline_level = "T0"
#'   )
#' }
#' @export
compute_permanova <- function(dist_obj,
                              ps,
                              group_col = NULL,
                              time_col = NULL,
                              subject_col = "subject_id",
                              permutations = 999,
                              contrasts = "none",
                              baseline_level = NULL) {
  # ----------------------------------------------------------------------------
  # input validation (chk)
  # ----------------------------------------------------------------------------
  chk::chk_s3_class(dist_obj, "dist")
  if (!is.null(group_col)) chk::chk_string(group_col)
  if (!is.null(time_col)) chk::chk_string(time_col)
  chk::chk_string(subject_col)
  chk::chk_count(permutations)
  chk::chk_gt(permutations, 0)
  chk::chk_character(contrasts)
  if (!is.null(baseline_level)) chk::chk_string(baseline_level)

  # if ps is <phip_data>, overwrite ps with ps$data_long;
  # otherwise treat ps as data_long
  if (inherits(ps, "phip_data")) {
    chk::chk_not_null(ps$data_long)
    ps <- ps$data_long
  }

  # ---------------------------------------------------------------------------
  # prepare result collector
  # ---------------------------------------------------------------------------
  results_list <- list()
  add_result <- function(scope, contrast, term = NA, p_value = NA,
                         F_stat = NA, R2 = NA) {
    results_list[[length(results_list) + 1]] <<- tibble::tibble(
      scope   = scope,
      contrast = contrast,
      term    = term,
      F_stat  = F_stat,
      R2      = R2,
      p_value = p_value,
      n_perm  = permutations
    )
  }

  # ---------------------------------------------------------------------------
  # start from distances + labels
  # ---------------------------------------------------------------------------
  d_full  <- dist_obj
  labels_full <- attr(d_full, "Labels")
  n_full  <- attr(d_full, "Size")

  # ---------------------------------------------------------------------------
  # build metadata from ps and align to dist labels
  # ---------------------------------------------------------------------------
  dat <- ps
  if (is.null(dat)) {
    .ph_abort("`ps` is missing. Cannot construct metadata.",
              step = "compute_permanova")
  }

  dat_cols <- dplyr::tbl_vars(dat)
  if (!"sample_id" %in% dat_cols) {
    .ph_abort("`ps` must contain a `sample_id` column.",
              step = "compute_permanova")
  }

  has_group   <- !is.null(group_col)   && group_col   %in% dat_cols
  has_time    <- !is.null(time_col)    && time_col    %in% dat_cols
  has_subject <- !is.null(subject_col) && subject_col %in% dat_cols

  if (!is.null(group_col) && !has_group) {
    .ph_abort(
      paste0("Column `", group_col, "` not found in `ps`."),
      step = "compute_permanova"
    )
  }
  if (!is.null(time_col) && !has_time) {
    .ph_abort(
      paste0("Column `", time_col, "` not found in `ps`."),
      step = "compute_permanova"
    )
  }

  cols_needed <- c(
    "sample_id",
    if (has_group) group_col else character(0L),
    if (has_time)  time_col  else character(0L),
    if (has_subject) subject_col else character(0L)
  )
  cols_needed <- unique(cols_needed)

  # build metadata from ps
  .ph_log_info(
    "building metadata from `ps`.",
    step = "compute_permanova"
  )

  meta_all <- dat |>
    dplyr::select(dplyr::all_of(cols_needed)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect() |>
    as.data.frame()

  if (nrow(meta_all) == 0L) {
    .ph_abort(
      "constructed metadata has zero rows. check that `ps` is not empty.",
      step = "compute_permanova"
    )
  }
  rownames(meta_all) <- meta_all$sample_id

  # align metadata to distance order
  if (!is.null(labels_full)) {
    idx_align <- match(labels_full, meta_all$sample_id)
    missing_samples <- labels_full[is.na(idx_align)]
    if (length(missing_samples) > 0L) {
      .ph_abort(
        paste0(
          "the following samples from `dist_obj` are missing in `ps`: ",
          paste(missing_samples, collapse = ", ")
        ),
        step = "compute_permanova"
      )
    }
    meta_sub <- meta_all[idx_align, , drop = FALSE]
    rownames(meta_sub) <- labels_full
  } else {
    meta_sub <- meta_all
    labels_full <- rownames(meta_sub)
    n_full <- length(labels_full)
    .ph_warn(
      "no labels found in `dist_obj`; assuming metadata row order matches
      the distance order.",
      step = "compute_permanova"
    )
  }

  # coerce group/time to factor for safety
  if (has_group) {
    meta_sub[[group_col]] <- as.factor(meta_sub[[group_col]])
  }
  if (has_time) {
    meta_sub[[time_col]] <- as.factor(meta_sub[[time_col]])
  }

  # ---------------------------------------------------------------------------
  # drop samples with na in constrained variables (+ subject if used)
  # ---------------------------------------------------------------------------
  vars_for_na <- c(
    if (has_group) group_col else character(0L),
    if (has_time)  time_col  else character(0L),
    if (has_time && has_subject) subject_col else character(0L)  # subject used for strata in longitudinal
  )
  vars_for_na <- unique(vars_for_na)

  if (length(vars_for_na) > 0L) {
    keep <- stats::complete.cases(meta_sub[, vars_for_na, drop = FALSE])
  } else {
    keep <- rep(TRUE, nrow(meta_sub))
  }

  if (!all(keep)) {
    dropped <- sum(!keep)
    .ph_log_info(
      paste0("dropping ", dropped,
             " samples with missing values in constrained/strata variables."),
      step = "compute_permanova"
    )
  }

  meta_df <- meta_sub[keep, , drop = FALSE]
  if (nrow(meta_df) == 0L) {
    .ph_abort(
      "all samples have missing values in constrained/strata variables; cannot run permanova.",
      step = "compute_permanova"
    )
  }
  keep_labels <- rownames(meta_df)

  # subset distance matrix to complete-case samples
  mat_d_full <- as.matrix(d_full)
  mat_d_sub  <- mat_d_full[keep_labels, keep_labels, drop = FALSE]
  d <- stats::as.dist(mat_d_sub)
  labels <- attr(d, "Labels")
  n      <- attr(d, "Size")

  # update factor presence after na-drop
  has_group <- has_group && length(unique(meta_df[[group_col]])) > 1L
  has_time  <- has_time  && length(unique(meta_df[[time_col]]))  > 1L

  # ---------------------------------------------------------------------------
  # global permanova
  # ---------------------------------------------------------------------------
  rhs_terms <- character(0L)
  if (has_group) rhs_terms <- c(rhs_terms, group_col)
  if (has_time)  rhs_terms <- c(rhs_terms, time_col)
  if (has_group && has_time) {
    rhs_terms <- c(rhs_terms, paste(group_col, "*", time_col))
  }

  if (length(rhs_terms) == 0L) {
    .ph_log_info(
      "global permanova skipped (insufficient number of factor levels).",
      step = "compute_permanova"
    )
  } else {
    formula_str <- paste("d_resp ~", paste(rhs_terms, collapse = " + "))
    d_resp <- d
    form <- stats::as.formula(formula_str, env = environment())

    df <- meta_df
    rownames(df) <- labels

    # strata for global test: subject if time present and repeated measures
    strata_var <- NULL
    if (has_time && has_subject && subject_col %in% names(meta_df)) {
      subj <- meta_df[[subject_col]]
      subj_counts <- table(subj)
      if (any(subj_counts > 1L)) {
        strata_var <- subj
      }
    }

    .ph_log_info(
      "running global permanova",
      step = "compute_permanova",
      bullets = c(
        paste("model:", formula_str),
        if (!is.null(strata_var)) "permutations stratified by subject"
      )
    )

    adonis_res <- try(
      vegan::adonis2(
        form,
        data        = df,
        permutations = permutations,
        by          = "terms",
        strata      = strata_var,
        parallel    = 1
      ),
      silent = TRUE
    )

    if (inherits(adonis_res, "try-error")) {
      .ph_warn("global permanova failed; no global results.",
               step = "compute_permanova")
    } else {
      res_df <- as.data.frame(adonis_res)
      res_df$term <- rownames(res_df)
      # check terms: group, time, and interaction
      for (term in c(
        if (!is.null(group_col)) group_col else NULL,
        if (!is.null(time_col))  time_col  else NULL,
        if (!is.null(group_col) && !is.null(time_col))
          paste(group_col, ":", time_col, sep = "") else NULL
      )) {
        if (!is.null(term) && term %in% res_df$term) {
          row <- res_df[res_df$term == term, , drop = FALSE]
          add_result(
            "global",
            "<global>",
            term    = term,
            F_stat  = row$F[1],
            R2      = row$R2[1],
            p_value = row$`Pr(>F)`[1]
          )
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # post-hoc contrasts (pairwise / each_vs_rest / baseline)
  # ---------------------------------------------------------------------------
  contrasts <- tolower(unique(contrasts))
  contrasts[contrasts == "group_vs_rest"] <- "each_vs_rest"

  # determine baseline if needed
  bl <- NULL
  if ("baseline" %in% contrasts && !is.null(baseline_level)) {
    bl <- baseline_level
  }
  if ("baseline" %in% contrasts && is.null(baseline_level) && has_time) {
    bl <- sort(unique(meta_df[[time_col]]))[1]
    .ph_log_info(
      paste("baseline_level not provided; using", bl, "as baseline for time."),
      step = "compute_permanova"
    )
  }

  # helper: run two-level adonis on subset idx of current meta_df/d
  run_two_level_adonis <- function(idx, fac, covar = NULL, strata = NULL,
                                   scope_label, contrast_label, term_label) {
    if (length(unique(fac)) < 2L || min(table(fac)) < 2L) {
      .ph_log_info(
        paste("skipping test", contrast_label,
              "- not enough samples in one or both groups."),
        step = "compute_permanova"
      )
      return(NULL)
    }
    df_sub <- meta_df[idx, , drop = FALSE]
    labels_sub <- labels[idx]
    rownames(df_sub) <- labels_sub

    form_rhs <- "fac"
    if (!is.null(covar)) form_rhs <- paste(form_rhs, "+", covar)

    fml <- stats::as.formula(paste("d_sub ~", form_rhs))
    d_mat <- as.matrix(d)
    d_sub <- stats::as.dist(d_mat[idx, idx, drop = FALSE])

    ad <- vegan::adonis2(
      fml,
      data        = cbind(df_sub, fac = fac),
      permutations = permutations,
      strata      = strata,
      by          = if (!is.null(covar)) "margin" else NULL
    )

    res <- as.data.frame(ad)
    res$term <- rownames(res)
    if ("fac" %in% res$term) {
      row <- res[res$term == "fac", , drop = FALSE]
      add_result(
        scope_label,
        contrast_label,
        term    = term_label,
        F_stat  = row$F[1],
        R2      = row$R2[1] %||% NA_real_,
        p_value = row$`Pr(>F)`[1]
      )
    }
  }

  # pairwise
  if ("pairwise" %in% contrasts) {
    # pairwise group comparisons
    if (has_group) {
      groups <- na.omit(unique(meta_df[[group_col]]))
      if (length(groups) > 1L) {
        pairs <- utils::combn(groups, 2, simplify = FALSE)
        for (p in pairs) {
          sel <- which(meta_df[[group_col]] %in% p)

          strata_use <- NULL
          if (has_subject && subject_col %in% names(meta_df)) {
            sub_sel <- meta_df[sel, subject_col]
            # if any subject appears in both groups -> stratify
            multi <- any(tapply(meta_df[sel, group_col], sub_sel,
                                function(x) length(unique(x)) > 1L))
            if (multi) strata_use <- sub_sel
          }

          fac_pair   <- factor(meta_df[sel, group_col], levels = p)
          covar_term <- if (has_time) time_col else NULL

          run_two_level_adonis(
            idx            = sel,
            fac            = fac_pair,
            covar          = covar_term,
            strata         = strata_use,
            scope_label    = "group_pairwise",
            contrast_label = paste(p, collapse = " vs "),
            term_label     = group_col
          )
        }
      }
    }

    # pairwise time comparisons
    if (has_time) {
      times <- na.omit(unique(meta_df[[time_col]]))
      if (length(times) > 1L) {
        pairs <- utils::combn(times, 2, simplify = FALSE)
        for (p in pairs) {
          sel <- which(meta_df[[time_col]] %in% p)

          strata_use <- NULL
          if (has_subject && subject_col %in% names(meta_df)) {
            sub_sel <- meta_df[sel, subject_col]
            if (any(table(sub_sel) > 1L)) {
              strata_use <- sub_sel
            }
          }

          fac_pair   <- factor(meta_df[sel, time_col], levels = p)
          covar_term <- if (has_group) group_col else NULL

          run_two_level_adonis(
            idx            = sel,
            fac            = fac_pair,
            covar          = covar_term,
            strata         = strata_use,
            scope_label    = "time_pairwise",
            contrast_label = paste(p, collapse = " vs "),
            term_label     = time_col
          )
        }
      }
    }
  }

  # each_vs_rest
  if ("each_vs_rest" %in% contrasts) {
    # each group vs rest
    if (has_group && length(unique(meta_df[[group_col]])) > 1L) {
      for (lvl in unique(meta_df[[group_col]])) {
        fac_vec <- factor(
          ifelse(meta_df[[group_col]] == lvl, lvl, "other"),
          levels = c(lvl, "other")
        )
        if (length(unique(fac_vec)) < 2L || min(table(fac_vec)) < 2L) next
        sel <- which(!is.na(fac_vec))
        covar_term <- if (has_time) time_col else NULL

        run_two_level_adonis(
          idx            = sel,
          fac            = fac_vec[sel],
          covar          = covar_term,
          strata         = NULL,
          scope_label    = "each_vs_rest",
          contrast_label = paste(lvl, "vs other"),
          term_label     = group_col
        )
      }
    }

    # each time vs rest
    if (has_time && length(unique(meta_df[[time_col]])) > 1L) {
      for (lvl in unique(meta_df[[time_col]])) {
        fac_vec <- factor(
          ifelse(meta_df[[time_col]] == lvl, lvl, "other"),
          levels = c(lvl, "other")
        )
        if (length(unique(fac_vec)) < 2L || min(table(fac_vec)) < 2L) next
        sel <- which(!is.na(fac_vec))

        strata_use <- NULL
        if (has_subject && subject_col %in% names(meta_df)) {
          sub_sel <- meta_df[sel, subject_col]
          if (any(table(sub_sel) > 1L)) strata_use <- sub_sel
        }

        covar_term <- if (has_group) group_col else NULL

        run_two_level_adonis(
          idx            = sel,
          fac            = fac_vec[sel],
          covar          = covar_term,
          strata         = strata_use,
          scope_label    = "each_vs_rest",
          contrast_label = paste(lvl, "vs other"),
          term_label     = time_col
        )
      }
    }
  }

  # baseline
  if ("baseline" %in% contrasts) {
    if (is.null(bl)) {
      .ph_warn(
        "Baseline contrast requested but `baseline_level` not properly specified; skipping baseline tests.",
        step = "compute_permanova"
      )
    } else {
      # baseline for time
      if (has_time && bl %in% meta_df[[time_col]]) {
        fac_vec <- factor(
          ifelse(meta_df[[time_col]] == bl, bl, "not_baseline"),
          levels = c(bl, "not_baseline")
        )
        sel <- which(!is.na(fac_vec))

        strata_use <- NULL
        if (has_subject && subject_col %in% names(meta_df)) {
          sub_sel <- meta_df[sel, subject_col]
          if (any(table(sub_sel) > 1L)) strata_use <- sub_sel
        }

        covar_term <- if (has_group) group_col else NULL

        run_two_level_adonis(
          idx            = sel,
          fac            = fac_vec[sel],
          covar          = covar_term,
          strata         = strata_use,
          scope_label    = "baseline",
          contrast_label = paste(bl, "vs others"),
          term_label     = time_col
        )
      }

      # baseline for group
      if (has_group && bl %in% meta_df[[group_col]]) {
        fac_vec <- factor(
          ifelse(meta_df[[group_col]] == bl, bl, "not_baseline"),
          levels = c(bl, "not_baseline")
        )
        sel <- which(!is.na(fac_vec))

        covar_term <- if (has_time) time_col else NULL

        run_two_level_adonis(
          idx            = sel,
          fac            = fac_vec[sel],
          covar          = covar_term,
          strata         = NULL,
          scope_label    = "baseline",
          contrast_label = paste(bl, "vs others"),
          term_label     = group_col
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # combine and return
  # ---------------------------------------------------------------------------
  result_df <- if (length(results_list) > 0L) {
    dplyr::bind_rows(results_list)
  } else {
    tibble::tibble()
  }

  result_df
}

#' @title Test Homogeneity of Dispersion (Beta Dispersion)
#' @description Computes distances of samples to group centroids (using
#'   \code{vegan::betadisper}) and tests for differences in dispersion among
#'   groups or time levels. Optionally performs pairwise and other post-hoc
#'   tests on dispersion.
#'
#' @param dist_obj A \code{dist} object of sample distances (e.g. from
#'   \code{compute_distance()}).
#' @param ps A \code{phip_data} object or a table providing sample-level
#'   metadata. This table must contain \code{sample_id} and the
#'   columns specified in \code{group_col} and/or \code{time_col}.
#' @param group_col Name of the group factor column in \code{ps}
#'   (between-subjects). Use \code{NULL} if no group factor.
#' @param time_col Name of the time factor column in \code{ps}
#'   (within-subjects, categorical only). Use \code{NULL} if not applicable.
#' @param subject_col Name of subject identifier column (for reference only;
#'   not used directly in dispersion test calculations, but kept for API
#'   consistency). Default \code{"subject_id"}.
#' @param permutations Number of permutations for significance testing in
#'   \code{vegan::permutest}. Default 999.
#' @param contrasts Which dispersion contrasts to perform. Options:
#'   \code{"none"} (default), \code{"pairwise"}, \code{"each_vs_rest"},
#'   \code{"baseline"}. Interpretation analogous to \code{compute_permanova},
#'   but applied to dispersion.
#' @param baseline_level If \code{contrasts} includes \code{"baseline"},
#'   specify the baseline level of group or time to compare against others.
#'
#' @return A list of class \code{"beta_dispersion"} with:
#' \item{distances}{Tibble of per-sample distances to centroid. Columns:
#'   \code{sample_id}, \code{distance}, \code{level} (group/time level for
#'   a given scope), \code{scope} (e.g. \code{"group"}, \code{"time"},
#'   \code{"group:time"}), \code{contrast} (e.g. \code{"<global>"},
#'   \code{"A vs B"}).}
#' \item{tests}{Tibble of dispersion test results. Columns: \code{scope},
#'   \code{contrast}, \code{term = "dispersion"}, \code{p_value}, \code{n_perm}.}
#'
#' @examples
#' \dontrun{
#'   dispersion_res <- compute_dispersion(
#'     dist_bc,
#'     ps        = ps,
#'     group_col = "type_person",
#'     time_col  = "timepoint",
#'     contrasts = "pairwise"
#'   )
#'   dispersion_res$tests
#'   head(dispersion_res$distances)
#' }
#' @export
compute_dispersion <- function(dist_obj,
                               ps,
                               group_col = NULL,
                               time_col = NULL,
                               subject_col = "subject_id",
                               permutations = 999,
                               contrasts = "none",
                               baseline_level = NULL) {
  # ----------------------------------------------------------------------------
  # input validation (chk)
  # ----------------------------------------------------------------------------
  chk::chk_s3_class(dist_obj, "dist")
  if (!is.null(group_col)) chk::chk_string(group_col)
  if (!is.null(time_col)) chk::chk_string(time_col)
  chk::chk_string(subject_col)
  chk::chk_count(permutations)
  chk::chk_gt(permutations, 0)
  chk::chk_character(contrasts)
  if (!is.null(baseline_level)) chk::chk_string(baseline_level)

  # if ps is <phip_data>, overwrite ps with ps$data_long;
  # otherwise treat ps as data_long
  if (inherits(ps, "phip_data")) {
    chk::chk_not_null(ps$data_long)
    ps <- ps$data_long
  }

  # ---------------------------------------------------------------------------
  # prepare metadata from ps and align to dist labels
  # ---------------------------------------------------------------------------
  d_full      <- dist_obj
  labels_full <- attr(d_full, "Labels")

  dat <- ps
  if (is.null(dat)) {
    .ph_abort("`ps` is missing. Cannot construct metadata.",
              step = "compute_dispersion")
  }

  dat_cols <- dplyr::tbl_vars(dat)
  if (!"sample_id" %in% dat_cols) {
    .ph_abort("`ps` must contain a `sample_id` column.",
              step = "compute_dispersion")
  }

  has_group <- !is.null(group_col) && group_col %in% dat_cols
  has_time  <- !is.null(time_col)  && time_col  %in% dat_cols

  if (!is.null(group_col) && !has_group) {
    .ph_abort(
      paste0("Column `", group_col, "` not found in `ps`."),
      step = "compute_dispersion"
    )
  }
  if (!is.null(time_col) && !has_time) {
    .ph_abort(
      paste0("Column `", time_col, "` not found in `ps`."),
      step = "compute_dispersion"
    )
  }

  cols_needed <- c(
    "sample_id",
    if (has_group) group_col else character(0L),
    if (has_time)  time_col  else character(0L)
  )
  cols_needed <- unique(cols_needed)

  # build metadata from ps
  .ph_log_info(
    "building metadata from `ps`.",
    step = "compute_dispersion"
  )

  meta_all <- dat |>
    dplyr::select(dplyr::all_of(cols_needed)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect() |>
    as.data.frame()

  if (nrow(meta_all) == 0L) {
    .ph_abort(
      "constructed metadata has zero rows. check that `ps` is not empty.",
      step = "compute_dispersion"
    )
  }
  rownames(meta_all) <- meta_all$sample_id

  # align metadata to distance labels
  if (!is.null(labels_full)) {
    idx_align <- match(labels_full, meta_all$sample_id)
    missing_samples <- labels_full[is.na(idx_align)]
    if (length(missing_samples) > 0L) {
      .ph_abort(
        paste0(
          "the following samples from `dist_obj` are missing in `ps`: ",
          paste(missing_samples, collapse = ", ")
        ),
        step = "compute_dispersion"
      )
    }
    meta_sub <- meta_all[idx_align, , drop = FALSE]
    rownames(meta_sub) <- labels_full
  } else {
    meta_sub   <- meta_all
    labels_full <- rownames(meta_sub)
    .ph_warn(
      "no labels found in `dist_obj`; assuming metadata row order matches the distance order.",
      step = "compute_dispersion"
    )
  }

  # coerce grouping vars to factor
  if (has_group) {
    meta_sub[[group_col]] <- as.factor(meta_sub[[group_col]])
  }
  if (has_time) {
    meta_sub[[time_col]] <- as.factor(meta_sub[[time_col]])
  }

  # ---------------------------------------------------------------------------
  # drop samples with na in group/time and subset distance matrix
  # ---------------------------------------------------------------------------
  vars_for_na <- c(
    if (has_group) group_col else character(0L),
    if (has_time)  time_col  else character(0L)
  )
  vars_for_na <- unique(vars_for_na)

  if (length(vars_for_na) > 0L) {
    keep <- stats::complete.cases(meta_sub[, vars_for_na, drop = FALSE])
  } else {
    keep <- rep(TRUE, nrow(meta_sub))
  }

  if (!all(keep)) {
    dropped <- sum(!keep)
    .ph_log_info(
      paste0("dropping ", dropped,
             " samples with missing values in dispersion grouping variables."),
      step = "compute_dispersion"
    )
  }

  meta_df <- meta_sub[keep, , drop = FALSE]
  if (nrow(meta_df) == 0L) {
    .ph_abort(
      "all samples have missing values in grouping variables; cannot run dispersion tests.",
      step = "compute_dispersion"
    )
  }
  keep_labels <- rownames(meta_df)

  d_mat_full <- as.matrix(d_full)
  d_mat_sub  <- d_mat_full[keep_labels, keep_labels, drop = FALSE]
  d          <- stats::as.dist(d_mat_sub)
  labels     <- attr(d, "Labels")

  # re-evaluate factor presence after na drop
  has_group <- has_group && length(unique(meta_df[[group_col]])) > 1L
  has_time  <- has_time  && length(unique(meta_df[[time_col]]))  > 1L

  # ---------------------------------------------------------------------------
  # prepare collectors for distances and tests
  # ---------------------------------------------------------------------------
  distances_list <- list()
  tests_list     <- list()

  add_distance_rows <- function(sample_ids, dists, level_vals,
                                scope_val, contrast_val) {
    distances_list[[length(distances_list) + 1L]] <<- tibble::tibble(
      sample_id = sample_ids,
      distance  = dists,
      level     = level_vals,
      scope     = scope_val,
      contrast  = contrast_val
    )
  }

  add_test_row <- function(scope_val, contrast_val, p_val) {
    tests_list[[length(tests_list) + 1L]] <<- tibble::tibble(
      scope    = scope_val,
      contrast = contrast_val,
      term     = "dispersion",
      p_value  = p_val,
      n_perm   = permutations
    )
  }

  # ---------------------------------------------------------------------------
  # global dispersion tests (group, time, group:time)
  # ---------------------------------------------------------------------------
  if (has_group) {
    fac <- factor(meta_df[[group_col]])
    bd <- try(vegan::betadisper(d, fac), silent = TRUE)
    if (inherits(bd, "try-error")) {
      .ph_warn(
        "betadisper failed for group factor; skipping group dispersion test.",
        step = "compute_dispersion"
      )
    } else {
      add_distance_rows(
        sample_ids = names(bd$distances),
        dists      = as.numeric(bd$distances),
        level_vals = as.character(bd$group),
        scope_val  = "group",
        contrast_val = "<global>"
      )
      pt <- vegan::permutest(bd, permutations = permutations, parallel = 1)
      pval <- tryCatch(pt$tab[1, "Pr(>F)"], error = function(e) NA_real_)
      add_test_row("group", "<global>", pval)
    }
  }

  if (has_time) {
    if (!is.numeric(meta_df[[time_col]])) {
      fac <- factor(meta_df[[time_col]])
      bd <- try(vegan::betadisper(d, fac), silent = TRUE)
      if (inherits(bd, "try-error")) {
        .ph_warn(
          "betadisper failed for time factor; skipping time dispersion test.",
          step = "compute_dispersion"
        )
      } else {
        add_distance_rows(
          sample_ids = names(bd$distances),
          dists      = as.numeric(bd$distances),
          level_vals = as.character(bd$group),
          scope_val  = "time",
          contrast_val = "<global>"
        )
        pt <- vegan::permutest(bd, permutations = permutations, parallel = 1)
        pval <- tryCatch(pt$tab[1, "Pr(>F)"], error = function(e) NA_real_)
        add_test_row("time", "<global>", pval)
      }
    } else {
      .ph_warn(
        "`time_col` is numeric; continuous dispersion by time not supported. skipping time dispersion test.",
        step = "compute_dispersion"
      )
    }
  }

  if (has_group && has_time && !is.numeric(meta_df[[time_col]])) {
    inter_factor <- factor(
      paste(meta_df[[group_col]], meta_df[[time_col]], sep = " * ")
    )
    if (length(unique(inter_factor)) > 1L) {
      bd <- try(vegan::betadisper(d, inter_factor), silent = TRUE)
      if (!inherits(bd, "try-error")) {
        add_distance_rows(
          sample_ids = names(bd$distances),
          dists      = as.numeric(bd$distances),
          level_vals = as.character(bd$group),
          scope_val  = "group:time",
          contrast_val = "<global>"
        )
        pt <- vegan::permutest(bd, permutations = permutations, parallel = 1)
        pval <- tryCatch(pt$tab[1, "Pr(>F)"], error = function(e) NA_real_)
        add_test_row("group:time", "<global>", pval)
      }
    }
  }

  # ---------------------------------------------------------------------------
  # helper for subset dispersion tests
  # ---------------------------------------------------------------------------
  run_disp_test <- function(idx, fac_vec, scope_lab, contrast_lab) {
    if (length(unique(fac_vec)) < 2L || min(table(fac_vec)) < 2L) {
      .ph_log_info(
        paste("skipping dispersion test for", contrast_lab, "- not enough data."),
        step = "compute_dispersion"
      )
      return(NULL)
    }
    d_mat <- as.matrix(d)
    d_sub <- stats::as.dist(d_mat[idx, idx, drop = FALSE])
    bd <- try(vegan::betadisper(d_sub, fac_vec), silent = TRUE)
    if (inherits(bd, "try-error")) return(NULL)

    add_distance_rows(
      sample_ids = names(bd$distances),
      dists      = as.numeric(bd$distances),
      level_vals = as.character(bd$group),
      scope_val  = scope_lab,
      contrast_val = contrast_lab
    )
    pt <- vegan::permutest(bd, permutations = permutations, parallel = 1)
    pval <- tryCatch(pt$tab[1, "Pr(>F)"], error = function(e) NA_real_)
    add_test_row(scope_lab, contrast_lab, pval)
  }

  # ---------------------------------------------------------------------------
  # post-hoc contrasts
  # ---------------------------------------------------------------------------
  contrasts <- tolower(unique(contrasts))
  contrasts[contrasts == "group_vs_rest"] <- "each_vs_rest"

  # pairwise
  if ("pairwise" %in% contrasts) {
    # groups
    if (has_group) {
      grps <- na.omit(unique(meta_df[[group_col]]))
      if (length(grps) > 1L) {
        for (pair in utils::combn(grps, 2, simplify = FALSE)) {
          sel <- which(meta_df[[group_col]] %in% pair)
          fac <- factor(meta_df[sel, group_col], levels = pair)
          run_disp_test(sel, fac_vec = fac, scope_lab = "group",
                        contrast_lab = paste(pair, collapse = " vs "))
        }
      }
    }
    # time
    if (has_time && !is.numeric(meta_df[[time_col]])) {
      times <- na.omit(unique(meta_df[[time_col]]))
      if (length(times) > 1L) {
        for (pair in utils::combn(times, 2, simplify = FALSE)) {
          sel <- which(meta_df[[time_col]] %in% pair)
          fac <- factor(meta_df[sel, time_col], levels = pair)
          run_disp_test(sel, fac_vec = fac, scope_lab = "time",
                        contrast_lab = paste(pair, collapse = " vs "))
        }
      }
    }
  }

  # each vs rest
  if ("each_vs_rest" %in% contrasts) {
    if (has_group && length(unique(meta_df[[group_col]])) > 1L) {
      for (lvl in unique(meta_df[[group_col]])) {
        fac_vec <- factor(
          ifelse(meta_df[[group_col]] == lvl, lvl, "other"),
          levels = c(lvl, "other")
        )
        sel <- which(!is.na(fac_vec))
        run_disp_test(sel, fac_vec = fac_vec[sel], scope_lab = "group",
                      contrast_lab = paste(lvl, "vs rest"))
      }
    }
    if (has_time && !is.numeric(meta_df[[time_col]]) &&
        length(unique(meta_df[[time_col]])) > 1L) {
      for (lvl in unique(meta_df[[time_col]])) {
        fac_vec <- factor(
          ifelse(meta_df[[time_col]] == lvl, lvl, "other"),
          levels = c(lvl, "other")
        )
        sel <- which(!is.na(fac_vec))
        run_disp_test(sel, fac_vec = fac_vec[sel], scope_lab = "time",
                      contrast_lab = paste(lvl, "vs rest"))
      }
    }
  }

  # baseline
  if ("baseline" %in% contrasts) {
    if (is.null(baseline_level)) {
      .ph_warn(
        "Baseline contrast requested but no `baseline_level` provided; skipping.",
        step = "compute_dispersion"
      )
    } else {
      # baseline for group
      if (has_group && baseline_level %in% meta_df[[group_col]]) {
        fac_vec <- factor(
          ifelse(meta_df[[group_col]] == baseline_level,
                 baseline_level, "other"),
          levels = c(baseline_level, "other")
        )
        sel <- which(!is.na(fac_vec))
        run_disp_test(sel, fac_vec = fac_vec[sel], scope_lab = "group",
                      contrast_lab = paste(baseline_level, "vs others"))
      }
      # baseline for time
      if (has_time && !is.numeric(meta_df[[time_col]]) &&
          baseline_level %in% meta_df[[time_col]]) {
        fac_vec <- factor(
          ifelse(meta_df[[time_col]] == baseline_level,
                 baseline_level, "other"),
          levels = c(baseline_level, "other")
        )
        sel <- which(!is.na(fac_vec))
        run_disp_test(sel, fac_vec = fac_vec[sel], scope_lab = "time",
                      contrast_lab = paste(baseline_level, "vs others"))
      }
    }
  }

  # ---------------------------------------------------------------------------
  # combine and return
  # ---------------------------------------------------------------------------
  distances_tbl <- if (length(distances_list)) {
    dplyr::bind_rows(distances_list)
  } else {
    tibble::tibble()
  }

  tests_tbl <- if (length(tests_list)) {
    dplyr::bind_rows(tests_list)
  } else {
    tibble::tibble()
  }

  result <- list(
    distances = distances_tbl,
    tests     = tests_tbl
  )
  class(result) <- "beta_dispersion"

  result
}


# ------------------------------------------------------------------------------
# internals (shared)
# ------------------------------------------------------------------------------
# normalizer used before distance
.beta_normalize <- function(m, # wide samples x ranks matrix
                            how, # the normalization method
                            used_exist_rule) {
  # when auto, do nothing; binary data actually encode the absence/presence well
  # and there is no need to scale anything
  if (identical(how, "auto")) {
    how <- if (isTRUE(used_exist_rule)) "none" else "relative"
  }
  if (identical(how, "none")) {
    return(m)
  } # return m unchanged

  # compute rowsums (per sample_id) and guard against 0 division; we actually
  # filter the 0s out before applyihg the normalization, but as a fallback
  rs <- rowSums(m, na.rm = TRUE)
  rs[rs == 0] <- 1
  rel <- m / rs # per sample_id prevalence; rows sum up to 1
  switch(how,
    relative  = rel, # each row sums up to 1
    hellinger = sqrt(rel),
    log       = log1p(m),
    rel
  )

  # the logic behind it:
  # * rel -> distances compare compositional profiles only (size-invariant)
  # * hellinger -> stabilize variance, downweight dominant features, and make
  #                euclidean geometry behave better on compositional data; rare
  #                features get a bit more voice, very abundant ones a bit less
  # * log -> big values get pulled closer to small ones; zeros stay zero;
  #          you don’t want a few huge counts to dominate
}

# supported by parallelDist; others --> vegan::vegdist fallback
.pd_supported <- function(method) {
  method <- tolower(method)
  method %in% c(
    "euclidean", "minkowski", "manhattan", "canberra",
    "binary", "maximum", "cosine", "chebyshev"
  )
}

# distance with parallelDist when possible; threaded Bray via L1 identity
.beta_dist <- function(mat, # wide, transformed matrix
                       method, # which distance
                       n_threads) { # how many threads to use
  method <- tolower(method) # case-insensitive

  # custom threaded Bray–Curtis via L1 identity to make it multithreaded; the
  # default vegan bray is single-threaded and painfully slow
  if (method == "bray" && rlang::is_installed("parallelDist")) {
    # firstly compute manhattan
    dL1 <- parallelDist::parDist(mat,
      method = "manhattan",
      threads = n_threads
    )
    # then transform to bray
    s <- rowSums(mat)
    n <- length(s)
    den <- vector("numeric", n * (n - 1) / 2L)
    k <- 1L
    for (i in seq_len(n - 1L)) {
      ni <- n - i
      den[k:(k + ni - 1L)] <- s[i] + s[(i + 1L):n]
      k <- k + ni
    }
    d <- dL1
    d[] <- as.numeric(dL1) / den
    return(d)
  }

  # other methods that have support in the parallelDist
  if (.pd_supported(method) && rlang::is_installed("parallelDist")) {
    return(parallelDist::parDist(mat, method = method, threads = n_threads))
  }

  vegan::vegdist(mat, method = method)
}

# long --> wide matrix transformer (fill=0); data.table fast path if installed
.long_to_wide_matrix <- function(agg_df) {
  # data.table path is suggested --> fast and efficient for bigger data
  if (rlang::is_installed("data.table")) {
    DT <- data.table::as.data.table(agg_df) # convert the input

    # force sample_id to character
    data.table::set(DT,
      j = "sample_id",
      value = as.character(DT[["sample_id"]])
    )

    # pivot to wide with one row per sample_id and one column per rank_val;
    # fill with the "abund" if present, if not --> fill with 0s
    wide_dt <- data.table::dcast(DT, sample_id ~ rank_val,
      value.var = "abund",
      fill = 0L
    )
    wide_df <- as.data.frame(wide_dt) # convert back to data.frame
    rn <- wide_df[["sample_id"]]

    # get uniques and return
    mat <- as.matrix(wide_df[, setdiff(names(wide_df), "sample_id"),
      drop = FALSE
    ])
    rownames(mat) <- rn
    return(mat)
  } else {
    # the same logic but with tidyr
    wide <- tidyr::pivot_wider(agg_df,
      names_from = "rank_val",
      values_from = "abund", values_fill = 0
    )
    rn <- as.character(wide$sample_id)
    mat <- as.matrix(wide[, setdiff(names(wide), "sample_id"), drop = FALSE])
    rownames(mat) <- rn
    return(mat)
  }
}

#' Compute t-SNE embeddings for samples
#'
#' @description
#' Runs t-distributed stochastic neighbour embedding (t-SNE) on a
#' sample-wise distance matrix (typically returned by [compute_distance()])
#' and returns a tibble with t-SNE coordinates and selected sample-level
#' metadata.
#'
#' The function expects that rows/columns (or labels) of the distance
#' object correspond to `sample_id`s present in `ps$data_long`.
#'
#' @param ps A `phip_data` object.
#' @param dist_obj A sample-wise distance object. Can be either:
#'   * a [`dist`] object (e.g. from [compute_distance()]), or
#'   * a numeric, symmetric matrix with row/column names giving `sample_id`s.
#' @param dims Integer number of t-SNE dimensions to compute (2 or 3).
#'   Defaults to `3L` so that both 2D and 3D plots can be generated from
#'   the same result.
#' @param perplexity Numeric perplexity parameter passed to [Rtsne::Rtsne()].
#'   Must be smaller than the number of samples; if it is too large, it is
#'   automatically reduced with a warning.
#' @param theta Speed/accuracy tradeoff for the Barnes-Hut approximation,
#'   passed to [Rtsne::Rtsne()]. Defaults to `0.5`.
#' @param max_iter Maximum number of iterations for t-SNE, passed to
#'   [Rtsne::Rtsne()]. Defaults to `1000L`.
#' @param meta_cols Optional character vector of column names in
#'   `ps$data_long` to attach as sample-level metadata. If `NULL`
#'   (default), the function tries to use `ps$meta$extra_cols` and keeps
#'   the intersection with columns available in `ps$data_long`.
#' @param seed Optional integer random seed for reproducibility. If
#'   provided, the function temporarily sets the seed for the duration of
#'   the t-SNE computation and restores the previous RNG state afterwards.
#' @param check_duplicates Logical; passed to [Rtsne::Rtsne()]. For
#'   distance input, duplicates are not expected, so the default is
#'   `FALSE`.
#' @param ... Additional arguments passed on to [Rtsne::Rtsne()].
#'
#' @return
#' A tibble with class `c("phip_tsne", "tbl_df", "tbl", "data.frame")`
#' containing at least:
#' \describe{
#'   \item{sample_id}{Sample identifier (from distance labels or matrix row names).}
#'   \item{tSNE1, tSNE2}{t-SNE coordinates for the first two dimensions.}
#'   \item{tSNE3}{Third t-SNE dimension if `dims >= 3`, otherwise `NA`.}
#' }
#'
#' Additional columns specified in `meta_cols` are attached by a left join
#' on `sample_id` (one unique row per `sample_id` is enforced).
#'
#' Attributes:
#' * `"distance"`: the original `dist_obj` as supplied.
#' * `"tsne_params"`: a list of key t-SNE parameters and the original call.
#' * `"meta_cols"`: the character vector of metadata columns actually used.
#'
#' @details
#' This function runs t-SNE in *distance mode* (`is_distance = TRUE`), using
#' the supplied distance object directly. Distance computation itself is
#' handled elsewhere, typically via [compute_distance()].
#'
#' @examples
#' \dontrun{
#' # Assume ps is a <phip_data> and D is a distance from compute_distance():
#' D <- compute_distance(ps, value_col = "fold_change")
#'
#' tsne_res <- compute_tsne(
#'   ps       = ps,
#'   dist_obj = D,
#'   dims     = 3L,
#'   perplexity = 30,
#'   meta_cols  = c("type_person", "sex", "age")
#' )
#'
#' # 2D plot, coloured by type_person
#' plot_tsne(tsne_res, view = "2d", colour = "type_person")
#'
#' # 3D interactive plot (requires plotly)
#' plot_tsne(tsne_res, view = "3d", colour = "type_person")
#' }
#'
#' @export
compute_tsne <- function(ps,
                         dist_obj,
                         dims = 3L,
                         perplexity = 30,
                         theta = 0.5,
                         max_iter = 1000L,
                         meta_cols = NULL,
                         seed = NULL,
                         check_duplicates = FALSE,
                         ...) {
  # ---------------------------------------------------------------------------
  # Basic checks
  # ---------------------------------------------------------------------------
  if (!inherits(ps, "phip_data")) {
    .ph_abort("Input `ps` must be a <phip_data> object.",
              step = "compute_tsne")
  }

  if (missing(dist_obj)) {
    .ph_abort("Argument `dist_obj` is required. Run `compute_distance()` first.",
              step = "compute_tsne")
  }

  if (!inherits(dist_obj, "dist") && !is.matrix(dist_obj)) {
    .ph_abort(
      "Argument `dist_obj` must be either a <dist> object or a numeric matrix.",
      step = "compute_tsne"
    )
  }

  if (!dims %in% c(2L, 3L)) {
    .ph_abort("Argument `dims` must be either 2 or 3.",
              step = "compute_tsne")
  }

  if (!rlang::is_installed("Rtsne")) {
    .ph_abort(
      "Package 'Rtsne' is required to compute t-SNE. Please install it.",
      step = "compute_tsne"
    )
  }

  # ---------------------------------------------------------------------------
  # Extract labels and size
  # ---------------------------------------------------------------------------
  if (inherits(dist_obj, "dist")) {
    n_samples <- attr(dist_obj, "Size")
    labels <- attr(dist_obj, "Labels")

    if (is.null(labels)) {
      # Fallback: integer sequence if labels are missing
      labels <- as.character(seq_len(n_samples))
      .ph_warn(
        "Distance object has no labels; using integer indices as `sample_id`.",
        step = "compute_tsne"
      )
    }
  } else {
    # matrix input
    if (!is.numeric(dist_obj)) {
      .ph_abort("If `dist_obj` is a matrix, it must be numeric.",
                step = "compute_tsne")
    }
    if (nrow(dist_obj) != ncol(dist_obj)) {
      .ph_abort("Distance matrix must be square (same number of rows and columns).",
                step = "compute_tsne")
    }
    n_samples <- nrow(dist_obj)
    labels <- rownames(dist_obj)
    if (is.null(labels)) {
      labels <- as.character(seq_len(n_samples))
      .ph_warn(
        "Distance matrix has no rownames; using integer indices as `sample_id`.",
        step = "compute_tsne"
      )
    }
  }

  if (n_samples < 3L) {
    .ph_abort(
      "t-SNE requires at least 3 samples.",
      step = "compute_tsne"
    )
  }

  if (length(labels) != n_samples) {
    .ph_abort(
      "Length of distance labels does not match the number of samples.",
      step = "compute_tsne"
    )
  }

  if (anyDuplicated(labels)) {
    .ph_abort(
      "Distance labels (sample IDs) must be unique.",
      step = "compute_tsne"
    )
  }

  # ---------------------------------------------------------------------------
  # Perplexity sanity check
  # ---------------------------------------------------------------------------
  max_perplexity <- floor((n_samples - 1L) / 3L)
  if (perplexity >= n_samples) {
    .ph_abort(
      paste0(
        "Perplexity (", perplexity, ") must be smaller than the number of samples (",
        n_samples, ")."
      ),
      step = "compute_tsne"
    )
  }
  if (perplexity > max_perplexity) {
    .ph_warn(
      paste0(
        "Perplexity (", perplexity, ") is high for n = ", n_samples,
        "; reducing to ", max_perplexity, "."
      ),
      step = "compute_tsne"
    )
    perplexity <- max_perplexity
  }

  # ---------------------------------------------------------------------------
  # Prepare metadata column selection
  # ---------------------------------------------------------------------------
  dat <- ps$data_long

  if (!is.null(dat)) {
    dat_cols <- dplyr::tbl_vars(dat)
  } else {
    dat_cols <- character()
    .ph_warn(
      "`ps$data_long` is NULL; no metadata columns can be attached.",
      step = "compute_tsne"
    )
  }

  if (is.null(meta_cols)) {
    # Start from meta$extra_cols if available, then intersect with data_long
    extra <- ps$meta$extra_cols %||% character()
    meta_cols <- intersect(extra, dat_cols)
  } else {
    # User-specified meta_cols: keep only those present in data_long
    missing_cols <- setdiff(meta_cols, dat_cols)
    if (length(missing_cols) > 0L) {
      .ph_warn(
        paste0(
          "The following `meta_cols` are not present in `ps$data_long` and ",
          "will be ignored: ",
          paste(missing_cols, collapse = ", ")
        ),
        step = "compute_tsne"
      )
      meta_cols <- intersect(meta_cols, dat_cols)
    }
  }

  # ---------------------------------------------------------------------------
  # Run Rtsne in distance mode
  # ---------------------------------------------------------------------------
  .ph_log_info(
    paste0("Running t-SNE with dims = ", dims,
           ", perplexity = ", perplexity,
           " on ", n_samples, " samples (distance input)."),
    step = "compute_tsne"
  )

  # Temporary RNG seed control
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit({
      if (exists("old_seed", inherits = FALSE)) {
        .Random.seed <<- old_seed
      }
    }, add = TRUE)
    set.seed(seed)
  }

  tsne_fit <- if (inherits(dist_obj, "dist")) {
    Rtsne::Rtsne(
      dist_obj,
      is_distance      = TRUE,
      dims             = dims,
      perplexity       = perplexity,
      theta            = theta,
      max_iter         = max_iter,
      check_duplicates = check_duplicates,
      ...
    )
  } else {
    # matrix input
    Rtsne::Rtsne(
      as.matrix(dist_obj),
      is_distance      = TRUE,
      dims             = dims,
      perplexity       = perplexity,
      theta            = theta,
      max_iter         = max_iter,
      check_duplicates = check_duplicates,
      ...
    )
  }

  coords <- tsne_fit$Y

  if (!is.matrix(coords) || nrow(coords) != n_samples || ncol(coords) != dims) {
    .ph_abort(
      "Unexpected output from Rtsne::Rtsne(): coordinate matrix has wrong shape.",
      step = "compute_tsne"
    )
  }

  # ---------------------------------------------------------------------------
  # Build base result tibble
  # ---------------------------------------------------------------------------
  tsne_df <- tibble::tibble(
    sample_id = labels,
    tSNE1     = coords[, 1],
    tSNE2     = if (dims >= 2L) coords[, 2L] else NA_real_,
    tSNE3     = if (dims >= 3L) coords[, 3L] else NA_real_
  )

  # ---------------------------------------------------------------------------
  # Attach metadata (if requested / available)
  # ---------------------------------------------------------------------------
  if (!is.null(dat) && length(meta_cols) > 0L) {
    .ph_log_info(
      paste0(
        "Attaching metadata columns to t-SNE result: ",
        paste(meta_cols, collapse = ", ")
      ),
      step = "compute_tsne"
    )

    meta_tbl <- dat |>
      dplyr::select(sample_id, dplyr::all_of(meta_cols)) |>
      dplyr::distinct() |>
      dplyr::collect()

    if (anyDuplicated(meta_tbl$sample_id)) {
      .ph_warn(
        paste0(
          "Metadata rows are not unique per `sample_id`. ",
          "Keeping the first row per sample."
        ),
        step = "compute_tsne"
      )

      meta_tbl <- meta_tbl |>
        dplyr::group_by(sample_id) |>
        dplyr::slice(1L) |>
        dplyr::ungroup()
    }

    tsne_df <- tsne_df |>
      dplyr::left_join(meta_tbl, by = "sample_id")
  }

  .ph_log_info("t-SNE embedding computation finished.",
               step = "compute_tsne")

  # ---------------------------------------------------------------------------
  # Attach attributes and class
  # ---------------------------------------------------------------------------
  attr(tsne_df, "distance") <- dist_obj
  attr(tsne_df, "tsne_params") <- list(
    dims        = dims,
    perplexity  = perplexity,
    theta       = theta,
    max_iter    = max_iter,
    seed        = seed,
    check_dup   = check_duplicates,
    call        = match.call()
  )
  attr(tsne_df, "meta_cols") <- meta_cols

  class(tsne_df) <- c("phip_tsne", class(tsne_df))

  tsne_df
}
