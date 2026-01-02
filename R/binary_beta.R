#' @title Compute Pairwise Sample Distances
#'
#' @description This function builds a sample-by-feature abundance matrix from a
#' \code{phip_data} object (using \code{ps$data_long}), optionally normalizes
#' the matrix, and then computes pairwise distances between samples.
#'
#' The normalized abundance matrix used for distance calculation is attached
#' to the returned \code{dist} object as attribute \code{"abundances"}.
#'
#' Note: this function pivots to a wide matrix in the database (via dbplyr)
#' and then collects the result into memory. This can be large for big cohorts
#' and/or large peptide sets.
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
#' \code{exist},  \code{counts_hits}, \code{counts_control}, \code{fold_change}.
#'
#' @param method_normalization character scalar. Normalization applied to the
#' abundance matrix before distance computation. One of:
#' \itemize{
#'   \item \code{"auto"}: uses \code{"none"} for binary (0/1) data, otherwise
#'     uses \code{"relative"}.
#'   \item \code{"relative"}: divide each row by its row sum.
#'   \item \code{"hellinger"}: \code{sqrt(relative)}.
#'   \item \code{"log"}: \code{log1p} transform.
#'   \item \code{"none"}: no normalization.
#' }
#'
#' @param distance character scalar. Distance method name. The string is
#'   lowercased internally.
#'
#'   Supported methods depend on which packages are installed:
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
#'   If 'parallelDist' is installed but the requested method is not in the fast
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
#' # build an example <phip_data> object from the package example dataset
#' ps <- phip_load_example_data()
#'
#' # small subset for speed: 5 peptides at time t1
#' keep_pep <- c("16627", "5243", "24799", "16196", "18003")
#' dat_cols <- dplyr::tbl_vars(ps$data_long)
#' tp_col <- "time"
#'
#' ps_small <- ps |>
#'   dplyr::filter(
#'     peptide_id %in% keep_pep,
#'     !!rlang::sym(tp_col) == "T1"
#'   ) |>
#'   dplyr::collect()
#'
#' # compute distances (needs either 'parallelDist' or 'vegan')
#' val_col <- "exist"
#'
#' d <- compute_distance(
#'   ps_small,
#'   value_col = val_col,
#'   method_normalization = "hellinger",
#'   distance = "bray",
#'   n_threads = 2L
#' )
#'
#' a <- attr(d, "abundances")
#' a[1:min(5, nrow(a)), 1:min(5, ncol(a)), drop = FALSE]
#' @export
compute_distance <- function(ps,
                             value_col = NULL,
                             method_normalization = c(
                               "auto", "relative",
                               "hellinger", "log",
                               "none"
                             ),
                             distance = "bray",
                             n_threads = 1L) {

  # ----------------------------------------------------------------------------
  # input validation (chk)
  # ----------------------------------------------------------------------------
  if (!is.null(value_col)) {
    chk::chk_string(value_col)
  }
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
          "could not infer an abundance column in `ps`. ",
          "tried: ", paste(candidates, collapse = ", "),
          ". please specify `value_col` explicitly."
        )
      )
    }
    value_col <- hit[1L]
    .ph_log_info(
      paste0("auto-detected `value_col = \"", value_col, "\"` from `ps`.")
    )
  }

  if (!value_col %in% dat_cols) {
    .ph_abort(
      paste0("column `", value_col, "` not found in `ps`.")
    )
  }

  required_cols <- c("sample_id", "peptide_id", value_col)
  missing_cols <- setdiff(required_cols, dat_cols)
  if (length(missing_cols) > 0L) {
    .ph_abort(
      paste0(
        "missing required column(s) in `ps`: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  .ph_log_info(
    paste0("building abundance matrix from `ps` using `", value_col, "`.")
  )

  value_sym <- rlang::sym(value_col)

  # ----------------------------------------------------------------------------
  # 2) build pivot spec and pivot in dbplyr
  # ----------------------------------------------------------------------------
  .ph_log_info("building pivot spec (sample_id x peptide_id).")

  id_levels <- dat |>
    dplyr::distinct(sample_id, peptide_id) |>
    dplyr::collect()

  id_levels[[value_col]] <- 0

  # basic sanity checks on collected ids
  if (anyNA(id_levels$sample_id) || anyNA(id_levels$peptide_id)) {
    .ph_abort(
      "`ps` contains missing values in `sample_id` and/or
      `peptide_id`."
    )
  }

  # duplicates will break pivot_wider (or create list-cols); require uniqueness
  dup_check <- dat |>
    dplyr::count(sample_id, peptide_id, name = "n") |>
    dplyr::filter(.data$n > 1L) |>
    dplyr::collect()

  if (nrow(dup_check) > 0L) {
    .ph_abort(
      "found duplicated (sample_id, peptide_id) pairs in `ps`."
    )
  }

  pivot_spec <- tidyr::build_wider_spec(
    data = id_levels,
    names_from = peptide_id,
    values_from = !!value_sym
  )

  .ph_log_info("pivoting to wide abundance matrix in db.")

  wide_df <- dbplyr::dbplyr_pivot_wider_spec(
    data = dat,
    spec = pivot_spec,
    id_cols = sample_id,
    values_fill = 0
  ) |>
    dplyr::collect()

  if (!"sample_id" %in% names(wide_df)) {
    .ph_abort("failed to construct wide abundance table (no `sample_id`).")
  }

  value_cols <- setdiff(names(wide_df), "sample_id")
  if (!all(vapply(wide_df[value_cols], is.numeric, logical(1)))) {
    .ph_abort(
      paste0("`", value_col, "` must be numeric after collect().")
    )
  }

  mat <- wide_df |>
    tibble::column_to_rownames("sample_id") |>
    as.matrix()

  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    .ph_abort(
      "abundance matrix is empty after reshaping. check filters and
      `value_col`."
    )
  }

  .ph_log_info(
    paste0(
      "abundance matrix has ", nrow(mat), " samples and ",
      ncol(mat), " features after preprocessing."
    )
  )

  # ----------------------------------------------------------------------------
  # 3) normalization
  # ----------------------------------------------------------------------------
  normalize_abundance <- function(mat, method_normalization) {
    if (identical(method_normalization, "auto")) {
      vals <- mat[!is.na(mat)]
      is_binary_data <- length(vals) > 0L && all(vals == 0 | vals == 1)
      method_normalization <- if (is_binary_data) {
        "none"
      } else {
        "relative"
      }
      .ph_log_info(
        paste0("auto normalization selected -> using ", method_normalization)
      )
    }

    switch(method_normalization,
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
  }

  norm_mat <- normalize_abundance(mat, method_normalization)

  # ----------------------------------------------------------------------------
  # 4) distance computation
  # ----------------------------------------------------------------------------
  .ph_compute_distance <- function(norm_mat,
                                   dist_method,
                                   distance,
                                   n_threads) {
    if (dist_method == "chebyshev") {
      .ph_log_info(
        "Distance method 'chebyshev' mapped to 'maximum' for parallelDist."
      )
      dist_method <- "maximum"
    }

    .ph_log_info(
      paste0("computing distance: ", dist_method)
    )

    dist_obj <- NULL

    pd_methods <- c(
      "bray", "euclidean", "minkowski", "manhattan", "canberra",
      "binary", "maximum", "cosine"
    )

    veg_methods <- c(
      "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
      "jaccard", "gower", "altgower", "morisita", "horn", "mountford", "raup",
      "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger",
      "aitchison", "robust.aitchison"
    )

    use_pd <- dist_method %in% pd_methods
    is_veg_method <- !is.na(pmatch(dist_method, veg_methods))

    if (use_pd && rlang::is_installed("parallelDist")) {
      dist_obj <- parallelDist::parDist(
        norm_mat,
        method  = dist_method,
        threads = n_threads
      )
    } else if (is_veg_method) {
      if (use_pd) {
        .ph_log_info(
          paste0(
            "parallelDist not installed; using vegan::vegdist for method ",
            dist_method,
            "."
          )
        )
      }
      dist_obj <- vegan::vegdist(norm_mat, method = dist_method)
    } else if (use_pd) {
      .ph_abort(
        paste0(
          "`distance = \"", distance, "\"` requires the suggested package ",
          "'parallelDist'. Please install 'parallelDist' or choose a ",
          "vegan::vegdist() method."
        )
      )
    } else {
      dist_obj <- vegan::vegdist(norm_mat, method = dist_method)
    }

    .ph_log_info("distance matrix computation complete.")

    dist_obj
  }

  dist_obj <- .ph_compute_distance(
    norm_mat = norm_mat,
    dist_method = dist_method,
    distance = distance,
    n_threads = n_threads
  )

  # attach normalized abundance matrix as attribute
  attr(dist_obj, "abundances") <- norm_mat

  dist_obj
}

.ph_feature_associations <- function(coords,
                                     X,
                                     feature_assoc,
                                     top_features,
                                     max_axes = NULL,
                                     warn_missing_ids,
                                     warn_insufficient_overlap) {
  if (identical(feature_assoc, "none")) {
    return(tibble::tibble())
  }

  coords_ids <- rownames(coords)
  X_ids <- rownames(X)
  if (is.null(coords_ids) || is.null(X_ids)) {
    .ph_warn(warn_missing_ids)
    return(tibble::tibble())
  }

  common_ids <- intersect(coords_ids, X_ids)
  if (length(common_ids) < 2L) {
    .ph_warn(warn_insufficient_overlap)
    return(tibble::tibble())
  }

  ax_limit <- ncol(coords)
  if (!is.null(max_axes)) {
    ax_limit <- min(ax_limit, max_axes)
  }
  if (ax_limit < 1L) {
    return(tibble::tibble())
  }

  ax_idx <- seq_len(ax_limit)
  U <- coords[common_ids, ax_idx, drop = FALSE]
  Xsub <- X[common_ids, , drop = FALSE]

  if (identical(feature_assoc, "weighted_average")) {
    w <- colSums(Xsub, na.rm = TRUE)
    keep_feats <- which(w > 0)
    if (length(keep_feats) > 0L) {
      S <- t(Xsub[, keep_feats, drop = FALSE]) %*% U
      S <- sweep(S, 1, w[keep_feats], "/")
      rownames(S) <- colnames(Xsub)[keep_feats]
      colnames(S) <- colnames(U)
    } else {
      S <- matrix(0, nrow = 0, ncol = ncol(U))
    }
  } else {
    X_center <- scale(Xsub, center = TRUE, scale = FALSE)
    U_center <- scale(U, center = TRUE, scale = FALSE)
    n_common <- nrow(X_center)
    var_X <- colSums(X_center^2) / (n_common - 1)
    keep_feats <- which(var_X > 0)

    if (length(keep_feats) > 0L) {
      cov_XU <- crossprod(
        X_center[, keep_feats, drop = FALSE],
        U_center
      ) / (n_common - 1)

      if (identical(feature_assoc, "correlation")) {
        sd_X <- sqrt(var_X[keep_feats])
        sd_U <- sqrt(colSums(U_center^2) / (n_common - 1))
        S <- sweep(cov_XU, 1, sd_X, "/")
        S <- sweep(S, 2, sd_U, "/")
        if (any(sd_U == 0)) {
          S[, sd_U == 0] <- NA_real_
        }
      } else {
        S <- sweep(cov_XU, 1, var_X[keep_feats], "/")
      }

      rownames(S) <- colnames(Xsub)[keep_feats]
      colnames(S) <- colnames(U)
    } else {
      S <- matrix(0, nrow = 0, ncol = ncol(U))
    }
  }

  if (nrow(S) < 1L) {
    return(tibble::tibble())
  }

  assoc_tbl <- tibble::as_tibble(S, rownames = "feature")

  ax_names <- colnames(U)
  top_list <- unique(unlist(lapply(seq_along(ax_names), function(j) {
    ord <- order(abs(S[, j]), decreasing = TRUE, na.last = NA)
    head(rownames(S)[ord], top_features)
  })))

  dplyr::filter(
    assoc_tbl,
    .data$feature %in% top_list
  )
}

#' @title Principal Components Analysis (PCoA) on a Distance Matrix
#'
#' @description Performs PCoA on a distance matrix (typically from
#' \code{compute_distance()}), optionally correcting for negative eigenvalues,
#' and returns coordinates, eigenvalues, variance explained, and feature-axis
#' associations.
#'
#' @param dist_obj a \code{dist} object (for example returned by
#'   \code{compute_distance()}). The normalized abundance matrix used to compute
#'   the distances is attached as attribute \code{"abundances"} (numeric
#'   matrix with samples in rows and features in columns). If missing, feature
#'   associations are skipped.
#' @param neg_correction character scalar. Method for adjusting negative
#'   eigenvalues (if any). One of \code{"none"}, \code{"lingoes"}, or
#'   \code{"cailliez"}. Default is \code{"none"}.
#' @param n_axes integer scalar. Number of PCoA axes to return in the sample
#'   scores. Must be > 0. Internally, \code{k = min(n_axes, n_samples - 1)}.
#' @param top_features integer scalar. Number of features to keep per axis when
#'   reporting associations. Features are selected by taking the union of the
#'   top \code{top_features} features (by absolute association) for each
#'   returned axis. Must be > 0.
#' @param feature_assoc character scalar. Type of feature-axis association to
#'   return. \code{"weighted_average"} returns weighted-average feature scores (centroid of
#'   sample scores weighted by feature abundance). \code{"correlation"} returns
#'   feature-axis correlations. \code{"regression"} returns regression slopes
#'   for axis scores on feature abundance. \code{"none"} skips feature
#'   associations.
#'
#' @return a list of class \code{"beta_pcoa"} with elements:
#' \itemize{
#'   \item \code{sample_coords}: tibble with \code{sample_id} and columns
#'     \code{PCoA1}, \code{PCoA2}, ... up to \code{n_axes} (or fewer if
#'     \code{n_samples - 1} is smaller).
#'   \item \code{eigenvalues}: numeric vector of eigenvalues from the PCoA.
#'   \item \code{var_explained}: one-row tibble with percent variance explained
#'     by the returned axes and \code{%Other}. percentages are computed from the
#'     sum of positive eigenvalues.
#'   \item \code{eigen_diagnostics}: one-row tibble with eigenvalue diagnostics:
#'     \code{sum_negative}, \code{sum_positive}, their ratio, the minimum
#'     eigenvalue, and the number and fraction of negative eigenvalues.
#'   \item \code{correction_infos}: one-row tibble describing any negative
#'     eigenvalue correction applied, including \code{corrected},
#'     \code{correction_method}, \code{correction_add}, and
#'     \code{correction_note}.
#'   \item \code{feature_associations}: tibble of feature-axis associations for
#'     the returned axes (empty if \code{"abundances"} is missing, cannot be
#'     aligned, or \code{feature_assoc = "none"}).
#' }
#'
#' @details
#' Negative eigenvalues indicate that the distances are not perfectly euclidean.
#' If \code{neg_correction} is \code{"lingoes"} or \code{"cailliez"}, a
#' correction is applied via \code{vegan::wcmdscale(add = ...)}.
#'
#' Feature associations are post-hoc summaries of how features relate to PCoA
#' axes. Weighted-average scores (\code{feature_assoc = "weighted_average"}) compute
#' \code{t(X) %*% U / colSums(X)}, where \code{X} is the abundance matrix and
#' \code{U} are the sample coordinates. Correlation and regression associations
#' are computed between feature abundances and axis scores and are not "true"
#' PCA loadings unless distances are Euclidean and derived compatibly.
#'
#' @examples
#' \donttest{
#' # compute a distance matrix with an attached abundance matrix
#' # build an example <phip_data> object from the package example dataset
#' ps <- phip_load_example_data()
#'
#' # small subset for speed: 5 peptides at time t1
#' keep_pep <- c("16627", "5243", "24799", "16196", "18003")
#' dat_cols <- dplyr::tbl_vars(ps$data_long)
#' tp_col <- "time"
#'
#' ps_small <- ps |>
#'   dplyr::filter(
#'     peptide_id %in% keep_pep,
#'     !!rlang::sym(tp_col) == "T1"
#'   ) |>
#'   dplyr::collect()
#'
#' # compute distances (needs either 'parallelDist' or 'vegan')
#' val_col <- "exist"
#'
#' d <- compute_distance(
#'   ps_small,
#'   value_col = val_col,
#'   distance = "jaccard",
#'   n_threads = 2L
#' )
#'
#' pcoa_res <- compute_pcoa(d, neg_correction = "none", n_axes = 3L)
#' pcoa_res$sample_coords
#' pcoa_res$var_explained
#' pcoa_res$feature_associations
#' }
#' @export
compute_pcoa <- function(dist_obj,
                         neg_correction = c("none", "lingoes", "cailliez"),
                         n_axes = 5L,
                         top_features = 30L,
                         feature_assoc = c(
                           "weighted_average", "correlation",
                           "regression", "none"
                         )) {
  # ----------------------------------------------------------------------------
  # 1) input validation
  # ----------------------------------------------------------------------------
  chk::chk_s3_class(dist_obj, "dist")
  neg_correction <- match.arg(neg_correction)
  chk::chk_count(n_axes)
  chk::chk_gt(n_axes, 0)
  chk::chk_count(top_features)
  chk::chk_gt(top_features, 0)
  feature_assoc <- match.arg(feature_assoc)

  n <- attr(dist_obj, "Size")
  if (is.null(n) || n < 2L) {
    .ph_abort("`dist_obj` must contain at least 2 samples.")
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
    bullets = if (identical(neg_correction, "none")) {
      NULL
    } else {
      paste("using", neg_correction, "correction")
    }
  )

  # cmdscale requires k in [1, n - 1]
  k_cmd <- min(n_axes, n - 1L)

  is_corrected <- !identical(neg_correction, "none")
  correction_add <- NA_real_
  correction_note <- "No negative eigenvalue correction applied."
  if (is_corrected) {
    pcoa_fit <- vegan::wcmdscale(dist_obj, eig = TRUE,
                                 k = k_cmd, add = neg_correction)
    if (!is.null(pcoa_fit$add)) {
      correction_add <- as.numeric(pcoa_fit$add)
    }
    correction_note <- sprintf(
      "Distances modified using %s correction prior to PCoA.",
      neg_correction
    )
    .ph_warn(correction_note)
  } else {
    pcoa_fit <- stats::cmdscale(dist_obj, eig = TRUE, k = k_cmd)
  }

  eig_vals <- pcoa_fit$eig

  eig_vals <- as.numeric(eig_vals)

  coords <- as.matrix(pcoa_fit$points)

  # enforce stable axis and sample naming
  colnames(coords) <- paste0("PCoA", seq_len(ncol(coords)))

  rownames(coords) <- labels

  # ----------------------------------------------------------------------------
  # 3) sample coordinates (first n_axes, or fewer)
  # ----------------------------------------------------------------------------
  .ph_log_info("extracting sample coordinates.")
  k_use <- min(k_cmd, ncol(coords))
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
  .ph_log_info("summarizing eigenvalues and variance explained.")
  pos_eig <- pmax(eig_vals, 0)
  sum_pos <- sum(pos_eig, na.rm = TRUE)
  neg_eig <- pmin(eig_vals, 0)
  sum_neg <- sum(abs(neg_eig), na.rm = TRUE)

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
    c(as.list(pct_axes), `%Other` = pct_other)
  )

  min_eig <- if (length(eig_vals) > 0L) {
    min(eig_vals, na.rm = TRUE)
  } else {
    NA_real_
  }
  n_neg <- sum(eig_vals < 0, na.rm = TRUE)
  frac_neg <- if (length(eig_vals) > 0L) {
    n_neg / length(eig_vals)
  } else {
    NA_real_
  }
  ratio_neg_pos <- if (sum_pos > 0) {
    sum_neg / sum_pos
  } else {
    NA_real_
  }

  eigen_diagnostics <- tibble::as_tibble_row(list(
    sum_negative = sum_neg,
    sum_positive = sum_pos,
    ratio_negative_positive = ratio_neg_pos,
    min_eigenvalue = min_eig,
    n_negative = n_neg,
    frac_negative = frac_neg
  ))

  correction_infos <- tibble::as_tibble_row(list(
    corrected = is_corrected,
    correction_method = neg_correction,
    correction_add = correction_add,
    correction_note = correction_note
  ))

  # ----------------------------------------------------------------------------
  # 5) feature associations (requires abundances attribute)
  # ----------------------------------------------------------------------------
  .ph_log_info(
    paste0("computing feature associations: ", feature_assoc)
  )
  feature_associations <- tibble::tibble()

  X <- attr(dist_obj, "abundances")
  if (is.null(X)) {
    .ph_warn(
      "no 'abundances' attribute found on `dist_obj`;
      skipping feature associations."
    )
  } else if (identical(feature_assoc, "none") || k_use < 1L) {
    feature_associations <- tibble::tibble()
  } else {
    X <- as.matrix(X)
    feature_associations <- .ph_feature_associations(
      coords = coords,
      X = X,
      feature_assoc = feature_assoc,
      top_features = top_features,
      max_axes = k_use,
      warn_missing_ids = paste(
        "row names missing in coordinates or 'abundances'; cannot align",
        "samples for feature associations."
      ),
      warn_insufficient_overlap = paste(
        "insufficient overlap between distance labels and abundance rows;",
        "skipping feature associations."
      )
    )
  }

  # ---------------------------------------------------------------------------
  # 6) return
  # ---------------------------------------------------------------------------
  result <- list(
    sample_coords = sample_coords,
    eigenvalues = eig_vals,
    var_explained = var_explained,
    eigen_diagnostics = eigen_diagnostics,
    correction_infos = correction_infos,
    feature_associations = feature_associations
  )
  class(result) <- "beta_pcoa"

  .ph_log_info("pcoa analysis complete.")
  result
}

#' @title Constrained Ordination (db-rda / cap) on Distance Matrix
#'
#' @description Performs distance-based redundancy analysis (constrained pcoa,
#'   a.k.a. cap) on a distance matrix using \pkg{vegan}::\code{capscale}, with
#'   optional negative eigenvalue correction. Returns constrained sample scores,
#'   eigenvalues, variance partitioning, and feature associations.
#'
#' @param dist_obj A \code{dist} object returned by \code{compute_distance()}.
#'   The normalized abundance matrix used to compute the distances is expected
#'   to be attached as an attribute \code{"abundances"} (samples in rows,
#'   features in columns).
#' @param ps A \code{phip_data} object or a table providing sample-level
#'   metadata. This table must contain \code{sample_id} and all variables
#'   referenced on the right-hand side of \code{formula}. Variable detection
#'   uses \code{all.names(terms(formula), functions = FALSE)}, so transformed
#'   terms like \code{log(age)} are supported as long as \code{age} exists.
#' @param formula An R formula specifying the constraints (independent
#'   variables) for the ordination, e.g. \code{~ sex + age}. Do not include a
#'   response on the left-hand side; the distance matrix is provided via
#'   \code{dist_obj}.
#' @param neg_correction One of \code{"none"}, \code{"lingoes"},
#'   \code{"cailliez"}. Method for negative eigenvalue correction. Default is
#'   \code{"none"}. Passed to the \code{add} argument of
#'   \code{vegan::capscale()}.
#' @param top_features Integer scalar. Number of top features to return in
#'   associations (selected per constrained axis by absolute association, then
#'   unioned). Default is 30.
#' @param permutations Integer scalar. Number of permutations for per-term
#'   permutation tests via \code{vegan::anova.cca(by = "term")}. Default is 999.
#' @param feature_assoc character scalar. Type of feature-axis association to
#'   return. \code{"weighted_average"} returns weighted-average feature scores (centroid of
#'   sample scores weighted by feature abundance). \code{"correlation"} returns
#'   feature-axis correlations. \code{"regression"} returns regression slopes
#'   for axis scores on feature abundance. \code{"none"} skips feature
#'   associations.
#'
#' @return A list of class \code{"beta_capscale"} with elements:
#' \item{sample_coords}{Tibble of sample scores on constrained axes
#'   (\code{CAP1}, \code{CAP2}, ...). Contains \code{sample_id} and
#'   coordinates.}
#'   \item{eigenvalues}{Numeric vector of eigenvalues of the constrained axes.}
#' \item{variance_partition}{Tibble with total inertia and inertia partitioned
#'   into constrained and unconstrained components, with their proportion of
#'   total.}
#' \item{feature_associations}{Tibble of top feature-axis associations for
#'   constrained axes (possibly empty if the \code{"abundances"} attribute is
#'   missing or cannot be aligned). To limit runtime/memory, associations are
#'   computed for at most 10 constrained axes.}
#'   \item{r2}{Numeric scalar. Unadjusted R-squared from
#'   \code{vegan::RsquareAdj()}.} \item{r2_adj}{Numeric scalar. Adjusted
#'   R-squared from \code{vegan::RsquareAdj()}.}
#' \item{perm_terms}{Tibble of per-term permutation tests from
#'   \code{vegan::anova.cca(by = "term")}.}
#'   \item{cap_model}{The full \code{vegan::capscale} model object.}
#'
#' @examples
#' \donttest{
#' ps <- phip_load_example_data()
#'
#' # small subset for speed: 5 peptides at time t1
#' keep_pep <- c("16627", "5243", "24799", "16196", "18003")
#' dat_cols <- dplyr::tbl_vars(ps$data_long)
#' tp_col <- if ("time" %in% dat_cols) {
#' "time"
#' } else {
#' "timepoint_factor"
#' }
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
#'   ps = ps_small,
#'   formula = stats::as.formula(paste0("~ ", rhs_var)),
#'   neg_correction = "none",
#'   top_features = 30L
#' )
#'
#' cap_res$variance_partition
#' cap_res$sample_coords
#' cap_res$feature_associations
#' }
#' @export
compute_capscale <- function(dist_obj,
                             ps,
                             formula,
                             neg_correction = c("none", "lingoes", "cailliez"),
                             top_features = 30L,
                             permutations = 999L,
                             feature_assoc = c(
                               "weighted_average", "correlation",
                               "regression", "none"
                             )) {
  # ----------------------------------------------------------------------------
  # 0) input validation
  # ----------------------------------------------------------------------------
  chk::chk_s3_class(dist_obj, "dist")
  chk::chk_true(inherits(formula, "formula"))
  if (length(formula) == 3L) {
    .ph_abort("Do not supply a response; use `~ ...`.")
  }
  neg_correction <- match.arg(neg_correction)
  chk::chk_count(top_features)
  chk::chk_gt(top_features, 0)
  chk::chk_count(permutations)
  chk::chk_gt(permutations, 0)
  feature_assoc <- match.arg(feature_assoc)

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
  dat <- if ("phip_data" %in% class(ps)) {
    ps$data_long
  } else {
    ps
  }
  if (is.null(dat)) {
    .ph_abort("`ps` is missing. cannot construct metadata.")
  }

  dat_cols <- dplyr::tbl_vars(dat)
  if (!"sample_id" %in% dat_cols) {
    .ph_abort("`ps$data_long` must contain a `sample_id` column.")
  }

  # required variables from formula rhs (supports transformations like log(age))
  vars_needed <- all.names(stats::terms(formula), functions = FALSE)
  vars_needed <- setdiff(vars_needed, ".") # defensive

  if (length(vars_needed) == 0L) {
    .ph_abort(
      "No constraints provided in formula (rhs is empty).
      Use compute_pcoa() for unconstrained ordination.",
    )
  }

  missing_vars <- setdiff(vars_needed, dat_cols)
  if (length(missing_vars) > 0L) {
    .ph_abort(
      paste0(
        "the following variables from the formula are missing in `ps`: ",
        paste(missing_vars, collapse = ", ")
      )
    )
  }

  .ph_log_info("building metadata from `ps$data_long`.")

  # sample-level metadata (one row per sample_id)
  meta_all <- dat |>
    dplyr::select(sample_id, dplyr::all_of(vars_needed)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect() |>
    as.data.frame()

  if (nrow(meta_all) == 0L) {
    .ph_abort(
      "constructed metadata has zero rows. check that `ps` is not empty."
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
        )
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
      the distance order."
    )

    # critical: propagate assumed labels to the dist object so downstream
    # indexing works
    if (!is.null(attr(d, "Size")) && attr(d, "Size") != n) {
      .ph_abort(
        paste0(
          "cannot align distance and metadata without labels: dist size is ",
          attr(d, "Size"),
          " but metadata has ", n, " samples."
        )
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
      paste0(
        "dropping ", dropped,
        " samples with missing values in constrained variables."
      )
    )
  }

  meta_df <- meta_sub[keep, , drop = FALSE]
  if (nrow(meta_df) == 0L) {
    .ph_abort(
      "all samples have missing values in constrained variables; cannot fit
      cap."
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
        samples."
      )
    }
  }

  # ----------------------------------------------------------------------------
  # 4) build capscale formula: d_resp ~ rhs
  # ----------------------------------------------------------------------------
  d_resp <- d
  cap_formula <- stats::update.formula(formula, d_resp ~ .)
  cap_env <- new.env(parent = environment(formula))
  cap_env$d_resp <- d_resp
  environment(cap_formula) <- cap_env

  add_arg <- if (identical(neg_correction, "none")) {
    FALSE
  } else {
    neg_correction
  }

  .ph_log_info(
    "fitting constrained ordination (cap/db-rda)",
    bullets = c(
      paste0("formula: ", paste(deparse(formula), collapse = " ")),
      if (!identical(add_arg, FALSE)) {
        paste0("neg_correction: ", add_arg)
      }
    )
  )

  # ----------------------------------------------------------------------------
  # 5) fit capscale model
  # ----------------------------------------------------------------------------
  cap_fit <- vegan::capscale(cap_formula, data = meta_df, add = add_arg)

  # ----------------------------------------------------------------------------
  # 6) sample scores on constrained axes
  # ----------------------------------------------------------------------------
  .ph_log_info("extracting constrained sample scores.")
  rank_constrained <- cap_fit$CCA$rank
  if (is.null(rank_constrained)) {
    rank_constrained <- 0L
  }

  if (rank_constrained > 0L) {
    site_scores <- vegan::scores(
      cap_fit,
      display = "sites",
      choices = seq_len(rank_constrained)
    )
    pts <- as.matrix(site_scores)
  } else {
    pts <- matrix(0,
      nrow = n, ncol = 0L,
      dimnames = list(labels, character(0L))
    )
  }

  if (ncol(pts) > 0L && is.null(colnames(pts))) {
    colnames(pts) <- paste0("CAP", seq_len(ncol(pts)))
  }
  if (!is.null(labels) && nrow(pts) == length(labels)) {
    rownames(pts) <- labels
  }

  sample_coords <- tibble::as_tibble(pts, rownames = "sample_id")

  # ----------------------------------------------------------------------------
  # 7) eigenvalues, variance partition, and inference
  # ----------------------------------------------------------------------------
  .ph_log_info("computing variance partitioning and permutation tests.")
  eig_constrained <- cap_fit$CCA$eig
  if (is.null(eig_constrained)) {
    eig_constrained <- numeric()
  }

  tot_inertia <- cap_fit$tot.chi
  cons_inertia <- cap_fit$CCA$tot.chi
  if (is.null(cons_inertia)) {
    cons_inertia <- sum(cap_fit$CCA$eig %||% 0)
  }

  uncon_inertia <- cap_fit$CA$tot.chi
  if (is.null(uncon_inertia)) {
    uncon_inertia <- sum(cap_fit$CA$eig %||% 0)
  }

  variance_partition <- tibble::tibble(
    component = c("Total", "Constrained", "Unconstrained"),
    inertia = c(tot_inertia, cons_inertia, uncon_inertia),
    proportion = c(
      1,
      if (is.finite(tot_inertia) && tot_inertia > 0) {
        cons_inertia / tot_inertia
      } else {
        NA_real_
      },
      if (is.finite(tot_inertia) && tot_inertia > 0) {
        uncon_inertia / tot_inertia
      } else {
        NA_real_
      }
    )
  )

  r2_vals <- vegan::RsquareAdj(cap_fit)
  r2 <- as.numeric(r2_vals$r.squared)
  r2_adj <- as.numeric(r2_vals$adj.r.squared)

  perm_terms_raw <- vegan::anova.cca(
    cap_fit,
    by = "term",
    permutations = permutations
  )
  perm_terms <- tibble::as_tibble(perm_terms_raw, rownames = "term")

  # ----------------------------------------------------------------------------
  # 8) feature associations on constrained axes (using X_sub)
  # ----------------------------------------------------------------------------
  .ph_log_info(paste0("computing feature associations: ", feature_assoc, "."))
  feature_associations <- tibble::tibble()

  if (!is.null(X_sub) && rank_constrained > 0L) {
    feature_associations <- .ph_feature_associations(
      coords = pts,
      X = X_sub,
      feature_assoc = feature_assoc,
      top_features = top_features,
      max_axes = min(rank_constrained, 10L),
      warn_missing_ids = paste(
        "row names missing in sample scores or abundance matrix; cannot align",
        "for feature associations."
      ),
      warn_insufficient_overlap = paste(
        "insufficient overlap between distance labels and abundance rows;",
        "skipping feature associations."
      )
    )
  }

  # ----------------------------------------------------------------------------
  # 9) return
  # ----------------------------------------------------------------------------
  result <- list(
    sample_coords = sample_coords,
    eigenvalues = as.numeric(eig_constrained),
    variance_partition = variance_partition,
    feature_associations = feature_associations,
    r2 = r2,
    r2_adj = r2_adj,
    perm_terms = perm_terms,
    cap_model = cap_fit
  )
  class(result) <- "beta_capscale"

  .ph_log_info("cap analysis complete.")

  result
}

#' @title PERMANOVA with Global and Post-hoc Tests on Beta Diversity
#' @description Performs PERMANOVA (adonis2) on a distance matrix for overall
#'   group/time effects, and optionally conducts post-hoc pairwise or contrast
#'   tests (e.g., between each pair of groups, etc.). Supports stratified
#'   permutations for repeated measures.
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
#'   subject. This is a simplification and does not implement a full
#'   repeated-measures permutation design (see Details).
#' @param permutations Number of permutations for significance testing
#'   (default 999).
#' @param p_adjust P-value adjustment method applied within each contrast scope.
#'   Use \code{"none"} for raw p-values. Passed to \code{stats::p.adjust()}.
#'
#' @return A tibble with columns:
#' \item{scope}{The scope of the test (e.g., \code{"global"},
#'   \code{"group_pairwise"}, \code{"time_pairwise"}).}
#' \item{contrast}{Description of the contrast (e.g., \code{"<global>"} for
#' overall test,
#'   \code{"A vs B"} for pairwise group or time comparisons).}
#' \item{term}{The term being tested. For global tests, this will be the factor
#' name (or interaction).
#'   For post-hoc tests, it may be \code{"group"} or \code{"time"} indicating
#'   which factor is being contrasted.}
#'   \item{F_stat}{F statistic of the PERMANOVA test (for global tests and some
#'   contrasts where applicable).} \item{R2}{R-squared (variance explained) for
#'   the term (global tests).} \item{p_value}{Permutation p-value for the test.}
#' \item{p_adjust}{Adjusted p-value (within scope); equals \code{p_value} when
#'   \code{p_adjust = "none"}.}
#'   \item{n_perm}{Number of permutations used.}
#'
#' @details
#' Global PERMANOVA uses a model with main effects of \code{group_col} and
#' \code{time_col} and their interaction (if both are provided and have >1
#' level).
#' Samples with missing values in the constrained variables (group/time and,
#' when used for stratification, subject) are dropped \emph{before} fitting the
#' model, and the distance matrix is subset to the remaining samples so that
#' distances and metadata are always aligned.
#'
#' Pairwise post-hoc tests are always performed for \code{group_col} and
#' \code{time_col} (when present and with >1 level), using \code{adonis2} with
#' appropriate subsetting and, where applicable, subject stratification.
#' Stratification via \code{strata = subject} is a simplified approach and does
#' not replace a full repeated-measures permutation design (e.g., via
#' \code{permute::how()}), especially for multi-level time with missing
#' timepoints.
#'
#' P-values are optionally adjusted within each contrast scope (e.g.,
#' \code{"group_pairwise"}, \code{"time_pairwise"}) using
#' \code{stats::p.adjust()}.
#'
#' @examples
#' \donttest{
#' ps <- phip_load_example_data()
#'
#' # small subset for speed
#' keep_pep <- c("16627", "5243", "24799", "16196", "18003")
#' dat_cols <- dplyr::tbl_vars(ps$data_long)
#' tp_col <- "time"
#'
#' ps_small <- ps |>
#'   dplyr::filter(
#'     peptide_id %in% keep_pep,
#'     !!rlang::sym(tp_col) == "T1"
#'   ) |>
#'   dplyr::collect()
#'
#' # compute distance matrix
#' val_col <- "exist"
#'
#' dist_bc <- compute_distance(
#'   ps_small,
#'   value_col = val_col,
#'   distance = "jaccard",
#'   n_threads = 2L
#' )
#'
#' permanova_res <- compute_permanova(
#'   dist_bc,
#'   ps        = ps_small,
#'   group_col = "group",
#'   time_col  = "time"
#' )
#'
#' permanova_res2 <- compute_permanova(
#'   dist_bc,
#'   ps        = ps_small,
#'   group_col = "group",
#'   time_col  = "time",
#'   p_adjust  = "BH"
#' )
#' }
#' @export
compute_permanova <- function(dist_obj,
                              ps,
                              group_col = NULL,
                              time_col = NULL,
                              subject_col = "subject_id",
                              permutations = 999,
                              p_adjust = "none") {
  # ----------------------------------------------------------------------------
  # input validation (chk)
  # ----------------------------------------------------------------------------
  chk::chk_s3_class(dist_obj, "dist")
  if (!is.null(group_col)) {
    chk::chk_string(group_col)
  }
  if (!is.null(time_col)) {
    chk::chk_string(time_col)
  }
  chk::chk_string(subject_col)
  chk::chk_count(permutations)
  chk::chk_gt(permutations, 0)
  chk::chk_string(p_adjust)
  if (!p_adjust %in% c("none", stats::p.adjust.methods)) {
    .ph_abort(
      paste0(
        "`p_adjust` must be one of: ",
        paste(c("none", stats::p.adjust.methods), collapse = ", "),
        "."
      )
    )
  }

  # if ps is <phip_data>, overwrite ps with ps$data_long;
  # otherwise treat ps as data_long
  if (inherits(ps, "phip_data")) {
    chk::chk_not_null(ps$data_long)
    ps <- ps$data_long
  }

  # ----------------------------------------------------------------------------
  # prepare result collector
  # ----------------------------------------------------------------------------
  results_list <- list()
  add_result <- function(scope, contrast, term = NA, p_value = NA,
                         F_stat = NA, R2 = NA) {
    results_list[[length(results_list) + 1]] <<- tibble::tibble(
      scope = scope,
      contrast = contrast,
      term = term,
      F_stat = F_stat,
      R2 = R2,
      p_value = p_value,
      p_adjust = p_value,
      n_perm = permutations
    )
  }

  # ----------------------------------------------------------------------------
  # start from distances + labels
  # ----------------------------------------------------------------------------
  .ph_log_info("preparing distance labels and metadata.")
  d_full <- dist_obj
  labels_full <- attr(d_full, "Labels")
  n_full <- attr(d_full, "Size")

  # ----------------------------------------------------------------------------
  # build metadata from ps and align to dist labels
  # ----------------------------------------------------------------------------
  dat <- ps
  if (is.null(dat)) {
    .ph_abort("`ps` is missing. Cannot construct metadata.")
  }

  dat_cols <- dplyr::tbl_vars(dat)
  if (!"sample_id" %in% dat_cols) {
    .ph_abort("`ps` must contain a `sample_id` column.")
  }

  has_group <- !is.null(group_col) && group_col %in% dat_cols
  has_time <- !is.null(time_col) && time_col %in% dat_cols
  has_subject <- !is.null(subject_col) && subject_col %in% dat_cols

  if (!is.null(group_col) && !has_group) {
    .ph_abort(
      paste0("Column `", group_col, "` not found in `ps`."),
    )
  }
  if (!is.null(time_col) && !has_time) {
    .ph_abort(
      paste0("Column `", time_col, "` not found in `ps`."),
    )
  }

  cols_needed <- c(
    "sample_id",
    if (has_group) {
      group_col
    } else {
      character(0L)
    },
    if (has_time) {
      time_col
    } else {
      character(0L)
    },
    if (has_subject) {
      subject_col
    } else {
      character(0L)
    }
  )
  cols_needed <- unique(cols_needed)

  # build metadata from ps
  .ph_log_info("building metadata from `ps`.")

  meta_all <- dat |>
    dplyr::select(dplyr::all_of(cols_needed)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect() |>
    as.data.frame()

  if (nrow(meta_all) == 0L) {
    .ph_abort(
      "constructed metadata has zero rows. check that `ps` is not empty."
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
        )
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
      the distance order."
    )
  }

  # coerce group/time to factor for safety
  if (has_group) {
    meta_sub[[group_col]] <- as.factor(meta_sub[[group_col]])
  }
  if (has_time) {
    meta_sub[[time_col]] <- as.factor(meta_sub[[time_col]])
  }

  # ----------------------------------------------------------------------------
  # drop samples with na in constrained variables (+ subject if used)
  # ----------------------------------------------------------------------------
  .ph_log_info("filtering samples with missing grouping variables.")
  vars_for_na <- c(
    if (has_group) {
      group_col
    } else {
      character(0L)
    },
    if (has_time) {
      time_col
    } else {
      character(0L)
    },
    if (has_time && has_subject) {
      subject_col
    } else {
      character(0L)
    }
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
      paste0(
        "dropping ", dropped,
        " samples with missing values in constrained/strata variables."
      )
    )
  }

  meta_df <- meta_sub[keep, , drop = FALSE]
  if (nrow(meta_df) == 0L) {
    .ph_abort(
      "all samples have missing values in constrained/strata variables;
      cannot run permanova."
    )
  }
  keep_labels <- rownames(meta_df)

  # subset distance matrix to complete-case samples
  .ph_log_info("subsetting distance matrix to complete cases.")
  mat_d_full <- as.matrix(d_full)
  mat_d_sub <- mat_d_full[keep_labels, keep_labels, drop = FALSE]
  d <- stats::as.dist(mat_d_sub)
  labels <- attr(d, "Labels")
  n <- attr(d, "Size")

  # update factor presence after na-drop
  has_group <- has_group && length(unique(meta_df[[group_col]])) > 1L
  has_time <- has_time && length(unique(meta_df[[time_col]])) > 1L

  # ----------------------------------------------------------------------------
  # global permanova
  # ----------------------------------------------------------------------------
  .ph_log_info("preparing global permanova model.")
  rhs_terms <- character(0L)
  if (has_group) {
    rhs_terms <- c(rhs_terms, group_col)
  }
  if (has_time) {
    rhs_terms <- c(rhs_terms, time_col)
  }
  if (has_group && has_time) {
    rhs_terms <- c(rhs_terms, paste(group_col, "*", time_col))
  }

  if (length(rhs_terms) == 0L) {
    .ph_log_info(
      "global permanova skipped (insufficient number of factor levels)."
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
      bullets = c(
        paste("model:", formula_str),
        if (!is.null(strata_var)) {
          "permutations stratified by subject"
        }
      )
    )

    adonis_res <- try(
      vegan::adonis2(
        form,
        data = df,
        permutations = permutations,
        by = "terms",
        strata = strata_var,
        parallel = 1
      ),
      silent = TRUE
    )

    if (inherits(adonis_res, "try-error")) {
      .ph_warn("global permanova failed; no global results.")
    } else {
      res_df <- as.data.frame(adonis_res)
      res_df$term <- rownames(res_df)
      # check terms: group, time, and interaction
      for (term in c(
        if (!is.null(group_col)) {
          group_col
        } else {
          NULL
        },
        if (!is.null(time_col)) {
          time_col
        } else {
          NULL
        },
        if (!is.null(group_col) && !is.null(time_col)) {
          paste(group_col, ":", time_col, sep = "")
        } else {
          NULL
        }
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

  # ----------------------------------------------------------------------------
  # post-hoc contrasts (pairwise)
  # ----------------------------------------------------------------------------
  .ph_log_info("running pairwise permanova contrasts.")

  # helper: run two-level adonis on subset idx of current meta_df/d
  run_two_level_adonis <- function(idx, fac, covar = NULL, strata = NULL,
                                   scope_label, contrast_label, term_label) {
    if (length(unique(fac)) < 2L || min(table(fac)) < 2L) {
      .ph_log_info(
        paste(
          "skipping test", contrast_label,
          "- not enough samples in one or both groups."
        )
      )
      return(NULL)
    }
    df_sub <- meta_df[idx, , drop = FALSE]
    labels_sub <- labels[idx]
    rownames(df_sub) <- labels_sub

    form_rhs <- "fac"
    if (!is.null(covar)) {
      form_rhs <- paste(form_rhs, "+", covar)
    }

    fml <- stats::as.formula(paste("d_sub ~", form_rhs))
    d_mat <- as.matrix(d)
    d_sub <- stats::as.dist(d_mat[idx, idx, drop = FALSE])

    ad <- vegan::adonis2(
      fml,
      data = cbind(df_sub, fac = fac),
      permutations = permutations,
      strata = strata,
      by = if (!is.null(covar)) {
        "margin"
      } else {
        NULL
      }
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
          multi <- any(tapply(
            meta_df[sel, group_col], sub_sel,
            function(x) length(unique(x)) > 1L
          ))
          if (multi) {
            strata_use <- sub_sel
          }
        }

        fac_pair <- factor(meta_df[sel, group_col], levels = p)
        covar_term <- if (has_time) {
          time_col
        } else {
          NULL
        }

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

        fac_pair <- factor(meta_df[sel, time_col], levels = p)
        covar_term <- if (has_group) {
          group_col
        } else {
          NULL
        }

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

  # ----------------------------------------------------------------------------
  # combine and return
  # ----------------------------------------------------------------------------
  result_df <- if (length(results_list) > 0L) {
    dplyr::bind_rows(results_list)
  } else {
    tibble::tibble()
  }

  if (nrow(result_df) > 0L) {
    if (!identical(p_adjust, "none")) {
      result_df <- result_df |>
        dplyr::group_by(.data$scope) |>
        dplyr::mutate(
          p_adjust = stats::p.adjust(.data$p_value, method = .env$p_adjust)
        ) |>
        dplyr::ungroup()
    }
  }

  result_df
}

#' @title Test Homogeneity of Dispersion (Beta Dispersion)
#' @description Computes distances of samples to group centroids (using
#'   \code{vegan::betadisper}) and tests for differences in dispersion among
#'   groups or time levels. Pairwise post-hoc tests are always performed when
#'   \code{group_col} or \code{time_col} has >1 level.
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
#' @param p_adjust P-value adjustment method applied within each contrast scope.
#'   Use \code{"none"} for raw p-values. Passed to \code{stats::p.adjust()}.
#'
#' @return A list of class \code{"beta_dispersion"} with:
#' \item{distances}{Tibble of per-sample distances to centroid. Columns:
#'   \code{sample_id}, \code{distance}, \code{level} (group/time level for
#'   a given scope), \code{scope} (e.g. \code{"group"}, \code{"time"},
#'   \code{"group:time"}), \code{contrast} (e.g. \code{"<global>"},
#'   \code{"A vs B"}).}
#' \item{tests}{Tibble of dispersion test results. Columns: \code{scope},
#'   \code{contrast}, \code{term = "dispersion"}, \code{p_value},
#'   \code{p_adjust} (equals \code{p_value} when \code{p_adjust = "none"}),
#'   \code{n_perm}.}
#'
#' @examples
#' \donttest{
#' ps <- phip_load_example_data()
#'
#' # small subset for speed
#' keep_pep <- c("16627", "5243", "24799", "16196", "18003")
#' dat_cols <- dplyr::tbl_vars(ps$data_long)
#' tp_col <- "time"
#'
#' ps_small <- ps |>
#'   dplyr::filter(
#'     peptide_id %in% keep_pep,
#'     !!rlang::sym(tp_col) == "T1"
#'   ) |>
#'   dplyr::collect()
#'
#' # compute distance matrix
#' val_col <- "exist"
#'
#' dist_bc <- compute_distance(
#'   ps_small,
#'   value_col = val_col,
#'   distance = "jaccard",
#'   n_threads = 2L
#' )
#'
#' dispersion_res <- compute_dispersion(
#'   dist_bc,
#'   ps        = ps_small,
#'   group_col = "group",
#'   time_col  = "time",
#'   p_adjust  = "BH"
#' )
#' dispersion_res$tests
#' head(dispersion_res$distances)
#' }
#' @export
compute_dispersion <- function(dist_obj,
                               ps,
                               group_col = NULL,
                               time_col = NULL,
                               subject_col = "subject_id",
                               permutations = 999,
                               p_adjust = "none") {
  # ----------------------------------------------------------------------------
  # input validation (chk)
  # ----------------------------------------------------------------------------
  chk::chk_s3_class(dist_obj, "dist")
  if (!is.null(group_col)) {
    chk::chk_string(group_col)
  }
  if (!is.null(time_col)) {
    chk::chk_string(time_col)
  }
  chk::chk_string(subject_col)
  chk::chk_count(permutations)
  chk::chk_gt(permutations, 0)
  chk::chk_string(p_adjust)
  if (!p_adjust %in% c("none", stats::p.adjust.methods)) {
    .ph_abort(
      paste0(
        "`p_adjust` must be one of: ",
        paste(c("none", stats::p.adjust.methods), collapse = ", "),
        "."
      )
    )
  }

  # if ps is <phip_data>, overwrite ps with ps$data_long;
  # otherwise treat ps as data_long
  if (inherits(ps, "phip_data")) {
    chk::chk_not_null(ps$data_long)
    ps <- ps$data_long
  }

  # ----------------------------------------------------------------------------
  # prepare metadata from ps and align to dist labels
  # ----------------------------------------------------------------------------
  .ph_log_info("preparing distance labels and metadata.")
  d_full <- dist_obj
  labels_full <- attr(d_full, "Labels")

  dat <- ps
  if (is.null(dat)) {
    .ph_abort("`ps` is missing. Cannot construct metadata.")
  }

  dat_cols <- dplyr::tbl_vars(dat)
  if (!"sample_id" %in% dat_cols) {
    .ph_abort("`ps` must contain a `sample_id` column.")
  }

  has_group <- !is.null(group_col) && group_col %in% dat_cols
  has_time <- !is.null(time_col) && time_col %in% dat_cols

  if (!is.null(group_col) && !has_group) {
    .ph_abort(
      paste0("Column `", group_col, "` not found in `ps`.")
    )
  }
  if (!is.null(time_col) && !has_time) {
    .ph_abort(
      paste0("Column `", time_col, "` not found in `ps`.")
    )
  }

  cols_needed <- c(
    "sample_id",
    if (has_group) {
      group_col
    } else {
      character(0L)
    },
    if (has_time) {
      time_col
    } else {
      character(0L)
    }
  )
  cols_needed <- unique(cols_needed)

  # build metadata from ps
  .ph_log_info("building metadata from `ps`.")

  meta_all <- dat |>
    dplyr::select(dplyr::all_of(cols_needed)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect() |>
    as.data.frame()

  if (nrow(meta_all) == 0L) {
    .ph_abort(
      "constructed metadata has zero rows. check that `ps` is not empty."
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
        )
      )
    }
    meta_sub <- meta_all[idx_align, , drop = FALSE]
    rownames(meta_sub) <- labels_full
  } else {
    meta_sub <- meta_all
    labels_full <- rownames(meta_sub)
    .ph_warn(
      "no labels found in `dist_obj`; assuming metadata row order matches the
      distance order."
    )
  }

  # coerce grouping vars to factor
  if (has_group) {
    meta_sub[[group_col]] <- as.factor(meta_sub[[group_col]])
  }
  if (has_time) {
    meta_sub[[time_col]] <- as.factor(meta_sub[[time_col]])
  }

  # ----------------------------------------------------------------------------
  # drop samples with na in group/time and subset distance matrix
  # ----------------------------------------------------------------------------
  .ph_log_info("filtering samples with missing grouping variables.")
  vars_for_na <- c(
    if (has_group) {
      group_col
    } else {
      character(0L)
    },
    if (has_time) {
      time_col
    } else {
      character(0L)
    }
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
      paste0(
        "dropping ", dropped,
        " samples with missing values in dispersion grouping variables."
      )
    )
  }

  meta_df <- meta_sub[keep, , drop = FALSE]
  if (nrow(meta_df) == 0L) {
    .ph_abort(
      "all samples have missing values in grouping variables; cannot run
      dispersion tests."
    )
  }
  keep_labels <- rownames(meta_df)

  d_mat_full <- as.matrix(d_full)
  d_mat_sub <- d_mat_full[keep_labels, keep_labels, drop = FALSE]
  d <- stats::as.dist(d_mat_sub)
  labels <- attr(d, "Labels")

  # re-evaluate factor presence after na drop
  has_group <- has_group && length(unique(meta_df[[group_col]])) > 1L
  has_time <- has_time && length(unique(meta_df[[time_col]])) > 1L

  # ----------------------------------------------------------------------------
  # prepare collectors for distances and tests
  # ----------------------------------------------------------------------------
  distances_list <- list()
  tests_list <- list()

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
      p_adjust = p_val,
      n_perm   = permutations
    )
  }

  # ----------------------------------------------------------------------------
  # global dispersion tests (group, time, group:time)
  # ----------------------------------------------------------------------------
  .ph_log_info("computing global dispersion tests.")
  if (has_group) {
    fac <- factor(meta_df[[group_col]])
    bd <- try(vegan::betadisper(d, fac), silent = TRUE)
    if (inherits(bd, "try-error")) {
      .ph_warn(
        "betadisper failed for group factor; skipping group dispersion test.",
      )
    } else {
      add_distance_rows(
        sample_ids = names(bd$distances),
        dists = as.numeric(bd$distances),
        level_vals = as.character(bd$group),
        scope_val = "group",
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
        )
      } else {
        add_distance_rows(
          sample_ids = names(bd$distances),
          dists = as.numeric(bd$distances),
          level_vals = as.character(bd$group),
          scope_val = "time",
          contrast_val = "<global>"
        )
        pt <- vegan::permutest(bd, permutations = permutations, parallel = 1)
        pval <- tryCatch(pt$tab[1, "Pr(>F)"], error = function(e) NA_real_)
        add_test_row("time", "<global>", pval)
      }
    } else {
      .ph_warn(
        "`time_col` is numeric; continuous dispersion by time not supported.
        Skipping time dispersion test."
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
          dists = as.numeric(bd$distances),
          level_vals = as.character(bd$group),
          scope_val = "group:time",
          contrast_val = "<global>"
        )
        pt <- vegan::permutest(bd, permutations = permutations, parallel = 1)
        pval <- tryCatch(pt$tab[1, "Pr(>F)"], error = function(e) NA_real_)
        add_test_row("group:time", "<global>", pval)
      }
    }
  }

  # ----------------------------------------------------------------------------
  # helper for subset dispersion tests
  # ----------------------------------------------------------------------------
  run_disp_test <- function(idx, fac_vec, scope_lab, contrast_lab) {
    if (length(unique(fac_vec)) < 2L || min(table(fac_vec)) < 2L) {
      .ph_log_info(
        paste(
          "skipping dispersion test for", contrast_lab,
          "- not enough data."
        ),
      )
      return(NULL)
    }
    d_mat <- as.matrix(d)
    d_sub <- stats::as.dist(d_mat[idx, idx, drop = FALSE])
    bd <- try(vegan::betadisper(d_sub, fac_vec), silent = TRUE)
    if (inherits(bd, "try-error")) {
      return(NULL)
    }

    add_distance_rows(
      sample_ids = names(bd$distances),
      dists = as.numeric(bd$distances),
      level_vals = as.character(bd$group),
      scope_val = scope_lab,
      contrast_val = contrast_lab
    )
    pt <- vegan::permutest(bd, permutations = permutations, parallel = 1)
    pval <- tryCatch(pt$tab[1, "Pr(>F)"], error = function(e) NA_real_)
    add_test_row(scope_lab, contrast_lab, pval)
  }

  # ----------------------------------------------------------------------------
  # post-hoc contrasts (pairwise)
  # ----------------------------------------------------------------------------
  .ph_log_info("running pairwise dispersion contrasts.")
  # groups
  if (has_group) {
    grps <- na.omit(unique(meta_df[[group_col]]))
    if (length(grps) > 1L) {
      for (pair in utils::combn(grps, 2, simplify = FALSE)) {
        sel <- which(meta_df[[group_col]] %in% pair)
        fac <- factor(meta_df[sel, group_col], levels = pair)
        run_disp_test(sel,
          fac_vec = fac, scope_lab = "group",
          contrast_lab = paste(pair, collapse = " vs ")
        )
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
        run_disp_test(sel,
          fac_vec = fac, scope_lab = "time",
          contrast_lab = paste(pair, collapse = " vs ")
        )
      }
    }
  }

  # ----------------------------------------------------------------------------
  # combine and return
  # ----------------------------------------------------------------------------
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

  if (nrow(tests_tbl) > 0L) {
    if (!identical(p_adjust, "none")) {
      tests_tbl <- tests_tbl |>
        dplyr::group_by(.data$scope) |>
        dplyr::mutate(
          p_adjust = stats::p.adjust(.data$p_value, method = .env$p_adjust)
        ) |>
        dplyr::ungroup()
    }
  }

  result <- list(
    distances = distances_tbl,
    tests     = tests_tbl
  )
  class(result) <- "beta_dispersion"

  result
}

#' @title Compute t-SNE Embeddings for Sample Distances
#'
#' @description Performs t-distributed stochastic neighbor embedding (t-SNE)
#' on a sample distance matrix to create low-dimensional embeddings for
#' visualization. Returns t-SNE coordinates with optional sample metadata
#' attached.
#'
#' @param ps A \code{phip_data} object or a table providing sample-level
#'   metadata. This table must contain \code{sample_id} and any columns
#'   specified in \code{meta_cols}.
#' @param dist_obj A sample distance object. Either:
#'   \itemize{
#'     \item A \code{dist} object (e.g., from \code{compute_distance()}).
#'     \item A numeric, symmetric matrix with row/column names corresponding
#'       to \code{sample_id}s.
#'   }
#' @param dims Integer scalar. Number of t-SNE dimensions to compute (2 or 3).
#'   Default is 3 to enable both 2D and 3D visualizations.
#' @param perplexity Numeric scalar. Perplexity parameter for t-SNE. Must be
#'   smaller than the number of samples. If too large, it is automatically
#'   reduced with a warning.
#' @param theta Numeric scalar. Speed/accuracy tradeoff for Barnes-Hut
#'   approximation. Default is 0.5.
#' @param max_iter Integer scalar. Maximum number of t-SNE iterations.
#'   Default is 1000.
#' @param meta_cols Character vector. Column names from \code{ps} to attach
#'   as metadata. If \code{NULL}, uses \code{ps$meta$extra_cols} when
#'   available.
#' @param seed Integer scalar. Random seed for reproducibility. If provided,
#'   the seed is temporarily set and restored after computation.
#' @param check_duplicates Logical scalar. Whether to check for duplicate
#'   points. Default is \code{FALSE} for distance input.
#' @param ... Additional arguments passed to \code{Rtsne::Rtsne()}.
#'
#' @return A tibble of class \code{c("phip_tsne", "tbl_df", "tbl",
#' "data.frame")}
#' with columns:
#' \itemize{
#'   \item \code{sample_id}: Sample identifier from distance labels.
#'   \item \code{tSNE1}, \code{tSNE2}: First two t-SNE dimensions.
#'   \item \code{tSNE3}: Third dimension if \code{dims >= 3}, otherwise
#'   \code{NA}.
#'   \item Additional metadata columns as specified in \code{meta_cols}.
#' }
#'
#' Attributes:
#' \itemize{
#'   \item \code{"distance"}: The original \code{dist_obj}.
#'   \item \code{"tsne_params"}: List of t-SNE parameters and function call.
#'   \item \code{"meta_cols"}: Character vector of metadata columns used.
#' }
#'
#' @details
#' This function runs t-SNE in distance mode (\code{is_distance = TRUE}),
#' using the supplied distance matrix directly. The distance computation
#' should be performed separately using \code{compute_distance()}.
#'
#' t-SNE is a visualization method: it preserves local neighborhoods rather
#' than global distances, axes are not directly interpretable, and embeddings
#' can change with different seeds or perplexity settings.
#'
#' Samples with missing values in metadata columns are retained in the t-SNE
#' result but will have \code{NA} values for those metadata columns.
#'
#' @examples
#' \donttest{
#' # Build example phip_data object
#' ps <- phip_load_example_data()
#'
#' # Small subset for speed
#' keep_pep <- c("16627", "5243", "24799", "16196", "18003")
#' dat_cols <- dplyr::tbl_vars(ps$data_long)
#' tp_col <- "time"
#'
#' ps_small <- ps |>
#'   dplyr::filter(
#'     peptide_id %in% keep_pep,
#'     !!rlang::sym(tp_col) == "T1"
#'   ) |>
#'   dplyr::collect()
#'
#' # Compute distance matrix
#' val_col <- if ("fold_change" %in% dplyr::tbl_vars(ps_small)) {
#'   "fold_change"
#' } else {
#'   "exist"
#' }
#'
#' d <- compute_distance(
#'   ps_small,
#'   value_col = val_col,
#'   method_normalization = "hellinger",
#'   distance = "bray",
#'   n_threads = 2L
#' )
#'
#' # Compute t-SNE embeddings
#' tsne_res <- compute_tsne(
#'   ps = ps_small,
#'   dist_obj = d,
#'   dims = 3L,
#'   perplexity = 15,
#'   meta_cols = c("subject_id", "timepoint"),
#'   seed = 42
#' )
#'
#' # View results
#' head(tsne_res)
#' }
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
  # ----------------------------------------------------------------------------
  # input validation (chk)
  # ----------------------------------------------------------------------------
  if (missing(dist_obj)) {
    .ph_abort("Argument `dist_obj` is required. Run `compute_distance()`
              first.")
  }

  if (!inherits(dist_obj, "dist") && !is.matrix(dist_obj)) {
    .ph_abort(
      "Argument `dist_obj` must be either a <dist> object or a numeric matrix."
    )
  }

  chk::chk_count(dims)
  if (!dims %in% c(2L, 3L)) {
    .ph_abort("Argument `dims` must be either 2 or 3.")
  }

  chk::chk_number(perplexity)
  chk::chk_gt(perplexity, 0)

  chk::chk_number(theta)
  chk::chk_gte(theta, 0)

  chk::chk_count(max_iter)
  chk::chk_gt(max_iter, 0)

  if (!is.null(meta_cols)) {
    chk::chk_character(meta_cols)
  }
  if (!is.null(seed)) {
    chk::chk_count(seed)
  }
  chk::chk_flag(check_duplicates)

  if (!rlang::is_installed("Rtsne")) {
    .ph_abort(
      "Package 'Rtsne' is required to compute t-SNE. Please install it."
    )
  }

  # if ps is <phip_data>, validate it; otherwise treat as data table
  if (inherits(ps, "phip_data")) {
    chk::chk_not_null(ps$data_long)
  }

  # ----------------------------------------------------------------------------
  # extract labels and size from distance object
  # ----------------------------------------------------------------------------
  if (inherits(dist_obj, "dist")) {
    n_samples <- attr(dist_obj, "Size")
    labels <- attr(dist_obj, "Labels")

    if (is.null(labels)) {
      labels <- as.character(seq_len(n_samples))
      .ph_warn(
        "Distance object has no labels; using integer indices as `sample_id`."
      )
    }
  } else {
    # matrix input validation
    if (!is.numeric(dist_obj)) {
      .ph_abort("If `dist_obj` is a matrix, it must be numeric.")
    }
    if (nrow(dist_obj) != ncol(dist_obj)) {
      .ph_abort("Distance matrix must be square (same number of rows and
                columns).")
    }
    n_samples <- nrow(dist_obj)
    labels <- rownames(dist_obj)
    if (is.null(labels)) {
      labels <- as.character(seq_len(n_samples))
      .ph_warn(
        "Distance matrix has no rownames; using integer indices as `sample_id`."
      )
    }
  }

  if (n_samples < 3L) {
    .ph_abort(
      "t-SNE requires at least 3 samples."
    )
  }

  if (length(labels) != n_samples) {
    .ph_abort(
      "Length of distance labels does not match the number of samples."
    )
  }

  if (anyDuplicated(labels)) {
    .ph_abort(
      "Distance labels (sample IDs) must be unique."
    )
  }

  # ----------------------------------------------------------------------------
  # perplexity validation and adjustment
  # ----------------------------------------------------------------------------
  max_perplexity <- floor((n_samples - 1L) / 3L)
  if (perplexity >= n_samples) {
    .ph_abort(
      paste0(
        "Perplexity (", perplexity, ") must be smaller than the number of
        samples (", n_samples, ")."
      )
    )
  }
  if (perplexity > max_perplexity) {
    .ph_warn(
      paste0(
        "Perplexity (", perplexity, ") is high for n = ", n_samples,
        "; reducing to ", max_perplexity, "."
      )
    )
    perplexity <- max_perplexity
  }

  # ----------------------------------------------------------------------------
  # prepare metadata from ps
  # ----------------------------------------------------------------------------
  dat <- if (inherits(ps, "phip_data")) {
    ps$data_long
  } else {
    ps
  }

  if (!is.null(dat)) {
    dat_cols <- dplyr::tbl_vars(dat)
    if (!"sample_id" %in% dat_cols) {
      .ph_abort("`ps` must contain a `sample_id` column.")
    }
  } else {
    dat_cols <- character()
    .ph_warn(
      "`ps` is NULL or has no data; no metadata columns can be attached."
    )
  }

  if (is.null(meta_cols)) {
    # try to use meta$extra_cols if ps is phip_data
    if (inherits(ps, "phip_data")) {
      extra <- ps$meta$extra_cols %||% character()
      meta_cols <- intersect(extra, dat_cols)
    } else {
      meta_cols <- character()
    }
  } else {
    # user-specified meta_cols: keep only those present in data
    missing_cols <- setdiff(meta_cols, dat_cols)
    if (length(missing_cols) > 0L) {
      .ph_warn(
        paste0(
          "The following `meta_cols` are not present in `ps` and ",
          "will be ignored: ",
          paste(missing_cols, collapse = ", ")
        )
      )
      meta_cols <- intersect(meta_cols, dat_cols)
    }
  }

  # ----------------------------------------------------------------------------
  # run t-SNE computation
  # ----------------------------------------------------------------------------
  .ph_log_info(
    paste0(
      "Running t-SNE with dims = ", dims,
      ", perplexity = ", perplexity,
      " on ", n_samples, " samples (distance input)."
    )
  )

  run_tsne <- function() {
    if (inherits(dist_obj, "dist")) {
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
  }

  tsne_fit <- if (!is.null(seed)) {
    withr::with_seed(seed, run_tsne())
  } else {
    run_tsne()
  }

  coords <- tsne_fit$Y

  if (!is.matrix(coords) || nrow(coords) != n_samples || ncol(coords) != dims) {
    .ph_abort(
      "Unexpected output from Rtsne::Rtsne(): coordinate matrix has wrong
      shape."
    )
  }

  # ----------------------------------------------------------------------------
  # build base result tibble
  # ----------------------------------------------------------------------------
  tsne_df <- tibble::tibble(
    sample_id = labels,
    tSNE1     = coords[, 1],
    tSNE2     = if (dims >= 2L) {
      coords[, 2L]
    } else {
      NA_real_
    },
    tSNE3     = if (dims >= 3L) {
      coords[, 3L]
    } else {
      NA_real_
    }
  )

  # ----------------------------------------------------------------------------
  # attach metadata (if requested and available)
  # ----------------------------------------------------------------------------
  if (!is.null(dat) && length(meta_cols) > 0L) {
    .ph_log_info(
      paste0(
        "Attaching metadata columns to t-SNE result: ",
        paste(meta_cols, collapse = ", ")
      )
    )

    meta_raw <- dat |>
      dplyr::select(sample_id, dplyr::all_of(meta_cols)) |>
      dplyr::collect()

    meta_checks <- meta_raw |>
      dplyr::group_by(sample_id) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(meta_cols),
          ~ dplyr::n_distinct(.x[!is.na(.x)]),
          .names = "n_{.col}"
        ),
        .groups = "drop"
      )

    conflict_rows <- meta_checks |>
      dplyr::filter(dplyr::if_any(dplyr::starts_with("n_"), ~ .x > 1L))

    if (nrow(conflict_rows) > 0L) {
      conflict_cols <- names(meta_checks)[-1]
      conflict_cols <- conflict_cols[
        vapply(meta_checks[conflict_cols], function(x) any(x > 1L), logical(1))
      ]
      conflict_cols <- sub("^n_", "", conflict_cols)
      sample_preview <- paste(head(conflict_rows$sample_id, 5L),
        collapse = ", "
      )

      .ph_abort(
        paste0(
          "inconsistent metadata values per `sample_id` for: ",
          paste(conflict_cols, collapse = ", "),
          ". examples: ", sample_preview,
          ". resolve duplicates before running compute_tsne()."
        )
      )
    }

    meta_tbl <- meta_raw |>
      dplyr::distinct(sample_id, .keep_all = TRUE)

    tsne_df <- tsne_df |>
      dplyr::left_join(meta_tbl, by = "sample_id")
  }

  .ph_log_info("t-SNE embedding computation finished.")

  # ----------------------------------------------------------------------------
  # attach attributes and class
  # ----------------------------------------------------------------------------
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
