#' @return A \code{dist} object containing pairwise distances between samples.
#'   The normalized abundance matrix used to compute the distances is attached
#'   as an attribute \code{"abundances"}.
compute_distance <- function(ps,
                             value_col = NULL,
                             method_normalization = c("auto", "relative",
                                                      "hellinger", "log", "none"),
                             distance = "bray",
                             n_threads = 1L,
                             drop_all_zero_features = TRUE) {
  method_normalization <- match.arg(method_normalization)

  # ---------------------------------------------------------------------------
  # 1) Basic checks and column selection
  # ---------------------------------------------------------------------------
  if (!inherits(ps, "phip_data")) {
    .ph_abort("Input `ps` must be a <phip_data> object.",
              step = "compute_distance")
  }

  dat <- ps$data_long
  if (is.null(dat)) {
    .ph_abort("`ps$data_long` is missing. Cannot construct abundance matrix.",
              step = "compute_distance")
  }

  dat_cols <- dplyr::tbl_vars(dat)

  # Decide which abundance column to use -------------------------------------
  if (is.null(value_col)) {
    candidates <- c("counts_hit", "counts_input", "fold_change", "counts")
    hit <- candidates[candidates %in% dat_cols]
    if (length(hit) == 0L) {
      .ph_abort(
        paste0(
          "Could not infer an abundance column in `ps$data_long`. ",
          "Tried: ", paste(candidates, collapse = ", "),
          ". Please specify `value_col` explicitly."
        ),
        step = "compute_distance"
      )
    }
    value_col <- hit[1L]
    .ph_log_info(
      paste0("Auto-detected `value_col = \"", value_col, "\"` from `ps$data_long`."),
      step = "compute_distance"
    )
  }

  if (!value_col %in% dat_cols) {
    .ph_abort(
      paste0("Column `", value_col, "` not found in `ps$data_long`."),
      step = "compute_distance"
    )
  }

  required_cols <- c("sample_id", "peptide_id", value_col)
  missing_cols  <- setdiff(required_cols, dat_cols)
  if (length(missing_cols) > 0L) {
    .ph_abort(
      paste0("Missing required column(s) in `ps$data_long`: ",
             paste(missing_cols, collapse = ", ")),
      step = "compute_distance"
    )
  }

  .ph_log_info(
    paste0("Building abundance matrix from `ps$data_long` using `", value_col, "`."),
    step = "compute_distance"
  )

  value_sym <- rlang::sym(value_col)

  # ---------------------------------------------------------------------------
  # 2) Collect only needed columns and pivot in R (fast, avoids slow DB pivot)
  # ---------------------------------------------------------------------------
  .ph_log_info("Collecting long table (sample_id, peptide_id, value).",
               step = "compute_distance")

  dat_small <- dat |>
    dplyr::select(sample_id, peptide_id, !!value_sym) |>
    dplyr::collect()

  # Replace NAs with 0 in abundance column
  dat_small[[value_col]][is.na(dat_small[[value_col]])] <- 0

  .ph_log_info("Pivoting to wide abundance matrix in R.",
               step = "compute_distance")

  wide_df <- dat_small |>
    tidyr::pivot_wider(
      id_cols     = sample_id,
      names_from  = peptide_id,
      values_from = !!value_sym,
      values_fill = 0
    )

  if (!"sample_id" %in% names(wide_df)) {
    .ph_abort("Failed to construct wide abundance table (no `sample_id`).",
              step = "compute_distance")
  }

  mat <- wide_df |>
    tibble::column_to_rownames("sample_id") |>
    as.matrix()

  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    .ph_abort(
      "Abundance matrix is empty after reshaping. ",
      "Check filters and `value_col`.",
      step = "compute_distance"
    )
  }

  # Optional: drop all-zero features (they carry no information) ------------
  if (isTRUE(drop_all_zero_features)) {
    cs <- colSums(mat, na.rm = TRUE)
    keep <- cs != 0
    if (!all(keep)) {
      dropped <- sum(!keep)
      mat <- mat[, keep, drop = FALSE]
      .ph_log_info(
        paste0(
          "Dropped ", dropped,
          " all-zero features before distance computation."
        ),
        step = "compute_distance"
      )
    }
  }

  .ph_log_info(
    paste0("Abundance matrix has ", nrow(mat), " samples and ",
           ncol(mat), " features after preprocessing."),
    step = "compute_distance"
  )

  # ---------------------------------------------------------------------------
  # 3) Normalization
  # ---------------------------------------------------------------------------
  if (identical(method_normalization, "auto")) {
    vals <- mat[!is.na(mat)]
    is_binary_data <- length(vals) > 0L && all(vals == 0 | vals == 1)
    method_normalization <- if (is_binary_data) "none" else "relative"
    .ph_log_info(
      paste("Auto normalization selected -> using", method_normalization),
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
  # 4) Distance computation
  # ---------------------------------------------------------------------------
  dist_method <- tolower(distance)
  .ph_log_info(
    paste("Computing distance:", toupper(dist_method)),
    step = "compute_distance"
  )

  pd_supported <- c(
    "euclidean", "minkowski", "manhattan", "canberra",
    "binary", "maximum", "cosine", "chebyshev"
  )

  dist_obj <- NULL

  if (dist_method == "bray" && rlang::is_installed("parallelDist")) {
    # Bray-Curtis via threaded Manhattan
    d_L1 <- parallelDist::parDist(norm_mat, method = "manhattan",
                                  threads = n_threads)

    s <- rowSums(norm_mat)
    n <- length(s)
    denom <- numeric(n * (n - 1L) / 2L)
    k <- 1L
    for (i in seq_len(n - 1L)) {
      ni <- n - i
      denom[k:(k + ni - 1L)] <- s[i] + s[(i + 1L):n]
      k <- k + ni
    }

    bc_vals <- as.numeric(d_L1) / denom
    dist_obj <- d_L1
    dist_obj[] <- bc_vals

  } else if (rlang::is_installed("parallelDist") &&
             dist_method %in% pd_supported) {

    dist_obj <- parallelDist::parDist(
      norm_mat,
      method  = dist_method,
      threads = n_threads
    )

  } else {
    if (!rlang::is_installed("vegan")) {
      .ph_abort(
        "Requested distance method requires 'vegan'. Please install it.",
        step = "compute_distance"
      )
    }
    dist_obj <- vegan::vegdist(norm_mat, method = dist_method)
  }

  .ph_log_info("Distance matrix computation complete.",
               step = "compute_distance")

  ## --- NEW: attach normalized abundance matrix as attribute ----------------
  attr(dist_obj, "abundances") <- norm_mat

  dist_obj
}
#' @title Principal Coordinates Analysis (PCoA) on a Distance Matrix
#' @description
#' Performs PCoA on a distance matrix (typically from \code{compute_distance()}),
#' optionally correcting for negative eigenvalues, and returns coordinates,
#' eigenvalues, variance explained, and feature loadings.
#'
#' @param dist_obj A \code{dist} object returned by \code{compute_distance()}.
#'   The normalized abundance matrix used to compute the distances is expected
#'   to be attached as an attribute \code{"abundances"} (a numeric matrix with
#'   samples in rows and features in columns).
#' @param neg_correction Method for adjusting negative eigenvalues (if any).
#'   One of \code{"none"}, \code{"lingoes"}, or \code{"cailliez"}.
#'   Default is \code{"none"} (no correction). If set to \code{"lingoes"}
#'   or \code{"cailliez"}, the \pkg{vegan} package is required.
#' @param n_axes Number of principal coordinate axes to return in the
#'   sample scores. Default is 5.
#' @param top_features Integer number of top features to return in loadings
#'   (based on highest absolute scores on any axis). Default is 30.
#'
#' @return A list of class \code{"beta_pcoa"} with elements:
#' \item{sample_coords}{A tibble with sample coordinates on the first
#'   \code{n_axes} PCoA axes. Includes \code{sample_id} and one column per axis
#'   (\code{PCoA1}, \code{PCoA2}, ...).}
#' \item{eigenvalues}{Numeric vector of all eigenvalues from the PCoA.}
#' \item{var_explained}{A one-row tibble summarizing the percentage of
#'   variation explained by the first \code{n_axes} axes and the remainder
#'   as \code{\%Other}. Percentages are computed from the sum of positive
#'   eigenvalues.}
#' \item{feature_loadings}{A tibble of top feature loadings for the first
#'   \code{n_axes} axes. Each row is a feature, with columns for each axis
#'   showing the feature's loading. Only the top \code{top_features} features
#'   (by absolute loading on at least one axis) are included. Empty if the
#'   \code{"abundances"} attribute is missing or cannot be aligned.}
#'
#' @details
#' Negative eigenvalues indicate that the distances are not perfectly
#' Euclidean. If \code{neg_correction} is \code{"lingoes"} or \code{"cailliez"},
#' a correction is applied before PCoA using \code{vegan::wcmdscale()}.
#' The variance explained is computed using only positive eigenvalues.
#'
#' Feature loadings are computed as weighted averages of sample scores,
#' where weights are the feature abundances across samples (the normalized
#' abundance matrix attached as an attribute to \code{dist_obj}).
#'
#' @examples
#' \dontrun{
#'   dist_bc <- compute_distance(ps, value_col = "counts_hit",
#'                               method_normalization = "hellinger",
#'                               distance = "bray")
#'
#'   pcoa_res <- compute_pcoa(dist_bc, neg_correction = "none", n_axes = 3)
#'   pcoa_res$sample_coords
#'   pcoa_res$var_explained
#'   pcoa_res$feature_loadings
#' }
compute_pcoa <- function(dist_obj,
                         neg_correction = c("none", "lingoes", "cailliez"),
                         n_axes = 5,
                         top_features = 30) {
  neg_correction <- match.arg(neg_correction)

  # ---------------------------------------------------------------------------
  # 1) Check that dist_obj is a proper dist from compute_distance()
  # ---------------------------------------------------------------------------
  if (!inherits(dist_obj, "dist")) {
    .ph_abort(
      "`dist_obj` must be a `dist` object (e.g. returned by `compute_distance()`).",
      step = "compute_pcoa"
    )
  }

  d <- dist_obj
  n <- attr(d, "Size")
  labels <- attr(d, "Labels")

  # vegan requirement if correction requested
  if (!identical(neg_correction, "none") && !rlang::is_installed("vegan")) {
    .ph_abort(
      paste0(
        "Negative eigenvalue correction ('", neg_correction,
        "') requires the 'vegan' package."
      ),
      step = "compute_pcoa"
    )
  }

  # ---------------------------------------------------------------------------
  # 2) PCoA computation
  # ---------------------------------------------------------------------------
  .ph_log_info(
    "Performing Principal Coordinates Analysis",
    step = "compute_pcoa",
    bullets = if (identical(neg_correction, "none")) NULL else
      paste("using", neg_correction, "correction")
  )

  k_cmd <- min(max(1L, n_axes), n - 1L)

  pcoa_fit <- if (identical(neg_correction, "none")) {
    stats::cmdscale(d, eig = TRUE, k = k_cmd)
  } else {
    vegan::wcmdscale(d, eig = TRUE, k = k_cmd, add = neg_correction)
  }

  eig_vals <- as.numeric(pcoa_fit$eig %||% numeric(0L))

  coords <- as.matrix(pcoa_fit$points)
  if (is.null(coords)) {
    coords <- matrix(0, nrow = n, ncol = 0)
  }

  if (ncol(coords) > 0L) {
    colnames(coords) <- paste0("PCoA", seq_len(ncol(coords)))
  }

  if (!is.null(labels) && nrow(coords) == length(labels)) {
    rownames(coords) <- labels
  }

  # ---------------------------------------------------------------------------
  # 3) Sample coordinates (first n_axes)
  # ---------------------------------------------------------------------------
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

  # ---------------------------------------------------------------------------
  # 4) Variance explained table
  # ---------------------------------------------------------------------------
  pos_eig <- pmax(eig_vals, 0)
  sum_pos <- sum(pos_eig, na.rm = TRUE)

  if (sum_pos > 0) {
    pct_axes <- 100 * pos_eig[seq_len(k_use)] / sum_pos
    pct_other <- 100 * sum(pos_eig[-seq_len(k_use)], na.rm = TRUE) / sum_pos
  } else {
    pct_axes <- rep(NA_real_, k_use)
    pct_other <- NA_real_
  }

  names(pct_axes) <- paste0("%PCoA", seq_len(k_use))
  var_explained <- tibble::as_tibble_row(
    c(as.list(round(pct_axes, 3)), `%Other` = round(pct_other, 3))
  )

  # ---------------------------------------------------------------------------
  # 5) Feature loadings (use abundances attribute from dist_obj)
  # ---------------------------------------------------------------------------
  feature_loadings <- tibble::tibble()

  X <- attr(dist_obj, "abundances")
  if (is.null(X)) {
    .ph_warn(
      "No 'abundances' attribute found on `dist_obj`; skipping feature loadings.",
      step = "compute_pcoa"
    )
  } else {
    X <- as.matrix(X)

    coords_ids <- rownames(coords)
    X_ids <- rownames(X)

    if (is.null(coords_ids) || is.null(X_ids)) {
      .ph_warn(
        "Row names missing in coordinates or 'abundances'; cannot align samples for feature loadings.",
        step = "compute_pcoa"
      )
    } else {
      common_ids <- intersect(coords_ids, X_ids)

      if (length(common_ids) < 2L) {
        .ph_warn(
          "Insufficient overlap between distance labels and abundance rows; skipping feature loadings.",
          step = "compute_pcoa"
        )
      } else {
        # Use up to min(n_axes, ncol(coords), 10) axes for loadings
        ax_idx <- seq_len(min(n_axes, ncol(coords), 10L))
        U <- coords[common_ids, ax_idx, drop = FALSE]
        Xsub <- X[common_ids, , drop = FALSE]

        # weights: total abundance per feature
        w <- colSums(Xsub, na.rm = TRUE)
        keep_feats <- which(w > 0)

        if (length(keep_feats) > 0L) {
          # crossprod: features x axes
          S <- t(Xsub[, keep_feats, drop = FALSE]) %*% U
          S <- sweep(S, 1, w[keep_feats], "/")

          rownames(S) <- colnames(Xsub)[keep_feats]
          colnames(S) <- colnames(U)

          load_tbl <- tibble::as_tibble(S, rownames = "feature")

          if (!is.null(top_features) && is.finite(top_features)) {
            ax_names <- colnames(U)
            top_list <- unique(unlist(lapply(seq_along(ax_names), function(j) {
              ord <- order(abs(S[, j]), decreasing = TRUE)
              head(rownames(S)[ord], top_features)
            })))
            load_tbl <- dplyr::filter(load_tbl, .data$feature %in% top_list)
          }

          feature_loadings <- load_tbl
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 6) Return
  # ---------------------------------------------------------------------------
  result <- list(
    sample_coords    = sample_coords,
    eigenvalues      = eig_vals,
    var_explained    = var_explained,
    feature_loadings = feature_loadings
  )
  class(result) <- "beta_pcoa"

  .ph_log_info("PCoA analysis complete.", step = "compute_pcoa")

  result
}

#' @title Constrained Ordination (db-RDA / CAP) on Distance Matrix
#' @description
#' Performs distance-based redundancy analysis (constrained PCoA, a.k.a. CAP)
#' on a distance matrix using \pkg{vegan}::\code{capscale}, with optional
#' negative eigenvalue correction. Returns constrained sample scores, eigenvalues,
#' variance partitioning, and feature loadings.
#'
#' @param dist_obj A \code{dist} object returned by \code{compute_distance()}.
#'   The normalized abundance matrix used to compute the distances is expected
#'   to be attached as an attribute \code{"abundances"} (samples in rows,
#'   features in columns).
#' @param ps A \code{phip_data} object providing sample-level metadata in
#'   \code{ps$data_long}. This table must contain \code{sample_id} and all
#'   variables referenced on the right-hand side of \code{formula}.
#' @param formula An R formula specifying the constraints (independent variables)
#'   for the ordination, e.g. \code{~ type_person + age}. Do not include a
#'   response on the left-hand side; the distance matrix is provided via
#'   \code{dist_obj}.
#' @param neg_correction One of \code{"none"}, \code{"lingoes"}, \code{"cailliez"}.
#'   Method for negative eigenvalue correction. Default is \code{"none"}.
#'   This is passed to the \code{add} argument of \code{vegan::capscale()}.
#' @param top_features Integer number of top features to return in loadings
#'   (based on highest absolute scores on any constrained axis). Default is 30.
#'
#' @return A list of class \code{"beta_capscale"} with elements:
#' \item{sample_coords}{Tibble of sample scores on constrained axes
#'   (\code{CAP1}, \code{CAP2}, ...). Contains \code{sample_id} and coordinates.}
#' \item{eigenvalues}{Numeric vector of eigenvalues of the constrained axes.}
#' \item{variance_partition}{Tibble with total inertia and inertia partitioned
#'   into constrained and unconstrained components, with their proportion of total.}
#' \item{feature_loadings}{Tibble of top feature loadings for constrained axes
#'   (possibly empty if the \code{"abundances"} attribute is missing or cannot be aligned).}
#' \item{cap_model}{The full \code{vegan::capscale} model object for further
#'   inspection (e.g., \code{summary()}, \code{anova()}, etc.).}
#'
#' @examples
#' \dontrun{
#'   dist_bc <- compute_distance(ps, value_col = "counts_hit",
#'                               method_normalization = "hellinger",
#'                               distance = "bray")
#'
#'   cap_res <- compute_capscale(
#'     dist_bc,
#'     ps      = ps,
#'     formula = ~ type_person + age
#'   )
#'   cap_res$variance_partition
#'   cap_res$sample_coords
#'   cap_res$feature_loadings
#'   summary(cap_res$cap_model)
#' }
compute_capscale <- function(dist_obj,
                             ps,
                             formula,
                             neg_correction = c("none", "lingoes", "cailliez"),
                             top_features = 30) {
  neg_correction <- match.arg(neg_correction)

  if (!inherits(formula, "formula")) {
    .ph_abort("`formula` must be an R formula (e.g., ~ group + time).",
              step = "compute_capscale")
  }

  if (!inherits(dist_obj, "dist")) {
    .ph_abort(
      "`dist_obj` must be a `dist` object (e.g. returned by `compute_distance()`).",
      step = "compute_capscale"
    )
  }

  if (!inherits(ps, "phip_data")) {
    .ph_abort("`ps` must be a <phip_data> object.",
              step = "compute_capscale")
  }

  if (!rlang::is_installed("vegan")) {
    .ph_abort("`compute_capscale()` requires the 'vegan' package.",
              step = "compute_capscale")
  }

  # ---------------------------------------------------------------------------
  # 1) Distance + labels
  # ---------------------------------------------------------------------------
  d <- dist_obj
  labels <- attr(d, "Labels")
  n      <- attr(d, "Size")

  # Abundance matrix (normalized) from dist attributes (may be NULL)
  X_all <- attr(dist_obj, "abundances")

  # ---------------------------------------------------------------------------
  # 2) Metadata from ps$data_long + alignment + NA handling
  # ---------------------------------------------------------------------------
  dat <- ps$data_long
  if (is.null(dat)) {
    .ph_abort("`ps$data_long` is missing. Cannot construct metadata.",
              step = "compute_capscale")
  }

  dat_cols <- dplyr::tbl_vars(dat)
  if (!"sample_id" %in% dat_cols) {
    .ph_abort("`ps$data_long` must contain a `sample_id` column.",
              step = "compute_capscale")
  }

  # RHS terms from formula (constraints)
  rhs_terms <- attr(stats::terms(formula), "term.labels")
  if (length(rhs_terms) == 0L) {
    .ph_abort(
      "No constraints provided in formula (RHS is empty). Use compute_pcoa() for unconstrained ordination.",
      step = "compute_capscale"
    )
  }

  missing_vars <- setdiff(rhs_terms, dat_cols)
  if (length(missing_vars) > 0L) {
    .ph_abort(
      paste0(
        "The following variables from the formula are missing in `ps$data_long`: ",
        paste(missing_vars, collapse = ", ")
      ),
      step = "compute_capscale"
    )
  }

  .ph_log_info(
    "Building metadata from `ps$data_long`.",
    step = "compute_capscale"
  )

  # Wyciągamy unikalne metadane na poziomie sample
  meta_all <- dat |>
    dplyr::select(sample_id, dplyr::all_of(rhs_terms)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect() |>
    as.data.frame()

  if (nrow(meta_all) == 0L) {
    .ph_abort(
      "Constructed metadata has zero rows. Check that `ps$data_long` is not empty.",
      step = "compute_capscale"
    )
  }

  rownames(meta_all) <- meta_all$sample_id

  # Dopasowanie metadanych do kolejności w macierzy dystansów
  if (!is.null(labels)) {
    idx <- match(labels, meta_all$sample_id)
    missing_samples <- labels[is.na(idx)]
    if (length(missing_samples) > 0L) {
      .ph_abort(
        paste0(
          "The following samples from `dist_obj` are missing in `ps$data_long`: ",
          paste(missing_samples, collapse = ", ")
        ),
        step = "compute_capscale"
      )
    }
    meta_sub <- meta_all[idx, , drop = FALSE]
    rownames(meta_sub) <- labels
  } else {
    meta_sub <- meta_all
    labels   <- rownames(meta_sub)
    n        <- length(labels)
    .ph_warn(
      "No labels found in `dist_obj`; assuming metadata row order matches the distance order.",
      step = "compute_capscale"
    )
  }

  # Usuwamy próbki z NA w którejkolwiek zmiennej z formuły
  rhs_df <- meta_sub[, rhs_terms, drop = FALSE]
  keep   <- stats::complete.cases(rhs_df)

  if (!all(keep)) {
    dropped <- sum(!keep)
    .ph_log_info(
      paste0("Dropping ", dropped,
             " samples with missing values in constrained variables."),
      step = "compute_capscale"
    )
  }

  meta_df <- meta_sub[keep, , drop = FALSE]
  if (nrow(meta_df) == 0L) {
    .ph_abort(
      "All samples have missing values in constrained variables; cannot fit CAP.",
      step = "compute_capscale"
    )
  }

  keep_labels <- rownames(meta_df)

  # ---------------------------------------------------------------------------
  # 3) Subset distance and abundances to complete-case samples
  # ---------------------------------------------------------------------------
  # subset distance matrix
  mat_d    <- as.matrix(d)
  mat_d_sub <- mat_d[keep_labels, keep_labels, drop = FALSE]
  d        <- stats::as.dist(mat_d_sub)
  n        <- attr(d, "Size")
  labels   <- attr(d, "Labels")

  # subset abundances (if present)
  X_sub <- NULL
  if (!is.null(X_all)) {
    X_all <- as.matrix(X_all)
    if (!is.null(rownames(X_all))) {
      # align by sample names
      X_sub <- X_all[keep_labels, , drop = FALSE]
    } else {
      .ph_warn(
        "Abundance matrix has no row names; cannot align precisely with samples.",
        step = "compute_capscale"
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 4) Build formula for capscale: d_resp ~ RHS
  # ---------------------------------------------------------------------------
  d_resp <- d

  rhs_txt <- gsub("^~", "", paste(deparse(formula), collapse = " "))
  cap_formula <- stats::as.formula(
    paste("d_resp ~", rhs_txt),
    env = environment()
  )

  add_arg <- if (identical(neg_correction, "none")) NULL else neg_correction

  .ph_log_info(
    "Fitting constrained ordination (CAP/db-RDA)",
    step = "compute_capscale",
    bullets = c(
      paste("formula:", paste(deparse(formula), collapse = " ")),
      if (!is.null(add_arg)) paste("neg_correction:", add_arg)
    )
  )

  # ---------------------------------------------------------------------------
  # 5) Fit capscale model
  # ---------------------------------------------------------------------------
  cap_fit <- vegan::capscale(cap_formula, data = meta_df, add = add_arg)

  # ---------------------------------------------------------------------------
  # 6) Sample scores on constrained axes
  # ---------------------------------------------------------------------------
  rank_constrained <- cap_fit$CCA$rank %||% 0L

  if (rank_constrained > 0L) {
    site_scores <- vegan::scores(
      cap_fit,
      display = "sites",
      choices = seq_len(rank_constrained)
    )
    pts <- as.matrix(site_scores)
  } else {
    pts <- matrix(0, nrow = n, ncol = 0,
                  dimnames = list(labels, character(0L)))
  }

  if (ncol(pts) > 0L && is.null(colnames(pts))) {
    colnames(pts) <- paste0("CAP", seq_len(ncol(pts)))
  }

  if (!is.null(labels) && nrow(pts) == length(labels)) {
    rownames(pts) <- labels
  }

  sample_coords <- tibble::as_tibble(pts, rownames = "sample_id")

  # ---------------------------------------------------------------------------
  # 7) Eigenvalues and variance partition
  # ---------------------------------------------------------------------------
  eig_constrained <- cap_fit$CCA$eig %||% numeric()

  tot_inertia   <- cap_fit$tot.chi
  cons_inertia  <- cap_fit$CCA$tot.chi %||% sum(cap_fit$CCA$eig %||% 0)
  uncon_inertia <- cap_fit$CA$tot.chi  %||% sum(cap_fit$CA$eig %||% 0)

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
  # 8) Feature loadings on constrained axes (using X_sub)
  # ---------------------------------------------------------------------------
  feature_loadings <- tibble::tibble()

  if (!is.null(X_sub) && rank_constrained > 0L) {
    coords_ids <- rownames(pts)
    X_ids      <- rownames(X_sub)

    if (is.null(coords_ids) || is.null(X_ids)) {
      .ph_warn(
        "Row names missing in sample scores or abundance matrix; cannot align for feature loadings.",
        step = "compute_capscale"
      )
    } else {
      common_ids <- intersect(coords_ids, X_ids)

      if (length(common_ids) < 2L) {
        .ph_warn(
          "Insufficient overlap between distance labels and abundance rows; skipping feature loadings.",
          step = "compute_capscale"
        )
      } else {
        ax_idx <- seq_len(min(rank_constrained, 10L))
        U    <- pts[common_ids, ax_idx, drop = FALSE]
        Xsub <- X_sub[common_ids, , drop = FALSE]

        w <- colSums(Xsub, na.rm = TRUE)
        keep_feats <- which(w > 0)

        if (length(keep_feats) > 0L) {
          S <- t(Xsub[, keep_feats, drop = FALSE]) %*% U
          S <- sweep(S, 1, w[keep_feats], "/")

          rownames(S) <- colnames(Xsub)[keep_feats]
          colnames(S) <- colnames(U)

          load_tbl <- tibble::as_tibble(S, rownames = "feature")

          if (!is.null(top_features) && is.finite(top_features)) {
            ax_names <- colnames(U)
            top_list <- unique(unlist(lapply(seq_along(ax_names), function(j) {
              ord <- order(abs(S[, j]), decreasing = TRUE)
              head(rownames(S)[ord], top_features)
            })))
            load_tbl <- dplyr::filter(load_tbl, .data$feature %in% top_list)
          }

          feature_loadings <- load_tbl
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 9) Return
  # ---------------------------------------------------------------------------
  result <- list(
    sample_coords      = sample_coords,
    eigenvalues        = as.numeric(eig_constrained),
    variance_partition = variance_partition,
    feature_loadings   = feature_loadings,
    cap_model          = cap_fit
  )
  class(result) <- "beta_capscale"

  .ph_log_info("CAP analysis complete.", step = "compute_capscale")

  result
}

#' @title PERMANOVA with Global and Post-hoc Tests on Beta Diversity
#' @description Performs PERMANOVA (adonis2) on a distance matrix for overall group/time effects,
#' and optionally conducts post-hoc pairwise or contrast tests (e.g., between each pair of groups, etc.).
#' Supports stratified permutations for repeated measures.
#'
#' @param dist_obj A \code{dist} object of distances between samples
#'   (e.g., output of \code{compute_distance()}).
#' @param ps A \code{phip_data} object providing sample-level metadata in
#'   \code{ps$data_long}. This table must contain \code{sample_id} and the
#'   columns specified in \code{group_col}, \code{time_col}, and optionally
#'   \code{subject_col}.
#' @param group_col Name of the grouping column in \code{ps$data_long}
#'   (between-subject factor). Use \code{NULL} if no group factor.
#' @param time_col Name of the time factor column in \code{ps$data_long}
#'   (within-subject factor for longitudinal data). Use \code{NULL} if not
#'   applicable. This should be a \emph{categorical} factor for this function
#'   (continuous time not supported).
#' @param subject_col Name of the subject identifier column in \code{ps$data_long}
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
compute_permanova <- function(dist_obj,
                              ps,
                              group_col = NULL,
                              time_col = NULL,
                              subject_col = "subject_id",
                              permutations = 999,
                              contrasts = "none",
                              baseline_level = NULL) {

  if (!inherits(dist_obj, "dist")) {
    .ph_abort(
      "`dist_obj` must be a `dist` object (e.g., from `compute_distance()`).",
      step = "compute_permanova"
    )
  }
  if (!inherits(ps, "phip_data")) {
    .ph_abort("`ps` must be a <phip_data> object.",
              step = "compute_permanova")
  }
  if (!rlang::is_installed("vegan")) {
    .ph_abort("`compute_permanova()` requires the 'vegan' package.",
              step = "compute_permanova")
  }

  # ---------------------------------------------------------------------------
  # 0) Prepare result collector
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
  # 1) Start from distances + labels
  # ---------------------------------------------------------------------------
  d_full  <- dist_obj
  labels_full <- attr(d_full, "Labels")
  n_full  <- attr(d_full, "Size")

  # ---------------------------------------------------------------------------
  # 2) Build metadata from ps$data_long and align to dist labels
  # ---------------------------------------------------------------------------
  dat <- ps$data_long
  if (is.null(dat)) {
    .ph_abort("`ps$data_long` is missing. Cannot construct metadata.",
              step = "compute_permanova")
  }

  dat_cols <- dplyr::tbl_vars(dat)
  if (!"sample_id" %in% dat_cols) {
    .ph_abort("`ps$data_long` must contain a `sample_id` column.",
              step = "compute_permanova")
  }

  has_group   <- !is.null(group_col)   && group_col   %in% dat_cols
  has_time    <- !is.null(time_col)    && time_col    %in% dat_cols
  has_subject <- !is.null(subject_col) && subject_col %in% dat_cols

  if (!is.null(group_col) && !has_group) {
    .ph_abort(
      paste0("Column `", group_col, "` not found in `ps$data_long`."),
      step = "compute_permanova"
    )
  }
  if (!is.null(time_col) && !has_time) {
    .ph_abort(
      paste0("Column `", time_col, "` not found in `ps$data_long`."),
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

  .ph_log_info(
    "Building metadata from `ps$data_long`.",
    step = "compute_permanova"
  )

  meta_all <- dat |>
    dplyr::select(dplyr::all_of(cols_needed)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect() |>
    as.data.frame()

  if (nrow(meta_all) == 0L) {
    .ph_abort(
      "Constructed metadata has zero rows. Check that `ps$data_long` is not empty.",
      step = "compute_permanova"
    )
  }
  rownames(meta_all) <- meta_all$sample_id

  # Align metadata to distance labels
  if (!is.null(labels_full)) {
    idx_align <- match(labels_full, meta_all$sample_id)
    missing_samples <- labels_full[is.na(idx_align)]
    if (length(missing_samples) > 0L) {
      .ph_abort(
        paste0(
          "The following samples from `dist_obj` are missing in `ps$data_long`: ",
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
      "No labels found in `dist_obj`; assuming metadata row order matches the distance order.",
      step = "compute_permanova"
    )
  }

  # Coerce group/time to factor for safety
  if (has_group) {
    meta_sub[[group_col]] <- as.factor(meta_sub[[group_col]])
  }
  if (has_time) {
    meta_sub[[time_col]] <- as.factor(meta_sub[[time_col]])
  }

  # ---------------------------------------------------------------------------
  # 3) Drop samples with NA in constrained variables (+ subject if used)
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
      paste0("Dropping ", dropped,
             " samples with missing values in constrained/strata variables."),
      step = "compute_permanova"
    )
  }

  meta_df <- meta_sub[keep, , drop = FALSE]
  if (nrow(meta_df) == 0L) {
    .ph_abort(
      "All samples have missing values in constrained/strata variables; cannot run PERMANOVA.",
      step = "compute_permanova"
    )
  }
  keep_labels <- rownames(meta_df)

  # Subset distance matrix to kept samples
  mat_d_full <- as.matrix(d_full)
  mat_d_sub  <- mat_d_full[keep_labels, keep_labels, drop = FALSE]
  d <- stats::as.dist(mat_d_sub)
  labels <- attr(d, "Labels")
  n      <- attr(d, "Size")

  # Update factor presence after NA-drop
  has_group <- has_group && length(unique(meta_df[[group_col]])) > 1L
  has_time  <- has_time  && length(unique(meta_df[[time_col]]))  > 1L

  # ---------------------------------------------------------------------------
  # 4) Global PERMANOVA
  # ---------------------------------------------------------------------------
  rhs_terms <- character(0L)
  if (has_group) rhs_terms <- c(rhs_terms, group_col)
  if (has_time)  rhs_terms <- c(rhs_terms, time_col)
  if (has_group && has_time) {
    rhs_terms <- c(rhs_terms, paste(group_col, "*", time_col))
  }

  if (length(rhs_terms) == 0L) {
    .ph_log_info(
      "Global PERMANOVA skipped (insufficient number of factor levels).",
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
      "Running global PERMANOVA",
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
      .ph_warn("Global PERMANOVA failed; no global results.",
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
  # 5) Post-hoc contrasts (pairwise / each_vs_rest / baseline)
  # ---------------------------------------------------------------------------
  contrasts <- tolower(unique(contrasts))
  contrasts[contrasts == "group_vs_rest"] <- "each_vs_rest"

  # Determine baseline if needed
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
        paste("Skipping test", contrast_label,
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

  # ---------------- pairwise -----------------------------------------------
  if ("pairwise" %in% contrasts) {
    # Pairwise group comparisons
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

    # Pairwise time comparisons
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

  # ---------------- each_vs_rest -------------------------------------------
  if ("each_vs_rest" %in% contrasts) {
    # Each group vs rest
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

    # Each time vs rest
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

  # ---------------- baseline ------------------------------------------------
  if ("baseline" %in% contrasts) {
    if (is.null(bl)) {
      .ph_warn(
        "Baseline contrast requested but `baseline_level` not properly specified; skipping baseline tests.",
        step = "compute_permanova"
      )
    } else {
      # Baseline for time
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

      # Baseline for group
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
  # 6) Combine and return
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
#' @param ps A \code{phip_data} object providing sample-level metadata in
#'   \code{ps$data_long}. This table must contain \code{sample_id} and the
#'   columns specified in \code{group_col} and/or \code{time_col}.
#' @param group_col Name of the group factor column in \code{ps$data_long}
#'   (between-subjects). Use \code{NULL} if no group factor.
#' @param time_col Name of the time factor column in \code{ps$data_long}
#'   (within-subjects, categorical only). Use \code{NULL} if not applicable.
#' @param subject_col Name of subject identifier column (for reference only;
#'   not used directly in dispersion test calculations, but kept for API
#'   consistency). Default \code{"subject_id"}.
#' @param permutations Number of permutations for significance testing in
#'   \code{vegan::permutest}. Default 999.
#' @param contrasts Which dispersion contrasts to perform. Options:
#'   \code{"none"} (default), \code{"pairwise"}, \code{"each_vs_rest"},
#'   \code{"baseline"}. Interpretation analogiczna do \code{compute_permanova},
#'   ale zastosowana do dyspersji.
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
compute_dispersion <- function(dist_obj,
                               ps,
                               group_col = NULL,
                               time_col = NULL,
                               subject_col = "subject_id",
                               permutations = 999,
                               contrasts = "none",
                               baseline_level = NULL) {

  if (!inherits(dist_obj, "dist")) {
    .ph_abort(
      "`dist_obj` must be a `dist` object (e.g., from `compute_distance()`).",
      step = "compute_dispersion"
    )
  }
  if (!inherits(ps, "phip_data")) {
    .ph_abort("`ps` must be a <phip_data> object.",
              step = "compute_dispersion")
  }
  if (!rlang::is_installed("vegan")) {
    .ph_abort("`compute_dispersion()` requires the 'vegan' package.",
              step = "compute_dispersion"
    )
  }

  # ---------------------------------------------------------------------------
  # 1) Prepare metadata from ps$data_long and align to dist labels
  # ---------------------------------------------------------------------------
  d_full      <- dist_obj
  labels_full <- attr(d_full, "Labels")

  dat <- ps$data_long
  if (is.null(dat)) {
    .ph_abort("`ps$data_long` is missing. Cannot construct metadata.",
              step = "compute_dispersion")
  }

  dat_cols <- dplyr::tbl_vars(dat)
  if (!"sample_id" %in% dat_cols) {
    .ph_abort("`ps$data_long` must contain a `sample_id` column.",
              step = "compute_dispersion")
  }

  has_group <- !is.null(group_col) && group_col %in% dat_cols
  has_time  <- !is.null(time_col)  && time_col  %in% dat_cols

  if (!is.null(group_col) && !has_group) {
    .ph_abort(
      paste0("Column `", group_col, "` not found in `ps$data_long`."),
      step = "compute_dispersion"
    )
  }
  if (!is.null(time_col) && !has_time) {
    .ph_abort(
      paste0("Column `", time_col, "` not found in `ps$data_long`."),
      step = "compute_dispersion"
    )
  }

  cols_needed <- c(
    "sample_id",
    if (has_group) group_col else character(0L),
    if (has_time)  time_col  else character(0L)
  )
  cols_needed <- unique(cols_needed)

  .ph_log_info(
    "Building metadata from `ps$data_long`.",
    step = "compute_dispersion"
  )

  meta_all <- dat |>
    dplyr::select(dplyr::all_of(cols_needed)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect() |>
    as.data.frame()

  if (nrow(meta_all) == 0L) {
    .ph_abort(
      "Constructed metadata has zero rows. Check that `ps$data_long` is not empty.",
      step = "compute_dispersion"
    )
  }
  rownames(meta_all) <- meta_all$sample_id

  # Align metadata to distance labels
  if (!is.null(labels_full)) {
    idx_align <- match(labels_full, meta_all$sample_id)
    missing_samples <- labels_full[is.na(idx_align)]
    if (length(missing_samples) > 0L) {
      .ph_abort(
        paste0(
          "The following samples from `dist_obj` are missing in `ps$data_long`: ",
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
      "No labels found in `dist_obj`; assuming metadata row order matches the distance order.",
      step = "compute_dispersion"
    )
  }

  # Coerce grouping vars to factor
  if (has_group) {
    meta_sub[[group_col]] <- as.factor(meta_sub[[group_col]])
  }
  if (has_time) {
    meta_sub[[time_col]] <- as.factor(meta_sub[[time_col]])
  }

  # ---------------------------------------------------------------------------
  # 2) Drop samples with NA in group/time and subset distance matrix
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
      paste0("Dropping ", dropped,
             " samples with missing values in dispersion grouping variables."),
      step = "compute_dispersion"
    )
  }

  meta_df <- meta_sub[keep, , drop = FALSE]
  if (nrow(meta_df) == 0L) {
    .ph_abort(
      "All samples have missing values in grouping variables; cannot run dispersion tests.",
      step = "compute_dispersion"
    )
  }
  keep_labels <- rownames(meta_df)

  d_mat_full <- as.matrix(d_full)
  d_mat_sub  <- d_mat_full[keep_labels, keep_labels, drop = FALSE]
  d          <- stats::as.dist(d_mat_sub)
  labels     <- attr(d, "Labels")

  # Re-evaluate factor presence after NA drop
  has_group <- has_group && length(unique(meta_df[[group_col]])) > 1L
  has_time  <- has_time  && length(unique(meta_df[[time_col]]))  > 1L

  # ---------------------------------------------------------------------------
  # 3) Prepare collectors for distances and tests
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
  # 4) Global dispersion tests (group, time, group:time)
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
        "`time_col` is numeric; continuous dispersion by time not supported. Skipping time dispersion test.",
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
  # 5) Helper for subset dispersion tests
  # ---------------------------------------------------------------------------
  run_disp_test <- function(idx, fac_vec, scope_lab, contrast_lab) {
    if (length(unique(fac_vec)) < 2L || min(table(fac_vec)) < 2L) {
      .ph_log_info(
        paste("Skipping dispersion test for", contrast_lab, "- not enough data."),
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
  # 6) Post-hoc contrasts
  # ---------------------------------------------------------------------------
  contrasts <- tolower(unique(contrasts))
  contrasts[contrasts == "group_vs_rest"] <- "each_vs_rest"

  # Pairwise
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

  # Each vs rest
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

  # Baseline
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
  # 7) Combine and return
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

  # vegan fallback section (not-recommended; slow)
  if (!rlang::is_installed("vegan")) {
    .ph_abort("Requested distance method requires 'vegan'. Please install it.")
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
#' @param ps A [`phip_data`] object.
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


# get the phiper-package specific options; can be modified by the user
.ph_opt <- function(key, default) {
  op <- getOption(key, NULL)
  if (is.null(op)) default else op
}

# simple weighted-average "peptides/feature scores" on axes (approximate biplot)
# U: n x a sample scores; X: n x p (same samples), non-negative abundances used
# in distance (post-normalization)
.feature_loadings_wa <- function(U,
                                 X,
                                 feature_names = colnames(X),
                                 top = .ph_opt(
                                   "phiper.beta.loadings_top",
                                   30
                                 )) {
  # coerce to matrices
  U <- as.matrix(U)
  X <- as.matrix(X)

  # align rows by rownames if both provided
  if (!is.null(rownames(U)) && !is.null(rownames(X))) {
    common <- intersect(rownames(U), rownames(X))
    if (length(common) == 0L) {
      return(tibble::tibble())
    } # no overlap
    U <- U[common, , drop = FALSE]
    X <- X[common, , drop = FALSE]
  }

  # axis names
  ax_names <- colnames(U)
  if (is.null(ax_names) || anyNA(ax_names) || any(ax_names == "")) {
    ax_names <- paste0("PCoA", seq_len(ncol(U)))
    colnames(U) <- ax_names
  }

  # feature names
  if (is.null(feature_names)) feature_names <- paste0("feat_", seq_len(ncol(X)))

  # weights per feature (columns of X)
  w <- colSums(X, na.rm = TRUE)
  keep <- which(w > 0)
  if (!length(keep)) {
    return(tibble::tibble())
  }

  # weighted averages: S (p_kept x a) = WA biplot coords for each axis
  S <- crossprod(X[, keep, drop = FALSE], U) # sum_i X_ij * U_ia
  S <- sweep(S, 1L, w[keep], "/") # divide each feature-row by its weight
  rownames(S) <- feature_names[keep]
  colnames(S) <- ax_names

  # build wide tibble (feature + one column per axis)
  out <- tibble::as_tibble(S, rownames = "feature", .name_repair = "minimal")

  # keep only features that are in the top |loading| for at least one axis
  # (if top is set)
  if (!is.null(top) && is.finite(top)) {
    # union of top 'top' features per axis
    top_feats <- unique(unlist(lapply(seq_along(ax_names), function(j) {
      ord <- order(abs(S[, j]), decreasing = TRUE, na.last = TRUE)
      head(rownames(S)[ord], top)
    })))
    out <- dplyr::filter(out, .data$feature %in% top_feats)
  } else {
    top_feats <- out$feature
  }

  # ---- Axis-block ordering (default) ----
  idx_all <- match(top_feats, rownames(S))
  S_sub <- S[idx_all, , drop = FALSE]
  feats <- rownames(S_sub)

  picked <- setNames(rep(FALSE, length(feats)), feats)
  order_names <- character(0)

  per_axis_cap <- if (!is.null(top) && is.finite(top)) top else Inf
  for (ax in colnames(S_sub)) {
    ord <- order(abs(S_sub[, ax]), decreasing = TRUE, na.last = TRUE)
    cand <- feats[ord]
    cand <- cand[!picked[cand]]
    if (is.finite(per_axis_cap)) cand <- head(cand, per_axis_cap)
    order_names <- c(order_names, cand)
    picked[cand] <- TRUE
  }

  leftovers <- names(picked)[!picked]
  if (length(leftovers)) {
    maxabs <- do.call(pmax, lapply(
      colnames(S_sub),
      function(ax) abs(S_sub[leftovers, ax])
    ))
    leftovers <- leftovers[order(maxabs, decreasing = TRUE, na.last = TRUE)]
    order_names <- c(order_names, leftovers)
  }

  out[match(order_names, out$feature), , drop = FALSE]
}

# ---- contrast builders -------------------------------------------------------

# levels -> list of pairs for "pairwise"
.pairs_of <- function(levels) {
  if (length(levels) < 2) {
    return(list())
  } # only one group
  utils::combn(levels, 2, simplify = FALSE)
}

# "each_vs_rest" builds a two-level factor for each reference level
.each_vs_rest_factors <- function(vec, levels) {
  lapply(levels, function(lv) {
    factor(ifelse(vec == lv, lv, "other"),
      levels = c(lv, "other")
    )
  })
}

# "baseline" builds a factor comparing baseline_level vs others (two-level)
.baseline_factor <- function(vec, baseline_level) {
  if (is.null(baseline_level)) {
    return(NULL)
  }
  if (!(baseline_level %in% as.character(unique(vec)))) {
    return(NULL)
  }
  factor(
    ifelse(vec == baseline_level, baseline_level,
      paste0("not_", baseline_level)
    ),
    levels = c(baseline_level, paste0("not_", baseline_level))
  )
}

# ---- adonis2 / betadisper wrappers -------------------------------------------
# runs the PERMANOVA + convenience wrappers + phiper-style logging
.adonis_terms <- function(dist_obj,
                          meta_df,
                          rhs_txt,
                          permutations,
                          step_label,
                          strata = NULL,
                          n_threads,
                          log_results = TRUE) { # turn off if you ever need silence
  # ensure rownames for adonis2
  if (is.null(rownames(meta_df))) rownames(meta_df) <- meta_df$sample_id

  # build formula in the parent env (so dist_obj is visible)
  form <- stats::as.formula(paste("dist_obj ~", rhs_txt), env = environment())

  # run the test (adonis2 is single-threaded)
  out <- try(
    vegan::adonis2(form,
      data         = droplevels(meta_df),
      permutations = permutations,
      by           = "terms",
      strata       = strata,
      parallel     = n_threads
    ),
    silent = TRUE
  )

  if (inherits(out, "try-error")) {
    .ph_warn("PERMANOVA failed; skipping.",
      step = step_label,
      bullets = as.character(out)
    )
    return(NULL)
  }

  df <- as.data.frame(out)
  df$term <- rownames(df)
  res <- tibble::tibble(
    term       = df$term,
    p_value    = df$`Pr(>F)`,
    F_stat     = df$F,
    R2         = df$R2,
    n_perm     = permutations,
    p_at_floor = !is.na(df$`Pr(>F)`) & df$`Pr(>F)` <= 1 / (permutations + 1)
  )

  # ---- phiper-style log ------------------------------------------------------
  if (isTRUE(log_results)) {
    # header bullets about the run
    bullets_head <- c(
      sprintf("model: dist ~ %s  (by = terms)", rhs_txt),
      sprintf(
        "permutations: %d%s",
        permutations,
        if (is.null(strata)) {
          ""
        } else {
          sprintf(" (stratified: %d strata)", length(unique(strata)))
        }
      )
    )

    .ph_log_info(
      "PERMANOVA results",
      step    = step_label,
      bullets = c(bullets_head)
    )
  }

  res
}

# calculate the homogeneity of group dispersiosn (one of the PERMANOVA
# assumptions)
.betadisper_table <- function(dist_obj,
                              fac,
                              step_label) {
  # check if sufficient levels
  if (length(unique(fac)) < 2) {
    return(NULL)
  }

  # calculate the distances to centroids for each group separately
  d <- try(vegan::betadisper(dist_obj, fac), silent = TRUE)

  # soft warning if not succeeded
  if (inherits(d, "try-error")) {
    .ph_warn("betadisper failed; skipping.", step = step_label)
    return(NULL)
  }

  # return the distances itself (NO p-values!!!)
  tibble::tibble(
    sample_id = names(d$distances),
    distance = as.numeric(d$distances),
    group_disp = as.character(d$group)
  )
}

# actually test the dispersions from .betadisper_table
.betadisper_test <- function(dist_obj,
                             fac,
                             permutations,
                             step_label,
                             n_threads) {
  # early exit if only one group
  if (length(unique(fac)) < 2) {
    return(NULL)
  }

  # calculate the disper (a little overhead as i am calculating it
  # previously already --> optimize later)
  d <- try(vegan::betadisper(dist_obj, fac), silent = TRUE)
  if (inherits(d, "try-error")) {
    return(NULL)
  }

  # perform the test
  pt <- try(
    vegan::permutest(d,
      permutations = permutations,
      parallel = n_threads
    ),
    silent = TRUE
  )
  if (inherits(pt, "try-error")) {
    return(NULL)
  }

  # return the results
  p <- tryCatch(pt$tab[1, "Pr(>F)"], error = function(e) NA_real_)
  tibble::tibble(term = "dispersion", p_value = p, n_perm = permutations)
}

# ---------- helpers for permutation ANCOVA on dispersion (continuous time) ----

# Freedman–Lane permutation test for a single term in a linear model
# y: numeric response; X_full/X_red: model matrices; nperm: permutations
.FL_perm_F <- function(y,
                       X_full,
                       X_red,
                       nperm = 999,
                       seed = NULL,
                       n_threads = 1L,
                       block_size = 1000L) {
  stopifnot(is.numeric(y), is.matrix(X_full), is.matrix(X_red))
  n <- length(y)
  if (nrow(X_full) != n || nrow(X_red) != n) stop("X dims ≠ length(y)")

  if (!is.null(seed)) set.seed(seed)

  # --- Precompute QR once (big speedup vs refitting) ------------------------
  qr_full <- qr(X_full)
  r_full <- qr_full$rank
  qr_red <- qr(X_red)
  r_red <- qr_red$rank

  # Observed fits (also via QR — no model.frame overhead)
  coef_red <- qr.coef(qr_red, y)
  mu_red <- as.numeric(X_red %*% coef_red)
  res_red <- as.numeric(y - mu_red)
  SSE_red <- sum(res_red^2)

  coef_full <- qr.coef(qr_full, y)
  res_full <- as.numeric(y - X_full %*% coef_full)
  SSE_full <- sum(res_full^2)

  df_full <- n - r_full
  df_term <- r_full - r_red
  SS_term_obs <- SSE_red - SSE_full
  F_obs <- (SS_term_obs / df_term) / (SSE_full / df_full)

  # Worker to compute F for one permutation ordering
  .F_for_perm <- function(ord) {
    y_perm <- mu_red + res_red[ord]
    # reuse QR:
    coef_red_p <- qr.coef(qr_red, y_perm)
    coef_full_p <- qr.coef(qr_full, y_perm)
    SSE_red_p <- sum((y_perm - as.numeric(X_red %*% coef_red_p))^2)
    SSE_full_p <- sum((y_perm - as.numeric(X_full %*% coef_full_p))^2)
    SS_term_p <- SSE_red_p - SSE_full_p
    (SS_term_p / df_term) / (SSE_full_p / df_full)
  }

  # Count how many permuted F >= observed F (Freedman-Lane)
  p_ge <- 1L

  # --- Parallel over permutations in blocks ---------------------------------
  do_block <- function(K) {
    # draw K permutations on the master to keep RNG reproducible across OS/backends
    ords <- replicate(K, sample.int(n, n, replace = FALSE), simplify = FALSE)

    if (n_threads <= 1L) {
      Fs <- lapply(ords, .F_for_perm)
    } else if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_threads)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      # export only the small objects; big matrices go via closures/environment once
      parallel::clusterExport(cl,
        varlist = c(
          "X_full", "X_red", "qr_full", "qr_red", "mu_red", "res_red",
          "df_full", "df_term", ".F_for_perm"
        ),
        envir = environment()
      )
      Fs <- parallel::parLapplyLB(cl, ords, .F_for_perm)
    } else {
      Fs <- parallel::mclapply(ords, .F_for_perm, mc.cores = n_threads)
    }
    sum(vapply(Fs, function(z) !is.na(z) && z >= F_obs, logical(1)))
  }

  left <- nperm
  while (left > 0L) {
    K <- min(block_size, left)
    p_ge <- p_ge + do_block(K)
    left <- left - K
  }

  tibble::tibble(
    F_stat = F_obs,
    p_value = p_ge / (nperm + 1),
    n_perm = nperm
  )
}

# Build tidy tests for dispersion with continuous time for a given data frame
# df needs: dist (numeric), fac (2+ level factor), time_cont (numeric)
# returns 2 rows: dispersion|time (main effect) & dispersion:time_cont (interaction)
.dispersion_tests_continuous <- function(df, scope_label, contrast_label, nperm,
                                         seed = NULL, n_threads = 1L, block_size = 1000L) {
  fac <- droplevels(df$fac)
  if (nlevels(fac) < 2 || min(table(fac)) < 2) {
    return(NULL)
  }

  X_full <- model.matrix(~ fac * time_cont, df) # fac + time + interaction
  X_red_main <- model.matrix(~ time_cont + fac:time_cont, df) # test fac|time
  X_red_int <- model.matrix(~ fac + time_cont, df) # test fac:time

  y <- df$dist

  tm_main <- .FL_perm_F(y, X_full, X_red_main,
    nperm = nperm,
    seed = seed, n_threads = n_threads, block_size = block_size
  )
  tm_int <- .FL_perm_F(y, X_full, X_red_int,
    nperm = nperm,
    seed = seed, n_threads = n_threads, block_size = block_size
  )

  out_main <- dplyr::mutate(tm_main,
    term = "dispersion|time",
    scope = scope_label,
    contrast = contrast_label
  )
  out_int <- dplyr::mutate(tm_int,
    term = "dispersion:time_cont",
    scope = scope_label,
    contrast = contrast_label
  )
  dplyr::bind_rows(out_main, out_int)
}

# ------------------------------------------------------------------------------
# Core engine (per table subset & per rank), with method_pcoa and corrections
# ------------------------------------------------------------------------------
.compute_beta_block <- function(tbl,
                                view_name,
                                group_col,
                                ranks,
                                fc_threshold = NULL,
                                method_normalization = c(
                                  "auto",
                                  "relative",
                                  "hellinger",
                                  "log",
                                  "none"
                                ),
                                distance = "bray",
                                permutations = 999,
                                time_col = NULL,
                                carry_cols = NULL,
                                filter_rank = NULL,
                                baseline_level = NULL,
                                contrasts = "none",
                                mtp = NULL,
                                map_provider,
                                n_threads = 1L,
                                method_pcoa = c(
                                  "joint",
                                  "separate_group",
                                  "separate_time",
                                  "separate_all",
                                  "cap"
                                ),
                                neg_correction = c(
                                  "none",
                                  "lingoes",
                                  "cailliez"
                                ),
                                time_force_continuous = FALSE) {
  # arg handling and tidy eval
  .data <- rlang::.data
  method_normalization <- match.arg(method_normalization)
  method_pcoa <- match.arg(method_pcoa)
  neg_correction <- match.arg(neg_correction)

  # ---------- required columns ------------------------------------------------
  need <- c("sample_id", "peptide_id") # bare minimum
  if (!is.null(group_col)) need <- union(need, group_col)

  # depending on how you define existence
  if (!is.null(fc_threshold)) {
    need <- union(need, "fold_change")
  } else {
    need <- union(need, "exist")
  }

  # abort when required is mising
  miss <- setdiff(need, colnames(tbl))
  if (length(miss)) {
    .ph_abort(
      headline = "Missing required columns.",
      step = sprintf("beta-div (%s) input validation", view_name),
      bullets = sprintf("missing: %s", paste(add_quotes(miss, 1L),
        collapse = ", "
      ))
    )
  }

  # ---------- unify group -----------------------------------------------------
  # basically the whole code is based on grouping, so the grouping column can
  # not be absent, it would complicate things --> add a column with single group
  # instead --> problem solved
  tbl <- if (is.null(group_col)) {
    dplyr::mutate(tbl, group = "All samples")
  } else {
    dplyr::mutate(tbl, group = .data[[group_col]])
  }

  # ---------- presence rule ---------------------------------------------------
  # presence can be defined simply as exist, or by fold_change treshold --> when
  # the user wants a sensitivity analysis and filtering the peptides based on
  # the signal strength
  used_exist_rule <- is.null(fc_threshold)
  pres_tbl <- if (used_exist_rule) {
    dplyr::filter(tbl, .data$exist > 0)
  } else {
    dplyr::filter(tbl, .data$fold_change > !!fc_threshold)
  }

  pres_tbl <- dplyr::collect(pres_tbl) # small collect is ok

  # ---------- time handling ---------------------------------------------------
  # right now only categorical times are supported, idk how to do the same for
  # purely continuous time; an option would be to group the data into bins based
  # on the continuous time, but the question is, how to actually do this? what
  # are the best cutoffs for bins etc.? i dont really have time for this now,
  # but i feel like it would require complete analysis/paradigm switch in how we
  # actually perform the beta_diveristy --> to think about later; for now i
  # generalized this approach prolly as far as meaningfully possible
  # ---------- time handling (auto + CAP-only continuous) ----------------------
  has_time <- !is.null(time_col)
  time_levels <- NULL
  time_is_continuous <- FALSE

  if (has_time) {
    if (!(time_col %in% colnames(tbl))) {
      .ph_abort(
        headline = "time_col not found in input.",
        step = "time handling",
        bullets = add_quotes(time_col, 1L)
      )
    }

    raw_time <- pres_tbl[[time_col]]
    thr <- .ph_opt("phiper.time.numeric_frac_threshold", 1.00)

    suppressWarnings(time_num_try <- as.numeric(as.character(raw_time)))
    frac_numeric <- mean(!is.na(time_num_try))
    cont_candidate <- is.numeric(raw_time) || (frac_numeric >= thr)

    if (isTRUE(time_force_continuous)) cont_candidate <- TRUE

    if (cont_candidate && identical(method_pcoa, "cap")) {
      if (all(!is.na(time_num_try))) {
        time_is_continuous <- TRUE
        pres_tbl[["time_cont"]] <- time_num_try
        .ph_log_info(
          "Time handling: continuous detected",
          step = "time handling",
          bullets = c(
            sprintf("time_col: %s", add_quotes(time_col, 1L)),
            sprintf(
              "numeric-conversion success: %.1f%% (threshold: %.0f%%)",
              100 * frac_numeric, 100 * thr
            ),
            if (isTRUE(time_force_continuous)) {
              "forced: TRUE"
            } else {
              "forced: FALSE"
            },
            "method_pcoa == 'cap' >> using numeric time (column: time_cont).",
            "Note: time pairwise contrasts & dispersion-by-time are disabled
            (continuous)."
          )
        )
      } else {
        .ph_warn(
          headline = "Forced/auto continuous time requested, but conversion
          has NAs.",
          step = "time handling",
          bullets = c(
            sprintf(
              "numeric-conversion success: %.1f%% (threshold: %.0f%%)",
              100 * frac_numeric, 100 * thr
            ),
            "Falling back to categorical time."
          )
        )
      }
    }

    # Categorical fallback (non-CAP, or failed numeric, or not forced)
    if (!time_is_continuous) {
      time_levels <- pres_tbl |>
        dplyr::distinct(!!rlang::sym(time_col)) |>
        dplyr::arrange(!!rlang::sym(time_col)) |>
        dplyr::pull(1) |>
        as.character()

      n_tlev <- sum(!is.na(time_levels))

      .ph_log_info(
        if (!identical(method_pcoa, "cap") && isTRUE(time_force_continuous)) {
          "Time handling: numeric input forced but non-CAP mode >>
          coerced to factor"
        } else if (cont_candidate && !identical(method_pcoa, "cap")) {
          "Time handling: numeric input but non-CAP mode >> coerced to factor"
        } else {
          "Time handling: categorical detected"
        },
        step = "time handling",
        bullets = c(
          sprintf("time_col: %s", add_quotes(time_col, 1L)),
          sprintf(
            "numeric-conversion success: %.1f%% (threshold: %.0f%%)",
            100 * frac_numeric, 100 * thr
          ),
          sprintf("n_levels: %d%s", n_tlev, if (n_tlev > 10L) {
            " (many levels; plots/tests may be slow)"
          } else {
            ""
          })
        )
      )
    }
  }


  # ---------- ranks -----------------------------------------------------------
  # as rank we understand the level of peptide metadata we want to aggregate -->
  # they can be raw peptides or species or another taxa levels
  ranks <- unique(ranks)
  if (!length(ranks) || !all(vapply(ranks, is.character, logical(1)))) {
    .ph_abort("`ranks` must be a non-empty character vector of exact
              column names.")
  }

  # ---------- helper to compute one RANK on a given subset --------------------
  # a pretty long one...
  # by the ranks in the previous version of the pipeline, we meant mainly the
  # peptide metadata/taxonomy; it is a convenient wrapping point for everything,
  # cause we can do all the stuff exactly the same but with different ranks;
  # so i built a wrapper which is actually the main workhorse of the beta
  # computing and laaplied it later to all the ranks defined by the user in the
  # call
  run_one_rank <- function(rank_name,
                           tbl_subset,
                           sublabel) {
    # 1) map rank --------------------------------------------------------------
    if (identical(rank_name, "peptide_id")) {
      ranked <- tbl_subset |>
        dplyr::transmute(
          sample_id = .data$sample_id,
          group = .data$group,
          time_fac = if (has_time && !time_is_continuous) {
            .data[[time_col]]
          } else {
            NA
          },
          time_cont = if (has_time && time_is_continuous) {
            .data[["time_cont"]]
          } else {
            NA_real_
          },
          rank_val = .data$peptide_id,
          value = if (used_exist_rule) {
            1L
          } else {
            .data$fold_change
          }
        )
    } else {
      # map provider is a convenient safe wrapper, to check if the rank defined
      # by the user is in the library; if not return NULL
      map_tbl <- map_provider(rank_name)
      if (is.null(map_tbl)) {
        return(NULL)
      } # end exec

      # build clean data for the subsequent analyses with custom rank; join on
      # peptides to take only those peptides, which actually ARE in the data,
      # not all peptides
      ranked <- tbl_subset |>
        dplyr::inner_join(map_tbl, by = "peptide_id") |>
        dplyr::filter(!is.na(.data$rank_val)) |>
        dplyr::transmute(
          sample_id = .data$sample_id,
          group = .data$group,
          time_fac = if (has_time && !time_is_continuous) {
            .data[[time_col]]
          } else {
            NA
          },
          time_cont = if (has_time && time_is_continuous) {
            .data[["time_cont"]]
          } else {
            NA_real_
          },
          rank_val = .data$rank_val,
          value = if (used_exist_rule) {
            1L
          } else {
            .data$fold_change
          }
        )
    }

    # optional rank filter --> filter only those ranks which are in the
    # filter_rank; eg when the rank is peptide_id and the filter_rank is
    # c("agilent12345", "twist12345"), only those two peptides will be taken
    # into account in the subsequent analyses
    if (is.function(filter_rank)) {
      # user can pass custom function
      allowed <- ranked |>
        dplyr::distinct(.data$rank_val) |>
        dplyr::collect() |>
        dplyr::pull()
      ok <- try(as.logical(filter_rank(allowed)), silent = TRUE)

      # filter out the not-allowed
      if (!inherits(ok, "try-error") && length(ok) == length(allowed)) {
        ranked <- ranked |>
          dplyr::filter(.data$rank_val %in% !!allowed[ok])
      }
    } else if (!is.null(filter_rank)) {
      # is a character vector
      ranked <- ranked |> dplyr::filter(.data$rank_val %in% !!filter_rank)
    }
    # in the future maybe delete entirely the function branch of the
    # if-statement

    # 2) aggregate -------------------------------------------------------------
    # we calculate the beta diversity per-rank, so we have to aggregate the data
    # per-rank; if the rank is peptide_id, this actually doesn't change anything
    agg <- ranked |>
      dplyr::group_by(sample_id, group, time_fac, rank_val) |>
      dplyr::summarise(abund = sum(.data$value), .groups = "drop") |>
      dplyr::collect()

    if (!nrow(agg)) {
      return(NULL)
    }

    # 3) metadata (subject_id/carry from tbl_subset; time_fac from ranked) -----
    # keep the carry_cols in the analysis (maybe the user would like to plot
    # it after the analysis?)
    meta_base <- tbl_subset |>
      dplyr::distinct(
        sample_id, group,
        dplyr::across(dplyr::all_of(intersect(
          c("subject_id", carry_cols),
          colnames(tbl_subset)
        )))
      ) |>
      dplyr::collect()

    time_map <- ranked |>
      dplyr::distinct(sample_id, time_fac, time_cont) |>
      dplyr::collect()

    meta <- dplyr::left_join(meta_base, time_map, by = "sample_id")
    meta$sample_id <- as.character(meta$sample_id)
    if (has_time && !time_is_continuous && "time_fac" %in% names(meta)) {
      meta$time_fac <- factor(as.character(meta$time_fac), levels = time_levels)
    }
    if (has_time && time_is_continuous && "time_cont" %in% names(meta)) {
      meta$time_cont <- as.numeric(meta$time_cont)
    }
    meta$.permutations <- permutations

    # 4) wide + align --> using previous helper (necessary for the pcoas) ------
    mat <- .long_to_wide_matrix(agg)
    meta <- meta[match(rownames(mat), meta$sample_id), , drop = FALSE]
    meta_df <- as.data.frame(meta, stringsAsFactors = FALSE)

    # drop all-zero rows --> now as the data is wide, we can drop all rows,
    # which are 0; technically it doesn't happen at all, or at least shouldn't
    # happen, cause this would mean, that a pearson has 0 enriched peptides;
    # nevertheless i keep it as fallback, these rows do not contribute anything
    # to the results of the analysis and can be safely dropped
    keep <- rowSums(mat) > 0
    if (!all(keep)) {
      mat <- mat[keep, , drop = FALSE]
      meta_df <- meta_df[keep, , drop = FALSE]
    }

    # if only one subject/sample is present in the wide matrix, there is nothing
    # to compute the distances from
    if (nrow(mat) < 2) {
      .ph_warn("Not enough rows after filtering to compute distances.",
        step = sprintf("beta-div (%s)", view_name)
      )
      return(NULL)
    }

    # 5) distance --------------------------------------------------------------
    # here i use the fast normalizer and calculate the distance using multi-
    # threaded engines where possible; from my understanding, the Carlos's
    # previous script used only the vegan pacakge which is single-threaded and
    # can slow down the computations a lot (several times)
    mat_norm_full <- .beta_normalize(
      mat,
      method_normalization,
      used_exist_rule
    )
    dist_full <- .ph_with_timing(
      headline = sprintf("Distance (%s)", distance),
      step = sprintf("beta-div (%s)", view_name),
      expr = {
        t0 <- proc.time()
        out <- .beta_dist(mat_norm_full, distance, n_threads)
        el <- (proc.time() - t0)[["elapsed"]]

        out
      },
      verbose = .ph_opt("verbose", TRUE)
    )

    # 6) ordination: PCoA or CAP -----------------------------------------------
    is_cap <- identical(method_pcoa, "cap")
    if (is_cap && !rlang::is_installed("vegan")) {
      .ph_abort("CAP mode requires 'vegan'. Please install it.")
    }

    # for the summary table i show only k_use axes (e.g. 10), but we will
    # normalize % by ALL positive eigenvalues --> REALLY IMPORTANT
    k_eig_req <- .ph_opt("phiper.beta.eig_axes", 10)
    n_eff <- max(1L, nrow(mat_norm_full) - 1L) # max axes available
    k_use <- max(2L, min(as.integer(k_eig_req), n_eff))
    k_points_full <- n_eff # for pcoa_full we request ALL axes

    # phiper-style log for the sizes
    .ph_log_info(
      "PCoA: choosing number of axes",
      step = sprintf("beta-div (%s)", view_name),
      bullets = c(
        sprintf("rank: %s", rank_name),
        sprintf("subset: %s", sublabel),
        sprintf("n_samples: %d", nrow(mat_norm_full)),
        sprintf("requested eig_axes option: %d", as.integer(k_eig_req)),
        sprintf("n-1 limit: %d", n_eff),
        sprintf("k used: %d%s", k_use, if (k_use < as.integer(k_eig_req)) {
          " (clipped by n-1 / min k=2)"
        } else {
          ""
        })
      )
    )

    if (!is_cap) {
      # ----- PCoA path --------------------------------------------------------
      ## no correction for the negative eigenvalues
      # ----- PCoA path --------------------------------------------------------
      fit <- .ph_with_timing(
        headline = sprintf(
          "PCoA (%s)",
          if (identical(neg_correction, "none")) "cmdscale" else paste0("wcmdscale + ", neg_correction)
        ),
        step = sprintf("beta-div (%s)", view_name),
        expr = {
          t0 <- proc.time()
          out <- if (identical(neg_correction, "none")) {
            # fast base cmdscale
            stats::cmdscale(dist_full, eig = TRUE, k = k_points_full)
          } else {
            # vegan with negative-eigen correction
            if (!rlang::is_installed("vegan")) {
              .ph_abort("Negative-eigen correction requires 'vegan'.")
            }
            vegan::wcmdscale(dist_full, eig = TRUE, k = k_points_full, add = neg_correction)
          }
          el <- (proc.time() - t0)[["elapsed"]]
          pts_nc <- tryCatch(ncol(as.matrix(out$points)), error = function(e) NA_integer_)
          ev <- as.numeric(out$eig %||% numeric())
          n_pos <- sum(ev > 0, na.rm = TRUE)
          n_neg <- sum(ev < 0, na.rm = TRUE)

          out
        },
        verbose = .ph_opt("verbose", TRUE)
      )

      # get the points/sample coordinates for ALL pcoa axes
      pts <- as.matrix(fit$points)

      # name rows and axes if not named
      if (is.null(rownames(pts))) {
        rownames(pts) <- rownames(as.matrix(dist_full))
      }
      colnames(pts) <- paste0("PCoA", seq_len(ncol(pts)))

      # --- eigenvalues (full + summary) ---------------------------------------
      # extract the data from the fit
      eig_full_vec <- as.numeric(fit$eig %||% numeric())

      # robust 0-row/2-col tibble even if no eigs (very rare, should not happen)
      eig_full <- tibble::tibble(
        axis = paste0(
          "PCoA",
          seq_len(length(eig_full_vec))
        ),
        eigenvalue = as.numeric(eig_full_vec)
      )

      # positive/negative bookkeeping --> for the count/energy ratios later and
      # percentages rn
      n_pos_all <- sum(eig_full_vec > 0, na.rm = TRUE)
      n_neg_all <- sum(eig_full_vec < 0, na.rm = TRUE)
      sum_pos_all <- sum(pmax(eig_full_vec, 0), na.rm = TRUE)
      sum_abs_neg <- sum(abs(eig_full_vec[eig_full_vec < 0]), na.rm = TRUE)

      # summarizing the "other" into a common tail of the eigen_summary table
      idx_head <- seq_len(k_use)
      idx_tail <- if (length(eig_full_vec) > k_use) {
        (k_use + 1L):length(eig_full_vec)
      } else {
        integer(0)
      }

      # vectorized head percentages (normal and cumulative)
      eig_head <- eig_full_vec[idx_head]
      pct_head <- if (sum_pos_all > 0) {
        100 * pmax(eig_head, 0) / sum_pos_all
      } else {
        rep(NA_real_, length(idx_head))
      }

      cum_head <- if (sum_pos_all > 0) {
        cumsum(pct_head)
      } else {
        rep(NA_real_, length(idx_head))
      }

      # "other" = raw sum of all remaining eigenvalues (pos + neg),
      # and % relative to the sum of positive eigenvalues (like the head rows)
      eig_tail_sum_raw <- if (length(idx_tail)) {
        sum(eig_full_vec[idx_tail])
      } else {
        0
      }
      pct_tail <- if (sum_pos_all > 0) {
        100 * sum(pmax(eig_full_vec[idx_tail], 0)) / sum_pos_all
      } else {
        NA_real_
      }

      # combining everything into a clean summary output table
      eig_tbl <- tibble::tibble(
        axis            = c(paste0("PCoA", idx_head), "Other"),
        eigenvalue      = c(eig_head, eig_tail_sum_raw),
        pct_of_pos      = c(pct_head, pct_tail),
        cum_pct_of_pos  = c(cum_head, if (!is.na(pct_tail)) 100 else NA_real_),
        n_pos           = n_pos_all,
        n_neg           = n_neg_all,
        rank            = rank_name,
        view            = view_name
      )

      # --- var_explained: %PCoA1..%PCoA{k_use} + %Other -----------------------
      var_tbl <- tibble::as_tibble_row(
        c(stats::setNames(as.list(pct_head), paste0("%PCoA", idx_head)),
          `%Other` = pct_tail
        )
      ) |>
        dplyr::mutate(rank = rank_name, view = view_name)

      # --- spectrum diagnostics & phiper-style logging/warning ----------------
      # ratio pos/neg + energy
      ratio_counts <- if (n_pos_all > 0) {
        n_neg_all / n_pos_all
      } else {
        NA_real_
      }
      ratio_energy <- if (sum_pos_all > 0) {
        sum_abs_neg / sum_pos_all
      } else {
        NA_real_
      }

      # thresholds (user-overridable) --> default 5% of axes negative and 5%
      # of energy negative
      thr_counts <- .ph_opt("phiper.beta.neg_ratio_warn_counts", 0.05)
      thr_energy <- .ph_opt("phiper.beta.neg_energy_warn", 0.05)

      # Info log
      .ph_log_info(
        "PCoA spectrum diagnostics",
        step = sprintf("beta-div (%s)", view_name),
        bullets = c(
          sprintf("rank: %s", rank_name),
          sprintf(
            "n_pos: %d; n_neg: %d; neg/pos (count ratio): %s",
            n_pos_all, n_neg_all, if (is.na(ratio_counts)) {
              "NA"
            } else {
              sprintf("%.3f", ratio_counts)
            }
          ),
          sprintf(
            "sum_pos: %.6g; sum|neg|: %.6g; |neg|/pos (energy ratio): %s",
            sum_pos_all, sum_abs_neg, if (is.na(ratio_energy)) {
              "NA"
            } else {
              sprintf("%.3f", ratio_energy)
            }
          ),
          sprintf("neg_eigen_correction: %s", neg_correction)
        )
      )

      # warn the user if noticeable non-Euclidean signal
      if (isTRUE(!is.na(ratio_counts) && ratio_counts > thr_counts) ||
        isTRUE(!is.na(ratio_energy) && ratio_energy > thr_energy)) {
        .ph_warn(
          headline = "Noticeable non-Euclidean signal in distances.",
          step = sprintf("beta-div (%s)", view_name),
          bullets = c(
            sprintf(
              "neg/pos (count) = %s  [threshold %.3f]",
              if (is.na(ratio_counts)) {
                "NA"
              } else {
                sprintf("%.3f", ratio_counts)
              }, thr_counts
            ),
            sprintf(
              "|neg|/pos (energy) = %s  [threshold %.3f]",
              if (is.na(ratio_energy)) {
                "NA"
              } else {
                sprintf("%.3f", ratio_energy)
              }, thr_energy
            ),
            "Consider applying a Cailliez or Lingoes correction if you need a
            strictly Euclidean embedding.",
            "Percentages are reported relative to the sum of positive
            eigenvalues."
          )
        )
      }

      # --- principal coordinate scores (full + summary) -----------------------
      # pcoa_full consists of ALL axes --> for diagnostic purposes only, but
      # nice to have in the output
      pcoa_full <- tibble::as_tibble(pts, rownames = "sample_id") |>
        dplyr::left_join(meta_df, by = "sample_id") |>
        dplyr::mutate(rank = rank_name, view = view_name)

      # pcoa_tbl summarizes the first k_use axes (usually the most important),
      # for convenient plotting (pad with zeros if needed); it also attaches
      # the metadata for the plotting/subsequent analyses
      ptsk <- if (ncol(pts) < k_use) {
        cbind(pts, matrix(0, nrow(pts), k_use - ncol(pts)))
      } else {
        pts[, 1:k_use, drop = FALSE]
      }

      colnames(ptsk) <- paste0("PCoA", seq_len(ncol(ptsk)))

      # nice summamary table for plotting/analyses
      pcoa_tbl <- tibble::as_tibble(ptsk, rownames = "sample_id") |>
        dplyr::left_join(meta_df, by = "sample_id") |>
        dplyr::mutate(rank = rank_name, view = view_name)

      # top-N feature/peptide loadings onto PCoA axes (approx. biplot via WA)
      load_tbl <- .feature_loadings_wa(
        U = pts[, seq_len(min(10, ncol(pts))), drop = FALSE],
        X = mat_norm_full,
        feature_names = colnames(mat_norm_full),
        top = .ph_opt("phiper.beta.loadings_top", 100)
      ) |>
        dplyr::mutate(rank = rank_name, view = view_name)

      ord_out <- list(
        pcoa = pcoa_tbl,
        pcoa_full = pcoa_full,
        var_explained = var_tbl,
        eigen_summary = eig_tbl,
        eig_full = eig_full,
        feature_loadings = load_tbl,
        pcoa_fit = fit
      )
    } else {
      # ----- CAP/dbRDA path ---------------------------------------------------
      # dbRDA - distance based redundancy analysis
      # the capscale/dbRDA is a method, which encompases all the previous
      # previous scenarios covered with the 3 pcoa methods - thus i will treat
      # them in the future as legacy/backward-compatibility options for the
      # previous analyses, but not develop them further, as there is no point
      # in further granulation/testing of multiple, smaller and smaller sub-
      # groups, when we can do all at once using capscale
      #
      #
      # so the capscale is a redundancy analysis, it starts with raw distances,
      # just like the pcoa; firstly it actually performs the PCoA to place the
      # samples on an euclidean space, it is also not a problem when the eigen-
      # values are negative --> we can easily apply a correction to deal with
      # them
      #
      # then it does a linear regression/ordination (RDA) of the PCoA
      # coordinates against the predictors, which the user defined and finds the
      # axes, that are best linear combination of the predictors
      #
      # the output is divided into a "constrained" and "unconstrained" part; the
      # constrained one are our predictors --> how much variance overall do they
      # explain??? the unconstrained part is everything else -> the variance not
      # explained by the predictors; we can calculate R2 on this
      #
      # then we can partition the "constrained" explained variance into multiple
      # principal components; it is a lot better than the regular PCoA, because
      # you can actually meaningfully interpret the components, instead of
      # doing the pcoa, hoping that there are any differences and dont interpret
      # them, cause the axes are not meaningful
      #
      # eg in the babies study, i performed the capscale on time and group
      # (mom_serum, kid_serum and mom_milk) and the first two components
      # explained ~90% variance and corresponded cleanly to the person (mom vs.
      # kid) + time. I assume the third component could correspond to the
      # material type (milk vs serum) but i am not sure --> to be tested
      #
      # The big win in this is: the time CAN BE CONTINUOUS --> it was completely
      # not doable with the old framework
      # ------------------------------------------------------------------------

      # create the rhs of the formula based on the user-input (within and
      # between-subject effects)
      rhs <- NULL
      if (!is.null(group_col)) rhs <- "group"
      if (has_time) {
        if (time_is_continuous) {
          rhs <- if (is.null(rhs)) {
            "time_cont"
          } else {
            paste(rhs, "+ time_cont + group:time_cont")
          }
        } else {
          rhs <- if (is.null(rhs)) {
            "time_fac"
          } else {
            paste(rhs, "+ time_fac + group:time_fac")
          }
        }
      }

      # if after that the rhs is still null, this means that no grouping factor
      # or time_col were provided --> there is no point in doing CAP when no
      # grouping
      if (is.null(rhs)) {
        .ph_abort(
          headline = "CAP mode requires a constraining factor.",
          step = sprintf("beta-div (%s) CAP", view_name),
          bullets = c(
            "No 'group' or 'time' available to constrain the ordination.",
            "Provide 'group_col' and/or 'time_col', or use method_pcoa =
            'joint'."
          )
        )
      }

      # process the negative eigenvalues correction
      add_arg <- if (identical(neg_correction, "none")) {
        NULL
      } else {
        neg_correction
      }

      # prepare the data.frame for the analysis --> only necessary data
      design_df <- data.frame(
        row.names = meta_df$sample_id,
        stringsAsFactors = FALSE
      )
      if (!is.null(group_col)) design_df$group <- meta_df$group

      # depending on the time being continuous/factor we have different naming
      # conventions
      if (has_time && time_is_continuous) {
        design_df$time_cont <- meta_df$time_cont
      }
      if (has_time && !time_is_continuous) {
        design_df$time_fac <- meta_df$time_fac
      }

      # pasting the formula together
      form <- stats::as.formula(paste("dist_full ~", rhs), env = environment())

      # the workhorse of the whole section --> performs the capscale;
      # maybe think about optimization in the future, as the capscale is single-
      # threaded and a little bit slow on larger dat
      cap <- .ph_with_timing(
        headline = "CAP/dbRDA fit (vegan::capscale)",
        step = sprintf("beta-div (%s) CAP", view_name),
        expr = {
          t0 <- proc.time()
          fit <- vegan::capscale(form, data = design_df, add = add_arg)
          el <- (proc.time() - t0)[["elapsed"]]

          # Quick diagnostics for the log
          eig_con <- tryCatch(fit$CCA$eig, error = function(e) numeric())
          n_axes <- length(eig_con)
          n_pos <- sum(eig_con > 0, na.rm = TRUE)

          tot_iner <- tryCatch(fit$tot.chi, error = function(e) NA_real_)
          con_iner <- tryCatch(fit$CCA$tot.chi, error = function(e) NA_real_)
          unc_iner <- tryCatch(fit$CA$tot.chi, error = function(e) NA_real_)
          prop_con <- if (is.finite(tot_iner) && tot_iner > 0) con_iner / tot_iner else NA_real_

          fit
        },
        verbose = .ph_opt("verbose", TRUE)
      )

      # ---- processing the principal component scores -------------------------
      # extract the scores and convert to a matrix
      scr <- vegan::scores(cap, display = "sites")
      pts <- as.matrix(scr)

      # safe fallback and colnames/rownames handling
      if (is.null(pts) || !nrow(pts)) {
        pts <- matrix(0, nrow = nrow(meta_df), ncol = 2)
      }
      colnames(pts) <- paste0("CAP", seq_len(ncol(pts)))
      rownames(pts) <- meta_df$sample_id

      # the logic is actuall similar/almost the same as with the normal PCoA, so
      # i will refrain from commenting it
      eig_con <- cap$CCA$eig %||% numeric()
      sum_pos_all <- sum(pmax(eig_con, 0))

      k_use <- min(length(eig_con), .ph_opt("phiper.beta.eig_axes", 10))
      head_eig <- eig_con[seq_len(k_use)]

      pct_head <- if (sum_pos_all > 0) {
        100 * pmax(head_eig, 0) / sum_pos_all
      } else {
        rep(NA_real_, length(head_eig))
      }

      cum_head <- if (sum_pos_all > 0) {
        cumsum(pct_head)
      } else {
        rep(NA_real_, length(head_eig))
      }

      tail_idx <- if (length(eig_con) > k_use) {
        (k_use + 1L):length(eig_con)
      } else {
        integer(0)
      }

      tail_sum <- if (length(tail_idx)) sum(eig_con[tail_idx]) else 0
      pct_tail <- if (sum_pos_all > 0) {
        100 * sum(pmax(eig_con[tail_idx], 0)) / sum_pos_all
      } else {
        NA_real_
      }

      eig_full <- tibble::tibble(
        axis = paste0("CAP", seq_len(length(eig_con))),
        eigenvalue = as.numeric(eig_con)
      )

      eig_tbl <- tibble::tibble(
        axis           = c(paste0("CAP", seq_len(k_use)), "Other"),
        eigenvalue     = c(head_eig, tail_sum),
        pct_of_pos     = c(pct_head, pct_tail),
        cum_pct_of_pos = c(cum_head, if (!is.na(pct_tail)) 100 else NA_real_),
        n_pos          = sum(eig_con > 0),
        n_neg          = sum(eig_con < 0),
        rank           = rank_name,
        view           = view_name
      )

      # ---- processing the variance of principal components -------------------
      var_tbl <- tibble::as_tibble_row(
        c(stats::setNames(as.list(pct_head), paste0("%CAP", seq_len(k_use))),
          `%Other` = pct_tail
        )
      ) |>
        dplyr::mutate(rank = rank_name, view = view_name)

      # ---- processing the eigenvalues ----------------------------------------
      pcoa_full <- tibble::as_tibble(pts, rownames = "sample_id") |>
        dplyr::left_join(meta_df, by = "sample_id") |>
        dplyr::mutate(rank = rank_name, view = view_name)
      ptsk <- pts[, seq_len(min(ncol(pts), max(2L, k_use))), drop = FALSE]
      pcoa_tbl <- tibble::as_tibble(ptsk, rownames = "sample_id") |>
        dplyr::left_join(meta_df, by = "sample_id") |>
        dplyr::mutate(rank = rank_name, view = view_name)

      # ---- extracting the loadings  ------------------------------------------
      load_tbl <- .feature_loadings_wa(
        U = pts[, seq_len(min(10, ncol(pts))), drop = FALSE],
        X = mat_norm_full,
        feature_names = colnames(mat_norm_full),
        top = .ph_opt("phiper.beta.loadings_top", 100)
      ) |>
        dplyr::mutate(rank = rank_name, view = view_name)

      # --- inertia partitioning (Total / Constrained / Unconstrained) ---------
      # this section is unique for the capscale; it processes the table with
      # constrained/unconstrained/total variance and sums of positve eigens, so
      # we can add it to the final output; it is important for converting the
      # scaled percentages to the % of total variance
      tot <- cap$tot.chi %||% # total variance
        (sum(cap$CCA$eig %||% 0) + sum(cap$CA$eig %||% 0))
      con <- cap$CCA$tot.chi %||% sum(cap$CCA$eig %||% 0) # constrained
      uncon <- cap$CA$tot.chi %||% sum(cap$CA$eig %||% 0) # unconstrained

      # combining them into a tibble
      part_tbl <- tibble::tibble(
        component = c("Total", "Constrained", "Unconstrained"),
        inertia = c(tot, con, uncon),
        proportion = c(
          1, ifelse(tot > 0, con / tot, NA_real_),
          ifelse(tot > 0, uncon / tot, NA_real_)
        ),
        rank = rank_name, view = view_name
      )

      # --- the final capscale output ------------------------------------------
      ord_out <- list(
        pcoa = pcoa_tbl, pcoa_full = pcoa_full,
        var_explained = var_tbl, eigen_summary = eig_tbl,
        eig_full = eig_full, feature_loadings = load_tbl,
        cap_fit = cap,
        cap_partitioning = part_tbl
      )
    }

    # 7) PERMANOVA: ALWAYS global model (+ optional contrasts)
    ## placeholder for the results
    tests_out <- list()

    # add results to the tests_out; superassignment in-place
    add_tests <- function(tt) {
      if (is.null(tt) || !nrow(tt)) {
        return()
      }
      tt$rank <- rank_name
      tt$view <- view_name
      tests_out[[length(tests_out) + 1]] <<- tt
    }

    # the tests depend on the time if it is specified --> we perform another set
    # for the continuous time and factor time actually supports the usual
    # options ("pairwise", "each_vs_rest" etc.)
    time_var <- if (has_time && time_is_continuous) {
      "time_cont"
    } else if (has_time) {
      "time_fac"
    } else {
      NULL
    }

    # building the formula for the tests; analogical to the upper one for
    # capscale BUT necessary when the path is simple pcoa; DO NOT DELETE
    rhs_global <- NULL

    # defining the grouping terms
    if (length(unique(meta_df$group)) > 1) {
      rhs_global <- "group"
    }
    if (!is.null(time_var) && length(unique(meta_df[[time_var]])) > 1) {
      rhs_global <- if (is.null(rhs_global)) {
        time_var
      } else {
        paste("group +", time_var, "+ group:", time_var)
      }
    }

    # ---- global test ---------------------------------------------------------
    if (!is.null(rhs_global)) {
      # previously defined helper
      tt <- .adonis_terms(
        dist_full,
        meta_df,
        rhs_global,
        permutations,
        step_label = sprintf("global_PERMANOVA (%s)", view_name),
        strata = if (has_time && ("subject_id" %in% names(meta_df))) {
          meta_df$subject_id
        } else {
          NULL
        },
        n_threads = n_threads
      )

      ## add to the tt list storing the tests
      if (!is.null(tt) && nrow(tt)) {
        tt$scope <- "global_2way"
        tt$contrast <- "<global>"
        add_tests(tt)
      }
    } else {
      # when no test can be performed --> log the info
      .ph_log_info("Global PERMANOVA skipped (insufficient number of factors).",
        step = sprintf("beta-div (%s)", view_name)
      )
    }

    # ---- optional contrasts --------------------------------------------------
    # argument processing
    contrasts <- unique(tolower(contrasts))
    contrasts[contrasts == "group_vs_rest"] <- "each_vs_rest"

    # if the method for contrasts is defined, perform or at least try to perform
    # them
    if (!("none" %in% contrasts)) {
      # -------- small helper to tidy an adonis2 row --------
      .tidy_row <- function(fit, term, permutations, scope, contrast) {
        if (inherits(fit, "try-error")) {
          return(NULL)
        }
        df <- as.data.frame(fit)
        df$term <- rownames(df)
        df <- df[df$term == term, c("term", "Pr(>F)", "F", "R2")]
        if (!nrow(df)) {
          return(NULL)
        }
        colnames(df) <- c("term", "p_value", "F_stat", "R2")
        df$n_perm <- permutations
        df$p_at_floor <- !is.na(df$p_value) &&
          df$p_value <= 1 / (permutations + 1)
        df$scope <- scope
        df$contrast <- contrast
        tibble::as_tibble(df)
      }

      # -------- ALWAYS: pairwise GROUP comparisons controlling for time -------
      # (a) time categorical -> test group while adjusting for time_fac
      run_group_pairwise_cat <- function() {
        # generate the pairs for comparisons
        prs <- .pairs_of(levels(factor(meta_df$group)))
        if (!length(prs)) {
          return(NULL)
        } # fallback

        # placaeholder list
        out <- list()
        for (p in prs) {
          # keep only samples from the pair
          sel <- which(meta_df$group %in% p)

          # we need at least 2 samples per group to run PERMANOVA
          if (length(sel) < 4) next

          # subset the distance matrix to select only the pairs
          d_sub <- stats::as.dist(as.matrix(dist_full)[sel, sel])
          meta_s <- droplevels(meta_df[sel, , drop = FALSE])
          rownames(meta_s) <- meta_s$sample_id

          # 2-level factor of the two groups (keep order of 'p')
          meta_s$fac <- droplevels(factor(as.character(meta_s$group),
            levels = p
          ))
          if (nlevels(meta_s$fac) < 2 || min(table(meta_s$fac)) < 2) next

          # decide strata: only if the tested factor (fac) varies within
          # subjects

          # if the tested factor (fac) varies within a subject (rare but
          # possible), we stratify permutations by subject to respect
          # pairing/repeated measures
          #
          # if every subject belongs to only one of the groups (usual
          # between-subject case) we dont stratify (and log that choice)

          strata_use <- NULL
          if ("subject_id" %in% names(meta_s)) {
            si <- meta_s$subject_id
            has_within <- any(tapply(
              meta_s$fac,
              si,
              function(v) length(unique(v)) > 1
            ))

            if (has_within) {
              strata_use <- si
            } else {
              .ph_log_info("No within-subject var for this group contrast;
                           using unstratified perms.",
                step = "pairwise PERMANOVA (group|time,
                           categorical)"
              )
            }
          }

          # group main effect controlling for time_fac
          # (Type III via by='margin')
          fit <- vegan::adonis2(d_sub ~ fac + time_fac,
            data = meta_s,
            by = "margin",
            permutations = permutations,
            strata = strata_use,
            n_threads = n_threads
          )

          # extract the results
          tt <- .tidy_row(fit, "fac", permutations,
            scope = "group_pairwise|time_cat",
            contrast = paste(p, collapse = " vs ")
          )

          # append them to the output
          if (!is.null(tt)) out[[length(out) + 1]] <- tt
        }
        if (!length(out)) NULL else dplyr::bind_rows(out)
      }

      # (b) time continuous -> two tests per pair: group|time and group:time
      # in one we test the group conditional on time, in the other we test the
      # time trends between the groups
      run_group_pairwise_ct <- function() {
        if (!("time_cont" %in% names(meta_df))) {
          return(NULL)
        }

        # similar logic as in the a)
        prs <- .pairs_of(levels(factor(meta_df$group)))
        if (!length(prs)) {
          return(NULL)
        }

        # output list
        out <- list()
        for (p in prs) {
          # select the groups; at least 3 observations (cause we are testing
          # trends!!!)
          sel <- which(meta_df$group %in% p)
          if (length(sel) < 6) next

          # subset the data to extract only the desired pairs
          d_sub <- stats::as.dist(as.matrix(dist_full)[sel, sel])
          meta_s <- droplevels(meta_df[sel, , drop = FALSE])
          rownames(meta_s) <- meta_s$sample_id

          # 2-level factor of the two groups
          meta_s$fac <- droplevels(factor(as.character(meta_s$group),
            levels = p
          ))
          if (nlevels(meta_s$fac) < 2 || min(table(meta_s$fac)) < 2) next

          # center time within group (helps with interaction interpretation)
          meta_s$time_c <- ave(meta_s$time_cont, meta_s$fac,
            FUN = function(z) z - mean(z, na.rm = TRUE)
          )

          # enough distinct time points?
          if (length(unique(meta_s$time_c)) < 3) {
            .ph_log_info("Too few distinct time points for slope test;
                         skipping group:time_cont.",
              step = "pairwise PERMANOVA (group*time, continuous)"
            )
          }

          # strata logic; same as up
          strata_use <- NULL
          if ("subject_id" %in% names(meta_s)) {
            si <- meta_s$subject_id
            has_within_group <- any(tapply(
              meta_s$fac,
              si,
              function(v) length(unique(v)) > 1
            ))
            has_within_time <- any(tapply(
              meta_s$time_c,
              si,
              function(v) length(unique(v)) > 1
            ))
            if (has_within_group || has_within_time) {
              strata_use <- si
            } else {
              .ph_log_info("No within-subject var for this pair;
                           using unstratified perms.",
                step = "pairwise PERMANOVA (group & group:time,
                           continuous)"
              )
            }
          }

          # (1) group main effect | time
          fit1 <- vegan::adonis2(d_sub ~ fac + time_c,
            data = meta_s, by = "margin",
            permutations = permutations,
            strata = strata_use,
            n_threads = n_threads
          )

          tt1 <- .tidy_row(fit1, "fac", permutations,
            scope = "group_pairwise|time_cont",
            contrast = paste(p, collapse = " vs ")
          )
          if (!is.null(tt1)) {
            tt1$term <- "group|time"
            out[[length(out) + 1]] <- tt1
          }

          # (2) difference in trends
          if (length(unique(meta_s$time_c)) >= 3) {
            fit2 <- vegan::adonis2(d_sub ~ fac * time_c,
              data = meta_s, by = "margin",
              permutations = permutations,
              strata = strata_use,
              n_threads = n_threads
            )
            tt2 <- .tidy_row(fit2, "fac:time_c", permutations,
              scope = "group_pairwise|time_cont",
              contrast = paste(p, collapse = " vs ")
            )
            if (!is.null(tt2)) {
              tt2$term <- "group:time_cont"
              out[[length(out) + 1]] <- tt2
            }
          }
        }
        if (!length(out)) NULL else dplyr::bind_rows(out)
      }

      # c) simple pairwise --> the difference between a) and this is, that this
      # one can run also for categorical time, a) tests only the group_col
      # effects
      run_pairwise <- function(idx,
                               col_name,
                               scope_label,
                               ctx_cols) {
        # list all pairs of levels
        fac_here <- factor(as.character(meta_df[[col_name]][idx]))
        prs <- .pairs_of(levels(fac_here))
        if (!length(prs)) {
          return(NULL)
        }

        # empty placeholder; then loop over all defined pairs and perform
        # pairwise PERMANOVAs
        res <- list()
        for (p in prs) {
          # select samples of those two levels in p
          sel <- idx[which(fac_here %in% p)]

          # at least two samples per group necessary (2x2 = 4)
          if (length(sel) < 4) next
          meta_s <- meta_df[sel, , drop = FALSE]
          rownames(meta_s) <- meta_s$sample_id

          # build a 2-level factor fac and check balance
          fac <- droplevels(factor(as.character(meta_s[[col_name]]),
            levels = p
          ))
          if (nlevels(fac) < 2 || min(table(fac)) < 2) next

          # subset the distance matrix to those samples
          d_sub <- stats::as.dist(as.matrix(dist_full)[sel, sel])

          meta_s$fac <- fac

          # decide strata correctly; similar logic as in the helpers a) and b)
          strata_use <- NULL
          if ("subject_id" %in% names(meta_s)) {
            si <- meta_s$subject_id
            has_within <- any(tapply(
              meta_s$fac,
              si,
              function(v) length(unique(v)) > 1
            ))

            if (has_within) {
              strata_use <- si
            } else {
              .ph_log_info("No within-subject variation for this contrast;
                           using unstratified perms.",
                step = "pairwise PERMANOVA"
              )
            }
          }

          # workhorse
          tt <- .adonis_terms(d_sub,
            meta_s,
            "fac",
            permutations,
            scope_label,
            strata = strata_use,
            n_threads = n_threads
          )

          # append output
          if (!is.null(tt) && nrow(tt)) {
            tt <- dplyr::filter(tt, .data$term == "fac")
            tt$scope <- scope_label
            tt$contrast <- paste(p, collapse = " vs ")

            if (!is.null(ctx_cols$time_level)) {
              tt$time_level <- ctx_cols$time_level
            }

            if (!is.null(ctx_cols$group_level)) {
              tt$group_level <- ctx_cols$group_level
            }

            res[[length(res) + 1]] <- tt
          }
        }
        if (!length(res)) NULL else dplyr::bind_rows(res)
      }

      # d) simple each vs the rest --> when only two levels present it actually
      # is the same as pairwise; when more, then it takes each level and
      # compares it against all other levels combined as one --> pretty self-
      # explanatory; the interface is the same as in all previous helpers
      run_each_vs_rest <- function(idx,
                                   col_name,
                                   scope_label,
                                   ctx_cols) {
        # collect levels and build “each vs rest” factors
        vec_chr <- as.character(meta_df[[col_name]][idx])
        lv <- unique(vec_chr)
        lv <- lv[!is.na(lv)]
        if (length(lv) < 2) {
          return(NULL)
        }

        # empty placeholder for the results + build the pairs for comparisons
        res <- list()
        evr <- .each_vs_rest_factors(vec_chr, lv)
        for (i in seq_along(lv)) {
          # keep rows with a defined binary label
          keep <- !is.na(evr[[i]])
          sel <- idx[keep]
          if (length(sel) < 4) next

          # make a 2-level factor and check balance
          meta_s <- meta_df[sel, , drop = FALSE]
          rownames(meta_s) <- meta_s$sample_id
          fac <- droplevels(evr[[i]][keep])

          # at least two levels; otherwise nothin to compare
          if (nlevels(fac) < 2 || min(table(fac)) < 2) next

          # subset the distance matrix
          d_sub <- stats::as.dist(as.matrix(dist_full)[sel, sel])

          meta_s$fac <- fac

          # the same strata logic as usual
          strata_use <- NULL
          if ("subject_id" %in% names(meta_s)) {
            si <- meta_s$subject_id
            has_within <- any(tapply(
              meta_s$fac,
              si,
              function(v) length(unique(v)) > 1
            ))
            if (has_within) {
              strata_use <- si
            } else {
              .ph_log_info("No within-subject variation for this contrast;
                           using unstratified perms.",
                step = "each_vs_rest PERMANOVA"
              )
            }
          }

          # workhorse
          tt <- .adonis_terms(d_sub,
            meta_s,
            "fac",
            permutations,
            scope_label,
            strata = strata_use,
            n_threads = n_threads
          )

          # appending to the output
          if (!is.null(tt) && nrow(tt)) {
            tt <- dplyr::filter(tt, .data$term == "fac")
            tt$scope <- scope_label
            tt$contrast <- paste0(lv[i], " vs other")

            if (!is.null(ctx_cols$time_level)) {
              tt$time_level <- ctx_cols$time_level
            }

            if (!is.null(ctx_cols$group_level)) {
              tt$group_level <- ctx_cols$group_level
            }
            res[[length(res) + 1]] <- tt
          }
        }
        if (!length(res)) NULL else dplyr::bind_rows(res)
      }

      # e) each vs baseline helper --> comes in handy especially for time/
      # longitudinal analyses, where we want to picka a baseline level (usually
      # the earliest timepoint) and compare everything else to it; the interface
      # is the same as in the rest
      run_baseline <- function(idx,
                               col_name,
                               baseline,
                               scope_label,
                               ctx_cols) {
        # build a binary factor “baseline vs other”
        vec_chr <- as.character(meta_df[[col_name]][idx])
        fac <- .baseline_factor(vec_chr, baseline)
        if (is.null(fac)) {
          return(NULL)
        }

        # keep valid rows and basic size checks
        keep <- !is.na(fac)
        sel <- idx[keep]
        if (length(sel) < 4) {
          return(NULL)
        }

        # subset meta and distances to those rows
        meta_s <- meta_df[sel, , drop = FALSE]
        rownames(meta_s) <- meta_s$sample_id
        fac <- droplevels(fac[keep])
        if (nlevels(fac) < 2 || min(table(fac)) < 2) {
          return(NULL)
        }
        d_sub <- stats::as.dist(as.matrix(dist_full)[sel, sel])

        meta_s$fac <- fac

        # the same strata logic
        strata_use <- NULL
        if ("subject_id" %in% names(meta_s)) {
          si <- meta_s$subject_id
          has_within <- any(tapply(
            meta_s$fac,
            si,
            function(v) length(unique(v)) > 1
          ))
          if (has_within) {
            strata_use <- si
          } else {
            .ph_log_info("No within-subject variation for this contrast;
                         using unstratified perms.",
              step = "baseline PERMANOVA"
            )
          }
        }

        # workhorse
        tt <- .adonis_terms(d_sub,
          meta_s,
          "fac",
          permutations,
          scope_label,
          strata = strata_use,
          n_threads = n_threads
        )

        # append output
        if (!is.null(tt) && nrow(tt)) {
          tt <- dplyr::filter(tt, .data$term == "fac")
          tt$scope <- scope_label
          tt$contrast <- paste0(baseline, " vs other")
          if (!is.null(ctx_cols$time_level)) {
            tt$time_level <- ctx_cols$time_level
          }
          if (!is.null(ctx_cols$group_level)) {
            tt$group_level <- ctx_cols$group_level
          }
          return(tt)
        }
        NULL
      }

      # -------- dispatcher: ALWAYS do group comparisons -----------------------
      if (length(unique(meta_df$group)) > 1) {
        if (isTRUE(time_is_continuous)) {
          add_tests(run_group_pairwise_ct())
        } else {
          add_tests(run_group_pairwise_cat())
        }
      }

      # -------- time contrasts ONLY when categorical (unchanged) --------------
      # dispatch to the appropriate helper, depending on the method
      if (!time_is_continuous &&
        has_time &&
        length(unique(meta_df$time_fac)) > 1) {
        if ("pairwise" %in% contrasts) {
          add_tests(run_pairwise(
            seq_len(nrow(meta_df)),
            "time_fac",
            "time",
            list()
          ))
        } else if ("each_vs_rest" %in% contrasts) {
          add_tests(run_each_vs_rest(
            seq_len(nrow(meta_df)),
            "time_fac",
            "time",
            list()
          ))
        } else if ("baseline" %in% contrasts) {
          t_levels <- levels(factor(meta_df$time_fac))
          t_base <- (baseline_level %||% t_levels[1])
          add_tests(run_baseline(
            seq_len(nrow(meta_df)),
            "time_fac",
            t_base,
            "time",
            list()
          ))
        }
      }

      # -------- nested scopes using time ONLY when categorical ----------------
      if (!time_is_continuous &&
        has_time &&
        length(unique(meta_df$group)) > 1 &&
        length(unique(meta_df$time_fac)) > 1) {
        # compare the groups WITHIN each timepoint if possible (when there is
        # no overlap between the groups, the comparisons will be NULL --> no
        # error)
        for (t in levels(factor(meta_df$time_fac))) {
          idx_t <- which(meta_df$time_fac == t)
          g_here <- droplevels(factor(meta_df$group[idx_t]))
          if (nlevels(g_here) < 2) next

          # limited possibilities here --> maybe add baseline later; for now it
          # suffices
          if ("pairwise" %in% contrasts) {
            add_tests(run_pairwise(
              idx_t,
              "group",
              sprintf("group_within_time[%s]", t),
              list(time_level = t)
            ))
          } else if ("each_vs_rest" %in% contrasts) {
            add_tests(run_each_vs_rest(
              idx_t,
              "group",
              sprintf("group_within_time[%s]", t),
              list(time_level = t)
            ))
          }
        }

        # compare the times WITHIN each group
        for (g in levels(factor(meta_df$group))) {
          idx_g <- which(meta_df$group == g)
          t_here <- droplevels(factor(as.character(meta_df$time_fac[idx_g])))
          if (nlevels(t_here) < 2) next

          # also limited possibilities --> expand in the future
          if ("pairwise" %in% contrasts) {
            add_tests(run_pairwise(
              idx_g,
              "time_fac",
              sprintf("time_within_group[%s]", g),
              list(group_level = g)
            ))
          } else if ("each_vs_rest" %in% contrasts) {
            add_tests(run_each_vs_rest(
              idx_g,
              "time_fac",
              sprintf("time_within_group[%s]", g),
              list(group_level = g)
            ))
          }
        }
      }
    }

    # ---- finalize tests table (always define tests_tbl) ----------------------
    tests_tbl <- if (!length(tests_out)) {
      tibble::tibble()
    } else {
      out <- dplyr::bind_rows(tests_out)
      if (!identical(mtp, "none") && nrow(out)) {
        out <- out |>
          dplyr::group_by(rank, view) |>
          dplyr::mutate(p_adj = stats::p.adjust(p_value, method = mtp)) |>
          dplyr::ungroup()
      }
      out
    }
    # 8) DISPERSION TESTING: ALWAYS global model (+ optional contrasts)
    # contrasts for dispersion (same scopes as PERMANOVA) --> helpers
    #
    # a) pairwise helper
    disp_pairwise <- function(idx,
                              col_name,
                              scope_label,
                              ctx_cols,
                              n_threads) {
      # pull the factor levels for this slice
      fac_here <- factor(as.character(meta_df[[col_name]][idx]))
      prs <- .pairs_of(levels(fac_here))
      if (!length(prs)) {
        return(NULL)
      }

      # initialize collectors
      res <- list()
      res_t <- list()
      for (p in prs) {
        # select just those samples in this pair and require enough data
        sel <- idx[which(fac_here %in% p)]
        if (length(sel) < 3) next

        # subset distance
        d_sub <- stats::as.dist(as.matrix(dist_full)[sel, sel])

        # build a 2-level factor aligned to the sliced data
        fac <- droplevels(factor(as.character(meta_df[[col_name]][sel]),
          levels = p
        ))
        if (nlevels(fac) < 2 || min(table(fac)) < 2) next

        # compute distances to centroids (one per sample)
        dd <- .betadisper_table(d_sub, fac, scope_label)

        # test dispersion difference between the two groups:
        dt <- .betadisper_test(
          d_sub,
          fac,
          permutations,
          scope_label,
          n_threads
        )

        # report the results
        if (!is.null(dd)) {
          dd$scope <- scope_label
          dd$contrast <- paste(p, collapse = " vs ")
          if (!is.null(ctx_cols$time_level)) {
            dd$time_level <- ctx_cols$time_level
          }
          if (!is.null(ctx_cols$group_level)) {
            dd$group_level <- ctx_cols$group_level
          }
          res[[length(res) + 1]] <- dd
        }
        if (!is.null(dt)) {
          dt$scope <- scope_label
          dt$contrast <- paste(p, collapse = " vs ")
          res_t[[length(res_t) + 1]] <- dt
        }
      }
      list(
        dist = if (!length(res)) NULL else dplyr::bind_rows(res),
        test = if (!length(res_t)) NULL else dplyr::bind_rows(res_t)
      )
    }

    # b) each_vs_rest helper
    disp_each_vs_rest <- function(idx,
                                  col_name,
                                  scope_label,
                                  ctx_cols,
                                  n_threads) {
      # collect levels present in idx, drop NA, ensure >=2 levels
      vec_chr <- as.character(meta_df[[col_name]][idx])
      lv <- unique(vec_chr)
      lv <- lv[!is.na(lv)]
      if (length(lv) < 2) {
        return(NULL)
      }

      # build lists and “each-vs-rest” factors
      res <- list()
      res_t <- list()
      evr <- .each_vs_rest_factors(vec_chr, lv)
      for (i in seq_along(lv)) {
        # keep samples where the 2-level factor is not NA and require enough
        # data
        keep <- !is.na(evr[[i]])
        sel <- idx[keep]
        if (length(sel) < 3) next

        # subset matrix
        d_sub <- stats::as.dist(as.matrix(dist_full)[sel, sel])

        # 2-level factor (target vs other) and guard
        fac <- droplevels(evr[[i]][keep])
        if (nlevels(fac) < 2 || min(table(fac)) < 2) next
        # compute per-sample centroid distances
        dd <- .betadisper_table(d_sub, fac, scope_label)

        # permutation test for dispersion difference
        dt <- .betadisper_test(d_sub, fac, permutations, scope_label, n_threads)

        # output
        if (!is.null(dd)) {
          dd$scope <- scope_label
          dd$contrast <- paste0(lv[i], " vs other")
          if (!is.null(ctx_cols$time_level)) {
            dd$time_level <- ctx_cols$time_level
          }
          if (!is.null(ctx_cols$group_level)) {
            dd$group_level <- ctx_cols$group_level
          }
          res[[length(res) + 1]] <- dd
        }
        if (!is.null(dt)) {
          dt$scope <- scope_label
          dt$contrast <- paste0(lv[i], " vs other")
          res_t[[length(res_t) + 1]] <- dt
        }
      }
      list(
        dist = if (!length(res)) NULL else dplyr::bind_rows(res),
        test = if (!length(res_t)) NULL else dplyr::bind_rows(res_t)
      )
    }

    # helper to add the tests and dispersion to the output; modifies in-place
    disp_out <- list()
    disp_tests_out <- list()
    add_disp <- function(dd) {
      if (!is.null(dd) && nrow(dd)) {
        dd$rank <- rank_name
        dd$view <- view_name
        disp_out[[length(disp_out) + 1]] <<- dd
      }
    }

    add_dpt <- function(tt) {
      if (!is.null(tt) && nrow(tt)) {
        tt$rank <- rank_name
        tt$view <- view_name
        disp_tests_out[[length(disp_tests_out) + 1]] <<- tt
      }
    }

    # -------- TESTING GLOBAL --------------------------------------------------
    # group: always test (betadisper), and if time is continuous, also ANCOVA
    # on distances
    if (length(unique(meta_df$group)) > 1) {
      # classic betadisper global by group
      dd <- .betadisper_table(
        dist_full,
        factor(meta_df$group),
        "dispersion[group]"
      )
      dt <- .betadisper_test(
        dist_full,
        factor(meta_df$group),
        permutations,
        "dispersion[group]",
        n_threads
      )

      ## append to the output
      if (!is.null(dd)) {
        dd$scope <- "group"
        dd$contrast <- "<global>"
        add_disp(dd)
      }

      if (!is.null(dt)) {
        dt$scope <- "group"
        dt$contrast <- "<global>"
        add_dpt(dt)
      }

      # when time is continuous --> dispersion ANCOVA: dist ~ group * time_cont
      if (isTRUE(time_is_continuous)) {
        bd <- try(vegan::betadisper(
          dist_full,
          factor(meta_df$group)
        ), silent = TRUE)
        if (!inherits(bd, "try-error")) {
          ddf <- tibble::tibble(
            sample_id = names(bd$distances),
            dist = as.numeric(bd$distances)
          ) |>
            dplyr::left_join(meta_df, by = "sample_id")
          # build the fac for the global multi-level test
          ddf$fac <- factor(ddf$group)

          tt_ct <- .dispersion_tests_continuous(ddf,
            scope_label    = "group|time_cont",
            contrast_label = "<global>",
            nperm          = permutations,
            n_threads      = n_threads
          )
          add_dpt(tt_ct)
        }
      }
    }

    # Dispatcher for the pairwise comp --> exact the same logic as with
    # PERMANOVA
    #
    # time: ONLY if categorical (betadisper requires a factor)
    if (!time_is_continuous &&
      has_time &&
      length(unique(meta_df$time_fac)) > 1) {
      # disper and tests
      dd <- .betadisper_table(
        dist_full,
        factor(meta_df$time_fac),
        "dispersion[time]"
      )
      dt <- .betadisper_test(
        dist_full,
        factor(meta_df$time_fac),
        permutations,
        "dispersion[time]",
        n_threads
      )
      # add output
      if (!is.null(dd)) {
        dd$scope <- "time"
        dd$contrast <- "<global>"
        add_disp(dd)
      }

      if (!is.null(dt)) {
        dt$scope <- "time"
        dt$contrast <- "<global>"
        add_dpt(dt)
      }
    }

    # group:time: again only when time is categorical
    if (!time_is_continuous &&
      has_time &&
      length(unique(meta_df$group)) > 1 &&
      length(unique(meta_df$time_fac)) > 1) {
      # defining the interaction
      inter <- factor(paste(meta_df$group, meta_df$time_fac, sep = " * "))

      # disper and test
      dd <- .betadisper_table(
        dist_full,
        inter,
        "dispersion[group:time]"
      )
      dt <- .betadisper_test(
        dist_full,
        inter,
        permutations,
        "dispersion[group:time]",
        n_threads
      )

      if (!is.null(dd)) {
        dd$scope <- "group:time"
        dd$contrast <- "<global>"
        add_disp(dd)
      }

      if (!is.null(dt)) {
        dt$scope <- "group:time"
        dt$contrast <- "<global>"
        add_dpt(dt)
      }
    }

    # -------- TESTING CONTRASTS -----------------------------------------------
    if (!("none" %in% contrasts)) {
      # existing helpers (pairwise / each_vs_rest) for categorical factors
      # stay as-is...

      # dispatch dispersion contrasts across GROUPS (always possible)
      if ("pairwise" %in% contrasts && length(unique(meta_df$group)) > 1) {
        out <- disp_pairwise(
          seq_len(nrow(meta_df)),
          "group",
          "group",
          list(),
          n_threads
        )
        # output
        if (!is.null(out$dist)) add_disp(out$dist)
        if (!is.null(out$test)) add_dpt(out$test)
      } else if ("each_vs_rest" %in% contrasts &&
        length(unique(meta_df$group)) > 1) {
        out <- disp_each_vs_rest(
          seq_len(nrow(meta_df)),
          "group",
          "group",
          list(),
          n_threads
        )
        # output
        if (!is.null(out$dist)) add_disp(out$dist)
        if (!is.null(out$test)) add_dpt(out$test)
      }

      # pairwise GROUP dispersion with continuous time (ANCOVA on betadisper
      # distances)
      if (isTRUE(time_is_continuous) &&
        "pairwise" %in% contrasts &&
        length(unique(meta_df$group)) > 1) {
        # building the pairs
        fac_here <- factor(as.character(meta_df$group))
        prs <- .pairs_of(levels(fac_here))
        if (length(prs)) {
          for (p in prs) {
            sel <- which(fac_here %in% p)
            if (length(sel) < 4) next

            # sub-distance + pairwise betadisper
            d_sub <- stats::as.dist(as.matrix(dist_full)[sel, sel])
            fac <- droplevels(factor(as.character(meta_df$group[sel]),
              levels = p
            ))
            if (nlevels(fac) < 2 || min(table(fac)) < 2) next

            # calculating dispersion
            bd2 <- try(vegan::betadisper(d_sub, fac), silent = TRUE)
            if (inherits(bd2, "try-error")) next

            df2 <- tibble::tibble(
              sample_id = names(bd2$distances),
              dist = as.numeric(bd2$distances)
            ) |>
              dplyr::left_join(meta_df[sel, , drop = FALSE], by = "sample_id")
            df2$fac <- droplevels(factor(df2$group, levels = p))
            if (!all(c("time_cont", "fac", "dist") %in% names(df2))) next

            # testing dispersion
            tt2 <- .dispersion_tests_continuous(df2,
              scope_label = "group_pairwise|time_cont",
              contrast_label = paste(p, collapse = " vs "),
              nperm = permutations,
              n_threads = n_threads
            )
            add_dpt(tt2)
          }
        }
      }

      # time contrasts ONLY when categorical (your existing code)
      if (!time_is_continuous &&
        has_time &&
        length(unique(meta_df$time_fac)) > 1) {
        if ("pairwise" %in% contrasts) {
          out <- disp_pairwise(
            seq_len(nrow(meta_df)),
            "time_fac",
            "time",
            list(),
            n_threads
          )

          # out
          if (!is.null(out$dist)) add_disp(out$dist)
          if (!is.null(out$test)) add_dpt(out$test)
        }
        if ("each_vs_rest" %in% contrasts) {
          out <- disp_each_vs_rest(
            seq_len(nrow(meta_df)),
            "time_fac",
            "time",
            list(),
            n_threads
          )

          # out
          if (!is.null(out$dist)) add_disp(out$dist)
          if (!is.null(out$test)) add_dpt(out$test)
        }
      }
    }

    # finalize tables
    disp_tbl <- if (!length(disp_out)) {
      tibble::tibble()
    } else {
      dplyr::bind_rows(disp_out)
    }

    disp_tests_tbl <- if (!length(disp_tests_out)) {
      tibble::tibble()
    } else {
      out <- dplyr::bind_rows(disp_tests_out)
      if (!identical(mtp, "none") && nrow(out)) {
        out <- out |>
          dplyr::group_by(rank, view) |>
          dplyr::mutate(p_adj = stats::p.adjust(p_value, method = mtp)) |>
          dplyr::ungroup()
      }
      out
    }


    # return
    c(
      ord_out,
      list(
        tests = tests_tbl,
        dispersion = disp_tbl,
        dispersion_tests =
          disp_tests_tbl
      )
    )
  }

  # ---------- dispatch by method_pcoa over subsets ----------------------------
  subsets <- list(list(
    name = "all",
    idx = rep(TRUE, nrow(pres_tbl))
  ))

  if (identical(method_pcoa, "separate_group")) {
    if (is.null(group_col)) {
      .ph_abort("separate_group mode requires group_col.")
    }

    levels_g <- pres_tbl |>
      dplyr::distinct(group = .data[[group_col]]) |>
      dplyr::pull(1) |>
      as.character()

    subsets <- lapply(
      levels_g,
      function(g) {
        list(
          name = paste0("group=", g),
          idx = pres_tbl[[group_col]] == g
        )
      }
    )
  } else if (identical(method_pcoa, "separate_time")) {
    if (!has_time) {
      .ph_abort("separate_time mode requires time_col.")
    }

    levels_t <- time_levels %||% (pres_tbl |>
      dplyr::distinct(.data[[time_col]]) |>
      dplyr::pull(1) |>
      as.character())

    subsets <- lapply(levels_t, function(t) {
      list(
        name = paste0("time=", t),
        idx = pres_tbl[[time_col]] == t
      )
    })
  } else if (identical(method_pcoa, "separate_all")) {
    if (is.null(group_col) || !has_time) {
      .ph_abort("separate_all mode requires both group_col and time_col.")
    }
    levels_g <- pres_tbl |>
      dplyr::distinct(group = .data[[group_col]]) |>
      dplyr::pull(1) |>
      as.character()

    levels_t <- time_levels %||% (pres_tbl |>
      dplyr::distinct(.data[[time_col]]) |>
      dplyr::pull(1) |>
      as.character())
    subsets <- unlist(lapply(levels_g, function(g) {
      lapply(levels_t, function(t) {
        list(
          name = paste0("group=", g, ";time=", t),
          idx = pres_tbl[[group_col]] == g &
            pres_tbl[[time_col]] == t
        )
      })
    }), recursive = FALSE)
  } else if (identical(method_pcoa, "cap")) {
    subsets <- list(list(name = "CAP_global", idx = rep(TRUE, nrow(pres_tbl))))
  }

  # run per subset and per rank
  outs <- list()
  for (ss in subsets) {
    tbl_sub <- pres_tbl[ss$idx, , drop = FALSE]
    if (!nrow(tbl_sub)) next
    res_ranks <- lapply(ranks, function(rn) {
      run_one_rank(rn, tbl_sub, sublabel = ss$name)
    })
    keep_idx <- !vapply(res_ranks, is.null, logical(1))
    res_ranks <- res_ranks[keep_idx]
    if (!length(res_ranks)) next
    names(res_ranks) <- ranks[keep_idx]
    outs[[ss$name]] <- list(
      pcoa = dplyr::bind_rows(lapply(res_ranks, `[[`, "pcoa")),
      pcoa_full = dplyr::bind_rows(lapply(res_ranks, `[[`, "pcoa_full")),
      var_explained = dplyr::bind_rows(lapply(
        res_ranks, `[[`,
        "var_explained"
      )),
      eigen_summary = dplyr::bind_rows(lapply(
        res_ranks, `[[`,
        "eigen_summary"
      )),
      eig_full = dplyr::bind_rows(lapply(res_ranks, `[[`, "eig_full")),
      feature_loadings = dplyr::bind_rows(lapply(
        res_ranks, `[[`,
        "feature_loadings"
      )),
      tests = dplyr::bind_rows(lapply(res_ranks, `[[`, "tests")),
      dispersion = dplyr::bind_rows(lapply(
        res_ranks, `[[`,
        "dispersion"
      )),
      dispersion_tests = dplyr::bind_rows(lapply(
        res_ranks, `[[`,
        "dispersion_tests"
      )),
      pcoa_fits = lapply(res_ranks, `[[`, "pcoa_fit"),
      cap_fits = lapply(res_ranks, `[[`, "cap_fit"),
      cap_partitioning = dplyr::bind_rows(lapply(
        res_ranks, `[[`,
        "cap_partitioning"
      ))
    )
  }

  outs
}

# ==============================================================================
# Beta diversity - generic + methods  (global tests + optional contrasts only)
# ==============================================================================

#' @title Compute beta diversity (PCoA / CAP, PERMANOVA, dispersion)
#'
#' @description Computes between-sample diversity for one or more **ranks**
#' via a dissimilarity matrix, **PCoA** scores, optional **CAP/dbRDA**,
#' **PERMANOVA** (adonis2; always global maximal model) and **dispersion**
#' (betadisper + permutest). Optional post-hoc contrasts can be requested.
#'
#' @details
#' ## Ranks
#' *The peptide identities or characteristics you aggregate by*. Must be **exact
#' column names**:
#' - For `<phip_data>`: columns from the Vogl Lab peptide library (e.g.,
#'   `peptide_id`, lineage/taxa fields).
#' - For `data.frame`: columns present in your long table.
#'
#' ## Presence rule
#' - Default (`fc_threshold = NULL`): presence is `exist > 0`.
#' - If `fc_threshold` is numeric, presence is `fold_change > fc_threshold`.
#'
#' ## Full-cross / auto-expanded zeros
#' For `<phip_data>` created with a full cross (synthetic zero rows), those are
#' pruned **upfront** (keep only present rows) to reduce compute, matching the
#' alpha-diversity behavior.
#'
#' ## Normalization
#' `method_normalization` controls how the wide abundance matrix is transformed
#' *before* distances:
#' - `"auto"`: if presence-by-`exist`, keep counts (binary) as-is; else use
#'   `"relative"`.
#' - `"relative"`: divide each row by its row sum.
#' - `"hellinger"`: sqrt of relative.
#' - `"log"`: `log1p`.
#' - `"none"`: leave counts as-is.
#'
#' ## Distance engine
#' Uses **parallelDist** when supported by the requested `distance` method.
#' For methods not implemented in parallelDist (e.g., `"bray"`), it uses a
#' **threaded Bray** via the identity `BC(x,y) = L1(x,y) / (sum(x)+sum(y))`.
#' Otherwise falls back to `vegan::vegdist`.
#'
#' ## Ordination mode
#' `method_pcoa` controls *how* ordinations are computed:
#' `"joint"`, `"separate_group"`, `"separate_time"`, `"separate_all"`, `"cap"`.
#' - `"cap"` uses `vegan::capscale` with constrained axes determined by
#'   available predictors; others use (w)cmdscale on the distance.
#'
#' ## Negative eigenvalues
#' `neg_correction`: `"none"`, `"lingoes"`, or `"cailliez"`. Applied via
#' `vegan::wcmdscale(add=)` (PCoA) or `capscale(add=)` (CAP).
#'
#' ## Testing (ALWAYS global, with optional contrasts)
#' - PERMANOVA always runs a **global maximal model**:
#'   `dist ~ group + time + group:time` (dropping terms not available).
#'   If `time_col` is present and `subject_id` exists, permutations are
#'   **stratified by subject_id**.
#' - Dispersion (betadisper + permutest) is run **globally** for each available
#'   factor (`group`, `time`, `group:time`) and can also run **contrasts**.
#' - `contrasts`: `"none"`, `"pairwise"`, `"each_vs_rest"` (alias `"group_vs_rest"`),
#'   or `"baseline"` (provide `baseline_level`).
#' - Multiple-testing correction is applied **within each (view × rank)** using
#'   `mtp` (default `.ph_opt("beta.mtp","BH")`). Use `"none"` to skip.
#'
#' @param x A `<phip_data>` object or a long `data.frame`.
#' @param group_cols Character vector of grouping columns, or `NULL` for a
#'   single **non-facetted** view called `"all_samples"`. Columns must exist.
#' @param ranks Character vector of **exact column names** to aggregate by.
#' @param fc_threshold Numeric or `NULL`. Presence rule (see Details).
#' @param method_normalization One of `"auto"`, `"relative"`, `"hellinger"`,
#'   `"log"`, `"none"`.
#' @param distance Distance method (e.g., `"bray"`, `"jaccard"`, `"euclidean"`,
#'   `"manhattan"`, `"canberra"`, `"cosine"`). (Backward-compat: `method`.)
#' @param permutations Number of permutations for adonis2 / permutest
#'   (default `999`). Use `0` to skip permutation tests.
#' @param time_col Optional **categorical** time column (e.g., `"timepoint"`).
#' @param carry_cols Optional character vector of extra columns to carry into
#'   outputs (joined by `sample_id`).
#' @param filter_rank Optional vector or function to limit levels of each `rank`.
#' @param baseline_level Optional baseline level name for `"baseline"` contrasts.
#'   For time, if `NULL`, the earliest level is used.
#' @param contrasts One or more of `"none"`, `"pairwise"`, `"each_vs_rest"`
#'   (alias `"group_vs_rest"`), `"baseline"`. Default `"none"`.
#' @param mtp Multiple-testing correction method passed to `p.adjust`
#'   (default `.ph_opt("beta.mtp","BH")`).
#' @param group_interaction Whether to also compute beta diversity on the interaction of all given groups (one extra view).
#'   (default FALSE)
#' @param interaction_sep Interaction separator to use for the group interaction. Only used when `group_interaction` is TRUE.
#'   (default `"*"``)
#' @param n_threads Integer; threads for parallelDist (default: all cores - 1).
#' @param method_pcoa One of `"joint"`, `"separate_group"`, `"separate_time"`,
#'   `"separate_all"`, `"cap"`. Controls ordination splitting / CAP.
#' @param neg_correction One of `"none"`, `"lingoes"`, `"cailliez"`.
#' @param time_force_continuous Wether to force time to being continuous (even if can be categorical).
#'   Currently only relevant in the case when `method_pcoa` is `"cap"`.
#' @inheritDotParams
#'
#' @return A **named list** (one element per view or subview) with class
#'   `"phip_beta_diversity"`. Each element contains:
#'   - `pcoa`: tibble with first k axes (see option `phiper.beta.eig_axes`).
#'   - `pcoa_full`: **all** axes (diagnostics).
#'   - `var_explained`: `%PCoA1…%PCoAk` plus `%Other` (normalized by sum of
#'     **positive** eigenvalues from the **full** spectrum).
#'   - `eigen_summary`: first k rows + **Other**; includes `n_pos`, `n_neg`.
#'   - `eig_full`: all raw eigenvalues.
#'   - `feature_loadings`: weighted-average loadings (wide, axis-block order).
#'   - `tests`: global PERMANOVA + requested contrasts (tidy).
#'   - `dispersion`: per-sample centroid distances (+ contrast tests).
#'
#' @export
compute_beta_diversity <- function(x,
                                             group_cols = NULL,
                                             ranks = "peptide_id",
                                             fc_threshold = NULL,
                                             method_normalization = c("auto", "relative", "hellinger", "log", "none"),
                                             distance = "bray",
                                             permutations = 999,
                                             time_col = NULL,
                                             carry_cols = NULL,
                                             filter_rank = NULL,
                                             baseline_level = NULL,
                                             contrasts = "none",
                                             mtp = NULL,
                                             group_interaction = FALSE,
                                             interaction_sep = " * ",
                                             n_threads = max(1L, (if (rlang::is_installed("parallel")) parallel::detectCores() else 1L) - 1L),
                                             method_pcoa = c("joint", "separate_group", "separate_time", "separate_all", "cap"),
                                             neg_correction = c("none", "lingoes", "cailliez"),
                                             time_force_continuous = FALSE,
                                             ...) {
  stopifnot(inherits(x, "phip_data"))
  .data <- rlang::.data
  method_pcoa <- match.arg(method_pcoa)
  neg_correction <- match.arg(neg_correction)
  method_normalization <- match.arg(method_normalization)

  .ph_with_timing(
    headline = "Computing beta diversity (<phip_data>)",
    step = paste0(
      "group_cols: ", if (is.null(group_cols)) "<none>" else paste(add_quotes(group_cols, 1L), collapse = ", "),
      "; ranks: ", paste(add_quotes(ranks, 1L), collapse = ", "),
      "; distance: ", distance,
      "; method_normalization: ", method_normalization,
      "; permutations: ", permutations,
      "; contrasts: ", if (length(contrasts)) paste(contrasts, collapse = ",") else "<none>",
      "; method_pcoa: ", method_pcoa,
      "; neg_correction: ", neg_correction
    ),
    expr = {
      tbl <- x$data_long
      # prune full-cross zeros (match alpha behavior)
      if (isTRUE(x$meta$full_cross) && ("exist" %in% colnames(tbl)) && is.null(fc_threshold)) {
        .ph_log_info("Full-cross detected; pruning non-existing rows before beta calc",
          bullets = c("rule: keep exist == 1")
        )
        tbl <- dplyr::filter(tbl, .data$exist == 1L)
      }

      # validate group columns
      if (!is.null(group_cols) && length(group_cols)) {
        miss_gc <- setdiff(group_cols, colnames(tbl))
        if (length(miss_gc)) {
          .ph_abort(
            headline = "Grouping columns not found in data_long.",
            step = "input validation",
            bullets = sprintf("missing: %s", paste(add_quotes(miss_gc, 1L), collapse = ", "))
          )
        }
      }

      # peptide library mapping on main connection
      peplib_main <- .ensure_peplib_on_main(x)
      peplib_cols <- colnames(peplib_main)
      .ph_log_info("Peptide library attached on main connection",
        bullets = c(sprintf(
          "available columns: %s%s",
          paste(utils::head(peplib_cols, 8), collapse = ", "),
          if (length(peplib_cols) > 8) sprintf(" …(+%d)", length(peplib_cols) - 8) else ""
        ))
      )

      map_provider <- function(rank_name) {
        if (!(rank_name %in% peplib_cols)) {
          .ph_warn(
            headline = "Rank not found in peptide_library (skipping).",
            step = "rank mapping",
            bullets = sprintf("rank: %s", add_quotes(rank_name, 1L))
          )
          return(NULL)
        }
        peplib_main |>
          dplyr::select(peptide_id, rank_val = .data[[rank_name]]) |>
          dplyr::distinct() |>
          dplyr::collect() # <- make it local so joins with pres_tbl work
      }

      # dispatcher across views
      views <- list()
      if (is.null(group_cols) || !length(group_cols)) {
        views[["all_samples"]] <- .compute_beta_block(
          tbl,
          view_name = "all_samples", group_col = NULL, ranks = ranks,
          fc_threshold = fc_threshold, method_normalization = method_normalization, distance = distance,
          permutations = permutations, time_col = time_col, carry_cols = carry_cols,
          filter_rank = filter_rank, baseline_level = baseline_level, contrasts = contrasts,
          mtp = mtp, map_provider = map_provider, n_threads = n_threads,
          method_pcoa = method_pcoa, neg_correction = neg_correction,
          time_force_continuous = time_force_continuous
        )
      } else {
        for (gc in group_cols) {
          views[[gc]] <- .compute_beta_block(
            tbl,
            view_name = gc, group_col = gc, ranks = ranks,
            fc_threshold = fc_threshold, method_normalization = method_normalization, distance = distance,
            permutations = permutations, time_col = time_col, carry_cols = carry_cols,
            filter_rank = filter_rank, baseline_level = baseline_level, contrasts = contrasts,
            mtp = mtp, map_provider = map_provider, n_threads = n_threads,
            method_pcoa = method_pcoa, neg_correction = neg_correction,
            time_force_continuous = time_force_continuous
          )
        }
        if (isTRUE(group_interaction) && length(group_cols) >= 2L) {
          inter_col <- "..phip_interaction.."
          combo_nm <- paste(group_cols, collapse = interaction_sep)
          tbl_inter <- dplyr::mutate(tbl, !!rlang::sym(inter_col) := paste(!!!rlang::syms(group_cols), sep = interaction_sep))
          views[[combo_nm]] <- .compute_beta_block(
            tbl_inter,
            view_name = combo_nm, group_col = inter_col, ranks = ranks,
            fc_threshold = fc_threshold, method_normalization = method_normalization, distance = distance,
            permutations = permutations, time_col = time_col, carry_cols = carry_cols,
            filter_rank = filter_rank, baseline_level = baseline_level, contrasts = contrasts,
            mtp = mtp, map_provider = map_provider, n_threads = n_threads,
            method_pcoa = method_pcoa, neg_correction = neg_correction,
            time_force_continuous = time_force_continuous
          )
        }
      }

      class(views) <- c("phip_beta_diversity", class(views))
      attr(views, "group_cols") <- group_cols
      attr(views, "ranks") <- ranks
      attr(views, "fc_threshold") <- fc_threshold
      attr(views, "method_normalization") <- method_normalization
      attr(views, "distance") <- distance
      attr(views, "permutations") <- permutations
      attr(views, "time_col") <- time_col
      attr(views, "contrasts") <- contrasts
      attr(views, "mtp") <- mtp %||% .ph_opt("beta.mtp", "BH")
      attr(views, "method_pcoa") <- method_pcoa
      attr(views, "neg_correction") <- neg_correction

      views
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

#' @export
print.phip_beta_diversity <- function(x,
                                      n_rows = getOption("phiper.print.rows", 8),
                                      n_terms = getOption("phiper.print.terms", 8),
                                      ...) {
  has_cli <- requireNamespace("cli", quietly = TRUE)
  has_knitr <- requireNamespace("knitr", quietly = TRUE)

  # simple helpers ------------------------------------------------------------
  rul <- function(left) {
    if (has_cli) cli::rule(left = left) else paste0("---- ", left, " ----")
  }
  bullet <- function(key, val) {
    line <- sprintf("%s: %s", key, val)
    if (has_cli) cli::cat_bullet(line) else cat(paste0(" - ", line, "\n"))
  }
  h1 <- function(txt) if (has_cli) cli::cat_line(cli::col_cyan(txt)) else cat(txt, "\n", sep = "")
  h2 <- function(txt) if (has_cli) cli::cat_line(cli::col_blue(txt)) else cat(txt, "\n", sep = "")
  fmt_num <- function(v, d = 3) if (is.numeric(v)) formatC(v, format = "f", digits = d) else v
  .kable <- function(df, n = NULL, align = NULL) {
    if (!is.null(n)) df <- utils::head(df, n)

    if (requireNamespace("knitr", quietly = TRUE)) {
      # capture printed table
      lines <- utils::capture.output(
        print(knitr::kable(df, format = "simple", align = align))
      )
      # remove leading blank lines that cause the “pause”
      while (length(lines) && grepl("^\\s*$", lines[1])) lines <- lines[-1]
      cat(paste(lines, collapse = "\n"), "\n\n", sep = "")
    } else {
      print(df)
      cat("\n")
    }
  }

  `%||%` <- function(a, b) if (is.null(a)) b else a
  .axis_cols <- function(nms) grep("^(PCoA|CAP)\\d+$", nms, value = TRUE)
  .time_mode <- function(pcoa_tbl) {
    if (is.null(pcoa_tbl) || !NROW(pcoa_tbl)) {
      return("absent")
    }
    has_cont <- ("time_cont" %in% names(pcoa_tbl)) && any(!is.na(pcoa_tbl$time_cont))
    has_fac <- ("time_fac" %in% names(pcoa_tbl)) && any(!is.na(pcoa_tbl$time_fac))
    if (has_cont) "continuous (time_cont)" else if (has_fac) "categorical (time_fac)" else "absent"
  }
  .split_vs <- function(x) {
    out <- t(vapply(strsplit(x, " vs ", fixed = TRUE), function(z) {
      c(z[1] %||% NA_character_, z[2] %||% NA_character_)
    }, character(2L)))
    colnames(out) <- c("group1", "group2")
    as.data.frame(out, stringsAsFactors = FALSE)
  }
  .pretty_term <- function(term, scope) {
    if (identical(term, "fac") && grepl("^group_pairwise\\|time_", scope)) "group|time" else term
  }

  # header --------------------------------------------------------------------
  cat(rul("<phip_beta_diversity>"), "\n\n", sep = "")
  bullet("method_pcoa", as.character(attr(x, "method_pcoa")))
  bullet("distance", as.character(attr(x, "distance")))
  bullet("neg_correction", as.character(attr(x, "neg_correction")))
  bullet("permutations", as.character(attr(x, "permutations")))
  bullet("views", length(x))
  bullet("ranks", paste(attr(x, "ranks"), collapse = ", "))
  bullet("time_col", attr(x, "time_col") %||% "<none>")
  bullet("contrasts", paste(attr(x, "contrasts"), collapse = ", "))
  bullet("mtp", attr(x, "mtp") %||% "BH")
  cat("\n") # single blank line after header bullets

  # body ----------------------------------------------------------------------
  for (view_name in names(x)) {
    view <- x[[view_name]]
    h1(rul(paste0("[view] ", view_name)))
    if (!length(view)) {
      cat("<empty view>\n\n")
      next
    }

    for (sub_name in names(view)) {
      block <- view[[sub_name]]
      h2(paste0("• subview: ", sub_name))

      # summary line + exactly one blank line
      pcoa_tbl <- block$pcoa %||% tibble::tibble()
      axis_cols <- .axis_cols(names(pcoa_tbl))
      cat(sprintf(
        "  samples: %s | groups: %s | time: %s\n",
        format(dplyr::n_distinct(pcoa_tbl$sample_id), big.mark = ","),
        if ("group" %in% names(pcoa_tbl)) dplyr::n_distinct(pcoa_tbl$group) else 1L,
        .time_mode(pcoa_tbl)
      ))
      cat("\n")

      # PCoA/CAP scores (raw tibble), then one newline ------------------------
      if (length(axis_cols)) {
        h2("  pcoa/cap scores (first rows):")
        keep_cols <- unique(c("sample_id", axis_cols, "group", "subject_id", "time_fac", "time_cont", "rank", "view"))
        keep_cols <- intersect(keep_cols, names(pcoa_tbl))
        if (inherits(pcoa_tbl, "tbl_df")) {
          print(dplyr::select(pcoa_tbl, dplyr::all_of(keep_cols)), n = n_rows, width = Inf)
        } else {
          print(utils::head(dplyr::select(pcoa_tbl, dplyr::all_of(keep_cols)), n_rows))
        }
        cat("\n") # exactly one newline after pcoa block
      } else {
        h2("  pcoa/cap scores (first rows):")
        cat("<no ordination scores found>\n\n")
      }

      # variance explained -----------------------------------------------------
      if (!is.null(block$var_explained) && NROW(block$var_explained)) {
        h2("  variance explained (% of positive constrained eigenvalues):")
        ve <- block$var_explained
        num_cols <- vapply(ve, is.numeric, logical(1))
        ve[num_cols] <- lapply(ve[num_cols], fmt_num, d = 2)
        .kable(ve, n = 1)
      }

      # eigen summary ----------------------------------------------------------
      if (!is.null(block$eigen_summary) && NROW(block$eigen_summary)) {
        h2("  eigen summary:")
        es <- dplyr::select(
          block$eigen_summary,
          axis, eigenvalue, pct_of_pos, cum_pct_of_pos,
          n_pos, n_neg, rank, view
        )
        num_cols <- vapply(es, is.numeric, logical(1))
        es[num_cols] <- lapply(es[num_cols], fmt_num, d = 3)
        .kable(utils::head(es, n_terms))
      }

      # CAP partitioning -------------------------------------------------------
      if (!is.null(block$cap_partitioning) && NROW(block$cap_partitioning)) {
        h2("  cap partitioning (total / constrained / unconstrained):")
        cp <- block$cap_partitioning
        cp$inertia <- fmt_num(cp$inertia, d = 3)
        cp$proportion <- fmt_num(cp$proportion, d = 4)
        .kable(cp, n = 3)
      }

      # PERMANOVA tests --------------------------------------------------------
      if (!is.null(block$tests) && NROW(block$tests)) {
        tests <- block$tests

        # Global
        glob <- dplyr::filter(tests, .data$contrast == "<global>")
        if (NROW(glob)) {
          h2("  tests (PERMANOVA, global):")
          gtbl <- dplyr::select(glob, term, F_stat, R2, p_value, p_adj, n_perm, scope)
          gtbl$F_stat <- fmt_num(gtbl$F_stat, 3)
          gtbl$R2 <- fmt_num(gtbl$R2, 3)
          gtbl$p_value <- fmt_num(gtbl$p_value, 4)
          if ("p_adj" %in% names(gtbl)) gtbl$p_adj <- fmt_num(gtbl$p_adj, 4)
          .kable(utils::head(gtbl, n_terms))
        }

        # Contrasts
        contr <- dplyr::filter(tests, .data$contrast != "<global>")
        if (NROW(contr)) {
          h2("  contrasts (PERMANOVA, pairwise & others):")
          contr$variable <- mapply(.pretty_term, contr$term, contr$scope)
          spl <- .split_vs(contr$contrast)
          ctbl <- dplyr::bind_cols(
            tibble::tibble(comparison = contr$contrast),
            spl[, c("group1", "group2"), drop = FALSE],
            tibble::tibble(
              variable = contr$variable,
              F_stat   = fmt_num(contr$F_stat, 3),
              R2       = if ("R2" %in% names(contr)) fmt_num(contr$R2, 3) else NA_character_,
              p_value  = fmt_num(contr$p_value, 4),
              p_adj    = if ("p_adj" %in% names(contr)) fmt_num(contr$p_adj, 4) else NA_character_,
              scope    = contr$scope
            )
          )
          ctbl <- dplyr::select(ctbl, comparison, group1, group2, variable, F_stat, R2, p_value, p_adj, scope)
          .kable(utils::head(ctbl, n_terms))
        }
      } else {
        h2("  tests (PERMANOVA):")
        cat("<none>\n\n")
      }

      # Dispersion tests -------------------------------------------------------
      if (!is.null(block$dispersion_tests) && NROW(block$dispersion_tests)) {
        dt <- block$dispersion_tests

        # Global
        glob_d <- dplyr::filter(dt, .data$contrast == "<global>")
        if (NROW(glob_d)) {
          h2("  dispersion tests (global):")
          gdt <- dplyr::select(glob_d, term, F_stat, p_value, n_perm, scope)
          if (!all(is.na(gdt$F_stat))) gdt$F_stat <- fmt_num(gdt$F_stat, 3)
          gdt$p_value <- fmt_num(gdt$p_value, 4)
          .kable(utils::head(gdt, n_terms))
        }

        # Contrasts
        contr_d <- dplyr::filter(dt, .data$contrast != "<global>")
        if (NROW(contr_d)) {
          h2("  dispersion contrasts:")
          spl2 <- .split_vs(contr_d$contrast)
          dct <- dplyr::bind_cols(
            tibble::tibble(comparison = contr_d$contrast),
            spl2[, c("group1", "group2"), drop = FALSE],
            tibble::tibble(
              variable = contr_d$term,
              F_stat   = if (!all(is.na(contr_d$F_stat))) fmt_num(contr_d$F_stat, 3) else contr_d$F_stat,
              p_value  = fmt_num(contr_d$p_value, 4),
              p_adj    = if ("p_adj" %in% names(contr_d)) fmt_num(contr_d$p_adj, 4) else NA_character_,
              scope    = contr_d$scope
            )
          )
          dct <- dplyr::select(dct, comparison, group1, group2, variable, F_stat, p_value, p_adj, scope)
          .kable(utils::head(dct, n_terms))
        }
      } else {
        h2("  dispersion tests:")
        cat("<none>\n\n")
      }
    }
  }

  invisible(x)
}
