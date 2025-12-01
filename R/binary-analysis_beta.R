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
compute_beta_diversity <- function(x, ...) UseMethod("compute_beta_diversity")

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

#' @rdname compute_beta_diversity
#' @export
compute_beta_diversity.phip_data <- function(x,
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


#' @rdname compute_beta_diversity
#' @export
compute_beta_diversity.data.frame <- function(x,
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
  tbl <- x
  .data <- rlang::.data
  method_pcoa <- match.arg(method_pcoa)
  neg_correction <- match.arg(neg_correction)
  method_normalization <- match.arg(method_normalization)

  .ph_with_timing(
    headline = "Computing beta diversity (data.frame)",
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
      if (!is.null(group_cols) && length(group_cols)) {
        miss_gc <- setdiff(group_cols, colnames(tbl))
        if (length(miss_gc)) {
          .ph_abort(
            headline = "Grouping columns not found in data.frame.",
            step = "input validation",
            bullets = sprintf("missing: %s", paste(add_quotes(miss_gc, 1L), collapse = ", "))
          )
        }
      }

      map_provider <- function(rank_name) {
        if (!(rank_name %in% colnames(tbl))) {
          .ph_warn(
            headline = "Rank not found in data.frame (skipping).",
            step = "rank mapping",
            bullets = sprintf("rank: %s", add_quotes(rank_name, 1L))
          )
          return(NULL)
        }
        tibble::tibble(peptide_id = tbl$peptide_id, rank_val = tbl[[rank_name]]) |> dplyr::distinct()
      }

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
