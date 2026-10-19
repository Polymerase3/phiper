# ==============================================================================
# compute_delta2(): resampling schemes and propensity models (phase 3)
# ==============================================================================
# Two generators serve the three terms: subject labels resampled between
# subjects (group term and interaction; a group main effect cancels out of
# within-subject deviations) and time labels resampled within subjects (time
# term). Each draw is one object shared by all peptides.

#' @title Internal helper: .ph_resample_group
#' @description Draw `B` subject labellings with the observed group sizes from
#'   the conditional law given the propensities.
#'
#' @param labels Observed labels, integers `1..K`, one per subject.
#' @param propensity `NULL` for constant propensities (a uniform random
#'   split, v1's null), a vector of `P(level 2)` for `K = 2`, or an `N x K`
#'   matrix.
#' @param B Number of draws.
#' @param seed Integer seed.
#' @param n_swap_steps Steps of the swap chain per draw; default `5 * N`.
#' @param max_tries Rejection attempts per draw for `K = 2` before falling
#'   back to the swap chain.
#' @return `list(labels = N x B integer matrix, method, accept_rate)`, with
#'   `method` one of `"uniform"`, `"rejection"`, `"swap"`.
#' @keywords internal
.ph_resample_group <- function(labels, propensity = NULL, B, seed,
                               n_swap_steps = NULL, max_tries = 1e4L) {
  chk::chk_whole_numeric(labels)
  chk::chk_whole_number(B)
  chk::chk_whole_number(seed)
  N <- length(labels)
  K <- max(labels)
  if (is.null(n_swap_steps)) n_swap_steps <- 5L * N
  chk::chk_whole_number(n_swap_steps)
  chk::chk_whole_number(max_tries)

  if (is.null(propensity)) {
    prob <- matrix(numeric(0), 0, 0)
  } else if (is.null(dim(propensity))) {
    if (K != 2L) .ph_abort("A propensity vector needs two groups.")
    chk::chk_true(length(propensity) == N)
    prob <- cbind(1 - propensity, propensity)
  } else {
    prob <- as.matrix(propensity)
    if (nrow(prob) != N || ncol(prob) != K) {
      .ph_abort("propensity must be an N x K matrix.")
    }
  }
  if (length(prob) && (anyNA(prob) || any(prob < 0))) {
    .ph_abort("propensity must be non-negative and not NA.")
  }

  res <- cpp_delta2_resample_group(
    as.integer(labels), prob, as.integer(K), as.integer(B), as.integer(seed),
    as.integer(n_swap_steps), as.integer(max_tries)
  )
  if (K == 2L && res$method == "swap") {
    .ph_warn(sprintf(
      "Rejection sampling found no labelling with the observed group sizes
      in %d tries; fell back to the swap chain. Propensity overlap is poor.",
      as.integer(max_tries)
    ))
  }
  accept_rate <- if (res$method == "swap") {
    res$n_swaps_accepted / (as.numeric(B) * n_swap_steps)
  } else {
    NA_real_
  }
  list(labels = res$labels, method = res$method, accept_rate = accept_rate)
}

#' @title Internal helper: .ph_resample_time
#' @description Draw `B` within-subject relabellings of the timepoints and
#'   their `F[t][u]` masks.
#'
#' @param sample_subject,time_code Per sample: subject index and timepoint
#'   code `1..n_time_levels`.
#' @param keep_prob Optional per-subject probability of keeping the observed
#'   orientation of a two-sample subject (see
#'   [.ph_delta2_time_propensity()]); `NULL` or `NA` is a fair coin.
#' @param sample_rows Bitset rows of the samples (default: in order).
#' @param n_words Words per bitset column.
#' @param masks Whether to build the `F[t][u]` masks.
#' @return `list(time_code = n_samples x B, masks)`; `masks` holds `B`
#'   blocks of `T * T` bitset columns, `F[t][u]` at column `(t - 1) * T + u`.
#' @keywords internal
.ph_resample_time <- function(sample_subject, time_code, n_time_levels, B,
                              seed, keep_prob = NULL,
                              sample_rows = seq_along(sample_subject),
                              n_words = (max(sample_rows) + 63L) %/% 64L,
                              masks = TRUE) {
  chk::chk_whole_number(B)
  chk::chk_whole_number(seed)
  if (is.null(keep_prob)) keep_prob <- numeric(0)
  if (length(keep_prob) && any(keep_prob < 0 | keep_prob > 1, na.rm = TRUE)) {
    .ph_abort("keep_prob must be in [0, 1].")
  }
  cpp_delta2_resample_time(
    as.integer(sample_rows), as.integer(sample_subject),
    as.integer(time_code), as.numeric(keep_prob),
    as.integer(n_time_levels), as.integer(n_words),
    as.integer(B), as.integer(seed), isTRUE(masks)
  )
}

#' @title Internal helper: .ph_delta2_propensity
#' @description Fit the group propensity model of one contrast and report
#'   overlap.
#'
#' @param subjects_meta Data frame with one row per subject and the
#'   covariates in `propensity_formula`.
#' @param group_code Group level per subject (`NA` = not in the contrast).
#' @param propensity_formula One-sided formula, e.g.
#'   `~ splines::ns(age, 3) + sex`. `NULL` means constant propensities.
#' @param n_ps_strata Optional number of quantile strata of the fitted
#'   propensity; each subject then gets its stratum's observed group share
#'   instead of its fitted value (two groups only).
#' @param trim_propensity Optional `c(lo, hi)`: subjects with a propensity
#'   outside it are dropped. Off by default; such subjects otherwise stay and
#'   contribute nothing.
#' @param min_overlap With a propensity model, warn when fewer than this
#'   share of subjects have every propensity in (0.05, 0.95).
#' @return `list(prob, keep, diagnostics)`: `prob` is `N x K` aligned to
#'   `subjects_meta` (`NULL` for constant propensities), `keep` the subjects
#'   left in the contrast, `diagnostics` a one-row tibble with `n`,
#'   `n_overlap`, `min_p`, `max_p`, `pseudo_r2`, `n_params`, `n_trimmed`.
#' @keywords internal
.ph_delta2_propensity <- function(subjects_meta, group_code,
                                  propensity_formula = NULL,
                                  n_ps_strata = NULL,
                                  trim_propensity = NULL,
                                  min_overlap = 0.8) {
  chk::chk_data(subjects_meta)
  chk::chk_true(length(group_code) == nrow(subjects_meta))
  if (!is.null(n_ps_strata)) {
    chk::chk_whole_number(n_ps_strata)
    chk::chk_true(n_ps_strata >= 2)
    if (is.null(propensity_formula)) {
      .ph_abort("n_ps_strata needs a propensity_formula to stratify on.")
    }
  }
  if (!is.null(trim_propensity)) {
    chk::chk_numeric(trim_propensity)
    chk::chk_true(length(trim_propensity) == 2L)
    chk::chk_true(trim_propensity[1] >= 0 &&
      trim_propensity[1] < trim_propensity[2] && trim_propensity[2] <= 1)
  }

  keep <- !is.na(group_code)
  code <- group_code[keep]
  K <- max(code)
  n <- length(code)
  n_k <- tabulate(code, nbins = K)
  ll0 <- sum(n_k[n_k > 0] * log(n_k[n_k > 0] / n))

  if (is.null(propensity_formula)) {
    prob <- NULL
    p_fit <- matrix(n_k / n, n, K, byrow = TRUE)
    pseudo_r2 <- 0
    n_params <- 0L
  } else {
    if (!inherits(propensity_formula, "formula") ||
      length(propensity_formula) != 2L) {
      .ph_abort("propensity_formula must be a one-sided formula, e.g. ~ age.")
    }
    dat <- as.data.frame(subjects_meta)[keep, , drop = FALSE]
    dat$.ph_y <- factor(code, levels = seq_len(K))
    f <- propensity_formula
    f[[3]] <- f[[2]]
    f[[2]] <- quote(.ph_y)

    if (K == 2L) {
      fit <- withCallingHandlers(
        stats::glm(f, family = stats::binomial(), data = dat),
        warning = function(w) {
          if (grepl("fitted probabilities", conditionMessage(w))) {
            invokeRestart("muffleWarning")
          }
        }
      )
      p2 <- unname(stats::fitted(fit))
      p_fit <- cbind(1 - p2, p2)
      n_params <- length(stats::coef(fit)) - 1L
    } else {
      if (!requireNamespace("nnet", quietly = TRUE)) {
        .ph_abort("A propensity model for 3 or more groups needs 'nnet'.")
      }
      fit <- nnet::multinom(f, data = dat, trace = FALSE)
      p_fit <- unname(as.matrix(stats::fitted(fit)))
      n_params <- length(stats::coef(fit)) - (K - 1L)
    }
    ll <- sum(log(pmax(p_fit[cbind(seq_len(n), code)], 1e-300)))
    pseudo_r2 <- if (ll0 < 0) 1 - ll / ll0 else 0

    if (!is.null(n_ps_strata)) {
      if (K != 2L) .ph_abort("n_ps_strata supports two groups only.")
      brks <- unique(stats::quantile(p_fit[, 2],
        probs = seq(0, 1, length.out = n_ps_strata + 1L), names = FALSE
      ))
      stratum <- cut(p_fit[, 2], brks, include.lowest = TRUE)
      p2 <- stats::ave(as.numeric(code == 2L), stratum)
      p_fit <- cbind(1 - p2, p2)
    }
    prob <- p_fit
  }

  # trimming (opt-in)
  n_trimmed <- 0L
  if (!is.null(trim_propensity) && !is.null(prob)) {
    lo <- trim_propensity[1]
    hi <- trim_propensity[2]
    in_range <- if (K == 2L) {
      prob[, 2] >= lo & prob[, 2] <= hi
    } else {
      apply(prob, 1, function(q) min(q) >= lo && max(q) <= hi)
    }
    n_trimmed <- sum(!in_range)
    keep[which(keep)[!in_range]] <- FALSE
  }

  prob_all <- NULL
  if (!is.null(prob)) {
    prob_all <- matrix(NA_real_, nrow(subjects_meta), K)
    prob_all[!is.na(group_code), ] <- prob
    prob_all[!keep, ] <- NA_real_
  }

  p_range <- if (K == 2L) p_fit[, 2, drop = FALSE] else p_fit
  overlap <- apply(p_range, 1, function(q) all(q > 0.05 & q < 0.95))
  diagnostics <- tibble::tibble(
    n = n, n_overlap = sum(overlap),
    min_p = min(p_range), max_p = max(p_range),
    pseudo_r2 = pseudo_r2, n_params = n_params, n_trimmed = n_trimmed
  )

  if (!is.null(propensity_formula) && sum(overlap) < min_overlap * n) {
    .ph_warn(sprintf(
      "Poor propensity overlap: %d of %d subjects have a propensity in
      (0.05, 0.95).",
      sum(overlap), n
    ))
  }
  # over-fitting: fewer than 10 subjects of the smallest group per parameter
  # of one logit equation, or a high pseudo-R2
  if (n_params > 0L && (min(n_k) / (n_params / (K - 1)) < 10 ||
    pseudo_r2 > 0.5)) {
    .ph_warn(sprintf(
      "Propensity model may be over-fitted (pseudo-R2 = %.2f, %d parameters,
      smallest group %d); propensities pinned near 0 or 1 cost power.",
      pseudo_r2, n_params, min(n_k)
    ))
  }

  list(prob = prob_all, keep = keep, diagnostics = diagnostics)
}

#' @title Internal helper: .ph_delta2_time_propensity
#' @description Biased coins for the within-subject time draw: a conditional
#'   logistic model on the within-subject covariate difference of each
#'   subject with two samples at two timepoints.
#'
#' @param samples_meta Data frame with one row per sample and the covariates
#'   in `time_propensity_formula`.
#' @param sample_subject,time_code Per sample: subject index and timepoint.
#' @param time_propensity_formula One-sided formula on sample-level
#'   covariates; `NULL` gives fair coins.
#' @return `list(keep_prob, diagnostics)`: `keep_prob` per subject is the
#'   probability of the observed orientation (`NA` for subjects the coin does
#'   not apply to, or when there is no formula).
#' @keywords internal
.ph_delta2_time_propensity <- function(samples_meta, sample_subject,
                                       time_code,
                                       time_propensity_formula = NULL) {
  chk::chk_data(samples_meta)
  n_sub <- max(sample_subject)
  keep_prob <- rep(NA_real_, n_sub)
  if (is.null(time_propensity_formula)) {
    return(list(keep_prob = keep_prob, diagnostics = NULL))
  }
  if (!inherits(time_propensity_formula, "formula") ||
    length(time_propensity_formula) != 2L) {
    .ph_abort("time_propensity_formula must be a one-sided formula.")
  }

  # two samples at two different timepoints; order each pair by time
  rows <- split(seq_along(sample_subject), sample_subject)
  pair <- vapply(rows, function(r) {
    length(r) == 2L && time_code[r[1]] != time_code[r[2]]
  }, logical(1))
  if (!any(pair)) {
    return(list(keep_prob = keep_prob, diagnostics = NULL))
  }
  rows <- lapply(rows[pair], function(r) r[order(time_code[r])])
  subj <- as.integer(names(rows))

  X <- stats::model.matrix(time_propensity_formula,
    data = as.data.frame(samples_meta)
  )
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  d <- t(vapply(rows, function(r) X[r[2], ] - X[r[1], ], numeric(ncol(X))))
  if (ncol(X) == 1L) d <- t(d)

  # the observed orientation always has y = 1; flipping every other pair to
  # (-d, y = 0) leaves the conditional likelihood unchanged and lets glm fit it
  flip <- seq_len(nrow(d)) %% 2L == 0L
  y <- as.numeric(!flip)
  dd <- d * ifelse(flip, -1, 1)
  fit <- withCallingHandlers(
    stats::glm(y ~ 0 + dd, family = stats::binomial()),
    warning = function(w) {
      if (grepl("fitted probabilities", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  beta <- stats::coef(fit)
  beta[is.na(beta)] <- 0
  p_keep <- as.numeric(stats::plogis(d %*% beta))
  keep_prob[subj] <- p_keep

  n <- length(p_keep)
  pseudo_r2 <- 1 - sum(log(pmax(p_keep, 1e-300))) / (n * log(0.5))
  diagnostics <- tibble::tibble(
    n = n, n_overlap = sum(p_keep > 0.05 & p_keep < 0.95),
    min_p = min(p_keep), max_p = max(p_keep),
    pseudo_r2 = pseudo_r2, n_params = ncol(X)
  )
  if (diagnostics$n_overlap < 0.8 * n) {
    .ph_warn(sprintf(
      "Poor time-propensity overlap: %d of %d subjects have a coin in
      (0.05, 0.95).",
      diagnostics$n_overlap, n
    ))
  }
  list(keep_prob = keep_prob, diagnostics = diagnostics)
}
