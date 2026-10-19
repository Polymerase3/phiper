# ==============================================================================
# compute_delta2(): design, rho and per-term pooling (phase 2)
# ==============================================================================
# The design says which terms a contrast tests and how each subject and sample
# is coded; cpp_delta2_score() turns it into a score z per peptide per term.

#' @title Internal helper: .ph_delta2_design
#' @description Code subjects and samples of one contrast for
#'   `cpp_delta2_score()` and list the testable terms.
#'
#' @param subjects_meta Data frame with one row per subject: `subject_id` and
#'   `group_col`. Row order fixes the subject index used by `group_code`.
#' @param samples_meta Data frame with one row per sample: `sample_id`,
#'   `subject_id` and, if given, `time_col`. Row order fixes `time_code`.
#' @param group_col Name of the between-subject group column.
#' @param group_type `"binary"` (exactly two levels) or `"multi"`.
#' @param time_col Optional name of the timepoint column.
#' @param time_type `"none"`, `"binary"` or `"multi"`; inferred from the
#'   number of timepoints when not given.
#' @param group_scores,time_scores Optional numeric scores for ordered levels,
#'   either named by level or in level order. A scored factor is tested as a
#'   1-df linear trend instead of an omnibus contrast.
#' @param min_samples_per_subject Subjects with fewer samples get an `NA`
#'   `group_code` and are left out of the contrast.
#'
#' @return A list with `terms`, `group_code` (per subject), `time_code` (per
#'   sample), `n_group_levels`, `n_time_levels`, `df` and `layer` (named per
#'   term), plus `group_levels`, `time_levels` and the aligned
#'   `group_scores` / `time_scores` (`numeric(0)` for omnibus).
#' @keywords internal
.ph_delta2_design <- function(subjects_meta, samples_meta,
                              group_col, group_type = c("binary", "multi"),
                              time_col = NULL,
                              time_type = c("none", "binary", "multi"),
                              group_scores = NULL, time_scores = NULL,
                              min_samples_per_subject = 1L) {
  chk::chk_data(subjects_meta)
  chk::chk_data(samples_meta)
  chk::chk_string(group_col)
  if (!is.null(time_col)) chk::chk_string(time_col)
  chk::chk_whole_number(min_samples_per_subject)
  chk::chk_true(min_samples_per_subject >= 1)

  miss <- setdiff(c("subject_id", group_col), colnames(subjects_meta))
  miss <- c(miss, setdiff(
    c("sample_id", "subject_id", time_col), colnames(samples_meta)
  ))
  if (length(miss)) {
    .ph_abort("Missing required columns", bullets = paste("-", miss))
  }
  if (anyDuplicated(subjects_meta$subject_id)) {
    .ph_abort("subjects_meta must have one row per subject_id.")
  }

  # --- group ------------------------------------------------------------------
  gvals <- subjects_meta[[group_col]]
  .ph_delta2_chk_categorical(gvals, group_col)
  group_levels <- .ph_delta2_levels(gvals)
  K <- length(group_levels)
  group_code <- match(gvals, group_levels)

  subj_idx <- match(samples_meta$subject_id, subjects_meta$subject_id)
  if (anyNA(subj_idx)) {
    .ph_abort("Every samples_meta subject_id must be in subjects_meta.")
  }
  n_i <- tabulate(subj_idx, nbins = nrow(subjects_meta))
  group_code[n_i < min_samples_per_subject] <- NA_integer_

  group_type <- if (missing(group_type) && K > 2L) "multi" else
    match.arg(group_type)
  if (group_type == "binary" && K > 2L && is.null(group_scores)) {
    .ph_abort(sprintf(
      "group_type = \"binary\" but %s has %d levels; use \"multi\".",
      group_col, K
    ))
  }
  group_scores <- .ph_delta2_scores(group_scores, group_levels, "group_scores")

  # --- time -------------------------------------------------------------------
  if (is.null(time_col)) {
    if (!missing(time_type) && match.arg(time_type) != "none") {
      .ph_abort("time_type needs a time_col.")
    }
    time_type <- "none"
    time_levels <- 1L
    time_code <- rep(1L, nrow(samples_meta))
  } else {
    tvals <- samples_meta[[time_col]]
    .ph_delta2_chk_categorical(tvals, time_col)
    time_levels <- .ph_delta2_levels(tvals)
    time_code <- match(tvals, time_levels)
    n_t <- length(time_levels)
    time_type <- if (missing(time_type)) {
      if (n_t > 2L) "multi" else "binary"
    } else {
      match.arg(time_type)
    }
    if (time_type == "none") {
      .ph_abort("time_type = \"none\" but a time_col was given.")
    }
    if (time_type == "binary" && n_t > 2L && is.null(time_scores)) {
      .ph_abort(sprintf(
        "time_type = \"binary\" but %s has %d levels; use \"multi\".",
        time_col, n_t
      ))
    }
  }
  T_lv <- length(time_levels)
  time_scores <- .ph_delta2_scores(time_scores, time_levels, "time_scores")

  # --- terms ------------------------------------------------------------------
  df <- c(
    group = if (length(group_scores)) 1L else K - 1L,
    time = if (length(time_scores)) 1L else T_lv - 1L
  )
  df["interaction"] <- df[["group"]] * df[["time"]]
  terms <- as.character(c(
    if (K >= 2L) "group",
    if (T_lv >= 2L) "time",
    if (K >= 2L && T_lv >= 2L) "interaction"
  ))
  layer <- c(group = "between", time = "within", interaction = "within")

  list(
    terms          = terms,
    group_code     = as.integer(group_code),
    time_code      = as.integer(time_code),
    n_group_levels = K,
    n_time_levels  = T_lv,
    df             = df[terms],
    layer          = layer[terms],
    group_levels   = group_levels,
    time_levels    = time_levels,
    group_scores   = group_scores,
    time_scores    = time_scores
  )
}

# Continuous covariates are deferred: only whole-numbered or categorical
# columns are accepted as group/time.
.ph_delta2_chk_categorical <- function(x, col) {
  if (anyNA(x)) .ph_abort(sprintf("%s must not contain NA.", col))
  if (is.double(x) && any(x != round(x))) {
    .ph_abort(sprintf(
      "%s looks continuous; continuous group/time columns are not supported
      yet. Bin it or pass it as a factor.",
      col
    ))
  }
  invisible(TRUE)
}

.ph_delta2_levels <- function(x) {
  if (is.factor(x)) {
    lv <- levels(droplevels(x))
  } else {
    lv <- sort(unique(x), method = "radix")
  }
  lv
}

# Scores aligned to levels, numeric(0) for none
.ph_delta2_scores <- function(scores, levels, arg) {
  if (is.null(scores)) {
    return(numeric(0))
  }
  chk::chk_numeric(scores)
  if (!is.null(names(scores))) {
    miss <- setdiff(as.character(levels), names(scores))
    if (length(miss)) {
      .ph_abort(sprintf("%s has no score for some levels", arg),
        bullets = paste("-", miss)
      )
    }
    scores <- scores[as.character(levels)]
  } else if (length(scores) != length(levels)) {
    .ph_abort(sprintf(
      "%s must have one score per level (%d).", arg, length(levels)
    ))
  }
  if (anyNA(scores) || length(unique(scores)) < 2L) {
    .ph_abort(sprintf("%s must be non-missing and not all equal.", arg))
  }
  unname(as.numeric(scores))
}

#' @title Internal helper: .ph_delta2_rho
#' @description Resolve the intraclass correlation used to weight subject
#'   means in the between-subject term.
#' @param rho `"pooled"` (estimate from all peptides), `"one"` (equal weight
#'   per subject), `"zero"` (equal weight per sample) or a number in `[0, 1]`.
#' @return `list(rho, rho_hat, subject_hits)`; `rho_hat` is `NA` unless
#'   estimated, and is `NA` with `rho = 0` when no subject has repeat samples.
#' @keywords internal
.ph_delta2_rho <- function(rho = "pooled", bitset_raw, n_words, pep_cols,
                           sample_rows, sample_subject) {
  if (is.numeric(rho)) {
    chk::chk_number(rho)
    chk::chk_range(rho, c(0, 1))
    return(list(rho = rho, rho_hat = NA_real_))
  }
  rho <- match.arg(rho, c("pooled", "one", "zero"))
  if (rho == "one") {
    return(list(rho = 1, rho_hat = NA_real_))
  }
  if (rho == "zero") {
    return(list(rho = 0, rho_hat = NA_real_))
  }
  est <- cpp_delta2_rho(
    bitset_raw, n_words, as.integer(pep_cols),
    as.integer(sample_rows), as.integer(sample_subject)
  )
  list(
    rho = if (is.na(est$rho_hat)) 0 else est$rho_hat,
    rho_hat = est$rho_hat,
    subject_hits = est$subject_hits
  )
}

#' @title Internal helper: .ph_delta2_pool
#' @description Pool a per-peptide term (rows = peptides, columns = term
#'   components) into the profile `E` (one weighted Stouffer per component)
#'   and `T = sum(E^2)`. For a 1-df term `T` is the squared Stouffer z.
#' @keywords internal
.ph_delta2_pool <- function(r, w = NULL) {
  r <- as.matrix(r)
  if (is.null(w)) w <- rep(1, nrow(r))
  chk::chk_numeric(w)
  chk::chk_true(length(w) == nrow(r))
  den <- sqrt(sum(w^2))
  E <- if (den > 0) colSums(w * r) / den else rep(0, ncol(r))
  list(E = E, T = sum(E^2))
}

#' @title Internal helper: .ph_delta2_dropout
#' @description Dropout diagnostics for one contrast: distribution of samples
#'   per subject, and its association with group and with breadth (mean
#'   positives per sample). Warns when either association exceeds
#'   `warn_cor` in absolute value.
#' @param n_i Samples per subject.
#' @param group_code Group level per subject.
#' @param subject_hits Total positives per subject.
#' @keywords internal
.ph_delta2_dropout <- function(n_i, group_code, subject_hits,
                               warn_cor = 0.3) {
  keep <- !is.na(group_code) & n_i > 0
  n_i <- n_i[keep]
  group_code <- group_code[keep]
  breadth <- subject_hits[keep] / n_i

  varies <- length(unique(n_i)) > 1L
  # n_i vs group: correlation ratio (eta), |point-biserial| for two groups
  cor_group <- NA_real_
  if (varies && length(unique(group_code)) > 1L) {
    fitted <- stats::ave(n_i, group_code)
    cor_group <- sqrt(sum((fitted - mean(n_i))^2) / sum((n_i - mean(n_i))^2))
  }
  cor_breadth <- NA_real_
  if (varies && length(unique(breadth)) > 1L) {
    cor_breadth <- stats::cor(n_i, breadth, method = "spearman")
  }

  if (!is.na(cor_group) && cor_group > warn_cor) {
    .ph_warn(sprintf(
      "Samples per subject differ by group (eta = %.2f); the time and
      interaction terms assume visits are missing at random.",
      cor_group
    ))
  }
  if (!is.na(cor_breadth) && abs(cor_breadth) > warn_cor) {
    .ph_warn(sprintf(
      "Samples per subject track breadth (Spearman = %.2f); dropout may be
      informative.",
      cor_breadth
    ))
  }

  list(
    n_dist = table(n_i),
    cor_group = cor_group,
    cor_breadth = cor_breadth
  )
}
