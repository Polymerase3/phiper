# ==============================================================================
# compute_delta2(): sample-level global prevalence shift
# ==============================================================================
# compute_delta() collapses a subject's samples into one bitset row
# (distinct(subject_id, peptide_id)), so positivity is an OR over all of a
# subject's samples and subjects with more samples get more chances to look
# positive. compute_delta2() indexes samples instead and tests a group term
# (between subjects), a time term (within subjects) and their interaction.
# The data structure is built by .ph_delta2_prepare() (below), the design and
# per-peptide statistic in R/delta2_design.R / src/delta2_score.cpp, the
# resampling in R/delta2_resample.R, and the pooled permutation test in
# src/delta2_contrast.cpp. compute_delta() itself is left as it is.

#' @title Sample-level Global Shift in Peptide Prevalence: Group, Time and
#'   Interaction
#'
#' @description
#' Tests, within each `(rank, feature, group_col)` stratum, whether peptide
#' prevalence differs between groups of subjects (**group** term), changes
#' over timepoints within subjects (**time** term), or changes differently
#' over time between groups (**interaction**). Unlike [compute_delta()],
#' every sample is its own row: a subject's samples are never collapsed, so
#' the number of samples per subject does not act as a hidden covariate.
#'
#' @details
#' **Per-peptide statistic.** Each term is a score statistic standardised by
#' its exact permutation moments:
#' - *group*: each subject contributes its mean positivity minus the
#'   weighted mean over subjects, times its weight
#'   \eqn{w_i = n_i / (1 + (n_i - 1)\rho)} (`rho` below); group labels are
#'   resampled across subjects with the group sizes fixed. Centring at the
#'   weighted mean keeps unequal sample counts from shifting every peptide
#'   together. With one sample
#'   per subject and two groups this is `compute_delta(stat_mode = "score")`
#'   up to the constant factor \eqn{\sqrt{(N-1)/N}}.
#' - *time*: each subject's deviations from its own mean (the fixed-effects
#'   "within" estimator), with timepoints reshuffled inside each subject. With
#'   two complete timepoints this is the McNemar z; with `T` it is Cochran's Q.
#'   Subjects with one sample contribute nothing.
#' - *interaction*: the within-subject deviations, residualised against the
#'   pooled time main effect, compared across group labels and standardised
#'   by a Welch-type variance built from each group's own covariance of the
#'   deviations, recomputed for every draw. The interaction's null is about
#'   mean changes, so this keeps it calibrated when groups differ in
#'   turnover or visit patterns.
#'
#' A term with one degree of freedom (two levels, or `group_scores` /
#' `time_scores` for a linear trend) gives a signed z per peptide; an omnibus
#' term (more levels) gives a vector whose squared length is the term's
#' quadratic form.
#'
#' **Pooling.** Per-peptide values are winsorised at `winsor_z` and pooled per
#' component as in [compute_delta()] (`weight_mode`, `aggregate_stat`,
#' strata of pooled prevalence from `strat_bins`) into the profile `E`. The
#' term statistic is `T = E` for one degree of freedom (signed, like
#' [compute_delta()]'s `T_obs`) and `T = sum(E^2)` otherwise; p-values compare
#' `|T|` with its permutation null. `n_eff_sqrt` weights by the square root
#' of a peptide's positive samples in the contrast; `se_invvar` by
#' \eqn{1/\sqrt{p_0(1 - p_0)}} of its pooled prevalence.
#'
#' **Resampling.** The group term, the interaction and the pairwise tests
#' share one set of `B_permutations` subject relabellings per `group_col`:
#' a uniform split of the subjects with the observed group sizes, or, with
#' `propensity_formula`, draws from the conditional law of the labels given
#' the covariates (rejection sampling for two groups, the pairwise-swap
#' chain of Berrett et al. 2020 for more). With `propensity_formula` the
#' group term and the pairwise tests are also adjusted for its covariates:
#' the subject means are residualised on them by weighted least squares
#' before the groups are compared (the score test for group in a model with
#' the covariates). The draws make the p-value valid; the adjustment keeps
#' the statistic from being dominated by the covariates, without which the
#' test is conservative and loses its power. The time term uses within-subject
#' relabellings: a coin per two-sample subject (biased by
#' `time_propensity_formula`) or a uniform permutation of the subject's
#' timepoints.
#'
#' **Type inference** (`"auto"`): two levels are binary; more levels of a
#' factor or character column are multi; a numeric column is multi with at
#' most 10 distinct values and otherwise an error, since continuous group and
#' time columns are not supported yet. A `group_col` must be constant within
#' subject; a column that varies within subject belongs in `time_col`.
#'
#' **Pairwise tests and max-T.** With `pairwise = TRUE` every pair of group
#' levels gets the binary group statistic, computed from the same subject
#' relabellings. `p_perm` is the pair's own permutation p-value; `p_maxT`
#' compares the pair's `|T|` with the maximum over all pairs in each draw,
#' which controls the family-wise error over pairs. Its null is "no group
#' differs": a claim about one pair while another group truly differs needs
#' that pair's own null.
#'
#' **Effect sizes.** `max_delta_*` is the largest absolute per-peptide effect:
#' the difference in subject-mean prevalence between the top and bottom
#' (scored) group for the group term, the within-subject change between the
#' extreme timepoints for the time term, and the difference of those changes
#' between groups for the interaction; for omnibus terms, the largest range
#' across levels (and the `frac_delta_pos_*` columns are `NA`). The `_adj`
#' version weights subjects by their inverse propensity (Hajek-normalised
#' within group), so it is adjusted for the covariates in
#' `propensity_formula`; without a propensity model it equals the raw one.
#' The time term's effect is a within-subject estimate and is not
#' re-weighted.
#'
#' @param x A `phip_data` object or a data frame in long format with the
#'   subject, sample, peptide, presence, group and time columns.
#' @param rank_cols Character vector of rank columns (e.g. `"peptide_id"`,
#'   `"species"`). Non-`"peptide_id"` ranks are joined from a peptide library
#'   as in [compute_delta()].
#' @param group_cols Character vector of between-subject grouping columns; each
#'   is tested separately.
#' @param exist_col Name of the 0/1 presence column. Default `"exist"`.
#' @param subject_col,sample_col Names of the subject and sample ID columns.
#' @param time_col Optional timepoint column. `NULL` tests the group term only.
#' @param group_type,time_type `"auto"` (see Details), `"binary"`, `"multi"`;
#'   `time_type` also takes `"none"`.
#' @param group_scores,time_scores Optional numeric scores of ordered levels
#'   (named by level), turning the omnibus term into a 1-df linear trend.
#'   `group_scores` may be a list named by `group_cols`.
#' @param test_terms `"auto"` (every term the design supports) or a subset of
#'   `c("group", "time", "interaction")`.
#' @param pairwise Logical; add the pairwise group tests (`$pairwise`).
#' @param propensity_formula Optional one-sided formula of subject-level
#'   covariates for the group labels, e.g. `~ splines::ns(age, 3) + sex`.
#' @param time_propensity_formula Optional one-sided formula of sample-level
#'   covariates for the within-subject time coins.
#' @param n_ps_strata Optional number of quantile strata of the fitted
#'   propensity (two groups only).
#' @param trim_propensity Optional `c(lo, hi)`: drop subjects whose propensity
#'   falls outside. Off by default.
#' @param n_swap_steps Steps of the swap chain per draw; default `5 * N`.
#' @param rho Intraclass correlation weighting subject means in the group
#'   term: `"pooled"` (estimated once per `group_col` from all peptides),
#'   `"one"` (equal weight per subject), `"zero"` (equal weight per sample)
#'   or a number in `[0, 1]`.
#' @param rho_by_stratum Logical; with `rho = "pooled"`, estimate it within
#'   each `strat_bins` prevalence stratum instead of globally.
#' @param min_samples_per_subject Subjects with fewer samples are dropped.
#' @param B_permutations Number of permutation draws (at least 100).
#' @param weight_mode One of `c("n_eff_sqrt", "equal", "se_invvar")`.
#' @param aggregate_stat One of `c("stouffer", "maxmean", "af")`.
#' @param perm_method One of `c("mid_p", "standard")`, as in
#'   [compute_delta()].
#' @param strat_bins Pooled-prevalence cutpoints, or `0` for no strata.
#' @param winsor_z Winsorisation threshold of per-peptide values.
#' @param min_m_eff Minimum effective number of peptides for a stratum to be
#'   tested.
#' @param rank_feature_keep Optional named list restricting features per rank.
#' @param peptide_library Optional peptide library with the rank columns.
#' @param log Logical; log progress.
#' @param log_file Log file when `log = TRUE`.
#'
#' @return An object of class `delta2_result`: a list with
#' - `terms`: one row per `(rank, feature, group_col)` with `design`
#'   (`"between"`, `"within"` or `"mixed"`), `n_subjects`,
#'   `n_samples_used`, `n_peptides_used`, `m_eff`, `rho_hat`, `n_overlap`,
#'   and per term (suffix `_group`, `_time`, `_interaction`): `df`, `T`,
#'   `T_null_mean`, `T_null_sd`, `T_stand`, `p`, `Z_from_p`, `b` (draws
#'   with a larger `|T|`), `E` (list-column, the pooled profile),
#'   `max_delta`, `max_delta_adj`, `frac_delta_pos`, `frac_delta_pos_w`.
#'   Columns of a term that is never tested are omitted.
#' - `pairwise`: one row per `(rank, feature, group_col, group1, group2)` with
#'   `T_obs`, `p_perm`, `p_maxT`, `max_delta`, `max_delta_adj`,
#'   `frac_delta_pos`, `frac_delta_pos_w` (empty without pairwise tests).
#' - `meta`: the call, `B_permutations`, and per `group_col` the design, the
#'   propensity and time-propensity diagnostics, the dropout diagnostics, the
#'   sampler used and `rho_hat` (per stratum with `rho_by_stratum = TRUE`).
#'
#' @seealso [compute_delta()]
#'
#' @examples
#' set.seed(1)
#' n <- 40
#' samples <- data.frame(
#'   sample_id = c(paste0("s", 1:n, "_T1"), paste0("s", 1:n, "_T2")),
#'   subject_id = rep(paste0("s", 1:n), 2),
#'   timepoint = rep(c("T1", "T2"), each = n),
#'   group = rep(rep(c("ctrl", "case"), each = n / 2), 2)
#' )
#' x <- merge(samples, data.frame(peptide_id = paste0("p", 1:20)))
#' x$exist <- rbinom(nrow(x), 1, 0.3)
#'
#' res <- compute_delta2(
#'   x,
#'   rank_cols = "peptide_id", group_cols = "group",
#'   time_col = "timepoint", rank_feature_keep = list(peptide_id = "p1"),
#'   B_permutations = 200L
#' )
#' res
#' res$terms
#'
#' @export
compute_delta2 <- function(
  x, rank_cols, group_cols,
  exist_col = "exist",
  subject_col = "subject_id",
  sample_col = "sample_id",
  time_col = NULL,
  group_type = c("auto", "binary", "multi"),
  time_type = c("auto", "none", "binary", "multi"),
  group_scores = NULL, time_scores = NULL,
  test_terms = "auto",
  pairwise = TRUE,
  propensity_formula = NULL,
  time_propensity_formula = NULL,
  n_ps_strata = NULL,
  trim_propensity = NULL,
  n_swap_steps = NULL,
  rho = "pooled",
  rho_by_stratum = FALSE,
  min_samples_per_subject = 1L,
  B_permutations = 2000L,
  weight_mode = c("n_eff_sqrt", "equal", "se_invvar"),
  aggregate_stat = c("stouffer", "maxmean", "af"),
  perm_method = c("mid_p", "standard"),
  strat_bins = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50),
  winsor_z = 4.0,
  min_m_eff = 0,
  rank_feature_keep = NULL,
  peptide_library = NULL,
  log = FALSE,
  log_file = "compute_delta2.log"
) {
  call <- match.call()

  # --- 0) Argument validation -------------------------------------------------
  chk::chk_character(rank_cols)
  chk::chk_true(length(rank_cols) >= 1)
  chk::chk_character(group_cols)
  chk::chk_true(length(group_cols) >= 1)
  chk::chk_string(exist_col)
  chk::chk_string(subject_col)
  chk::chk_string(sample_col)
  if (!is.null(time_col)) chk::chk_string(time_col)
  group_type <- match.arg(group_type)
  time_type <- match.arg(time_type)
  chk::chk_character(test_terms)
  all_terms <- c("group", "time", "interaction")
  if (!identical(test_terms, "auto")) {
    bad <- setdiff(test_terms, all_terms)
    if (length(bad)) {
      .ph_abort("Unknown test_terms", bullets = paste("-", bad))
    }
    if (is.null(time_col) && any(c("time", "interaction") %in% test_terms)) {
      .ph_abort("The time and interaction terms need a time_col.")
    }
  }
  chk::chk_flag(pairwise)
  chk::chk_flag(rho_by_stratum)
  chk::chk_whole_number(min_samples_per_subject)
  chk::chk_number(B_permutations)
  chk::chk_true(B_permutations >= 100)
  weight_mode <- match.arg(weight_mode)
  aggregate_stat <- match.arg(aggregate_stat)
  perm_method <- match.arg(perm_method)
  chk::chk_numeric(strat_bins)
  if (!(length(strat_bins) == 1L && isTRUE(all.equal(strat_bins, 0)))) {
    chk::chk_true(all(is.finite(strat_bins)))
    chk::chk_true(all(strat_bins >= 0 & strat_bins <= 1))
    strat_bins <- sort(unique(as.numeric(strat_bins)))
  }
  chk::chk_number(winsor_z)
  chk::chk_true(winsor_z > 0)
  chk::chk_number(min_m_eff)
  chk::chk_true(min_m_eff >= 0)
  if (is.list(group_scores)) {
    miss <- setdiff(names(group_scores), group_cols)
    if (is.null(names(group_scores)) || length(miss)) {
      .ph_abort("A group_scores list must be named by group_cols.")
    }
  } else if (!is.null(group_scores) && length(group_cols) > 1L) {
    .ph_abort("With several group_cols, pass group_scores as a named list.")
  }
  if (rho_by_stratum && !identical(rho, "pooled")) {
    .ph_warn("rho_by_stratum only applies to rho = \"pooled\"; ignored.")
    rho_by_stratum <- FALSE
  }
  B <- as.integer(B_permutations)

  # --- 1) Sample-level bitset -------------------------------------------------
  x_cols <- if (inherits(x, "phip_data")) {
    colnames(x$data_long)
  } else {
    colnames(x)
  }
  cov_vars <- function(f) {
    if (is.null(f)) character() else intersect(all.vars(f), x_cols)
  }
  subj_covs <- cov_vars(propensity_formula)
  samp_covs <- cov_vars(time_propensity_formula)
  extra_cols <- setdiff(
    unique(c(subj_covs, samp_covs)),
    c(subject_col, sample_col, time_col, group_cols, "peptide_id", exist_col)
  )
  prep <- .ph_delta2_prepare(x,
    group_cols = group_cols, exist_col = exist_col, time_col = time_col,
    subject_col = subject_col, sample_col = sample_col,
    extra_cols = if (length(extra_cols)) extra_cols
  )
  peptides_order <- prep$peptides_order
  samples_meta <- prep$samples_meta
  sample_subject <- prep$sample_subject

  # --- 2) Rank/feature map (shared with compute_delta) ------------------------
  rank_map_long <- .ph_delta_rank_map(
    x, rank_cols, peptides_order, peptide_library
  ) |>
    .ph_delta_keep_features(rank_feature_keep)
  rf_counts <- rank_map_long |>
    dplyr::count(rank, feature, name = "n_peptides")
  pep_col_map <- stats::setNames(seq_along(peptides_order), peptides_order)
  pep_cols_by_rf <- lapply(
    split(rank_map_long, rank_map_long$rank),
    function(d) {
      lapply(split(as.integer(pep_col_map[d$peptide_id]), d$feature), sort)
    }
  )

  # --- 3) One contrast per group_col: design, propensity, rho, draws ----------
  # propensity covariates must be complete (glm would drop those subjects),
  # and subject-level ones must not vary within subject
  for (v in unique(c(subj_covs, samp_covs))) {
    if (anyNA(samples_meta[[v]])) {
      .ph_abort(sprintf("Propensity covariate %s has missing values.", v))
    }
  }
  first_row <- !duplicated(sample_subject)
  for (v in subj_covs) {
    n_val <- tapply(samples_meta[[v]], sample_subject, function(z) {
      length(unique(z))
    })
    if (any(n_val > 1L)) {
      .ph_abort(sprintf(
        "propensity_formula covariate %s varies within subject; use
        time_propensity_formula for sample-level covariates.",
        v
      ))
    }
  }

  contrasts <- list()
  meta_gc <- list()
  for (gc in group_cols) {
    gvals <- prep$subject_group[[gc]]
    if (anyNA(gvals)) {
      .ph_abort(sprintf(
        "%s varies within subject; a group_col must be constant within
        subject. Pass a within-subject column as time_col.",
        gc
      ))
    }
    subjects_meta <- tibble::tibble(subject_id = prep$subjects_order)
    subjects_meta[[gc]] <- prep$group_levels[[gc]][gvals]
    for (v in subj_covs) {
      subjects_meta[[v]] <- samples_meta[[v]][first_row][
        order(sample_subject[first_row])
      ]
    }
    gs <- if (is.list(group_scores)) group_scores[[gc]] else group_scores
    design_args <- list(
      group_col = gc, group_type = group_type,
      time_col = time_col, time_type = time_type,
      group_scores = gs, time_scores = time_scores
    )
    des0 <- do.call(.ph_delta2_design, c(
      list(subjects_meta, samples_meta),
      design_args,
      list(min_samples_per_subject = min_samples_per_subject)
    ))

    # group propensity (and opt-in trimming) on the subjects left
    keep_subj <- !is.na(des0$group_code)
    ps <- NULL
    if (des0$n_group_levels >= 2L) {
      ps <- .ph_delta2_propensity(subjects_meta, des0$group_code,
        propensity_formula = propensity_formula, n_ps_strata = n_ps_strata,
        trim_propensity = trim_propensity
      )
      keep_subj <- ps$keep
    }

    # recode on the kept subjects and their samples
    subj_k <- which(keep_subj)
    samp_k <- which(sample_subject %in% subj_k)
    des <- do.call(.ph_delta2_design, c(
      list(subjects_meta[subj_k, ], samples_meta[samp_k, ]),
      design_args
    ))
    terms <- des$terms
    if (!identical(test_terms, "auto")) terms <- intersect(terms, test_terms)
    if (!length(terms)) {
      .ph_warn(sprintf("No testable terms for group_col %s; skipped.", gc))
      next
    }
    sub_k <- match(sample_subject[samp_k], subj_k)
    rows_k <- samp_k # bitset rows are the sample order
    K <- des$n_group_levels
    pairs <- if (pairwise && K >= 2L && "group" %in% terms) {
      t(utils::combn(K, 2L))
    } else {
      matrix(integer(0), 0, 2)
    }
    need_group <- any(c("group", "interaction") %in% terms) || nrow(pairs) > 0

    # rho: pooled over all peptides of the contrast (optionally per stratum)
    est <- cpp_delta2_rho(
      prep$bitset_raw, prep$bitset_words, seq_len(prep$bitset_m),
      rows_k, sub_k
    )
    rho_res <- if (identical(rho, "pooled")) {
      list(
        rho = if (is.na(est$rho_hat)) 0 else est$rho_hat,
        rho_hat = est$rho_hat
      )
    } else {
      .ph_delta2_rho(rho)
    }
    rho_vec <- rep(rho_res$rho, prep$bitset_m)
    rho_strata <- NULL
    if (rho_by_stratum) {
      bins <- .ph_delta2_prevalence_bins(
        est$pep_tot / est$n_samples, strat_bins
      )
      rho_strata <- tibble::tibble(stratum = sort(unique(bins)))
      rho_strata$rho_hat <- vapply(rho_strata$stratum, function(b) {
        cpp_delta2_rho(
          prep$bitset_raw, prep$bitset_words, which(bins == b), rows_k, sub_k
        )$rho_hat
      }, numeric(1))
      rho_strata$n_peptides <- as.integer(
        table(bins)[as.character(rho_strata$stratum)]
      )
      rb <- rho_strata$rho_hat
      rb[is.na(rb)] <- 0
      rho_vec <- rb[match(bins, rho_strata$stratum)]
    }

    n_i <- tabulate(sub_k, nbins = length(subj_k))
    dropout <- .ph_delta2_dropout(n_i, des$group_code, est$subject_hits[subj_k])

    # draws: seeds from R's RNG, so set.seed() reproduces the run
    prob_k <- if (!is.null(ps) && !is.null(ps$prob)) {
      ps$prob[subj_k, des0$group_levels %in% des$group_levels, drop = FALSE]
    }
    if (!is.null(prob_k)) prob_k <- prob_k / rowSums(prob_k)
    group_draws <- matrix(0L, 0, 0)
    sampler <- NULL
    if (need_group) {
      gd <- .ph_resample_group(des$group_code, prob_k,
        B = B, seed = sample.int(.Machine$integer.max, 1L),
        n_swap_steps = n_swap_steps
      )
      group_draws <- gd$labels
      sampler <- list(method = gd$method, accept_rate = gd$accept_rate)
    }
    time_draws <- matrix(0L, 0, 0)
    tp <- NULL
    if ("time" %in% terms) {
      tp <- .ph_delta2_time_propensity(samples_meta[samp_k, ], sub_k,
        des$time_code, time_propensity_formula
      )
      time_draws <- .ph_resample_time(sub_k, des$time_code,
        des$n_time_levels, B,
        seed = sample.int(.Machine$integer.max, 1L),
        keep_prob = tp$keep_prob, masks = FALSE
      )$time_code
    }
    adj_weight <- if (!is.null(prob_k)) {
      1 / prob_k[cbind(seq_len(nrow(prob_k)), des$group_code)]
    } else {
      numeric(0)
    }
    # the group term and pairs are adjusted for the propensity covariates
    covariates <- NULL
    if (!is.null(propensity_formula) && need_group) {
      mm <- stats::model.matrix(
        propensity_formula, as.data.frame(subjects_meta[subj_k, ])
      )
      mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
      if (ncol(mm)) covariates <- unname(mm)
    }

    contrasts[[gc]] <- list(
      design = des, terms = terms, pairs = pairs,
      sample_rows = rows_k, sample_subject = sub_k,
      group_draws = group_draws, time_draws = time_draws,
      adj_weight = adj_weight, covariates = covariates,
      rho = rho_vec, rho_hat = rho_res$rho_hat,
      n_overlap = if (!is.null(ps)) ps$diagnostics$n_overlap else NA_integer_
    )
    meta_gc[[gc]] <- list(
      design = des[c(
        "terms", "df", "layer", "group_levels", "time_levels",
        "group_scores", "time_scores"
      )],
      terms_tested = terms,
      propensity = if (!is.null(ps)) ps$diagnostics,
      time_propensity = if (!is.null(tp)) tp$diagnostics,
      dropout = dropout,
      sampler = sampler,
      rho = rho_res$rho,
      rho_hat = rho_res$rho_hat,
      rho_by_stratum = rho_strata
    )
  }

  # --- 4) Master plan: rank-feature strata x group_cols -----------------------
  master_plan <- tidyr::crossing(
    rf_counts,
    group_col = names(contrasts)
  ) |>
    dplyr::filter(.data$n_peptides > 0L)
  master_plan$work_weight <- master_plan$n_peptides *
    vapply(contrasts[master_plan$group_col], function(ct) {
      length(ct$terms) + nrow(ct$pairs)
    }, numeric(1))

  # --- 5) Logging header ------------------------------------------------------
  log_to_file <- isTRUE(log) && (!is.null(log_file) && nzchar(log_file))
  progress_path <- .progress_file(log_file)
  work_total <- sum(as.numeric(master_plan$work_weight)) * B
  if (isTRUE(log)) {
    .progress_init(progress_path)
    hdr_bullets <- c(
      sprintf("strata: %d", nrow(master_plan)),
      sprintf("B per stratum: %d", B),
      sprintf("peptides: %d", length(peptides_order)),
      sprintf("samples: %d", length(prep$samples_order))
    )
    if (log_to_file) {
      .ph_log_info_file(log_file, "Starting compute_delta2",
        bullets = hdr_bullets
      )
    } else {
      .ph_log_info("Starting compute_delta2", bullets = hdr_bullets)
    }
  }

  # --- 6) One stratum ---------------------------------------------------------
  .run_one <- function(i) {
    st <- master_plan[i, , drop = FALSE]
    pep_cols <- pep_cols_by_rf[[st$rank]][[st$feature]]
    if (length(pep_cols) == 0L) {
      return(NULL)
    }
    ct <- contrasts[[st$group_col]]
    des <- ct$design
    res <- cpp_delta2_contrast(
      bitset_raw = prep$bitset_raw, n_words = prep$bitset_words,
      pep_cols = pep_cols, sample_rows = ct$sample_rows,
      sample_subject = ct$sample_subject, group_code = des$group_code,
      time_code = des$time_code, group_scores = des$group_scores,
      time_scores = des$time_scores, terms = ct$terms,
      n_group_levels = des$n_group_levels, n_time_levels = des$n_time_levels,
      rho = ct$rho[pep_cols], group_draws = ct$group_draws,
      time_draws = ct$time_draws, pairs = ct$pairs,
      adj_weight = ct$adj_weight, weight_mode = weight_mode,
      aggregate_stat = aggregate_stat, perm_method = perm_method,
      strat_bins = strat_bins, winsor_z = winsor_z,
      min_m_eff = as.numeric(min_m_eff), covariates = ct$covariates
    )
    if (isTRUE(res$skipped)) {
      return(NULL)
    }
    out <- .ph_delta2_rows(st, ct, res)

    if (isTRUE(log)) {
      msg <- sprintf(
        "%s / %s / %s | n_pep=%d", st$rank, st$feature, st$group_col,
        as.integer(res$n_peptides_used)
      )
      done_w <- .progress_add(progress_path, st$work_weight * B)
      prog <- if (is.finite(done_w) && work_total > 0) {
        sprintf("%.1f%%", pmin(100, 100 * done_w / work_total))
      } else {
        "n/a"
      }
      msg <- paste0(msg, " | progress=", prog)
      if (log_to_file) .ph_log_info_file(log_file, msg) else .ph_log_info(msg)
    }
    out
  }

  # --- 7) Run (shared runner: sequential or future) ---------------------------
  run <- .ph_delta_run_tasks(master_plan$work_weight, .run_one)
  rr <- Filter(Negate(is.null), run$results)

  if (isTRUE(log)) {
    bullets <- c(
      sprintf("strata tested: %d", length(rr)),
      sprintf("B: %d", B),
      sprintf("workers: %d", run$n_workers)
    )
    if (log_to_file) {
      .ph_log_ok_file(log_file, "Done (compute_delta2)", bullets = bullets)
    } else {
      .ph_log_ok("Done (compute_delta2)", bullets = bullets)
    }
  }

  # --- 8) Assemble ------------------------------------------------------------
  terms_tbl <- dplyr::bind_rows(lapply(rr, `[[`, "terms"))
  pair_tbl <- dplyr::bind_rows(lapply(rr, `[[`, "pairwise"))
  if (!nrow(terms_tbl)) {
    terms_tbl <- tibble::tibble(
      rank = character(), feature = character(), group_col = character(),
      design = character(), n_subjects = integer(),
      n_samples_used = integer(), n_peptides_used = integer(),
      m_eff = numeric(), rho_hat = numeric(), n_overlap = integer()
    )
  }
  terms_tbl <- dplyr::arrange(terms_tbl, rank, feature, group_col)
  if (!nrow(pair_tbl)) {
    pair_tbl <- tibble::tibble(
      rank = character(), feature = character(), group_col = character(),
      group1 = character(), group2 = character(), T_obs = numeric(),
      p_perm = numeric(), p_maxT = numeric(), max_delta = numeric(),
      max_delta_adj = numeric(), frac_delta_pos = numeric(),
      frac_delta_pos_w = numeric()
    )
  }
  pair_tbl <- dplyr::arrange(pair_tbl, rank, feature, group_col, group1, group2)

  structure(
    list(
      terms = terms_tbl,
      pairwise = pair_tbl,
      meta = list(
        call = call,
        B_permutations = B,
        group_cols = meta_gc,
        rho_hat = vapply(meta_gc, function(m) m$rho_hat, numeric(1))
      )
    ),
    class = c("delta2_result", "list")
  )
}

# Result rows of one stratum: one $terms row (term columns suffixed) and the
# pairwise rows
.ph_delta2_rows <- function(st, ct, res) {
  des <- ct$design
  terms <- ct$terms
  design <- if (all(terms == "group")) {
    "between"
  } else if (!"group" %in% terms && !"interaction" %in% terms) {
    "within"
  } else {
    "mixed"
  }
  row <- list(
    rank = st$rank, feature = st$feature, group_col = st$group_col,
    design = design,
    n_subjects = as.integer(res$n_subjects),
    n_samples_used = as.integer(res$n_samples),
    n_peptides_used = as.integer(res$n_peptides_used),
    m_eff = as.numeric(res$m_eff),
    rho_hat = as.numeric(ct$rho_hat),
    n_overlap = as.integer(ct$n_overlap)
  )

  # E components: the trend, the second level of two, or every level
  lv_names <- function(levels, scores) {
    levels <- as.character(levels)
    if (length(scores)) {
      "trend"
    } else if (length(levels) == 2L) {
      levels[2]
    } else {
      levels
    }
  }
  g_names <- lv_names(des$group_levels, des$group_scores)
  t_names <- lv_names(des$time_levels, des$time_scores)
  comp_names <- list(
    group = g_names, time = t_names,
    interaction = as.vector(t(outer(g_names, t_names, paste, sep = ":")))
  )

  for (tm in terms) {
    r <- res$terms[[tm]]
    T_obs <- as.numeric(r$T_obs)
    mu <- as.numeric(r$T_null_mean)
    sdv <- as.numeric(r$T_null_sd)
    p <- as.numeric(r$p_perm)
    E <- stats::setNames(as.numeric(r$E), comp_names[[tm]])
    vals <- list(
      df = as.integer(des$df[[tm]]),
      T = T_obs,
      T_null_mean = mu,
      T_null_sd = sdv,
      T_stand = if (!is.na(sdv) && sdv > 0) (T_obs - mu) / sdv else NA_real_,
      p = p,
      Z_from_p = if (is.na(p)) {
        NA_real_
      } else {
        sign(T_obs) * stats::qnorm(p / 2, lower.tail = FALSE)
      },
      b = as.integer(r$b),
      E = list(E),
      max_delta = as.numeric(r$max_delta),
      max_delta_adj = as.numeric(r$max_delta_adj),
      frac_delta_pos = as.numeric(r$frac_delta_pos),
      frac_delta_pos_w = as.numeric(r$frac_delta_pos_w)
    )
    names(vals) <- paste0(names(vals), "_", tm)
    row <- c(row, vals)
  }

  pairwise <- NULL
  if (nrow(ct$pairs)) {
    pr <- res$pairs
    lv <- as.character(des$group_levels)
    pairwise <- tibble::tibble(
      rank = st$rank, feature = st$feature, group_col = st$group_col,
      group1 = lv[ct$pairs[, 1]], group2 = lv[ct$pairs[, 2]],
      T_obs = as.numeric(pr$T_obs), p_perm = as.numeric(pr$p_perm),
      p_maxT = as.numeric(pr$p_maxT), max_delta = as.numeric(pr$max_delta),
      max_delta_adj = as.numeric(pr$max_delta_adj),
      frac_delta_pos = as.numeric(pr$frac_delta_pos),
      frac_delta_pos_w = as.numeric(pr$frac_delta_pos_w)
    )
  }
  list(terms = tibble::as_tibble(row), pairwise = pairwise)
}

# Pooled-prevalence stratum of each peptide (0-based), with the breaks and
# bin rule of compute_delta()
.ph_delta2_prevalence_bins <- function(p0, strat_bins) {
  cuts <- strat_bins[is.finite(strat_bins) & strat_bins > 0 & strat_bins < 1]
  cuts <- sort(unique(cuts))
  if (!length(cuts) || (length(strat_bins) == 1L && strat_bins == 0)) {
    return(integer(length(p0)))
  }
  brks <- c(0, cuts, 1)
  vapply(p0, function(v) {
    for (b in seq_len(length(brks) - 1L)) {
      lo <- brks[b]
      hi <- brks[b + 1L]
      if ((if (b == 1L) v >= lo else v > lo) && v <= hi) {
        return(b - 1L)
      }
    }
    NA_integer_
  }, integer(1))
}

#' @export
print.delta2_result <- function(x, ...) {
  tt <- x$terms
  terms <- intersect(
    c("group", "time", "interaction"),
    sub("^p_", "", grep("^p_", names(tt), value = TRUE))
  )
  cat("<delta2_result>\n")
  cat(sprintf(
    "  %d stratum row(s) x %d group_col(s); B = %d\n",
    nrow(tt), length(unique(tt$group_col)), x$meta$B_permutations
  ))
  cat(sprintf(
    "  terms: %s\n",
    if (length(terms)) paste(terms, collapse = ", ") else "none"
  ))
  for (tm in terms) {
    p <- tt[[paste0("p_", tm)]]
    if (!length(p) || all(is.na(p))) next
    i <- which.min(p)
    cat(sprintf(
      "    %-11s smallest p = %.3g (%s / %s / %s)\n", tm, p[i],
      tt$rank[i], tt$feature[i], tt$group_col[i]
    ))
  }
  if (nrow(x$pairwise)) {
    cat(sprintf(
      "  pairwise: %d row(s); smallest p_maxT = %.3g\n",
      nrow(x$pairwise), min(x$pairwise$p_maxT, na.rm = TRUE)
    ))
  }
  rh <- x$meta$rho_hat
  if (length(rh) && any(!is.na(rh))) {
    cat(sprintf(
      "  rho_hat: %s\n",
      paste(sprintf("%s = %.3f", names(rh), rh), collapse = ", ")
    ))
  }
  invisible(x)
}

#' @title Internal helper: .ph_delta2_prepare
#' @description Build the sample-level presence bitset and the index vectors
#'   that map samples to subjects, timepoints and group levels.
#'
#' @param x A `phip_data` object or a data frame in long format with columns
#'   `sample_id`, `subject_id`, `peptide_id`, `exist_col`, `group_cols` and,
#'   if given, `time_col`.
#' @param group_cols Character vector of grouping columns.
#' @param exist_col Name of the 0/1 presence column. Default `"exist"`.
#' @param time_col Optional name of the timepoint column. If `NULL`, every
#'   sample gets the single timepoint level 1.
#' @param subject_col,sample_col Names of the subject and sample ID columns;
#'   they are returned as `subject_id` and `sample_id`.
#' @param extra_cols Further sample-level columns (e.g. propensity
#'   covariates) to carry into `samples_meta`; each must have one value per
#'   sample.
#'
#' @return A list with
#'   - `samples_order`: sample IDs ordered by `(subject_id, time, sample_id)`;
#'     the bitset rows.
#'   - `sample_subject`: 1-based index into `subjects_order`, per sample.
#'   - `sample_time`: 1-based index into `time_levels`, per sample.
#'   - `sample_group`: named list (one entry per `group_cols` element) of
#'     1-based indices into `group_levels[[gc]]`, per sample.
#'   - `subject_group`: named list as `sample_group`, per subject; `NA` where
#'     the group value varies within the subject (e.g. a timepoint column).
#'   - `subject_n`: number of samples per subject.
#'   - `subjects_order`, `time_levels`, `group_levels`: the level vectors the
#'     indices above point into.
#'   - `samples_meta`: one row per sample in `samples_order`, with
#'     `sample_id`, `subject_id`, `time_col`, `group_cols` and `extra_cols`.
#'   - `peptides_order`: peptide IDs with at least one positive sample; the
#'     bitset columns.
#'   - `bitset_raw`, `bitset_m`, `bitset_words`: packed 64-bit bitset of
#'     `bitset_m` peptide columns with `bitset_words` words each, as returned
#'     by [build_bitset_unpaired()].
#' @keywords internal
.ph_delta2_prepare <- function(x, group_cols, exist_col = "exist",
                               time_col = NULL, subject_col = "subject_id",
                               sample_col = "sample_id", extra_cols = NULL) {
  chk::chk_character(group_cols)
  chk::chk_true(length(group_cols) >= 1)
  chk::chk_string(exist_col)
  if (!is.null(time_col)) chk::chk_string(time_col)
  chk::chk_string(subject_col)
  chk::chk_string(sample_col)
  if (!is.null(extra_cols)) chk::chk_character(extra_cols)

  in_cols <- unique(c(
    sample_col, subject_col, time_col, group_cols, extra_cols,
    "peptide_id", exist_col
  ))
  sample_cols <- unique(c(
    "sample_id", "subject_id", time_col, group_cols, extra_cols
  ))

  if (inherits(x, "phip_data")) {
    df_long <- x$data_long |>
      dplyr::select(dplyr::any_of(in_cols))
  } else {
    chk::chk_data(x)
    miss <- setdiff(in_cols, colnames(x))
    if (length(miss)) {
      .ph_abort("Missing required columns",
        bullets = paste("-", miss)
      )
    }
    df_long <- tibble::as_tibble(x) |>
      dplyr::select(dplyr::any_of(in_cols))
  }
  if (sample_col != "sample_id") {
    df_long <- dplyr::rename(df_long, sample_id = !!rlang::sym(sample_col))
  }
  if (subject_col != "subject_id") {
    df_long <- dplyr::rename(df_long, subject_id = !!rlang::sym(subject_col))
  }

  # --- SAMPLE TABLE: one row per sample ---------------------------------------
  # No hits guard as in compute_delta(): rows are samples, so several positive
  # samples of one subject are kept apart instead of rejected.
  samples <- df_long |>
    dplyr::distinct(!!!rlang::syms(sample_cols)) |>
    dplyr::collect()

  if (anyNA(samples$sample_id) || anyNA(samples$subject_id)) {
    .ph_abort("sample_id and subject_id must not contain NA.")
  }
  dup_samples <- unique(samples$sample_id[duplicated(samples$sample_id)])
  if (length(dup_samples)) {
    .ph_abort(
      sprintf(
        "Each sample_id must have a single subject_id, %s value.",
        paste(c(time_col, group_cols, extra_cols), collapse = ", ")
      ),
      bullets = paste("- sample_id =", utils::head(dup_samples, 10))
    )
  }

  # --- FREEZE ORDERS ----------------------------------------------------------
  # Same arrange() as compute_delta(), so subjects (and, with one sample per
  # subject, samples) come out in the same order.
  subjects_order <- samples |>
    dplyr::distinct(subject_id) |>
    dplyr::arrange(subject_id) |>
    dplyr::pull(subject_id)

  if (is.null(time_col)) {
    time_levels <- 1L
    samples$.time <- 1L
  } else {
    time_levels <- samples |>
      dplyr::distinct(!!rlang::sym(time_col)) |>
      dplyr::arrange(!!rlang::sym(time_col)) |>
      dplyr::pull(!!rlang::sym(time_col))
    samples$.time <- match(samples[[time_col]], time_levels)
  }

  samples <- samples |>
    dplyr::arrange(subject_id, .data$.time, sample_id)

  samples_order <- samples$sample_id
  sample_subject <- match(samples$subject_id, subjects_order)
  sample_time <- as.integer(samples$.time)
  N_samples <- length(samples_order)
  N_subjects <- length(subjects_order)
  subject_n <- tabulate(sample_subject, nbins = N_subjects)

  group_levels <- list()
  sample_group <- list()
  subject_group <- list()
  for (gc in group_cols) {
    lv <- samples |>
      dplyr::distinct(!!rlang::sym(gc)) |>
      dplyr::arrange(!!rlang::sym(gc)) |>
      dplyr::pull(!!rlang::sym(gc))
    idx <- match(samples[[gc]], lv)
    # one level per subject, NA if the column varies within the subject
    lo <- vapply(split(idx, sample_subject), min, integer(1))
    hi <- vapply(split(idx, sample_subject), max, integer(1))
    group_levels[[gc]] <- lv
    sample_group[[gc]] <- idx
    subject_group[[gc]] <- unname(dplyr::if_else(lo == hi, lo, NA_integer_))
  }

  peptides_order <- df_long |>
    dplyr::distinct(peptide_id) |>
    dplyr::arrange(peptide_id) |>
    dplyr::collect() |>
    dplyr::pull(peptide_id)

  # --- BUILD HITS FROM DISTINCT POSITIVE SAMPLES ------------------------------
  pos_pairs <- df_long |>
    dplyr::filter(!!rlang::sym(exist_col) > 0L) |>
    dplyr::distinct(sample_id, peptide_id) |>
    dplyr::collect()

  samp_index <- match(pos_pairs$sample_id, samples_order)
  pep_index <- match(pos_pairs$peptide_id, peptides_order)

  hits_split <- split(samp_index, pep_index)

  nonempty_pep_ids <- as.integer(names(hits_split))
  if (length(nonempty_pep_ids) == 0L) {
    .ph_abort("All peptides are zero (no positives).")
  }
  peptides_order <- peptides_order[nonempty_pep_ids]
  hits_by_peptide <- lapply(unname(hits_split), function(i) sort(as.integer(i)))

  # --- BUILD 64-bit BITSET ONCE -----------------------------------------------
  bs <- build_bitset_unpaired(hits_by_peptide, N_samples)
  if (bs$m != length(peptides_order)) {
    .ph_abort("Bitset/peptide dimension mismatch.")
  }

  list(
    samples_order  = samples_order,
    sample_subject = sample_subject,
    sample_time    = sample_time,
    sample_group   = sample_group,
    subject_group  = subject_group,
    subject_n      = subject_n,
    subjects_order = subjects_order,
    time_levels    = time_levels,
    group_levels   = group_levels,
    peptides_order = peptides_order,
    samples_meta   = dplyr::select(samples, -dplyr::all_of(".time")),
    bitset_raw     = bs$data,
    bitset_m       = bs$m,
    bitset_words   = bs$n_words
  )
}
