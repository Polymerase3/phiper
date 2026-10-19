# Shared helpers for the compute_delta2() validation (see ../README.md).

# --- package and parallelism ---------------------------------------------------

# phiper from PHIPER_SOURCE (a source checkout, via pkgload) or the library
val_load_phiper <- function() {
  src <- Sys.getenv("PHIPER_SOURCE", "")
  if (nzchar(src)) {
    pkgload::load_all(src, quiet = TRUE, export_all = FALSE)
  } else {
    suppressPackageStartupMessages(library(phiper))
  }
  invisible(TRUE)
}

# one result per replicate; replicate r uses seed `seed + r`
val_reps <- function(n_reps, fun, seed, cores = val_cores()) {
  one <- function(r) {
    set.seed(seed + r)
    out <- fun(r)
    out$rep <- r
    out
  }
  res <- if (cores > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(seq_len(n_reps), one,
      mc.cores = cores,
      mc.preschedule = FALSE
    )
  } else {
    lapply(seq_len(n_reps), one)
  }
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop("replicate failed: ", as.character(res[[which(bad)[1]]]))
  dplyr::bind_rows(res)
}

val_cores <- function() {
  as.integer(Sys.getenv("VAL_CORES", max(1L, parallel::detectCores() - 0L)))
}

# Monte Carlo band of a rejection rate at level alpha over n replicates
val_band <- function(n, alpha = 0.05, z = 1.96) {
  half <- z * sqrt(alpha * (1 - alpha) / n)
  c(lo = alpha - half, hi = alpha + half)
}

# --- simulator -------------------------------------------------------------------

# baseline prevalences: the phiper example data (phip_mixture, 1950 peptides)
# or a sparse profile, log-uniform on [0.01, 0.3]
val_prevalence <- function(m, profile = c("example", "sparse")) {
  profile <- match.arg(profile)
  if (profile == "sparse") {
    return(exp(stats::runif(m, log(0.01), log(0.3))))
  }
  path <- file.path(val_dir(), "data", "phip_mixture_prevalence.csv")
  p <- utils::read.csv(path)$prevalence
  p <- p[p > 0 & p < 0.95]
  sample(p, m, replace = TRUE)
}

val_dir <- function() {
  d <- Sys.getenv("VAL_DIR", "")
  if (nzchar(d)) d else "."
}

#' Long PhIP-seq-like data for compute_delta2().
#'
#' Subjects carry a covariate `age` (standardised). Group labels come from
#' fixed `group_sizes`, or, with `confound != 0`, from a multinomial logit in
#' `age` (group k's logit is `confound * (k - 1) * age`). Each subject-peptide
#' has a persistent state copied into a sample with probability
#' `sqrt(rho)` (so the within-subject correlation of repeat samples is
#' `rho`); otherwise a sample is drawn afresh from
#' `logit p = logit(prev_j) + beta_age * age + effects`.
#'
#' @param visits Number of complete timepoints, or a function
#'   `(group, subject) -> visit indices`.
#' @param group_effect,time_effect,interaction Logit shifts: group effect on
#'   the affected peptides in the last group; time effect per visit step on
#'   all peptides; extra per-step change in the last group on the affected
#'   peptides.
#' @param dropout Optional informative dropout: subjects leave after the
#'   first visit whose number of positives is above this quantile of all
#'   samples' counts.
val_sim <- function(n_subjects = 60, K = 2, group_sizes = NULL,
                    confound = 0, beta_age = 0, visits = 1L, m = 100,
                    prev = NULL, profile = "example", rho = 0.3,
                    group_effect = 0, frac_affected = 0.2, time_effect = 0,
                    interaction = 0, dropout = NULL) {
  age <- stats::rnorm(n_subjects)
  if (confound != 0) {
    eta <- outer(age, confound * (seq_len(K) - 1))
    pr <- exp(eta) / rowSums(exp(eta))
    grp <- apply(pr, 1, function(q) sample.int(K, 1, prob = q))
  } else {
    if (is.null(group_sizes)) {
      group_sizes <- rep(n_subjects %/% K, K)
      group_sizes[K] <- n_subjects - sum(group_sizes[-K])
    }
    grp <- sample(rep(seq_len(K), group_sizes))
    n_subjects <- length(grp)
    age <- age[seq_len(n_subjects)]
  }
  if (is.null(prev)) prev <- val_prevalence(m, profile)
  affected <- seq_len(m) <= round(frac_affected * m)

  vis <- if (is.function(visits)) {
    lapply(seq_len(n_subjects), function(i) visits(grp[i], i))
  } else {
    rep(list(seq_len(visits)), n_subjects)
  }
  n_v <- lengths(vis)
  subj <- rep(seq_len(n_subjects), n_v)
  tp <- unlist(vis)
  last <- grp[subj] == K

  base <- outer(beta_age * age, rep(1, m)) +
    matrix(stats::qlogis(prev), n_subjects, m, byrow = TRUE)
  base <- base + outer(grp == K, affected * group_effect)
  state <- matrix(stats::rbinom(n_subjects * m, 1, stats::plogis(base)), n_subjects)

  eta <- base[subj, , drop = FALSE] +
    time_effect * (tp - 1) +
    outer(last * (tp - 1), affected * interaction)
  fresh <- matrix(stats::rbinom(length(eta), 1, stats::plogis(eta)), nrow(eta))
  copy <- matrix(stats::runif(length(eta)) < sqrt(rho), nrow(eta))
  Y <- ifelse(copy, state[subj, , drop = FALSE], fresh)

  keep <- rep(TRUE, length(subj))
  if (!is.null(dropout)) {
    pos <- rowSums(Y)
    thr <- stats::quantile(pos, dropout)
    for (i in seq_len(n_subjects)) {
      r <- which(subj == i)
      r <- r[order(tp[r])]
      hit <- which(pos[r] > thr)
      if (length(hit) && hit[1] < length(r)) keep[r[(hit[1] + 1):length(r)]] <- FALSE
    }
  }

  samples <- data.frame(
    sample_id = sprintf("s%03d_v%d", subj, tp),
    subject_id = sprintf("s%03d", subj),
    group = LETTERS[grp[subj]],
    timepoint = sprintf("T%d", tp),
    age = age[subj],
    stringsAsFactors = FALSE
  )[keep, ]
  Y <- Y[keep, , drop = FALSE]
  peps <- sprintf("p%03d", seq_len(m))
  out <- samples[rep(seq_len(nrow(samples)), times = m), ]
  out$peptide_id <- rep(peps, each = nrow(samples))
  out$exist <- as.integer(Y)
  rownames(out) <- NULL
  attr(out, "Y") <- Y
  attr(out, "samples") <- samples
  out
}

# all peptides in one stratum
val_lib <- function(df) {
  data.frame(peptide_id = unique(df$peptide_id), set = "all")
}

# compute_delta2() on one pooled stratum; returns one row with the p-value of
# each term (p_group, p_time, p_interaction) and any warnings raised
val_run <- function(df, B = 500L, time_col = "timepoint", ...) {
  warn <- character()
  res <- withCallingHandlers(
    compute_delta2(df,
      rank_cols = "set", group_cols = "group", time_col = time_col,
      peptide_library = val_lib(df), B_permutations = B, pairwise = FALSE,
      ...
    ),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  tt <- res$terms
  out <- tibble::as_tibble(tt[, intersect(
    c("p_group", "p_time", "p_interaction"),
    names(tt)
  )])
  out$warnings <- paste(unique(warn), collapse = " | ")
  out
}

# The interaction term of one pooled stratum straight from the kernel, with
# or without its time main-effect residualisation (compute_delta2() always
# residualises). Both settings see the same group draws.
val_run_interaction <- function(df, B = 500L) {
  ns <- asNamespace("phiper")
  prep <- ns$.ph_delta2_prepare(df, "group", time_col = "timepoint")
  sm <- tibble::tibble(
    subject_id = prep$subjects_order,
    group = prep$group_levels$group[prep$subject_group$group]
  )
  des <- ns$.ph_delta2_design(sm, prep$samples_meta, "group",
    time_col = "timepoint"
  )
  ss <- prep$sample_subject
  cols <- seq_len(prep$bitset_m)
  est <- ns$cpp_delta2_rho(
    prep$bitset_raw, prep$bitset_words, cols, seq_along(ss), ss
  )
  rho <- if (is.na(est$rho_hat)) 0 else est$rho_hat
  draws <- ns$.ph_resample_group(des$group_code, NULL,
    B = B, seed = sample.int(.Machine$integer.max, 1L)
  )$labels
  p <- vapply(c(TRUE, FALSE), function(resid) {
    ns$cpp_delta2_contrast(
      prep$bitset_raw, prep$bitset_words, cols, seq_along(ss), ss,
      des$group_code, des$time_code, des$group_scores, des$time_scores,
      "interaction", des$n_group_levels, des$n_time_levels, rho, draws,
      matrix(0L, 0, 0), matrix(0L, 0, 2), numeric(0), "n_eff_sqrt",
      "stouffer", "mid_p", c(0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50),
      4, 0, resid
    )$terms$interaction$p_perm
  }, numeric(1))
  tibble::tibble(p_residualised = p[1], p_raw = p[2])
}
