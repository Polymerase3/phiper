# tests/testthat/test-delta2_compute.R

# --- reference pooling, written from compute_delta()'s rules -----------------
.ref_pool <- function(z, w, bins, stat) {
  af_one <- function(p) {
    p <- sort(p)
    max(cumsum(-log(p)) / cumsum(1 / seq_along(p)))
  }
  vals <- vapply(sort(unique(bins)), function(b) {
    zb <- z[bins == b]
    wb <- w[bins == b]
    if (stat == "stouffer") {
      sum(wb * zb) / sqrt(max(sum(wb^2), 1e-300))
    } else if (stat == "maxmean") {
      sp <- if (any(zb > 0)) mean(zb[zb > 0]) else 0
      sn <- if (any(zb < 0)) mean(-zb[zb < 0]) else 0
      if (sp >= sn) sp else -sn
    } else {
      tp <- af_one(pmax(stats::pnorm(-zb), 1e-300))
      tn <- af_one(pmax(stats::pnorm(zb), 1e-300))
      if (tp >= tn) tp else -tn
    }
  }, numeric(1))
  mean(vals)
}

.ref_T <- function(r, w, bins, stat, winsor) {
  r <- as.matrix(r)
  r[] <- pmax(-winsor, pmin(winsor, r))
  E <- apply(r, 2, .ref_pool, w = w, bins = bins, stat = stat)
  list(E = E, T = if (length(E) == 1L) E else sum(E^2))
}

.ref_p <- function(Tb, T_obs, mid_p = TRUE) {
  tol <- 1e-9 * pmax(1, pmax(abs(Tb), abs(T_obs)))
  gt <- sum(abs(Tb) - abs(T_obs) > tol)
  eq <- sum(abs(abs(Tb) - abs(T_obs)) <= tol)
  B <- length(Tb)
  if (mid_p) max((gt + 0.5 * eq) / B, 0.5 / B) else (1 + gt + eq) / (1 + B)
}

# 3 groups x 3 timepoints, unbalanced visits, several samples per visit
.d2_panel <- function(seed = 1, n = 36, m = 12) {
  set.seed(seed)
  rows <- lapply(seq_len(n), function(i) {
    tp <- sort(sample(1:3, sample(1:3, 1)))
    tp <- c(tp, if (runif(1) < 0.2) tp[1])
    data.frame(
      sample_id = sprintf("s%02d_%d", i, seq_along(tp)),
      subject_id = sprintf("s%02d", i),
      timepoint = sprintf("T%d", tp),
      group = c("a", "b", "c")[(i - 1) %% 3 + 1],
      stringsAsFactors = FALSE
    )
  })
  samples <- do.call(rbind, rows)
  df <- merge(samples, data.frame(peptide_id = sprintf("p%02d", 1:m)))
  prev <- stats::runif(m, 0.05, 0.6)
  pep <- match(df$peptide_id, sprintf("p%02d", 1:m))
  df$exist <- stats::rbinom(nrow(df), 1, prev[pep])
  df
}

test_that("kernel null equals the engine on relabelled data, pooled", {
  df <- .d2_panel()
  prep <- .ph_delta2_prepare(df, "group", time_col = "timepoint")
  sm <- tibble::tibble(
    subject_id = prep$subjects_order,
    group = prep$group_levels$group[prep$subject_group$group]
  )
  ss <- prep$sample_subject
  cols <- seq_len(prep$bitset_m)
  B <- 12L
  sb <- c(0.1, 0.3)

  set.seed(9)
  age <- cbind(age = stats::rnorm(length(prep$subjects_order)))
  for (cfg in list(
    list(gs = NULL, ts = NULL, stat = "stouffer", wm = "n_eff_sqrt"),
    list(gs = NULL, ts = NULL, stat = "stouffer", wm = "equal", cov = age),
    list(
      gs = c(a = 0, b = 1, c = 3), ts = NULL, stat = "maxmean", wm = "equal"
    ),
    list(
      gs = NULL, ts = c(T1 = 0, T2 = 1, T3 = 2), stat = "af", wm = "se_invvar"
    )
  )) {
    des <- .ph_delta2_design(sm, prep$samples_meta, "group",
      time_col = "timepoint", group_scores = cfg$gs, time_scores = cfg$ts
    )
    gd <- .ph_resample_group(des$group_code, NULL, B = B, seed = 5)$labels
    td <- .ph_resample_time(ss, des$time_code, 3L, B,
      seed = 6, masks = FALSE
    )$time_code
    pairs <- t(utils::combn(3L, 2L))
    rho <- 0.3
    k <- cpp_delta2_contrast(
      prep$bitset_raw, prep$bitset_words, cols, seq_along(ss), ss,
      des$group_code, des$time_code, des$group_scores, des$time_scores,
      des$terms, 3L, 3L, rho, gd, td, pairs, numeric(0),
      cfg$wm, cfg$stat, "mid_p", sb, 3, 0, cfg$cov
    )

    engine <- function(gc, tc, terms) {
      cpp_delta2_score(
        prep$bitset_raw, prep$bitset_words, cols, seq_along(ss), ss,
        gc, tc, des$group_scores, des$time_scores, terms, 3L, 3L, rho, cfg$cov
      )
    }
    obs <- engine(des$group_code, des$time_code, des$terms)
    p0 <- obs$tot / obs$n_samples
    bins <- .ph_delta2_prevalence_bins(p0, sb)
    w <- switch(cfg$wm,
      equal = rep(1, length(p0)),
      n_eff_sqrt = sqrt(pmax(obs$tot, 1e-12)),
      se_invvar = ifelse(p0 > 0 & p0 < 1, 1 / sqrt(p0 * (1 - p0)), 0)
    )
    pair_z <- function(gc) {
      sapply(seq_len(nrow(pairs)), function(p) {
        code <- match(gc, pairs[p, ])
        z <- cpp_delta2_score(
          prep$bitset_raw, prep$bitset_words, cols, seq_along(ss), ss,
          code, des$time_code, numeric(0), numeric(0), "group", 2L, 3L, rho,
          cfg$cov
        )$group
        .ref_T(z, w, bins, cfg$stat, 3)$T
      })
    }

    for (tm in des$terms) {
      ref_obs <- .ref_T(obs[[tm]], w, bins, cfg$stat, 3)
      expect_equal(k$terms[[tm]]$T_obs, ref_obs$T, tolerance = 1e-10)
      expect_equal(k$terms[[tm]]$E, unname(ref_obs$E), tolerance = 1e-10)
      Tb <- vapply(seq_len(B), function(b) {
        e <- if (tm == "time") {
          engine(des$group_code, td[, b], "time")
        } else {
          engine(gd[, b], des$time_code, tm)
        }
        .ref_T(e[[tm]], w, bins, cfg$stat, 3)$T
      }, numeric(1))
      expect_equal(k$terms[[tm]]$T_null_mean, mean(Tb), tolerance = 1e-10)
      expect_equal(k$terms[[tm]]$T_null_sd, stats::sd(Tb), tolerance = 1e-10)
      expect_equal(k$terms[[tm]]$p_perm, .ref_p(Tb, ref_obs$T))
    }

    pz_obs <- pair_z(des$group_code)
    pz <- sapply(seq_len(B), function(b) pair_z(gd[, b]))
    expect_equal(k$pairs$T_obs, pz_obs, tolerance = 1e-10)
    for (p in seq_len(nrow(pairs))) {
      expect_equal(k$pairs$p_perm[p], .ref_p(pz[p, ], pz_obs[p]))
      expect_equal(k$pairs$p_maxT[p], .ref_p(apply(abs(pz), 2, max), pz_obs[p]))
    }
  }
})

.d2_cross <- function(seed, N = 60, m = 15) {
  set.seed(seed)
  samples <- data.frame(
    sample_id = sprintf("x%02d", 1:N), subject_id = sprintf("s%02d", 1:N),
    group = rep(c("A", "B"), c(26, N - 26)), stringsAsFactors = FALSE
  )
  df <- merge(samples, data.frame(peptide_id = sprintf("p%02d", 1:m)))
  hit <- df$group == "B" & df$peptide_id %in% sprintf("p%02d", 1:5)
  shift <- ifelse(hit, 0.25, 0)
  df$exist <- stats::rbinom(nrow(df), 1, 0.2 + shift)
  lib <- data.frame(
    peptide_id = sprintf("p%02d", 1:m),
    species = rep(c("sp1", "sp2", "sp3"), each = 5)
  )
  list(df = df, lib = lib, N = N)
}

test_that("one sample per subject, two groups: v1 up to sqrt((N-1)/N)", {
  skip_on_cran()
  d <- .d2_cross(31)
  for (wm in c("n_eff_sqrt", "equal")) {
    for (agg in c("stouffer", "maxmean")) {
      args <- list(
        x = d$df, rank_cols = "species", group_cols = "group",
        peptide_library = d$lib, B_permutations = 1000L,
        weight_mode = wm, aggregate_stat = agg, winsor_z = 1e6
      )
      set.seed(1)
      v1 <- do.call(compute_delta, c(args, stat_mode = "score"))
      set.seed(1)
      v2 <- do.call(compute_delta2, args)
      tt <- v2$terms

      expect_identical(tt$feature, v1$feature)
      f <- sqrt((d$N - 1) / d$N)
      expect_equal(tt$T_group, v1$T_obs * f, tolerance = 1e-10)
      expect_equal(tt$max_delta_group, v1$max_delta, tolerance = 1e-12)
      expect_equal(tt$frac_delta_pos_group, v1$frac_delta_pos)
      expect_equal(tt$frac_delta_pos_w_group, v1$frac_delta_pos_w)
      expect_equal(tt$m_eff, v1$m_eff, tolerance = 1e-10)
      # different draws, same null: p-values agree up to Monte Carlo error
      expect_lt(max(abs(tt$p_group - v1$p_perm)), 0.06)
    }
  }
  expect_identical(unique(tt$design), "between")
  expect_true(all(is.na(tt$rho_hat)))
  expect_false(any(grepl("_time$|_interaction$", names(tt))))
  # two groups: the pairwise test is the group term
  expect_equal(v2$pairwise$p_perm, tt$p_group)
  expect_equal(v2$pairwise$p_maxT, tt$p_group)
})

test_that("two complete timepoints: the time term is v1's paired McNemar", {
  skip_on_cran()
  set.seed(32)
  N <- 40
  samples <- data.frame(
    sample_id = c(sprintf("s%02d_T1", 1:N), sprintf("s%02d_T2", 1:N)),
    subject_id = rep(sprintf("s%02d", 1:N), 2),
    timepoint = rep(c("T1", "T2"), each = N), stringsAsFactors = FALSE
  )
  df <- merge(samples, data.frame(peptide_id = sprintf("p%02d", 1:12)))
  prev <- ifelse(df$timepoint == "T2", 0.35, 0.2)
  df$exist <- stats::rbinom(nrow(df), 1, prev)
  df$cohort <- "all"
  lib <- data.frame(peptide_id = sprintf("p%02d", 1:12), species = "sp")

  v1 <- compute_delta(df,
    rank_cols = "species", group_cols = "timepoint",
    paired_by = "subject_id", stat_mode = "mcnemar", winsor_z = 1e6,
    weight_mode = "n_eff_sqrt", peptide_library = lib, B_permutations = 1000L
  )
  v2 <- compute_delta2(df,
    rank_cols = "species", group_cols = "cohort",
    time_col = "timepoint", winsor_z = 1e6,
    peptide_library = lib, B_permutations = 1000L
  )
  tt <- v2$terms
  expect_identical(tt$design, "within")
  # v1's McNemar z counts g1-only minus g2-only; v2 counts gains over time
  expect_equal(tt$T_time, -v1$T_obs, tolerance = 1e-10)
  expect_equal(tt$max_delta_time, v1$max_delta, tolerance = 1e-12)
  expect_equal(tt$frac_delta_pos_time, v1$frac_delta_pos, tolerance = 1e-12)
  expect_lt(abs(tt$p_time - v1$p_perm), 0.06)
  expect_false(any(grepl("_group$|_interaction$", names(tt))))
  expect_identical(nrow(v2$pairwise), 0L)
})

test_that("three groups with time: all terms, profiles and pairwise max-T", {
  df <- .d2_panel(seed = 3, n = 45, m = 10)
  df$species <- "sp"
  lib <- data.frame(peptide_id = unique(df$peptide_id), species = "sp")
  set.seed(4)
  res <- compute_delta2(df,
    rank_cols = c("peptide_id", "species"), group_cols = "group",
    time_col = "timepoint", peptide_library = lib, B_permutations = 200L
  )
  expect_s3_class(res, "delta2_result")
  tt <- res$terms
  expect_identical(nrow(tt), 11L)
  expect_identical(unique(tt$df_group), 2L)
  expect_identical(unique(tt$df_time), 2L)
  expect_identical(unique(tt$df_interaction), 4L)
  E <- tt$E_group[[1]]
  expect_identical(names(E), c("a", "b", "c"))
  expect_equal(tt$T_group, vapply(tt$E_group, function(e) sum(e^2), 1))
  expect_length(tt$E_interaction[[1]], 9L)
  expect_true(all(is.na(tt$frac_delta_pos_group)))
  expect_true(all(tt$p_group > 0 & tt$p_group <= 1))

  pw <- res$pairwise
  expect_identical(nrow(pw), 33L)
  expect_true(all(pw$p_maxT >= pw$p_perm - 1e-12))
  expect_identical(unique(paste(pw$group1, pw$group2)), c("a b", "a c", "b c"))
  expect_false(is.na(res$meta$rho_hat[["group"]]))
  expect_output(print(res), "delta2_result")

  # same seed, same result
  set.seed(4)
  again <- compute_delta2(df,
    rank_cols = c("peptide_id", "species"), group_cols = "group",
    time_col = "timepoint", peptide_library = lib, B_permutations = 200L
  )
  expect_identical(again$terms, res$terms)

  # test_terms restricts the terms and their columns
  only <- compute_delta2(df,
    rank_cols = "species", group_cols = "group", time_col = "timepoint",
    peptide_library = lib, test_terms = "time", B_permutations = 100L
  )
  expect_true(any(grepl("_time$", names(only$terms))))
  expect_false(any(grepl("_group$|_interaction$", names(only$terms))))
  expect_identical(nrow(only$pairwise), 0L)
})

test_that("inverse-propensity effect sizes and propensity diagnostics", {
  set.seed(33)
  N <- 120
  age <- stats::rnorm(N, 50, 10)
  grp <- ifelse(stats::runif(N) < stats::plogis(0.08 * (age - 50)), "B", "A")
  samples <- data.frame(
    sample_id = sprintf("x%03d", 1:N), subject_id = sprintf("s%03d", 1:N),
    group = grp, age = age, stringsAsFactors = FALSE
  )
  df <- merge(samples, data.frame(peptide_id = sprintf("p%d", 1:6)))
  # positivity driven by age only: confounded raw effect, no adjusted one
  prev <- stats::plogis(-1 + 0.06 * (df$age - 50))
  df$exist <- stats::rbinom(nrow(df), 1, prev)

  set.seed(1)
  raw <- compute_delta2(df, "peptide_id", "group", B_permutations = 100L)
  expect_equal(raw$terms$max_delta_adj_group, raw$terms$max_delta_group)

  set.seed(1)
  adj <- compute_delta2(df, "peptide_id", "group",
    propensity_formula = ~age, B_permutations = 100L
  )
  ps <- adj$meta$group_cols$group$propensity
  expect_identical(ps$n, as.integer(N))
  expect_identical(adj$meta$group_cols$group$sampler$method, "rejection")
  expect_identical(adj$terms$n_overlap[1], ps$n_overlap)

  # adjusted group difference, by hand, for one peptide
  fit <- stats::glm(I(grp == "B") ~ age, family = stats::binomial())
  p <- stats::fitted(fit)
  p1 <- df[df$peptide_id == "p1", ]
  y <- p1$exist[match(samples$sample_id, p1$sample_id)]
  wB <- 1 / p[grp == "B"]
  wA <- 1 / (1 - p[grp == "A"])
  d_adj <- sum(wB * y[grp == "B"]) / sum(wB) - sum(wA * y[grp == "A"]) / sum(wA)
  d_raw <- mean(y[grp == "B"]) - mean(y[grp == "A"])
  row <- adj$terms[adj$terms$feature == "p1", ]
  expect_equal(row$max_delta_adj_group, abs(d_adj), tolerance = 1e-8)
  expect_equal(row$max_delta_group, abs(d_raw), tolerance = 1e-12)
  # adjustment pulls the confounded differences towards 0
  tt <- adj$terms
  expect_lt(mean(tt$max_delta_adj_group), mean(tt$max_delta_group))

  # the group statistic is adjusted for the propensity covariates: one
  # peptide per stratum, so T_group is the engine's z with age as covariate
  prep <- .ph_delta2_prepare(df, "group")
  sub <- match(prep$subjects_order, samples$subject_id)
  z <- cpp_delta2_score(
    prep$bitset_raw, prep$bitset_words, seq_len(prep$bitset_m),
    seq_along(prep$sample_subject), prep$sample_subject,
    match(grp[sub], c("A", "B")), rep(1L, length(prep$sample_subject)),
    numeric(0), numeric(0), "group", 2L, 1L,
    adj$meta$group_cols$group$rho, cbind(age = age[sub])
  )$group[, 1]
  z <- pmax(-4, pmin(4, z))[match(tt$feature, prep$peptides_order)]
  expect_equal(tt$T_group, z, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(tt$T_group, raw$terms$T_group)))
})

test_that("a timepoint missing from a group triggers a warning", {
  df <- .d2_panel(seed = 6, n = 30, m = 6)
  df <- df[!(df$group == "a" & df$timepoint == "T3"), ]
  expect_warning(
    compute_delta2(df, "peptide_id", "group",
      time_col = "timepoint", B_permutations = 100L, pairwise = FALSE,
      rank_feature_keep = list(peptide_id = "p01")
    ),
    "not estimable"
  )
  full <- .d2_panel(seed = 6, n = 30, m = 6)
  expect_no_warning(
    compute_delta2(full, "peptide_id", "group",
      time_col = "timepoint", B_permutations = 100L, pairwise = FALSE,
      rank_feature_keep = list(peptide_id = "p01")
    )
  )
})

test_that("rho options, custom ID columns, and input checks", {
  df <- .d2_panel(seed = 5, n = 30, m = 8)
  names(df)[names(df) == "subject_id"] <- "patient"
  names(df)[names(df) == "sample_id"] <- "specimen"
  base <- list(
    x = df, rank_cols = "peptide_id", group_cols = "group",
    subject_col = "patient", sample_col = "specimen",
    time_col = "timepoint", B_permutations = 100L
  )
  set.seed(1)
  r1 <- do.call(compute_delta2, c(base, rho = "one"))
  expect_true(all(is.na(r1$terms$rho_hat)))
  expect_identical(r1$meta$group_cols$group$rho, 1)
  set.seed(1)
  r0 <- do.call(compute_delta2, c(base, rho = 0.25))
  expect_identical(r0$meta$group_cols$group$rho, 0.25)
  # rho only weights the group term
  expect_equal(r1$terms$T_time, r0$terms$T_time)

  rs <- do.call(compute_delta2, c(base, rho_by_stratum = TRUE))
  strata <- rs$meta$group_cols$group$rho_by_stratum
  expect_s3_class(strata, "tbl_df")
  expect_identical(sum(strata$n_peptides), 8L)

  # min_samples_per_subject drops single-sample subjects
  ms <- do.call(compute_delta2, c(base, min_samples_per_subject = 2L))
  n_i <- table(df$patient[!duplicated(df$specimen)])
  expect_identical(ms$terms$n_subjects[1], sum(n_i >= 2))

  bad <- df
  bad$group <- paste0(bad$group, bad$timepoint)
  expect_error(do.call(compute_delta2, utils::modifyList(base, list(x = bad))),
    "constant"
  )
  expect_error(do.call(compute_delta2, c(base, B_permutations = 10L)))
  expect_error(
    compute_delta2(df, "peptide_id", "group",
      subject_col = "patient", sample_col = "specimen",
      test_terms = "time"
    ),
    "time_col"
  )

  # "auto": a many-valued numeric group reads as continuous
  num <- df
  num$group <- as.numeric(factor(num$patient))
  expect_error(do.call(compute_delta2, utils::modifyList(base, list(x = num))),
    "continuous"
  )
})
