# tests/testthat/test-delta2_standardize.R

.std_pool <- function(z, w, bins) {
  vals <- vapply(sort(unique(bins)), function(b) {
    sum(w[bins == b] * z[bins == b]) / sqrt(sum(w[bins == b]^2))
  }, numeric(1))
  mean(vals)
}

# 2 groups x 2 timepoints, some single-visit subjects
.std_panel <- function(seed = 1, n = 30, m = 10) {
  set.seed(seed)
  rows <- lapply(seq_len(n), function(i) {
    tp <- if (stats::runif(1) < 0.2) sample(1:2, 1) else 1:2
    data.frame(
      sample_id = sprintf("s%02d_%d", i, tp),
      subject_id = sprintf("s%02d", i),
      timepoint = sprintf("T%d", tp),
      group = c("a", "b")[(i - 1) %% 2 + 1],
      stringsAsFactors = FALSE
    )
  })
  samples <- do.call(rbind, rows)
  df <- merge(samples, data.frame(peptide_id = sprintf("p%02d", 1:m)))
  prev <- stats::runif(m, 0.1, 0.5)
  df$exist <- stats::rbinom(nrow(df), 1, prev[match(
    df$peptide_id, sprintf("p%02d", 1:m)
  )])
  df
}

test_that("two_pass standardises each peptide over the observed and draws", {
  df <- .std_panel()
  prep <- .ph_delta2_prepare(df, "group", time_col = "timepoint")
  sm <- tibble::tibble(
    subject_id = prep$subjects_order,
    group = prep$group_levels$group[prep$subject_group$group]
  )
  ss <- prep$sample_subject
  cols <- seq_len(prep$bitset_m)
  des <- .ph_delta2_design(sm, prep$samples_meta, "group",
    time_col = "timepoint"
  )
  B <- 15L
  # biased labels: propensities that are far from uniform
  prob <- stats::runif(length(des$group_code), 0.2, 0.8)
  gd <- .ph_resample_group(des$group_code, prob, B = B, seed = 2)$labels
  td <- .ph_resample_time(ss, des$time_code, 2L, B,
    seed = 3, masks = FALSE
  )$time_code
  pairs <- matrix(c(1L, 2L), 1)
  winsor <- 2.5
  run <- function(std) {
    cpp_delta2_contrast(
      prep$bitset_raw, prep$bitset_words, cols, seq_along(ss), ss,
      des$group_code, des$time_code, numeric(0), numeric(0), des$terms,
      2L, 2L, 0.2, gd, td, pairs, numeric(0), "equal", "stouffer",
      "mid_p", 0.3, winsor, 0,
      standardize = std
    )
  }
  k <- run("two_pass")

  engine <- function(gc, tc, tm) {
    cpp_delta2_score(
      prep$bitset_raw, prep$bitset_words, cols, seq_along(ss), ss,
      gc, tc, numeric(0), numeric(0), tm, 2L, 2L, 0.2
    )[[tm]]
  }
  obs <- cpp_delta2_score(
    prep$bitset_raw, prep$bitset_words, cols, seq_along(ss), ss,
    des$group_code, des$time_code, numeric(0), numeric(0), "group", 2L, 2L,
    0.2
  )
  bins <- .ph_delta2_prevalence_bins(obs$tot / obs$n_samples, 0.3)
  w <- rep(1, length(cols))

  for (tm in des$terms) {
    draws <- lapply(seq_len(B), function(b) {
      if (tm == "time") {
        engine(des$group_code, td[, b], tm)
      } else {
        engine(gd[, b], des$time_code, tm)
      }
    })
    z <- cbind(engine(des$group_code, des$time_code, tm), sapply(draws, c))
    mu <- rowMeans(z)
    sdv <- apply(z, 1, stats::sd)
    zs <- (z - mu) / ifelse(sdv > 1e-12, sdv, Inf)
    zs[] <- pmax(-winsor, pmin(winsor, zs))
    Tz <- apply(zs, 2, .std_pool, w = w, bins = bins)
    expect_equal(k$terms[[tm]]$T_obs, Tz[1], tolerance = 1e-10)
    expect_equal(k$terms[[tm]]$T_null_mean, mean(Tz[-1]), tolerance = 1e-10)
    expect_equal(k$terms[[tm]]$T_null_sd, stats::sd(Tz[-1]), tolerance = 1e-10)
  }

  # the pair is the binary group statistic, standardised the same way
  expect_equal(k$pairs$T_obs, k$terms$group$T_obs, tolerance = 1e-10)

  # without two_pass the kernel is unchanged: winsorised raw values
  k0 <- run("none")
  z <- engine(des$group_code, des$time_code, "group")
  expect_equal(k0$terms$group$T_obs,
    .std_pool(pmax(-winsor, pmin(winsor, z)), w, bins),
    tolerance = 1e-10
  )
  expect_error(run("three_pass"), "standardize")
})

test_that("two_pass recentres a null that propensity draws push off zero", {
  # positivity tracks age, groups track age: under the conditional draws a
  # peptide's z is off-centre and too narrow even after the covariate
  # adjustment, which two_pass removes
  set.seed(8)
  N <- 80
  age <- stats::rnorm(N)
  grp <- ifelse(stats::runif(N) < stats::plogis(1.5 * age), "B", "A")
  samples <- data.frame(
    sample_id = sprintf("x%02d", 1:N), subject_id = sprintf("s%02d", 1:N),
    group = grp, age = age, stringsAsFactors = FALSE
  )
  df <- merge(samples, data.frame(peptide_id = sprintf("p%02d", 1:20)))
  df$exist <- stats::rbinom(nrow(df), 1, stats::plogis(-1 + 1.5 * df$age))

  fit <- function(std) {
    set.seed(3)
    compute_delta2(df, "peptide_id", "group",
      propensity_formula = ~age, strat_bins = 0, B_permutations = 300L,
      standardize = std, return_peptides = FALSE
    )
  }
  none <- fit("none")
  two <- fit("two_pass")
  expect_identical(two$meta$standardize, "two_pass")
  expect_null(two$peptides)
  expect_output(print(two), "two_pass")
  # one peptide per stratum, so T is that peptide's value: its null is
  # off-centre and narrow without two_pass; with it, mean 0 and SD 1 over the
  # observed and the draws, so over the draws alone mean -T/B and SD ~ 1
  expect_gt(mean(abs(none$terms$T_null_mean_group)), 0.1)
  expect_lt(mean(none$terms$T_null_sd_group), 0.9)
  tt <- two$terms
  # (up to winsorising a rare draw beyond winsor_z)
  expect_lt(max(abs(tt$T_null_mean_group + tt$T_group / 300)), 0.01)
  expect_equal(tt$T_null_sd_group, rep(1, 20), tolerance = 0.01)
})
