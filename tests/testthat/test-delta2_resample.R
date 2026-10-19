# tests/testthat/test-delta2_resample.R

# all distinct labellings with the sizes of `labels`
.labellings <- function(labels) {
  if (length(labels) <= 1L) {
    return(list(labels))
  }
  out <- list()
  for (v in unique(labels)) {
    rest <- labels[-match(v, labels)]
    for (p in .labellings(rest)) out[[length(out) + 1L]] <- c(v, p)
  }
  out
}

.key <- function(L) paste(L, collapse = "")

# exact conditional law: pi(L) proportional to prod_i q[i, L_i]
.exact_law <- function(states, q) {
  w <- vapply(states, function(L) prod(q[cbind(seq_along(L), L)]), numeric(1))
  stats::setNames(w / sum(w), vapply(states, .key, ""))
}

# chi-square goodness of fit of drawn labellings (columns) to a law
.gof_p <- function(draws, law) {
  obs <- table(factor(apply(draws, 2, .key), levels = names(law)))
  suppressWarnings(stats::chisq.test(as.numeric(obs), p = law)$p.value)
}

.q6 <- function() {
  set.seed(21)
  q <- matrix(rgamma(18, 2), 6, 3)
  q / rowSums(q)
}

test_that("the swap step leaves the exact conditional law stationary", {
  obs <- c(1L, 1L, 2L, 2L, 3L, 3L)
  states <- .labellings(obs)
  expect_length(states, 90L)
  keys <- vapply(states, .key, "")
  N <- length(obs)

  # one MH step: ordered pair uniform over i != j, swap if labels differ
  transition <- function(q) {
    P <- matrix(0, 90, 90)
    for (a in seq_along(states)) {
      L <- states[[a]]
      for (i in 1:N) {
        for (j in setdiff(1:N, i)) {
          if (L[i] == L[j]) {
            P[a, a] <- P[a, a] + 1 / (N * (N - 1))
            next
          }
          acc <- min(1, q[i, L[j]] * q[j, L[i]] / (q[i, L[i]] * q[j, L[j]]))
          L2 <- L
          L2[c(i, j)] <- L[c(j, i)]
          b <- match(.key(L2), keys)
          P[a, b] <- P[a, b] + acc / (N * (N - 1))
          P[a, a] <- P[a, a] + (1 - acc) / (N * (N - 1))
        }
      }
    }
    P
  }

  q <- .q6()
  law <- .exact_law(states, q)
  stat <- as.numeric(law %*% transition(q))
  expect_equal(stat, unname(law), tolerance = 1e-12)

  # constant propensities: stationary law is uniform over the 90 labellings
  P1 <- transition(matrix(1 / 3, 6, 3))
  u <- rep(1 / 90, 90)
  expect_equal(as.numeric(u %*% P1), u, tolerance = 1e-12)
  expect_equal(rowSums(P1), rep(1, 90), tolerance = 1e-12)
})

test_that("swap sampler frequencies match the exact law over 90 labellings", {
  obs <- c(1L, 1L, 2L, 2L, 3L, 3L)
  states <- .labellings(obs)
  q <- .q6()
  law <- .exact_law(states, q)

  res <- .ph_resample_group(obs, q, B = 60000, seed = 1, n_swap_steps = 60)
  expect_identical(res$method, "swap")
  expect_true(all(apply(res$labels, 2, function(L) all(tabulate(L, 3) == 2))))
  expect_gt(.gof_p(res$labels, law), 1e-3)

  unif <- .ph_resample_group(obs, matrix(1 / 3, 6, 3),
    B = 60000, seed = 2, n_swap_steps = 60
  )
  expect_gt(.gof_p(unif$labels, rep(1 / 90, 90) |>
    stats::setNames(names(law))), 1e-3)
})

test_that("draws are exact even after a single swap step from the hub", {
  # the observed labelling drawn from the law: Besag-Clifford makes every
  # draw follow the same law however short the chain
  obs <- c(1L, 1L, 2L, 2L, 3L, 3L)
  states <- .labellings(obs)
  q <- .q6()
  law <- .exact_law(states, q)
  set.seed(3)
  start <- sample(length(states), 20000, replace = TRUE, prob = law)
  draws <- vapply(seq_along(start), function(r) {
    .ph_resample_group(states[[start[r]]], q,
      B = 1, seed = r, n_swap_steps = 1
    )$labels[, 1]
  }, integer(6))
  expect_gt(.gof_p(draws, law), 1e-3)
})

test_that("K = 2 rejection sampler matches the exact conditional law", {
  obs <- c(1L, 1L, 1L, 2L, 2L, 2L)
  states <- .labellings(obs)
  expect_length(states, 20L)
  p <- c(0.2, 0.35, 0.5, 0.6, 0.7, 0.9)
  law <- .exact_law(states, cbind(1 - p, p))

  res <- .ph_resample_group(obs, p, B = 40000, seed = 4)
  expect_identical(res$method, "rejection")
  expect_gt(.gof_p(res$labels, law), 1e-3)

  # constant propensities: uniform random split with the observed sizes
  unif <- .ph_resample_group(obs, NULL, B = 40000, seed = 5)
  expect_identical(unif$method, "uniform")
  u <- stats::setNames(rep(1 / 20, 20), names(law))
  expect_gt(.gof_p(unif$labels, u), 1e-3)
})

test_that("rejection falls back to the swap chain with a warning", {
  obs <- c(1L, 1L, 1L, 2L, 2L, 2L)
  p <- c(0.01, 0.02, 0.03, 0.97, 0.98, 0.99)
  expect_warning(
    res <- .ph_resample_group(obs, p, B = 50, seed = 6, max_tries = 1),
    "swap chain"
  )
  expect_identical(res$method, "swap")
  expect_true(all(colSums(res$labels == 2L) == 3L))
})

test_that("time coins and within-subject permutations", {
  # s1, s2: two samples at T1/T2; s3: three samples at T1, T2, T3
  sample_subject <- c(1L, 1L, 2L, 2L, 3L, 3L, 3L)
  time_code <- c(1L, 2L, 1L, 2L, 1L, 2L, 3L)
  B <- 30000
  res <- .ph_resample_time(sample_subject, time_code, 3L, B,
    seed = 7, keep_prob = c(0.8, NA, NA), masks = FALSE
  )
  tc <- res$time_code
  # each subject keeps its multiset of timepoints
  for (s in 1:3) {
    r <- which(sample_subject == s)
    tr <- apply(tc[r, , drop = FALSE], 2, sort)
    expect_true(all(tr == sort(time_code[r])))
  }
  kept1 <- mean(tc[1, ] == 1L)
  kept2 <- mean(tc[3, ] == 1L)
  expect_lt(abs(kept1 - 0.8), 4 * sqrt(0.8 * 0.2 / B))
  expect_lt(abs(kept2 - 0.5), 4 * sqrt(0.25 / B))
  perm3 <- table(apply(tc[5:7, ], 2, .key))
  expect_length(perm3, 6L)
  expect_gt(stats::chisq.test(as.numeric(perm3))$p.value, 1e-3)
})

test_that("F masks give the relabelled per-timepoint counts", {
  set.seed(8)
  vis <- c(2, 4, 3, 2, 3, 1, 4)
  samples <- do.call(rbind, lapply(seq_along(vis), function(i) {
    data.frame(
      sample_id = sprintf("s%d_T%d", i, seq_len(vis[i])),
      subject_id = sprintf("s%d", i),
      timepoint = sprintf("T%d", seq_len(vis[i])),
      group = "A", stringsAsFactors = FALSE
    )
  }))
  df <- samples[rep(seq_len(nrow(samples)), 5), ]
  df$peptide_id <- rep(sprintf("p%d", 1:5), each = nrow(samples))
  df$exist <- rbinom(nrow(df), 1, 0.4)
  prep <- .ph_delta2_prepare(df, group_cols = "group", time_col = "timepoint")
  ns <- length(prep$samples_order)
  Y <- matrix(0L, ns, prep$bitset_m)
  pos <- df[df$exist > 0, ]
  Y[cbind(
    match(pos$sample_id, prep$samples_order),
    match(pos$peptide_id, prep$peptides_order)
  )] <- 1L

  B <- 5
  res <- .ph_resample_time(prep$sample_subject, prep$sample_time, 4L, B,
    seed = 9, n_words = prep$bitset_words
  )
  for (b in seq_len(B)) {
    cnt <- cpp_delta2_time_counts(
      prep$bitset_raw, prep$bitset_words, seq_len(prep$bitset_m),
      res$masks, 4L, b
    )
    ref <- t(apply(Y, 2, function(y) tabulate(res$time_code[y == 1L, b], 4)))
    expect_equal(cnt, ref)
  }
})

test_that("mask counts reproduce the engine's resampled time z", {
  # binary time: U = sum_s x_s y_s - sum_i xbar_i c_i, and the second sum
  # and the variance are fixed under within-subject draws
  set.seed(10)
  N <- 25
  samples <- data.frame(
    sample_id = c(sprintf("s%02d_T1", 1:N), sprintf("s%02d_T2", 1:N)),
    subject_id = rep(sprintf("s%02d", 1:N), 2),
    timepoint = rep(c("T1", "T2"), each = N),
    group = "A", stringsAsFactors = FALSE
  )
  df <- samples[rep(1:(2 * N), 6), ]
  df$peptide_id <- rep(sprintf("p%d", 1:6), each = 2 * N)
  df$exist <- rbinom(nrow(df), 1, 0.4)
  prep <- .ph_delta2_prepare(df, group_cols = "group", time_col = "timepoint")
  ss <- prep$sample_subject
  score <- function(tc) {
    cpp_delta2_score(
      prep$bitset_raw, prep$bitset_words, seq_len(prep$bitset_m),
      seq_along(ss), ss, rep(1L, N), tc, numeric(0), numeric(0),
      "time", 1L, 2L, 0
    )$time[, 1]
  }
  z_obs <- score(prep$sample_time)
  obs_cnt <- cpp_delta2_score(
    prep$bitset_raw, prep$bitset_words, seq_len(prep$bitset_m),
    seq_along(ss), ss, rep(1L, N), prep$sample_time, numeric(0), numeric(0),
    "time", 1L, 2L, 0
  )$S_time
  # xbar_i = 1/2 for every subject, so U = count_T2 - tot / 2
  U_obs <- obs_cnt[, 2] - rowSums(obs_cnt) / 2

  res <- .ph_resample_time(ss, prep$sample_time, 2L, 3, seed = 11)
  for (b in 1:3) {
    cnt <- cpp_delta2_time_counts(
      prep$bitset_raw, prep$bitset_words, seq_len(prep$bitset_m),
      res$masks, 2L, b
    )
    U_b <- cnt[, 2] - rowSums(cnt) / 2
    nz <- U_obs != 0
    expect_equal(score(res$time_code[, b])[nz], (z_obs / U_obs * U_b)[nz])
  }
})

.ps_data <- function(n = 120, K = 2) {
  set.seed(12)
  age <- round(rnorm(n, 50, 10))
  sex <- sample(c("f", "m"), n, TRUE)
  eta <- 0.05 * (age - 50) + 0.4 * (sex == "m")
  if (K == 2) {
    code <- 1L + rbinom(n, 1, plogis(eta))
  } else {
    code <- vapply(eta, function(e) {
      sample(1:3, 1, prob = exp(c(0, e, -e)))
    }, integer(1))
  }
  meta <- data.frame(subject_id = seq_len(n), age = age, sex = sex)
  list(meta = meta, code = code)
}

test_that("two-group propensity is the glm fit, with diagnostics", {
  d <- .ps_data()
  ps <- .ph_delta2_propensity(d$meta, d$code, ~ splines::ns(age, 3) + sex)
  y <- d$code == 2L
  ref <- stats::glm(y ~ splines::ns(age, 3) + sex,
    family = stats::binomial(), data = d$meta
  )
  expect_equal(ps$prob[, 2], unname(stats::fitted(ref)), tolerance = 1e-10)
  expect_equal(rowSums(ps$prob), rep(1, nrow(d$meta)))
  expect_identical(ps$diagnostics$n_params, 4L)
  ll <- as.numeric(stats::logLik(ref))
  fit0 <- stats::glm(y ~ 1, family = stats::binomial())
  ll0 <- as.numeric(stats::logLik(fit0))
  expect_equal(ps$diagnostics$pseudo_r2, 1 - ll / ll0, tolerance = 1e-10)
  p <- ps$prob[, 2]
  expect_identical(ps$diagnostics$n_overlap, sum(p > 0.05 & p < 0.95))

  # constant propensities without a formula
  ps0 <- .ph_delta2_propensity(d$meta, d$code)
  expect_null(ps0$prob)
  expect_identical(ps0$diagnostics$pseudo_r2, 0)

  # strata: each subject gets its stratum's observed share
  pss <- .ph_delta2_propensity(d$meta, d$code, ~ age + sex, n_ps_strata = 4)
  expect_lte(length(unique(pss$prob[, 2])), 4L)
  expect_equal(mean(pss$prob[, 2]), mean(y), tolerance = 1e-12)

  # trimming drops subjects outside the range
  pst <- .ph_delta2_propensity(d$meta, d$code, ~ age + sex,
    trim_propensity = c(0.3, 0.7)
  )
  out <- pst$diagnostics$n_trimmed
  expect_gt(out, 0L)
  expect_identical(sum(!pst$keep), out)
  expect_true(all(is.na(pst$prob[!pst$keep, ])))
})

test_that("three-group propensity uses multinom", {
  skip_if_not_installed("nnet")
  d <- .ps_data(K = 3)
  ps <- .ph_delta2_propensity(d$meta, d$code, ~ age + sex)
  expect_identical(dim(ps$prob), c(120L, 3L))
  expect_equal(rowSums(ps$prob), rep(1, 120), tolerance = 1e-8)
  expect_identical(ps$diagnostics$n_params, 4L)

  # the sampler accepts it and keeps the group sizes
  res <- .ph_delta2_propensity(d$meta, d$code, ~ age)$prob |>
    (\(q) .ph_resample_group(d$code, q, B = 20, seed = 13))()
  expect_identical(res$method, "swap")
  expect_true(all(apply(res$labels, 2, tabulate, 3) == tabulate(d$code, 3)))
})

test_that("propensity warns on poor overlap and over-fitting", {
  d <- .ps_data(n = 40)
  sep <- d$meta
  sep$x <- ifelse(d$code == 2L, 1, -1) + rnorm(40, 0, 0.05)
  expect_warning(
    expect_warning(
      .ph_delta2_propensity(sep, d$code, ~x),
      "overlap"
    ),
    "over-fitted"
  )
})

test_that("time coins come from the within-subject covariate difference", {
  set.seed(14)
  N <- 200L
  dose <- rnorm(2 * N)
  sample_subject <- rep(1:N, each = 2)
  # the later visit tends to have the higher dose
  time_code <- unlist(lapply(1:N, function(i) {
    d <- dose[2 * i] - dose[2 * i - 1]
    if (runif(1) < plogis(0.8 * d)) c(1L, 2L) else c(2L, 1L)
  }))
  meta <- data.frame(dose = dose)
  tp <- .ph_delta2_time_propensity(meta, sample_subject, time_code, ~dose)

  late <- vapply(1:N, function(i) {
    r <- 2 * i - c(1, 0)
    dose[r[time_code[r] == 2L]] - dose[r[time_code[r] == 1L]]
  }, numeric(1))
  fit <- stats::glm(rep(1, N) ~ 0 + late, family = stats::binomial())
  expect_equal(tp$keep_prob, unname(stats::plogis(late * stats::coef(fit))),
    tolerance = 1e-6
  )
  expect_gt(stats::coef(fit), 0.3)
  expect_identical(tp$diagnostics$n, N)

  fair <- .ph_delta2_time_propensity(meta, sample_subject, time_code)
  expect_true(all(is.na(fair$keep_prob)))
})

test_that("constant propensities give the uniform-shuffle null (K = 2)", {
  # 8 subjects, 4 per group: 70 labellings, uniform under the null
  obs <- rep(1:2, each = 4)
  states <- .labellings(obs)
  expect_length(states, 70L)
  uniform <- stats::setNames(rep(1 / 70, 70), vapply(states, .key, ""))
  for (p in c(0.5, 0.2)) {
    res <- .ph_resample_group(obs, rep(p, 8), B = 35000, seed = 31)
    expect_identical(res$method, "rejection")
    expect_gt(.gof_p(res$labels, uniform), 1e-3)
  }

  # the pooled statistic's null: exact over all 70 labellings, then from
  # rejection draws with a constant propensity and from uniform shuffles
  set.seed(32)
  samples <- data.frame(
    sample_id = sprintf("x%d", 1:8), subject_id = sprintf("s%d", 1:8),
    stringsAsFactors = FALSE
  )
  samples$group <- c("a", "b")[obs]
  df <- merge(samples, data.frame(peptide_id = sprintf("p%d", 1:5)))
  df$exist <- stats::rbinom(nrow(df), 1, 0.4)
  prep <- .ph_delta2_prepare(df, "group")
  ss <- prep$sample_subject
  code <- obs[match(prep$subjects_order, samples$subject_id)]
  kernel <- function(draws) {
    cpp_delta2_contrast(
      prep$bitset_raw, prep$bitset_words, seq_len(prep$bitset_m),
      seq_along(ss), ss, code, rep(1L, length(ss)), numeric(0), numeric(0),
      "group", 2L, 1L, 0, draws, matrix(0L, 0, 0), matrix(0L, 0, 2),
      numeric(0), "n_eff_sqrt", "stouffer", "mid_p", 0, 4, 0
    )$terms$group
  }
  exact <- kernel(do.call(cbind, states))
  T_all <- vapply(states, function(L) {
    cpp_delta2_contrast(
      prep$bitset_raw, prep$bitset_words, seq_len(prep$bitset_m),
      seq_along(ss), ss, L, rep(1L, length(ss)), numeric(0), numeric(0),
      "group", 2L, 1L, 0, matrix(L), matrix(0L, 0, 0), matrix(0L, 0, 2),
      numeric(0), "n_eff_sqrt", "stouffer", "mid_p", 0, 4, 0
    )$terms$group$T_obs
  }, numeric(1))
  expect_equal(exact$T_null_mean, mean(T_all), tolerance = 1e-12)

  B <- 20000
  rej <- .ph_resample_group(code, rep(0.3, 8), B = B, seed = 33)$labels
  uni <- .ph_resample_group(code, NULL, B = B, seed = 34)$labels
  T_of <- function(draws) T_all[match(apply(draws, 2, .key), names(uniform))]
  t_rej <- T_of(rej)
  t_uni <- T_of(uni)
  sd_T <- sqrt(mean((T_all - mean(T_all))^2))
  for (t in list(t_rej, t_uni)) {
    expect_lt(abs(mean(t) - mean(T_all)), 4 * sd_T / sqrt(B))
  }
  ks <- suppressWarnings(stats::ks.test(t_rej, t_uni, exact = FALSE))
  expect_gt(ks$p.value, 1e-3)
  expect_equal(kernel(rej)$T_null_mean, mean(t_rej), tolerance = 1e-10)
})
