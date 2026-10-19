# tests/testthat/test-delta2_score.R

# long data from a sample table and a sample x peptide 0/1 matrix
.d2_long <- function(samples, Y) {
  peps <- sprintf("p%02d", seq_len(ncol(Y)))
  df <- samples[rep(seq_len(nrow(samples)), times = ncol(Y)), , drop = FALSE]
  df$peptide_id <- rep(peps, each = nrow(samples))
  df$exist <- as.integer(Y)
  rownames(df) <- NULL
  df
}

# prepare + design + score for one contrast
.d2_run <- function(df, group_col = "group", time_col = NULL, rho = 0,
                    group_scores = NULL, time_scores = NULL,
                    covariates = NULL, ...) {
  prep <- .ph_delta2_prepare(df, group_cols = group_col, time_col = time_col)
  gl <- prep$group_levels[[group_col]]
  subjects_meta <- data.frame(
    subject_id = prep$subjects_order,
    stringsAsFactors = FALSE
  )
  subjects_meta[[group_col]] <- gl[prep$subject_group[[group_col]]]
  samples_meta <- data.frame(
    sample_id = prep$samples_order,
    subject_id = prep$subjects_order[prep$sample_subject],
    stringsAsFactors = FALSE
  )
  if (!is.null(time_col)) {
    samples_meta[[time_col]] <- prep$time_levels[prep$sample_time]
  }
  des <- .ph_delta2_design(subjects_meta, samples_meta,
    group_col = group_col, time_col = time_col,
    group_scores = group_scores, time_scores = time_scores, ...
  )
  res <- cpp_delta2_score(
    prep$bitset_raw, prep$bitset_words, seq_len(prep$bitset_m),
    seq_along(prep$samples_order), prep$sample_subject,
    des$group_code, des$time_code, des$group_scores, des$time_scores,
    des$terms, des$n_group_levels, des$n_time_levels, rho, covariates
  )
  list(prep = prep, des = des, res = res)
}

# per-peptide sample bits in bitset row order, from the long data
.d2_bits <- function(run, df) {
  p <- run$prep
  Y <- matrix(0L, length(p$samples_order), length(p$peptides_order))
  pos <- df[df$exist > 0, ]
  Y[cbind(
    match(pos$sample_id, p$samples_order),
    match(pos$peptide_id, p$peptides_order)
  )] <- 1L
  Y
}

# all permutations of a vector (small n only)
.perms <- function(x) {
  if (length(x) <= 1L) {
    return(list(x))
  }
  out <- list()
  for (i in seq_along(x)) {
    for (p in .perms(x[-i])) out[[length(out) + 1L]] <- c(x[i], p)
  }
  out
}

# U' V^+ U for a covariance V
# Welch-type interaction, by hand: subject residuals R (rows) and group
# codes; S = sum_k c_k (x) RS_k and V = sum_k n_k c_k c_k' (x) Sigma_k, with
# c_k the centred group design row (G: levels x columns)
.welch_SV <- function(R, code, G) {
  R <- as.matrix(R)
  gbar <- colMeans(G[code, , drop = FALSE])
  S <- 0
  V <- 0
  for (k in sort(unique(code))) {
    Rk <- R[code == k, , drop = FALSE]
    ck <- G[k, ] - gbar
    S <- S + kronecker(ck, colSums(Rk))
    if (nrow(Rk) >= 2) {
      V <- V + nrow(Rk) * kronecker(tcrossprod(ck), stats::cov(Rk))
    }
  }
  list(S = S, V = V)
}

MASS_ginv <- function(M) {
  e <- eigen(M, symmetric = TRUE)
  keep <- e$values > 1e-10 * max(e$values)
  v <- e$vectors[, keep, drop = FALSE]
  v %*% (t(v) / e$values[keep])
}

.quad <- function(U, V) {
  e <- eigen(V, symmetric = TRUE)
  keep <- e$values > 1e-10 * max(e$values)
  sum((crossprod(e$vectors[, keep, drop = FALSE], U))^2 / e$values[keep])
}

test_that("two groups, one sample per subject: v1 score z up to (N-1)/N", {
  set.seed(11)
  N <- 40
  samples <- data.frame(
    sample_id = sprintf("x%02d", 1:N), subject_id = sprintf("s%02d", 1:N),
    group = rep(c("A", "B"), c(17, 23)), stringsAsFactors = FALSE
  )
  Y <- matrix(rbinom(N * 30, 1, 0.3), N)
  df <- .d2_long(samples, Y)
  run <- .d2_run(df)
  Yb <- .d2_bits(run, df)
  grp <- samples$group[match(run$prep$samples_order, samples$sample_id)]

  # v1 "score": (p2 - p1) / sqrt(p0 (1 - p0) (1/n1 + 1/n2)), 0 at boundary
  x1 <- colSums(Yb[grp == "A", ])
  x2 <- colSums(Yb[grp == "B", ])
  n1 <- sum(grp == "A")
  n2 <- sum(grp == "B")
  p0 <- (x1 + x2) / N
  z_v1 <- ifelse(p0 > 0 & p0 < 1,
    (x2 / n2 - x1 / n1) / sqrt(p0 * (1 - p0) * (1 / n1 + 1 / n2)), 0
  )

  # exact permutation variance has N - 1 where the binomial one has N
  expect_equal(
    run$res$group[, 1], z_v1 * sqrt((N - 1) / N),
    tolerance = 1e-10
  )
  expect_identical(run$des$terms, "group")
})

test_that("rho does not matter with one sample per subject", {
  set.seed(12)
  samples <- data.frame(
    sample_id = sprintf("x%02d", 1:20), subject_id = sprintf("s%02d", 1:20),
    group = rep(c("A", "B"), 10), stringsAsFactors = FALSE
  )
  df <- .d2_long(samples, matrix(rbinom(20 * 5, 1, 0.4), 20))
  expect_equal(
    .d2_run(df, rho = 0)$res$group,
    .d2_run(df, rho = 1)$res$group,
    tolerance = 1e-12
  )
})

test_that("two complete timepoints give the McNemar z", {
  set.seed(13)
  N <- 50
  samples <- data.frame(
    sample_id = c(sprintf("s%02d_T1", 1:N), sprintf("s%02d_T2", 1:N)),
    subject_id = rep(sprintf("s%02d", 1:N), 2),
    timepoint = rep(c("T1", "T2"), each = N),
    group = "A", stringsAsFactors = FALSE
  )
  Y <- matrix(rbinom(2 * N * 25, 1, 0.3), 2 * N)
  df <- .d2_long(samples, Y)
  run <- .d2_run(df, time_col = "timepoint")
  Yb <- .d2_bits(run, df)
  t2 <- run$prep$sample_time == 2L
  y1 <- Yb[!t2, ] # rows in subject order for both timepoints
  y2 <- Yb[t2, ]
  b <- colSums(y1 == 0 & y2 == 1)
  c <- colSums(y1 == 1 & y2 == 0)
  z_mcn <- ifelse(b + c > 0, (b - c) / sqrt(b + c), 0)

  expect_identical(run$des$terms, "time")
  expect_equal(run$res$time[, 1], z_mcn, tolerance = 1e-10)
})

test_that("T complete timepoints give Cochran's Q as a signed vector", {
  set.seed(14)
  N <- 30
  T_n <- 4
  samples <- data.frame(
    sample_id = sprintf("s%02d_T%d", rep(1:N, T_n), rep(1:T_n, each = N)),
    subject_id = rep(sprintf("s%02d", 1:N), T_n),
    timepoint = rep(sprintf("T%d", 1:T_n), each = N),
    group = "A", stringsAsFactors = FALSE
  )
  Y <- matrix(rbinom(N * T_n * 20, 1, 0.35), N * T_n)
  df <- .d2_long(samples, Y)
  run <- .d2_run(df, time_col = "timepoint")
  Yb <- .d2_bits(run, df)
  subj <- run$prep$sample_subject
  tim <- run$prep$sample_time

  Q <- apply(Yb, 2, function(y) {
    C_t <- tapply(y, tim, sum)
    R_i <- tapply(y, subj, sum)
    Ntot <- sum(y)
    den <- T_n * Ntot - sum(R_i^2)
    if (den == 0) 0 else (T_n - 1) * (T_n * sum(C_t^2) - Ntot^2) / den
  })
  r <- run$res$time

  expect_identical(ncol(r), 4L)
  expect_identical(run$des$df[["time"]], 3L)
  expect_equal(rowSums(r^2), unname(Q), tolerance = 1e-10)
  expect_equal(rowSums(r), rep(0, ncol(Yb)), tolerance = 1e-10)
  # balanced: r_t = (C_t - N_tot / T) / sigma has the sign of the excess
  C_t <- t(apply(Yb, 2, function(y) tapply(y, tim, sum)))
  excess <- C_t - rowMeans(C_t)
  nz <- abs(excess) > 1e-9
  expect_true(all(sign(r[nz]) == sign(excess[nz])))
})

test_that("2x2 interaction is the Welch t of gained - lost between groups", {
  set.seed(15)
  N <- 40
  grp <- rep(c("A", "B"), c(18, 22))
  samples <- data.frame(
    sample_id = c(sprintf("s%02d_T1", 1:N), sprintf("s%02d_T2", 1:N)),
    subject_id = rep(sprintf("s%02d", 1:N), 2),
    timepoint = rep(c("T1", "T2"), each = N),
    group = rep(grp, 2), stringsAsFactors = FALSE
  )
  Y <- matrix(rbinom(2 * N * 25, 1, 0.3), 2 * N)
  df <- .d2_long(samples, Y)
  run <- .d2_run(df, time_col = "timepoint")
  Yb <- .d2_bits(run, df)
  t2 <- run$prep$sample_time == 2L
  d <- Yb[t2, ] - Yb[!t2, ] # +1 gained, -1 lost, per subject
  inB <- run$des$group_code == 2L
  nA <- sum(!inB)
  nB <- sum(inB)

  # each group's own variance of the change (Welch), not the pooled one
  z_ref <- apply(d, 2, function(di) {
    se2 <- stats::var(di[!inB]) / nA + stats::var(di[inB]) / nB
    if (se2 == 0) {
      return(0)
    }
    (mean(di[inB]) - mean(di[!inB])) / sqrt(se2)
  })

  expect_identical(run$des$terms, c("group", "time", "interaction"))
  expect_equal(run$res$interaction[, 1], z_ref, tolerance = 1e-10)
})

test_that("ordered groups with scores give the linear-trend score z", {
  set.seed(16)
  N <- 45
  samples <- data.frame(
    sample_id = sprintf("x%02d", 1:N), subject_id = sprintf("s%02d", 1:N),
    group = rep(c("lo", "mid", "hi"), each = 15), stringsAsFactors = FALSE
  )
  df <- .d2_long(samples, matrix(rbinom(N * 20, 1, 0.3), N))
  sc <- c(lo = 0, mid = 1, hi = 3)
  run <- .d2_run(df, group_type = "multi", group_scores = sc)
  Yb <- .d2_bits(run, df)
  a <- sc[run$des$group_levels][run$des$group_code]

  z_ref <- apply(Yb, 2, function(y) {
    ss <- sum((y - mean(y))^2)
    if (ss == 0) {
      return(0)
    }
    sum((a - mean(a)) * y) / sqrt(sum((a - mean(a))^2) * ss / (N - 1))
  })
  expect_identical(run$des$df[["group"]], 1L)
  expect_equal(run$res$group[, 1], unname(z_ref), tolerance = 1e-10)
})

# unbalanced panel: 2/4/3/2 visits
.d2_unbalanced <- function() {
  vis <- list(s1 = c(1, 3), s2 = 1:4, s3 = c(1, 2, 4), s4 = c(2, 4))
  samples <- do.call(rbind, lapply(names(vis), function(s) {
    data.frame(
      sample_id = paste0(s, "_T", vis[[s]]), subject_id = s,
      timepoint = paste0("T", vis[[s]]), group = "A",
      stringsAsFactors = FALSE
    )
  }))
  Y <- cbind(
    c(0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1),
    c(1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1),
    c(0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1)
  )
  .d2_long(samples, Y)
}

test_that("within-subject variance matches all 576 permutations", {
  df <- .d2_unbalanced()
  sc <- c(T1 = 0, T2 = 1, T3 = 2, T4 = 3)

  for (scores in list(sc, NULL)) {
    run <- .d2_run(df, time_col = "timepoint", time_scores = scores)
    Yb <- .d2_bits(run, df)
    subj <- run$prep$sample_subject
    tim <- run$prep$sample_time
    X <- if (is.null(scores)) diag(4)[tim, ] else matrix(sc[tim])
    rows <- split(seq_along(subj), subj)
    # U = sum_i X_i' (y_i - ybar_i)
    U_of <- function(y) {
      colSums(do.call(rbind, lapply(rows, function(r) {
        crossprod(X[r, , drop = FALSE], y[r] - mean(y[r]))[, 1]
      })))
    }
    perm_rows <- lapply(rows, .perms)
    grid <- expand.grid(lapply(perm_rows, seq_along))
    expect_identical(nrow(grid), 576L)

    for (j in seq_len(ncol(Yb))) {
      y <- Yb[, j]
      Us <- t(apply(grid, 1, function(g) {
        yp <- y
        for (k in seq_along(rows)) yp[rows[[k]]] <- y[perm_rows[[k]][[g[k]]]]
        U_of(yp)
      }))
      if (ncol(X) == 1L) Us <- t(Us)
      expect_equal(colMeans(Us), rep(0, ncol(X)), tolerance = 1e-12)
      V <- crossprod(Us) / nrow(Us)
      U <- U_of(y)
      if (ncol(X) == 1L) {
        expect_equal(run$res$time[j, 1], U / sqrt(V[1, 1]), tolerance = 1e-10)
      } else {
        expect_equal(sum(run$res$time[j, ]^2), .quad(U, V), tolerance = 1e-10)
      }
    }
  }
})

test_that("a single-sample subject contributes 0 to the time term", {
  df <- .d2_unbalanced()
  extra <- df[df$sample_id == "s2_T1", ]
  extra$sample_id <- "s5_T2"
  extra$subject_id <- "s5"
  extra$timepoint <- "T2"
  extra$exist <- 1L
  a <- .d2_run(df, time_col = "timepoint")
  b <- .d2_run(rbind(df, extra), time_col = "timepoint")
  expect_equal(a$res$time, b$res$time, tolerance = 1e-12)
})

# 6 subjects with unequal sample counts, one timepoint
.d2_between <- function(groups) {
  n_i <- c(1, 3, 2, 4, 2, 1)
  samples <- data.frame(
    sample_id = sprintf("s%d_%d", rep(1:6, n_i), sequence(n_i)),
    subject_id = sprintf("s%d", rep(1:6, n_i)),
    group = rep(groups, n_i), stringsAsFactors = FALSE
  )
  Y <- cbind(
    c(1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0),
    c(0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1)
  )
  .d2_long(samples, Y)
}

test_that("between-subject z matches all label permutations for any rho", {
  designs <- list(
    binary = c("A", "A", "A", "B", "B", "B"),
    multi = c("a", "b", "c", "a", "b", "c")
  )
  for (type in names(designs)) {
    df <- .d2_between(designs[[type]])
    for (rho in c(0, 0.4, 1)) {
      run <- .d2_run(df, rho = rho, group_type = type)
      Yb <- .d2_bits(run, df)
      subj <- run$prep$sample_subject
      n_i <- tabulate(subj)
      w <- n_i / (1 + (n_i - 1) * rho)
      code <- run$des$group_code
      K <- run$des$n_group_levels
      G <- if (K == 2L) matrix(code - 1) else diag(K)[code, ]
      labs <- unique(.perms(code))

      for (j in seq_len(ncol(Yb))) {
        ybar <- as.numeric(tapply(Yb[, j], subj, mean))
        v <- w * (ybar - sum(w * ybar) / sum(w))
        Us <- t(sapply(labs, function(cd) {
          Gp <- if (K == 2L) matrix(cd - 1) else diag(K)[cd, ]
          colSums(Gp * v)
        }))
        if (K == 2L) Us <- t(Us)
        V <- stats::cov(Us) * (nrow(Us) - 1) / nrow(Us)
        U <- colSums(G * v) - colMeans(Us)
        r <- run$res$group[j, ]
        if (K == 2L) {
          expect_equal(r, U / sqrt(V[1, 1]), tolerance = 1e-10)
        } else {
          expect_equal(sum(r^2), .quad(U, V), tolerance = 1e-10)
        }
      }
    }
  }
})

test_that("unequal sample counts alone carry no group signal", {
  # every subject positive in half its samples, 2 to 6 samples each: subject
  # means are all 0.5, so the group term must be 0 whatever the weights
  n_i <- c(6, 6, 4, 2, 2, 4)
  samples <- data.frame(
    sample_id = sprintf("s%d_%d", rep(1:6, n_i), sequence(n_i)),
    subject_id = sprintf("s%d", rep(1:6, n_i)),
    group = rep(c("A", "A", "A", "B", "B", "B"), n_i),
    stringsAsFactors = FALSE
  )
  y <- unlist(lapply(n_i, function(n) rep(0:1, n / 2)))
  df <- .d2_long(samples, cbind(y, 1 - y))
  for (rho in c(0, 0.3, 1)) {
    expect_equal(.d2_run(df, rho = rho)$res$group[, 1], c(0, 0))
  }
})

test_that("covariates adjust the group term: weighted least-squares score", {
  set.seed(20)
  n_sub <- 36
  n_i <- sample(1:4, n_sub, TRUE)
  age <- stats::rnorm(n_sub)
  for (K in c(2, 3)) {
    grp <- letters[1:K][(seq_len(n_sub) - 1) %% K + 1]
    samples <- data.frame(
      sample_id = sprintf("s%02d_%d", rep(1:n_sub, n_i), sequence(n_i)),
      subject_id = sprintf("s%02d", rep(1:n_sub, n_i)),
      group = rep(grp, n_i), stringsAsFactors = FALSE
    )
    prob <- stats::plogis(-1 + rep(age, n_i))
    Y <- matrix(rbinom(nrow(samples) * 8, 1, prob), nrow(samples))
    df <- .d2_long(samples, Y)
    type <- if (K == 2) "binary" else "multi"
    for (rho in c(0, 0.5)) {
      run0 <- .d2_run(df, rho = rho, group_type = type)
      sub <- match(run0$prep$subjects_order, sprintf("s%02d", 1:n_sub))
      Zc <- cbind(age = age[sub])
      run <- .d2_run(df, rho = rho, covariates = Zc, group_type = type)
      Yb <- .d2_bits(run, df)
      subj <- run$prep$sample_subject
      ni <- tabulate(subj)
      w <- ni / (1 + (ni - 1) * rho)
      Z <- cbind(1, Zc)
      code <- run$des$group_code
      G <- if (K == 2) matrix(c(0, 1)) else diag(K)
      Gc <- sweep(G[code, , drop = FALSE], 2, colMeans(G[code, , drop = FALSE]))
      N <- length(code)
      for (j in seq_len(ncol(Yb))) {
        ybar <- as.numeric(tapply(Yb[, j], subj, mean))
        a <- solve(crossprod(Z, w * Z), crossprod(Z, w * ybar))
        u <- w * (ybar - Z %*% a)[, 1]
        U <- colSums(Gc * u)
        V <- sum(u^2) * crossprod(Gc) / (N - 1)
        if (K == 2) {
          z <- U / sqrt(V[1, 1])
          expect_equal(run$res$group[j, 1], z, tolerance = 1e-10)
        } else {
          q <- .quad(U, V)
          expect_equal(sum(run$res$group[j, ]^2), q, tolerance = 1e-10)
        }
      }
      # an intercept-only covariate matrix changes nothing
      ones <- .d2_run(df,
        rho = rho, covariates = matrix(1, N, 1), group_type = type
      )
      expect_equal(ones$res$group, run0$res$group, tolerance = 1e-10)
    }
  }
})

test_that("interaction is residualised against the time main effect", {
  # cases have 4 visits, controls 2; positivity depends on time only
  vis <- list(c1 = 1:4, c2 = 1:4, c3 = 1:4, k1 = 1:2, k2 = 1:2, k3 = 1:2)
  samples <- do.call(rbind, lapply(names(vis), function(s) {
    data.frame(
      sample_id = paste0(s, "_T", vis[[s]]), subject_id = s,
      timepoint = paste0("T", vis[[s]]),
      group = if (startsWith(s, "c")) "case" else "ctrl",
      stringsAsFactors = FALSE
    )
  }))
  f <- c(0, 0, 1, 1) # time main effect only
  Y <- cbind(f[as.integer(substring(samples$timepoint, 2))], c(1, rep(0, 17)))
  df <- .d2_long(samples, Y)

  # omnibus time absorbs any time main effect: interaction is exactly 0
  omni <- .d2_run(df, time_col = "timepoint")
  expect_equal(omni$res$interaction[1, ], rep(0, 4), tolerance = 1e-12)
  expect_true(sum(omni$res$time[1, ]^2) > 5)

  # a linear trend only absorbs the linear part of the time effect
  sc <- c(T1 = 1, T2 = 2, T3 = 3, T4 = 4)
  run <- .d2_run(df, time_col = "timepoint", time_scores = sc)

  # second peptide by hand: residualised D_i, Welch standardisation
  Yb <- .d2_bits(run, df)
  subj <- run$prep$sample_subject
  x <- sc[run$prep$sample_time]
  rows <- split(seq_along(subj), subj)
  y <- Yb[, 2]
  D <- sapply(rows, function(r) sum((x[r] - mean(x[r])) * (y[r] - mean(y[r]))))
  Sx <- sapply(rows, function(r) sum((x[r] - mean(x[r]))^2))
  R <- D - Sx * sum(D) / sum(Sx)
  sv <- .welch_SV(R, run$des$group_code, matrix(c(0, 1)))
  expect_equal(run$res$interaction[2, 1], as.numeric(sv$S) / sqrt(sv$V[1, 1]),
    tolerance = 1e-10
  )
})

test_that("omnibus interaction is the Welch quadratic form", {
  # 3 groups x 3 timepoints, unbalanced visits and unequal turnover
  set.seed(19)
  vis <- lapply(1:24, function(i) sort(sample(1:3, sample(2:3, 1))))
  samples <- do.call(rbind, lapply(seq_along(vis), function(i) {
    data.frame(
      sample_id = sprintf("s%02d_T%d", i, vis[[i]]),
      subject_id = sprintf("s%02d", i), timepoint = sprintf("T%d", vis[[i]]),
      group = c("a", "b", "c")[(i - 1) %% 3 + 1], stringsAsFactors = FALSE
    )
  }))
  prev <- ifelse(samples$group == "c", 0.5, 0.2)
  Y <- matrix(rbinom(nrow(samples) * 6, 1, prev), nrow(samples))
  df <- .d2_long(samples, Y)
  run <- .d2_run(df, time_col = "timepoint")
  Yb <- .d2_bits(run, df)
  subj <- run$prep$sample_subject
  X <- diag(3)[run$prep$sample_time, ]
  rows <- split(seq_along(subj), subj)
  A <- lapply(rows, function(r) {
    crossprod(scale(X[r, , drop = FALSE], scale = FALSE))
  })
  for (j in seq_len(ncol(Yb))) {
    y <- Yb[, j]
    D <- t(sapply(rows, function(r) {
      crossprod(X[r, , drop = FALSE], y[r] - mean(y[r]))[, 1]
    }))
    spread <- vapply(A, function(a) sum(diag(a)) > 1e-12, TRUE)
    sumA <- Reduce(`+`, A[spread])
    beta <- MASS_ginv(sumA) %*% colSums(D[spread, , drop = FALSE])
    R <- t(sapply(seq_along(rows), function(i) D[i, ] - A[[i]] %*% beta))
    R <- R[spread, ]
    sv <- .welch_SV(R, run$des$group_code[spread], diag(3))
    expect_equal(sum(run$res$interaction[j, ]^2), .quad(sv$S, sv$V),
      tolerance = 1e-8
    )
  }
})

test_that("pooled rho is the one-way ANOVA moment estimator", {
  df <- .d2_between(c("A", "A", "A", "B", "B", "B"))
  run <- .d2_run(df)
  Yb <- .d2_bits(run, df)
  subj <- run$prep$sample_subject
  n_i <- tabulate(subj)
  k <- length(n_i)
  Ns <- length(subj)

  ss <- sapply(seq_len(ncol(Yb)), function(j) {
    m_i <- tapply(Yb[, j], subj, mean)
    c(
      ssb = sum(n_i * (m_i - mean(Yb[, j]))^2),
      ssw = sum((Yb[, j] - m_i[subj])^2)
    )
  })
  MSB <- sum(ss["ssb", ]) / (k - 1)
  MSW <- sum(ss["ssw", ]) / (Ns - k)
  n0 <- (Ns - sum(n_i^2) / Ns) / (k - 1)
  ref <- min(1, max(0, (MSB - MSW) / (MSB + (n0 - 1) * MSW)))

  est <- .ph_delta2_rho("pooled",
    bitset_raw = run$prep$bitset_raw, n_words = run$prep$bitset_words,
    pep_cols = seq_len(run$prep$bitset_m),
    sample_rows = seq_along(subj), sample_subject = subj
  )
  expect_equal(est$rho_hat, ref, tolerance = 1e-12)
  expect_equal(est$rho, ref, tolerance = 1e-12)
  expect_equal(est$subject_hits, as.numeric(rowsum(rowSums(Yb), subj)))
  expect_identical(.ph_delta2_rho("one")$rho, 1)
  expect_identical(.ph_delta2_rho(0.25)$rho, 0.25)
  expect_error(.ph_delta2_rho(1.5))

  # no repeat samples: rho_hat is NA, weights are 1 anyway
  one <- .d2_run(df[!duplicated(df[c("subject_id", "peptide_id")]), ])
  est1 <- .ph_delta2_rho("pooled",
    bitset_raw = one$prep$bitset_raw, n_words = one$prep$bitset_words,
    pep_cols = seq_len(one$prep$bitset_m),
    sample_rows = seq_along(one$prep$samples_order),
    sample_subject = one$prep$sample_subject
  )
  expect_true(is.na(est1$rho_hat))
  expect_identical(est1$rho, 0)
})

test_that("pooling a 1-df term gives the squared Stouffer z", {
  set.seed(17)
  z <- matrix(rnorm(30))
  w <- runif(30)
  p <- .ph_delta2_pool(z, w)
  expect_equal(p$T, (sum(w * z) / sqrt(sum(w^2)))^2, tolerance = 1e-12)
  r <- matrix(rnorm(90), 30)
  p3 <- .ph_delta2_pool(r, w)
  expect_length(p3$E, 3L)
  expect_equal(p3$T, sum(p3$E^2))
})

test_that("design lists terms, df and layers", {
  subj <- data.frame(
    subject_id = sprintf("s%d", 1:6),
    group = rep(c("a", "b", "c"), 2), stringsAsFactors = FALSE
  )
  samp <- data.frame(
    sample_id = sprintf("x%d", 1:12),
    subject_id = rep(subj$subject_id, 2),
    tp = rep(c("T1", "T2"), each = 6), stringsAsFactors = FALSE
  )
  d <- .ph_delta2_design(subj, samp, "group", "multi", time_col = "tp")
  expect_identical(d$terms, c("group", "time", "interaction"))
  expect_identical(d$df, c(group = 2L, time = 1L, interaction = 2L))
  expect_identical(
    d$layer,
    c(group = "between", time = "within", interaction = "within")
  )
  expect_identical(d$time_code, rep(1:2, each = 6))

  d1 <- .ph_delta2_design(subj, samp, "group",
    group_scores = c(a = 1, b = 2, c = 3)
  )
  expect_identical(d1$terms, "group")
  expect_identical(d1$df, c(group = 1L))
  expect_identical(d1$group_scores, c(1, 2, 3))

  expect_error(.ph_delta2_design(subj, samp, "group", "binary"), "multi")
  cont <- subj
  cont$group <- c(0.5, 1.2, 3.3, 4.1, 5.9, 6.2)
  expect_error(.ph_delta2_design(cont, samp, "group"), "continuous")

  # min_samples_per_subject drops subjects
  d2 <- .ph_delta2_design(subj, samp[-(1:2), ], "group", "multi",
    min_samples_per_subject = 2
  )
  expect_identical(is.na(d2$group_code), rep(c(TRUE, FALSE), c(2, 4)))

  one <- subj
  one$group <- "a"
  expect_identical(.ph_delta2_design(one, samp, "group")$terms, character(0))
})

test_that("dropout diagnostics warn when visit counts follow group", {
  n_i <- c(4, 4, 4, 4, 2, 2, 2, 2)
  g <- rep(1:2, each = 4)
  expect_warning(
    out <- .ph_delta2_dropout(n_i, g, subject_hits = 3 * n_i),
    "differ by group"
  )
  expect_equal(out$cor_group, 1)
  quiet <- .ph_delta2_dropout(rep(3, 8), g, subject_hits = 1:8)
  expect_true(is.na(quiet$cor_group))
  expect_identical(as.integer(quiet$n_dist), 8L)
})

test_that("the group z is compute_delta()'s score z times sqrt((N-1)/N)", {
  skip_on_cran()
  set.seed(18)
  N <- 40
  samples <- data.frame(
    sample_id = sprintf("x%02d", 1:N), subject_id = sprintf("s%02d", 1:N),
    group = rep(c("A", "B"), c(17, 23)), stringsAsFactors = FALSE
  )
  df <- .d2_long(samples, matrix(rbinom(N * 4, 1, 0.3), N))
  run <- .d2_run(df)

  # one peptide per call: the Stouffer T_obs of a single peptide is its z
  z_v1 <- vapply(run$prep$peptides_order, function(p) {
    compute_delta(df[df$peptide_id == p, ],
      rank_cols = "peptide_id", group_cols = "group",
      stat_mode = "score", weight_mode = "equal",
      aggregate_stat = "stouffer", strat_bins = 0, winsor_z = 1e6,
      B_permutations = 100
    )$T_obs
  }, numeric(1))
  expect_equal(
    run$res$group[, 1], unname(z_v1) * sqrt((N - 1) / N),
    tolerance = 1e-10
  )
})
