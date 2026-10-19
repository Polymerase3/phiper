# 6.2 Pooled power of the per-peptide statistic: score z vs srlr,
# standardised srlr and the linear Wald t.
#
# Cross-sectional, two groups, sparse prevalence (log-uniform 0.01-0.3,
# 500 peptides, many with 1-5 positives). Every statistic is computed on the
# same data and the same label permutations, winsorised at 4 and pooled by
# an equal-weight Stouffer sum over all peptides, so the comparison isolates
# the per-peptide statistic. Report only: this justifies the score default,
# it is not a release gate.

# per-peptide statistics from group-2 counts (rows = labellings)
val_s2_stats <- function(X2, tot, n1, n2) {
  N <- n1 + n2
  Tot <- matrix(tot, nrow(X2), ncol(X2), byrow = TRUE)
  X1 <- Tot - X2
  p1 <- X1 / n1
  p2 <- X2 / n2
  p0 <- Tot / N
  edge <- Tot == 0 | Tot == N
  xlogy <- function(x, y) ifelse(x == 0, 0, x * log(y))
  ll <- function(x, n, p) xlogy(x, p) + xlogy(n - x, 1 - p)

  score <- (p2 - p1) / sqrt(p0 * (1 - p0) * (1 / n1 + 1 / n2))
  score[edge] <- 0
  LR <- 2 * (ll(X1, n1, p1) + ll(X2, n2, p2) - ll(X1, n1, p0) - ll(X2, n2, p0))
  srlr <- sign(p2 - p1) * sqrt(pmax(LR, 0))
  s2 <- (n1 * p1 * (1 - p1) + n2 * p2 * (1 - p2)) / (N - 2)
  wald <- ifelse(s2 > 0, (p2 - p1) / sqrt(s2 * (1 / n1 + 1 / n2)), 0)

  # srlr standardised by its exact permutation moments given the total
  mom <- val_s2_moments(n1, n2, "srlr")
  mu <- matrix(mom$mean[Tot + 1], nrow(X2))
  sdv <- matrix(mom$sd[Tot + 1], nrow(X2))
  srlr_std <- ifelse(sdv > 0, (srlr - mu) / sdv, 0)

  list(score = score, srlr = srlr, srlr_std = srlr_std, wald_t = wald)
}

# exact permutation mean and SD of a statistic for each total 0..N
val_s2_moments <- function(n1, n2, stat) {
  N <- n1 + n2
  mom <- lapply(0:N, function(k) {
    x2 <- max(0, k - n1):min(k, n2)
    w <- stats::dhyper(x2, n2, n1, k)
    z <- val_s2_stats_raw(x2, k, n1, n2)[[stat]]
    mu <- sum(w * z)
    c(mean = mu, sd = sqrt(sum(w * (z - mu)^2)))
  })
  mom <- do.call(rbind, mom)
  list(mean = mom[, "mean"], sd = mom[, "sd"])
}

# the raw statistics for a vector of group-2 counts at one total
val_s2_stats_raw <- function(x2, k, n1, n2) {
  N <- n1 + n2
  x1 <- k - x2
  p1 <- x1 / n1
  p2 <- x2 / n2
  p0 <- k / N
  xlogy <- function(x, y) ifelse(x == 0, 0, x * log(y))
  ll <- function(x, n, p) xlogy(x, p) + xlogy(n - x, 1 - p)
  edge <- k == 0 || k == N
  score <- if (edge) 0 * x2 else (p2 - p1) / sqrt(p0 * (1 - p0) * (1 / n1 + 1 / n2))
  LR <- 2 * (ll(x1, n1, p1) + ll(x2, n2, p2) - ll(x1, n1, p0) - ll(x2, n2, p0))
  s2 <- (n1 * p1 * (1 - p1) + n2 * p2 * (1 - p2)) / (N - 2)
  list(
    score = score,
    srlr = sign(p2 - p1) * sqrt(pmax(LR, 0)),
    wald_t = ifelse(s2 > 0, (p2 - p1) / sqrt(s2 * (1 / n1 + 1 / n2)), 0)
  )
}

# exact one-peptide null moments (the background table of 6.2)
val_s2_moment_table <- function(sizes = list(c(20, 20), c(10, 30)),
                                positives = c(1, 2, 3, 5, 20)) {
  rows <- list()
  for (sz in sizes) {
    for (stat in c("score", "srlr", "wald_t")) {
      mom <- val_s2_moments(sz[1], sz[2], stat)
      rows[[length(rows) + 1L]] <- tibble::tibble(
        n1 = sz[1], n2 = sz[2], stat = stat, positives = positives,
        null_mean = mom$mean[positives + 1], null_sd = mom$sd[positives + 1]
      )
    }
  }
  dplyr::bind_rows(rows)
}

val_s2_pooled_power <- function(n_reps = 500, B = 500L, m = 500,
                                effect = 0.6, seed = 620000) {
  configs <- expand.grid(
    sizes = c("20/20", "30/10", "10/30"),
    alt = c("null", "up", "down"), stringsAsFactors = FALSE
  )
  out <- list()
  for (cfg in seq_len(nrow(configs))) {
    sz <- as.integer(strsplit(configs$sizes[cfg], "/")[[1]])
    alt <- configs$alt[cfg]
    delta <- switch(alt, null = 0, up = effect, down = -effect)
    res <- val_reps(n_reps, function(r) {
      df <- val_sim(
        K = 2, group_sizes = sz, visits = 1L, m = m, profile = "sparse",
        group_effect = delta, frac_affected = 0.2
      )
      s <- attr(df, "samples")
      Y <- attr(df, "Y")
      g2 <- as.numeric(s$group == "B")
      L <- cbind(g2, replicate(B, sample(g2)))
      X2 <- crossprod(L, Y)
      st <- val_s2_stats(X2, colSums(Y), sz[1], sz[2])
      dplyr::bind_rows(lapply(names(st), function(nm) {
        Z <- st[[nm]]
        Z[] <- pmax(-4, pmin(4, Z))
        Tp <- rowSums(Z) / sqrt(m)
        Tb <- Tp[-1]
        gt <- sum(abs(Tb) > abs(Tp[1]) + 1e-9)
        eq <- sum(abs(abs(Tb) - abs(Tp[1])) <= 1e-9)
        tibble::tibble(
          stat = nm, T_obs = Tp[1], null_mean = mean(Tb), null_sd = stats::sd(Tb),
          p = max((gt + 0.5 * eq) / B, 0.5 / B)
        )
      }))
    }, seed = seed + 1000 * cfg)
    out[[cfg]] <- tibble::tibble(sizes = configs$sizes[cfg], alt = alt, res)
  }
  reps <- dplyr::bind_rows(out)
  power <- reps |>
    dplyr::group_by(.data$sizes, .data$alt, .data$stat) |>
    dplyr::summarise(
      n = dplyr::n(), rejection = mean(.data$p <= 0.05),
      pooled_null_mean = mean(.data$null_mean), .groups = "drop"
    )
  band <- val_band(n_reps)
  gates <- power |>
    dplyr::filter(.data$alt == "null") |>
    dplyr::transmute(
      section = "6.2",
      gate = sprintf("%s, %s: permutation test calibrated", .data$sizes, .data$stat),
      value = .data$rejection, lo = band[["lo"]], hi = band[["hi"]],
      pass = .data$rejection >= band[["lo"]] & .data$rejection <= band[["hi"]]
    )
  list(
    reps = reps, power = power, moments = val_s2_moment_table(),
    gates = gates
  )
}
