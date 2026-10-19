# 6.3 Unbalanced designs and term orthogonality.
#
# (a) Visit counts correlated with group: cases have all 4 visits, controls
#     2 of them at random; positivity rises over visits (a real time main
#     effect) but there is no interaction. Without the 3.5 residualisation
#     the interaction over-rejects (the time effect leaks into it); with it
#     the interaction must be calibrated. Acceptance test of that fix.
#     Reported alongside, not gated: controls seen at T1-T2 only, where the
#     interaction at T3-T4 is not estimable; the test is then valid but
#     conservative, and compute_delta2() must warn about the gap.
# (b) Informative dropout: subjects leave after a visit with many positives
#     (a "seroconversion"), under no time effect. Quantifies how far the time
#     term drifts and whether the dropout diagnostics flag it.

val_s3_unbalanced <- function(n_reps = 1000, B = 500L, seed = 630000) {
  band <- val_band(n_reps)

  # (a)
  random_two <- function(g, i) if (g == 2L) 1:4 else sort(sample(1:4, 2))
  first_two <- function(g, i) if (g == 2L) 1:4 else 1:2
  a <- val_reps(n_reps, function(r) {
    df <- val_sim(
      n_subjects = 60, K = 2, visits = random_two, m = 100,
      rho = 0.3, time_effect = 0.4
    )
    val_run_interaction(df, B = B)
  }, seed = seed)
  n_fixed <- max(100, n_reps %/% 5)
  a2 <- val_reps(n_fixed, function(r) {
    df <- val_sim(
      n_subjects = 60, K = 2, visits = first_two, m = 100,
      rho = 0.3, time_effect = 0.4
    )
    out <- val_run_interaction(df, B = B)
    out$warned <- grepl("not estimable", val_run(df, B = 100L)$warnings)
    out
  }, seed = seed + 50000)
  rates_a <- tibble::tibble(
    scenario = c(
      rep("controls miss 2 of 4 visits at random", 2),
      rep("controls seen at T1-T2 only", 2)
    ),
    statistic = rep(
      c("interaction, residualised", "interaction, not residualised"), 2
    ),
    rejection = c(
      mean(a$p_residualised <= 0.05), mean(a$p_raw <= 0.05),
      mean(a2$p_residualised <= 0.05), mean(a2$p_raw <= 0.05)
    ),
    n = c(n_reps, n_reps, n_fixed, n_fixed)
  )
  warned <- mean(a2$warned)

  # (b)
  b <- dplyr::bind_rows(lapply(c(FALSE, TRUE), function(drop) {
    res <- val_reps(n_reps, function(r) {
      df <- val_sim(
        n_subjects = 60, K = 2, visits = 4L, m = 100, rho = 0.3,
        dropout = if (drop) 0.75
      )
      val_run(df, B = B)
    }, seed = seed + 100000 + 50000 * drop)
    res$dropout <- drop
    res
  }))
  flagged <- grepl("track breadth|differ by group", b$warnings)
  rates_b <- b |>
    dplyr::mutate(flagged = flagged) |>
    dplyr::group_by(.data$dropout) |>
    dplyr::summarise(
      time_rejection = mean(.data$p_time <= 0.05),
      interaction_rejection = mean(.data$p_interaction <= 0.05),
      flagged = mean(.data$flagged), .groups = "drop"
    )

  in_band <- function(x) x >= band[["lo"]] & x <= band[["hi"]]
  nd <- rates_b[!rates_b$dropout, ]
  dr <- rates_b[rates_b$dropout, ]
  gates <- tibble::tibble(
    section = "6.3",
    gate = c(
      "(a) residualised interaction calibrated",
      "(a) without residualisation the interaction over-rejects",
      "(a) controls at T1-T2 only: gap warned in every run",
      "(b) no dropout: time term calibrated",
      "(b) informative dropout flagged by the diagnostics in >= 80% of runs",
      "(b) no dropout: diagnostics silent in >= 95% of runs"
    ),
    value = c(
      rates_a$rejection[1], rates_a$rejection[2], warned,
      nd$time_rejection, dr$flagged, 1 - nd$flagged
    ),
    lo = c(band[["lo"]], band[["hi"]], 1, band[["lo"]], 0.8, 0.95),
    hi = c(band[["hi"]], 1, 1, band[["hi"]], 1, 1)
  )
  gates$pass <- gates$value >= gates$lo & gates$value <= gates$hi
  gates$pass[2] <- gates$value[2] > band[["hi"]]
  list(
    reps_a = a, reps_a2 = a2, rates_a = rates_a, reps_b = b,
    rates_b = rates_b, gates = gates
  )
}

# 6.4 Interaction robustness: one group has more turnover (higher
# prevalence, so more gained and lost) and the groups are very unequal;
# no time effect, no interaction.
val_s4_turnover <- function(n_reps = 1000, B = 500L, seed = 640000) {
  band <- val_band(n_reps)
  sizes <- list("10 high-turnover vs 50" = c(50, 10), "50 high-turnover vs 10" = c(10, 50))
  res <- dplyr::bind_rows(lapply(seq_along(sizes), function(s) {
    out <- val_reps(n_reps, function(r) {
      df <- val_sim(
        K = 2, group_sizes = sizes[[s]], visits = 2L, m = 100, rho = 0.3,
        group_effect = 1.5, frac_affected = 1
      )
      val_run(df, B = B)
    }, seed = seed + 100000 * s)
    out$scenario <- names(sizes)[s]
    out
  }))
  rates <- res |>
    dplyr::group_by(.data$scenario) |>
    dplyr::summarise(
      interaction_rejection = mean(.data$p_interaction <= 0.05),
      time_rejection = mean(.data$p_time <= 0.05), .groups = "drop"
    )
  gates <- tibble::tibble(
    section = "6.4",
    gate = sprintf("%s: interaction calibrated", rates$scenario),
    value = rates$interaction_rejection, lo = band[["lo"]], hi = band[["hi"]]
  )
  gates$pass <- gates$value >= gates$lo & gates$value <= gates$hi
  list(reps = res, rates = rates, gates = gates)
}
