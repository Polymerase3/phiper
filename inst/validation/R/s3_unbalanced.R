# 6.3 Unbalanced designs and term orthogonality.
#
# (a) Visit counts correlated with group: cases have 4 visits, controls 2;
#     positivity rises over visits (a real time main effect) but there is no
#     interaction. Without the 3.5 residualisation the interaction
#     over-rejects (the time effect leaks into it); with it the interaction
#     must be calibrated. Acceptance test of that fix.
# (b) Informative dropout: subjects leave after a visit with many positives
#     (a "seroconversion"), under no time effect. Quantifies how far the time
#     term drifts and whether the dropout diagnostics flag it.

val_s3_unbalanced <- function(n_reps = 1000, B = 500L, seed = 630000) {
  band <- val_band(n_reps)

  # (a)
  visits_by_group <- function(g, i) if (g == 2L) 1:4 else 1:2
  a <- val_reps(n_reps, function(r) {
    df <- val_sim(
      n_subjects = 60, K = 2, visits = visits_by_group, m = 100,
      rho = 0.3, time_effect = 0.4
    )
    val_run_interaction(df, B = B)
  }, seed = seed)
  rates_a <- tibble::tibble(
    scenario = "visits 4 (cases) vs 2 (controls), time effect",
    statistic = c("interaction, residualised", "interaction, not residualised"),
    rejection = c(mean(a$p_residualised <= 0.05), mean(a$p_raw <= 0.05))
  )

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
      "(b) no dropout: time term calibrated",
      "(b) informative dropout flagged by the diagnostics in >= 80% of runs",
      "(b) no dropout: diagnostics silent in >= 95% of runs"
    ),
    value = c(
      rates_a$rejection[1], rates_a$rejection[2], nd$time_rejection,
      dr$flagged, 1 - nd$flagged
    ),
    lo = c(band[["lo"]], band[["hi"]], band[["lo"]], 0.8, 0.95),
    hi = c(band[["hi"]], 1, band[["hi"]], 1, 1)
  )
  gates$pass <- gates$value >= gates$lo & gates$value <= gates$hi
  gates$pass[2] <- gates$value[2] > band[["hi"]]
  list(reps_a = a, rates_a = rates_a, reps_b = b, rates_b = rates_b, gates = gates)
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
