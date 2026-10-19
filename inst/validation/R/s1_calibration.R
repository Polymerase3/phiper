# 6.1 Calibration per term under a confounded null.
#
# For {2, 3 groups} x {1, 2, 4 timepoints}: group labels depend on `age`,
# positivity depends on `age` but not on group, time or their interaction.
# Each design runs with and without propensity_formula = ~ age. With it,
# every term must reject at alpha = 0.05 within Monte Carlo error; without
# it, the group term is expected to over-reject (the confounded null), while
# time and interaction are within-subject and stay calibrated.

val_s1_calibration <- function(n_reps = 1000, B = 500L, seed = 610000) {
  designs <- expand.grid(K = c(2L, 3L), T = c(1L, 2L, 4L))
  out <- list()
  for (d in seq_len(nrow(designs))) {
    K <- designs$K[d]
    T_n <- designs$T[d]
    for (ps in c(FALSE, TRUE)) {
      res <- val_reps(n_reps, function(r) {
        df <- val_sim(
          n_subjects = 60, K = K, visits = T_n, confound = 0.7,
          beta_age = 0.8, m = 100, rho = 0.3
        )
        val_run(df,
          B = B, time_col = if (T_n > 1L) "timepoint",
          propensity_formula = if (ps) ~age
        )
      }, seed = seed + 10000 * d + 5000 * ps)
      out[[length(out) + 1L]] <- tibble::tibble(
        K = K, T = T_n, propensity = ps, res
      )
    }
  }
  reps <- dplyr::bind_rows(out)

  rates <- reps |>
    tidyr::pivot_longer(dplyr::starts_with("p_"),
      names_to = "term", values_to = "p", names_prefix = "p_"
    ) |>
    dplyr::filter(!is.na(.data$p)) |>
    dplyr::group_by(.data$K, .data$T, .data$propensity, .data$term) |>
    dplyr::summarise(
      n = dplyr::n(), rejection = mean(.data$p <= 0.05), .groups = "drop"
    )
  band <- val_band(n_reps)
  rates$lo <- band[["lo"]]
  rates$hi <- band[["hi"]]
  rates$within_band <- rates$rejection >= rates$lo & rates$rejection <= rates$hi

  # gate: with the propensity model every term is calibrated; without it,
  # the within-subject terms are
  gates <- rates |>
    dplyr::filter(.data$propensity | .data$term != "group") |>
    dplyr::transmute(
      section = "6.1",
      gate = sprintf(
        "K=%d, T=%d, %s: %s rejection in band",
        .data$K, .data$T,
        ifelse(.data$propensity, "with PS", "no PS"), .data$term
      ),
      value = .data$rejection, lo = .data$lo, hi = .data$hi,
      pass = .data$within_band
    )
  list(reps = reps, rates = rates, gates = gates)
}
