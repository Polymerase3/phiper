# 6.5 The intraclass correlation rho.
#
# (a) Recovery: simulate a known within-subject correlation with unequal
#     samples per subject (1-6) and estimate it with the pooled ANOVA
#     estimator; bias, spread and clamping at 0 and 1.
# (b) Power of rho = "pooled" vs "one" vs "zero" for a group effect, at low
#     and high true rho with unequal sample counts; same data and draws.
# (c) Real data: rho per prevalence stratum on the phiper example data
#     (phip_mixture), to decide the rho_by_stratum default.

val_unequal_visits <- function(g, i) seq_len(sample.int(6, 1))

val_s5_rho <- function(n_reps = 200, n_reps_power = 300, B = 500L,
                       data_path = NULL, seed = 650000) {
  ns <- asNamespace("phiper")

  # (a)
  rhos <- c(0, 0.05, 0.2, 0.5, 0.8, 0.95, 1)
  rec <- dplyr::bind_rows(lapply(seq_along(rhos), function(k) {
    out <- val_reps(n_reps, function(r) {
      df <- val_sim(
        n_subjects = 60, K = 1, visits = val_unequal_visits, m = 200,
        rho = rhos[k]
      )
      prep <- ns$.ph_delta2_prepare(df, "group")
      est <- ns$cpp_delta2_rho(
        prep$bitset_raw, prep$bitset_words, seq_len(prep$bitset_m),
        seq_along(prep$sample_subject), prep$sample_subject
      )
      tibble::tibble(rho_hat = est$rho_hat)
    }, seed = seed + 1000 * k)
    out$rho <- rhos[k]
    out
  }))
  recovery <- rec |>
    dplyr::group_by(.data$rho) |>
    dplyr::summarise(
      mean = mean(.data$rho_hat), bias = mean(.data$rho_hat) - .data$rho[1],
      sd = stats::sd(.data$rho_hat), at_0 = mean(.data$rho_hat == 0),
      at_1 = mean(.data$rho_hat == 1), .groups = "drop"
    )

  # (b)
  pw <- dplyr::bind_rows(lapply(c(0.2, 0.8), function(true_rho) {
    out <- val_reps(n_reps_power, function(r) {
      df <- val_sim(
        n_subjects = 60, K = 2, visits = val_unequal_visits, m = 100,
        rho = true_rho, group_effect = 0.6, frac_affected = 0.2
      )
      s <- sample.int(.Machine$integer.max, 1L)
      dplyr::bind_rows(lapply(c("pooled", "one", "zero"), function(opt) {
        set.seed(s) # same draws for every rho option
        p <- val_run(df, B = B, time_col = NULL, rho = opt)$p_group
        tibble::tibble(option = opt, p = p)
      }))
    }, seed = seed + 50000 + 10000 * true_rho)
    out$true_rho <- true_rho
    out
  }))
  power <- pw |>
    dplyr::group_by(.data$true_rho, .data$option) |>
    dplyr::summarise(power = mean(.data$p <= 0.05), .groups = "drop")

  # (c)
  strata <- NULL
  if (!is.null(data_path) && file.exists(data_path)) {
    real <- val_read_mixture(data_path)
    prep <- ns$.ph_delta2_prepare(real, "group", time_col = "timepoint")
    ss <- prep$sample_subject
    est <- ns$cpp_delta2_rho(
      prep$bitset_raw, prep$bitset_words, seq_len(prep$bitset_m),
      seq_along(ss), ss
    )
    sb <- c(0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50)
    bins <- ns$.ph_delta2_prevalence_bins(est$pep_tot / est$n_samples, sb)
    brks <- c(0, sb, 1)
    strata <- dplyr::bind_rows(lapply(sort(unique(bins)), function(b) {
      e <- ns$cpp_delta2_rho(
        prep$bitset_raw, prep$bitset_words, which(bins == b), seq_along(ss), ss
      )
      tibble::tibble(
        stratum = sprintf("(%.3g, %.3g]", brks[b + 1], brks[b + 2]),
        n_peptides = sum(bins == b), rho_hat = e$rho_hat
      )
    }))
    strata <- dplyr::bind_rows(
      strata,
      tibble::tibble(stratum = "pooled", n_peptides = length(bins), rho_hat = est$rho_hat)
    )
  }

  gates <- tibble::tibble(
    section = "6.5",
    gate = sprintf("(a) |bias| of rho_hat <= 0.05 at rho = %s", recovery$rho),
    value = recovery$bias, lo = -0.05, hi = 0.05
  ) |>
    dplyr::filter(recovery$rho > 0 & recovery$rho < 1)
  gates$pass <- gates$value >= gates$lo & gates$value <= gates$hi
  list(recovery = recovery, rec = rec, power = power, strata = strata, gates = gates)
}

# phip_mixture as long data for compute_delta2(): subjects are numbered within
# group, so the subject is group + id
val_read_mixture <- function(path) {
  real <- if (grepl("\\.parquet$", path)) {
    if (requireNamespace("arrow", quietly = TRUE)) {
      as.data.frame(arrow::read_parquet(path))
    } else {
      as.data.frame(nanoparquet::read_parquet(path))
    }
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE)
  }
  data.frame(
    sample_id = as.character(real$sample_id),
    subject_id = paste(real$group, real$subject_id, sep = "_"),
    group = real$group, timepoint = real$time,
    peptide_id = as.character(real$peptide_id), exist = real$exist,
    stringsAsFactors = FALSE
  )
}
