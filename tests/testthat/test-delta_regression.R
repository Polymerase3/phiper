# Regression snapshot of compute_delta() output across all option combinations.
# Later refactors must leave this snapshot unchanged.

test_that("compute_delta output is unchanged across all options", {
  skip_on_cran()

  fx <- readRDS(test_path("fixtures", "delta_regression.rds"))
  designs <- list(
    unpaired = list(
      x = fx$data[fx$data$timepoint == "T1", ],
      group_cols = "group",
      paired_by = NULL,
      stat_modes = c("srlr", "diff", "asin", "score")
    ),
    paired = list(
      x = fx$data,
      group_cols = "timepoint",
      paired_by = "subject_id",
      stat_modes = c("srlr", "diff", "asin", "score", "mcnemar", "srlr_paired")
    )
  )
  strat_bins <- list(
    none = 0,
    default = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50)
  )

  res <- list()
  for (d in names(designs)) {
    des <- designs[[d]]
    grid <- expand.grid(
      stat_mode = des$stat_modes,
      weight_mode = c("equal", "se_invvar", "n_eff_sqrt"),
      aggregate_stat = c("stouffer", "maxmean", "af"),
      strat_bins = names(strat_bins),
      perm_method = c("mid_p", "standard"),
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(grid))) {
      g <- grid[i, ]
      set.seed(1)
      out <- compute_delta(
        x               = des$x,
        rank_cols       = "species",
        group_cols      = des$group_cols,
        paired_by       = des$paired_by,
        peptide_library = fx$peptide_library,
        B_permutations  = 200L,
        stat_mode       = g$stat_mode,
        weight_mode     = g$weight_mode,
        aggregate_stat  = g$aggregate_stat,
        strat_bins      = strat_bins[[g$strat_bins]],
        perm_method     = g$perm_method
      )
      res[[length(res) + 1L]] <- dplyr::bind_cols(
        tibble::tibble(design_case = d, g),
        out
      )
    }
  }
  res <- dplyr::bind_rows(res)

  expect_snapshot_value(res, style = "deparse", tolerance = 1e-12)
})
