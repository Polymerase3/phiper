# tests/testthat/test-delta2_plots.R

.d2p_data <- function(seed = 1, n = 40, m = 12) {
  set.seed(seed)
  samples <- data.frame(
    sample_id = c(sprintf("s%02d_T1", 1:n), sprintf("s%02d_T2", 1:n)),
    subject_id = rep(sprintf("s%02d", 1:n), 2),
    timepoint = rep(c("T1", "T2"), each = n),
    group = rep(rep(c("ctrl", "case"), each = n / 2), 2),
    arm = rep(rep(c("x", "y", "z", "x"), length.out = n), 2),
    stringsAsFactors = FALSE
  )
  df <- merge(samples, data.frame(peptide_id = sprintf("p%02d", 1:m)))
  up <- df$group == "case" & df$peptide_id %in% sprintf("p%02d", 1:4)
  df$exist <- stats::rbinom(nrow(df), 1, 0.2 + 0.3 * up)
  df
}

test_that("$peptides holds the per-level prevalences", {
  df <- .d2p_data()
  set.seed(1)
  res <- compute_delta2(df, "peptide_id", c("group", "arm"),
    time_col = "timepoint", B_permutations = 100L
  )
  pep <- res$peptides
  expect_named(pep, c(
    "group_col", "peptide_id", "factor", "level", "prevalence",
    "prevalence_adj", "n"
  ))
  expect_setequal(unique(pep$group_col), c("group", "arm"))

  # group: mean over subjects of their share of positive samples
  sub <- stats::aggregate(exist ~ subject_id + group + peptide_id, df, mean)
  ref <- stats::aggregate(exist ~ group + peptide_id, sub, mean)
  got <- pep[pep$group_col == "group" & pep$factor == "group", ]
  i <- match(paste(ref$group, ref$peptide_id), paste(got$level, got$peptide_id))
  expect_equal(got$prevalence[i], ref$exist, tolerance = 1e-12)
  expect_equal(got$prevalence_adj, got$prevalence)
  expect_equal(unique(got$n), 20)

  # time: share of positive samples at the timepoint
  ref_t <- stats::aggregate(exist ~ timepoint + peptide_id, df, mean)
  got_t <- pep[pep$group_col == "group" & pep$factor == "time", ]
  i <- match(
    paste(ref_t$timepoint, ref_t$peptide_id),
    paste(got_t$level, got_t$peptide_id)
  )
  expect_equal(got_t$prevalence[i], ref_t$exist, tolerance = 1e-12)

  # restricted features: only their peptides
  set.seed(1)
  few <- compute_delta2(df, "peptide_id", "group",
    rank_feature_keep = list(peptide_id = c("p01", "p02")),
    B_permutations = 100L
  )
  expect_setequal(unique(few$peptides$peptide_id), c("p01", "p02"))
  expect_false(any(few$peptides$factor == "time"))
})

test_that("prevalence adapter gives compute_pop()'s shape", {
  df <- .d2p_data()
  set.seed(1)
  res <- compute_delta2(df, "peptide_id", c("group", "arm"),
    time_col = "timepoint", B_permutations = 100L
  )
  expect_error(.ph_delta2_prev_tbl(res), "several group_cols")
  g <- .ph_delta2_prev_tbl(res, group_col = "group")
  expect_identical(nrow(g), 12L)
  expect_identical(unique(paste(g$group1, g$group2)), "case ctrl")
  pep <- res$peptides
  sel <- pep$group_col == "group" & pep$factor == "group" & pep$level == "ctrl"
  expect_equal(g$prop2, pep$prevalence[sel][match(g$peptide_id, pep$peptide_id[sel])])

  a <- .ph_delta2_prev_tbl(res, group_col = "arm")
  expect_identical(nrow(a), 36L)
  t <- .ph_delta2_prev_tbl(res, "time", group_col = "arm")
  expect_identical(unique(paste(t$group1, t$group2)), "T1 T2")

  res$peptides <- NULL
  expect_error(.ph_delta2_prev_tbl(res, group_col = "group"), "return_peptides")
})

test_that("forest adapter maps terms and pairs, rejects omnibus terms", {
  df <- .d2p_data()
  set.seed(1)
  res <- compute_delta2(df, "peptide_id", c("group", "arm"),
    time_col = "timepoint", B_permutations = 100L
  )
  f <- .ph_delta2_forest_tbl(res, "group", group_col = "group")
  tt <- res$terms[res$terms$group_col == "group", ]
  expect_equal(f$T_obs, tt$T_group)
  expect_equal(f$p_perm, tt$p_group)
  expect_equal(f$T_obs_stand, tt$T_stand_group)
  expect_identical(unique(f$group1), "case")
  expect_identical(attr(f, "labels"), c("More in case", "More in ctrl"))

  ft <- .ph_delta2_forest_tbl(res, "time", group_col = "group")
  expect_identical(attr(ft, "labels"), c("Higher at T1", "Higher at T2"))
  fi <- .ph_delta2_forest_tbl(res, "interaction", group_col = "group")
  expect_equal(fi$T_obs, tt$T_interaction)
  expect_identical(
    attr(fi, "labels"), c("Rises more in case", "Rises more in ctrl")
  )

  expect_error(
    .ph_delta2_forest_tbl(res, "group", group_col = "arm"), "omnibus"
  )
  expect_error(
    .ph_delta2_forest_tbl(res, "pairwise", group_col = "arm"), "Several pairs"
  )
  fp <- .ph_delta2_forest_tbl(res, "pairwise",
    group_col = "arm", group_pair_values = c("x", "z")
  )
  pw <- res$pairwise[res$pairwise$group_col == "arm" &
    res$pairwise$group1 == "x" & res$pairwise$group2 == "z", ]
  expect_equal(fp$T_obs, pw$T_obs)
  expect_equal(fp$p_maxT, pw$p_maxT)
  expect_equal(
    fp$Z_from_p,
    sign(pw$T_obs) * stats::qnorm(pw$p_perm / 2, lower.tail = FALSE)
  )
})

test_that("plotters take a delta2_result", {
  df <- .d2p_data()
  set.seed(1)
  res <- compute_delta2(df, "peptide_id", "group",
    time_col = "timepoint", B_permutations = 100L
  )
  p <- deltaplot(res, add_smooth = FALSE)
  expect_s3_class(p, "ggplot")
  expect_identical(nrow(p$data), 12L)
  expect_s3_class(deltaplot(res, term = "time", add_smooth = FALSE), "ggplot")
  expect_s3_class(ecdf_plot(res, adjusted = TRUE), "ggplot")

  fp <- forestplot(res, rank_of_interest = "peptide_id", n_neg_each = 3,
    n_pos_each = 3
  )
  expect_s3_class(fp$plot, "ggplot")
  expect_match(fp$plot$labels$subtitle, "More in case vs More in ctrl")
  fp2 <- forestplot(res, "peptide_id", term = "time", left_label = "L",
    right_label = "R"
  )
  expect_match(fp2$plot$labels$subtitle, "L vs R")

  skip_if_not_installed("plotly")
  expect_s3_class(deltaplot_interactive(res, add_smooth = FALSE), "plotly")
  expect_s3_class(ecdf_plot_interactive(res), "plotly")
  expect_s3_class(
    forestplot_interactive(res, "peptide_id", term = "interaction")$plot,
    "plotly"
  )
})
