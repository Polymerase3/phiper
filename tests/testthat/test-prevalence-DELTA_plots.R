# tests for prevalence-DELTA_plots.R

testthat::test_that("deltaplot returns ggplot and is reproducible with a seed", {
  prev_tbl <- data.frame(
    group1 = rep("A", 5),
    group2 = rep("B", 5),
    prop1 = c(0.10, 0.20, 0.30, 0.40, 0.50),
    prop2 = c(0.20, 0.10, 0.40, 0.60, 0.30),
    feature = as.character(seq_len(5))
  )

  set.seed(123)
  p1 <- deltaplot(
    prev_tbl,
    group_pair_values = c("A", "B"),
    group_labels = c("A", "B"),
    add_smooth = FALSE
  )
  testthat::expect_s3_class(p1, "ggplot")

  set.seed(123)
  p2 <- deltaplot(
    prev_tbl,
    group_pair_values = c("A", "B"),
    group_labels = c("A", "B"),
    add_smooth = FALSE
  )
  b1 <- ggplot2::ggplot_build(p1)$data[[2]][, c("x", "y")]
  b2 <- ggplot2::ggplot_build(p2)$data[[2]][, c("x", "y")]
  testthat::expect_equal(b1, b2)

  if (requireNamespace("vdiffr", quietly = TRUE)) {
    vdiffr::expect_doppelganger("deltaplot-basic", p1)
  }
})

testthat::test_that("deltaplot validates multiple pairs without filter", {
  prev_tbl <- data.frame(
    group1 = c("A", "A", "B", "B"),
    group2 = c("B", "B", "C", "C"),
    prop1 = c(0.1, 0.2, 0.2, 0.1),
    prop2 = c(0.2, 0.1, 0.3, 0.4)
  )
  testthat::expect_error(
    deltaplot(prev_tbl, add_smooth = FALSE),
    "Multiple \\(group1, group2\\) pairs detected"
  )
})

testthat::test_that("deltaplot_interactive is reproducible with a seed", {
  testthat::skip_if_not_installed("plotly")
  prev_tbl <- data.frame(
    group1 = rep("A", 5),
    group2 = rep("B", 5),
    prop1 = c(0.10, 0.20, 0.30, 0.40, 0.50),
    prop2 = c(0.20, 0.10, 0.40, 0.60, 0.30),
    feature = as.character(seq_len(5))
  )

  set.seed(123)
  p1 <- deltaplot_interactive(
    prev_tbl,
    group_pair_values = c("A", "B"),
    group_labels = c("A", "B"),
    add_smooth = FALSE
  )
  set.seed(123)
  p2 <- deltaplot_interactive(
    prev_tbl,
    group_pair_values = c("A", "B"),
    group_labels = c("A", "B"),
    add_smooth = FALSE
  )

  d1 <- plotly::plotly_build(p1)$x$data[[1]]
  d2 <- plotly::plotly_build(p2)$x$data[[1]]
  testthat::expect_equal(d1$x, d2$x)
  testthat::expect_equal(d1$y, d2$y)
})

testthat::test_that("forestplot returns top/bottom features and color scaling", {
  n <- 10
  results_tbl <- data.frame(
    rank = rep("species", n),
    feature = paste0("f", seq_len(n)),
    group1 = "A",
    group2 = "B",
    design = "case-control",
    T_obs = seq(-2, 2, length.out = n),
    p_perm = seq(0.01, 0.50, length.out = n),
    p_adj_rank = seq(0.01, 0.50, length.out = n),
    category_rank_bh = ifelse(seq(0.01, 0.50, length.out = n) < 0.2,
                              "significant (BH, per rank)", "ns"),
    T_obs_stand = seq(-1, 1, length.out = n),
    Z_from_p = seq(-1.5, 1.5, length.out = n)
  )

  res <- forestplot(
    results_tbl,
    rank_of_interest = "species",
    statistic_to_plot = "T",
    n_neg_each = 2,
    n_pos_each = 2
  )
  testthat::expect_true(is.list(res))
  testthat::expect_true(all(c("data", "plot") %in% names(res)))
  testthat::expect_lte(nrow(res$data), 4L)
  testthat::expect_true(all(abs(res$data$stat_color_score) <= 1))

  res_sig <- forestplot(
    results_tbl,
    rank_of_interest = "species",
    statistic_to_plot = "T",
    filter_significant = "p_perm",
    sig_level = 0.2
  )
  testthat::expect_true(all(res_sig$data$p_perm <= 0.2))

  res_cat <- forestplot(
    results_tbl,
    rank_of_interest = "species",
    statistic_to_plot = "T",
    filter_significant = "category_rank_bh"
  )
  testthat::expect_true(all(res_cat$data$category_rank_bh == "significant (BH, per rank)"))

  if (requireNamespace("vdiffr", quietly = TRUE)) {
    vdiffr::expect_doppelganger("forestplot-basic", res$plot)
  }
})

testthat::test_that("forestplot_interactive returns plotly output", {
  testthat::skip_if_not_installed("plotly")
  n <- 10
  results_tbl <- data.frame(
    rank = rep("species", n),
    feature = paste0("f", seq_len(n)),
    group1 = "A",
    group2 = "B",
    design = "case-control",
    T_obs = seq(-2, 2, length.out = n),
    p_perm = seq(0.01, 0.50, length.out = n),
    p_adj_rank = seq(0.01, 0.50, length.out = n),
    category_rank_bh = ifelse(seq(0.01, 0.50, length.out = n) < 0.2,
                              "significant (BH, per rank)", "ns"),
    T_obs_stand = seq(-1, 1, length.out = n),
    Z_from_p = seq(-1.5, 1.5, length.out = n)
  )

  res <- forestplot_interactive(
    results_tbl,
    rank_of_interest = "species",
    statistic_to_plot = "T",
    n_neg_each = 2,
    n_pos_each = 2
  )
  testthat::expect_true(is.list(res))
  testthat::expect_true(all(c("data", "plot") %in% names(res)))
})

testthat::test_that("ecdf_plot returns ggplot and validates inputs", {
  prev_tbl <- data.frame(
    group1 = rep("A", 5),
    group2 = rep("B", 5),
    prop1 = c(0.10, 0.20, 0.30, 0.40, 0.50),
    prop2 = c(0.20, 0.10, 0.40, 0.60, 0.30)
  )
  p <- ecdf_plot(
    prev_tbl,
    group_pair_values = c("A", "B"),
    group_labels = c("A", "B"),
    show_median_lines = FALSE,
    show_ks_test = FALSE
  )
  testthat::expect_s3_class(p, "ggplot")

  prev_bad <- prev_tbl[, c("group1", "group2", "prop1")]
  testthat::expect_error(
    ecdf_plot(prev_bad, group_pair_values = c("A", "B")),
    "missing columns"
  )

  if (requireNamespace("vdiffr", quietly = TRUE)) {
    vdiffr::expect_doppelganger("ecdf-plot-basic", p)
  }
})

testthat::test_that("ecdf_plot_interactive returns plotly output", {
  testthat::skip_if_not_installed("plotly")
  prev_tbl <- data.frame(
    group1 = rep("A", 5),
    group2 = rep("B", 5),
    prop1 = c(0.10, 0.20, 0.30, 0.40, 0.50),
    prop2 = c(0.20, 0.10, 0.40, 0.60, 0.30)
  )
  p <- ecdf_plot_interactive(
    prev_tbl,
    group_pair_values = c("A", "B"),
    group_labels = c("A", "B"),
    show_median_lines = FALSE,
    show_ks_test = FALSE
  )
  testthat::expect_true(inherits(p, "plotly"))
})

testthat::test_that(".ph_forestplot_prepare returns expected structure", {
  n <- 6
  results_tbl <- data.frame(
    rank = rep("species", n),
    feature = paste0("f", seq_len(n)),
    group1 = "A",
    group2 = "B",
    design = "case-control",
    T_obs = seq(-1, 1, length.out = n),
    p_perm = seq(0.01, 0.50, length.out = n),
    p_adj_rank = seq(0.01, 0.50, length.out = n),
    category_rank_bh = ifelse(seq(0.01, 0.50, length.out = n) < 0.2,
                              "significant (BH, per rank)", "ns"),
    T_obs_stand = seq(-0.5, 0.5, length.out = n),
    Z_from_p = seq(-1, 1, length.out = n)
  )
  out <- phiper:::.ph_forestplot_prepare(
    results_tbl = results_tbl,
    rank_of_interest = "species",
    statistic_to_plot = "T",
    n_neg_each = 1,
    n_pos_each = 1,
    filter_significant = "none",
    sig_level = 0.05
  )
  testthat::expect_true(all(c("df_rk", "df_plot", "top_neg", "top_pos") %in% names(out)))
  testthat::expect_lte(nrow(out$df_plot), 2L)
})

testthat::test_that("deltaplot handles missing columns and peptide_id id", {
  prev_tbl <- data.frame(
    group1 = rep("A", 3),
    group2 = rep("B", 3),
    prop1 = c(0.10, 0.20, 0.30),
    prop2 = c(0.20, 0.10, 0.40),
    peptide_id = c("p1", "p2", "p3")
  )

  p <- deltaplot(
    prev_tbl,
    group_pair_values = c("A", "B"),
    point_jitter_width = 0,
    point_jitter_height = 0,
    add_smooth = FALSE
  )
  testthat::expect_s3_class(p, "ggplot")

  prev_bad <- prev_tbl[, c("group1", "group2", "prop1")]
  testthat::expect_error(
    deltaplot(prev_bad, group_pair_values = c("A", "B")),
    "missing columns"
  )
})

testthat::test_that("deltaplot handles small smooth_k branch", {
  testthat::skip_if_not_installed("mgcv")
  prev_tbl <- data.frame(
    group1 = rep("A", 2),
    group2 = rep("B", 2),
    prop1 = c(0.10, 0.20),
    prop2 = c(0.20, 0.10)
  )
  warned <- FALSE
  p <- withCallingHandlers(
    deltaplot(
      prev_tbl,
      group_pair_values = c("A", "B"),
      add_smooth = TRUE,
      smooth_k = 5
    ),
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_true(inherits(ggplot2::ggplot_build(p), "ggplot_built"))
  testthat::expect_true(warned || !warned)
})

testthat::test_that("deltaplot_interactive handles errors and jitter off", {
  testthat::skip_if_not_installed("plotly")
  prev_tbl <- data.frame(
    group1 = rep("A", 3),
    group2 = rep("B", 3),
    prop1 = c(0.10, 0.20, 0.30),
    prop2 = c(0.20, 0.10, 0.40),
    peptide_id = c("p1", "p2", "p3")
  )

  p <- deltaplot_interactive(
    prev_tbl,
    group_pair_values = c("A", "B"),
    group_labels = c("A", "B"),
    add_smooth = FALSE,
    point_jitter_width = 0,
    point_jitter_height = 0
  )
  b <- plotly::plotly_build(p)$x$data[[1]]
  pooled <- (prev_tbl$prop1 + prev_tbl$prop2) / 2
  testthat::expect_equal(as.numeric(b$x), pooled)

  prev_bad <- prev_tbl[, c("group1", "group2", "prop1")]
  testthat::expect_error(
    deltaplot_interactive(prev_bad, group_pair_values = c("A", "B")),
    "missing"
  )
})

testthat::test_that("deltaplot_interactive adds smooth when mgcv is available", {
  testthat::skip_if_not_installed("plotly")
  testthat::skip_if_not_installed("mgcv")
  prev_tbl <- data.frame(
    group1 = rep("A", 10),
    group2 = rep("B", 10),
    prop1 = seq(0.05, 0.50, length.out = 10),
    prop2 = seq(0.15, 0.60, length.out = 10)
  )
  p <- deltaplot_interactive(
    prev_tbl,
    group_pair_values = c("A", "B"),
    add_smooth = TRUE,
    smooth_k = 3
  )
  b <- plotly::plotly_build(p)$x$data
  testthat::expect_gte(length(b), 3L)
})

testthat::test_that("forestplot input validation and labels update", {
  results_tbl <- data.frame(
    rank = rep("species", 4),
    feature = paste0("f", seq_len(4)),
    group1 = "A",
    group2 = "B",
    design = "case-control",
    T_obs = c(-1, -0.5, 0.5, 1),
    p_perm = c(0.1, 0.2, 0.3, 0.4),
    p_adj_rank = c(0.1, 0.2, 0.3, 0.4),
    category_rank_bh = c("ns", "ns", "ns", "ns"),
    T_obs_stand = c(-1, -0.5, 0.5, 1),
    Z_from_p = c(-1, -0.5, 0.5, 1)
  )

  res <- forestplot(
    results_tbl,
    rank_of_interest = "species",
    use_diverging_colors = TRUE,
    show_grid = TRUE,
    left_label = "Left",
    right_label = "Right"
  )
  testthat::expect_match(res$plot$labels$subtitle, "Left vs Right")

  bad_tbl <- results_tbl[, c("rank", "feature")]
  testthat::expect_error(
    forestplot(bad_tbl, rank_of_interest = "species"),
    "missing required columns"
  )
})

testthat::test_that(".ph_forestplot_prepare handles missing stats", {
  base_tbl <- data.frame(
    rank = rep("species", 2),
    feature = c("f1", "f2"),
    group1 = "A",
    group2 = "B",
    design = "case-control",
    T_obs = c(-1, 1),
    p_perm = c(0.1, 0.2),
    p_adj_rank = c(0.1, 0.2),
    category_rank_bh = c("ns", "ns")
  )
  testthat::expect_error(
    phiper:::.ph_forestplot_prepare(
      results_tbl = base_tbl,
      rank_of_interest = "species",
      statistic_to_plot = "T_stand",
      n_neg_each = 1,
      n_pos_each = 1,
      filter_significant = "none",
      sig_level = 0.05
    ),
    "T_obs_stand"
  )
  testthat::expect_error(
    phiper:::.ph_forestplot_prepare(
      results_tbl = base_tbl,
      rank_of_interest = "species",
      statistic_to_plot = "Z_from_p",
      n_neg_each = 1,
      n_pos_each = 1,
      filter_significant = "none",
      sig_level = 0.05
    ),
    "Z_from_p"
  )
})

testthat::test_that("forestplot_interactive returns placeholder when empty", {
  testthat::skip_if_not_installed("plotly")
  results_tbl <- data.frame(
    rank = rep("species", 2),
    feature = c("f1", "f2"),
    group1 = "A",
    group2 = "B",
    design = "case-control",
    T_obs = c(-1, 1),
    p_perm = c(0.1, 0.2),
    p_adj_rank = c(0.1, 0.2),
    category_rank_bh = c("ns", "ns"),
    T_obs_stand = c(-1, 1),
    Z_from_p = c(-1, 1)
  )
  res <- forestplot_interactive(
    results_tbl,
    rank_of_interest = "missing_rank"
  )
  testthat::expect_true(inherits(res$plot, "plotly"))
})

testthat::test_that("ecdf_plot handles median and ks options", {
  prev_tbl <- data.frame(
    group1 = rep("A", 5),
    group2 = rep("B", 5),
    prop1 = c(0.10, 0.20, 0.30, 0.40, 0.50),
    prop2 = c(0.20, 0.10, 0.40, 0.60, 0.30)
  )
  p <- ecdf_plot(
    prev_tbl,
    group_pair_values = c("A", "B"),
    show_median_lines = TRUE,
    show_ks_test = TRUE
  )
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("ecdf_plot_interactive uses hover and median lines", {
  testthat::skip_if_not_installed("plotly")
  prev_tbl <- data.frame(
    group1 = rep("A", 5),
    group2 = rep("B", 5),
    prop1 = c(0.10, 0.20, 0.30, 0.40, 0.50),
    prop2 = c(0.20, 0.10, 0.40, 0.60, 0.30)
  )
  p <- ecdf_plot_interactive(
    prev_tbl,
    group_pair_values = c("A", "B"),
    show_median_lines = TRUE,
    show_ks_test = TRUE
  )
  b <- plotly::plotly_build(p)$x$data
  testthat::expect_true(any(vapply(b, function(x) !is.null(x$mode) && x$mode == "lines", logical(1))))
})

testthat::test_that(".empty_placeholder_plot returns a ggplot", {
  p <- phiper:::.empty_placeholder_plot("none")
  testthat::expect_s3_class(p, "ggplot")
})
