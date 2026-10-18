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
  seed_before <- .Random.seed
  p1 <- deltaplot(
    prev_tbl,
    group_pair_values = c("A", "B"),
    group_labels = c("A", "B"),
    add_smooth = FALSE
  )
  testthat::expect_identical(.Random.seed, seed_before)
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
  seed_before <- .Random.seed
  p1 <- deltaplot_interactive(
    prev_tbl,
    group_pair_values = c("A", "B"),
    group_labels = c("A", "B"),
    add_smooth = FALSE
  )
  testthat::expect_identical(.Random.seed, seed_before)
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
  warn_env <- new.env(parent = emptyenv())
  warn_env$count <- 0L
  p <- withCallingHandlers(
    deltaplot(
      prev_tbl,
      group_pair_values = c("A", "B"),
      add_smooth = TRUE,
      smooth_k = 5
    ),
    warning = function(w) {
      warn_env$count <- warn_env$count + 1L
      invokeRestart("muffleWarning")
    }
  )
  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_true(inherits(ggplot2::ggplot_build(p), "ggplot_built"))
  testthat::expect_true(warn_env$count >= 0L)
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
    p_perm = c(0.1, 0.2)
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

# ---- profile_plot ------------------------------------------------------------

.make_profile_peplib <- function() {
  data.frame(
    peptide_id     = paste0("pep", 1:6),
    protein_seq_id = rep(c("protA", "protB"), each = 3),
    pos            = c(0, 20, 40, 0, 25, 50),
    aa_seq         = strrep("A", 30),
    len_seq        = rep(c(70, 75), each = 3),
    taxon_with_protein_unique_annotation =
      rep(c("Species a: protein A", "Species b: protein B"), each = 3)
  )
}

.make_profile_prev <- function() {
  data.frame(
    rank    = "peptide_id",
    feature = rep(c("pep1", "pep2", "pep4", "pep5"), 2),
    group1  = "control",
    group2  = rep(c("sjogren", "sle"), each = 4),
    prop1   = rep(c(0.10, 0.20, 0.05, 0.00), 2),
    prop2   = c(0.40, 0.50, 0.10, 0.20, 0.30, 0.60, 0.00, 0.10),
    p_raw   = c(0.001, 0.002, 0.5, 0.5, 0.5, 0.001, 0.5, 0.5)
  )
}

testthat::test_that("profile_plot draws every library peptide per group", {
  p <- profile_plot(
    .make_profile_prev(),
    proteins = c("protA", "protB"),
    peplib   = .make_profile_peplib()
  )
  testthat::expect_s3_class(p, "ggplot")

  d <- p$data
  testthat::expect_equal(nrow(d), 6L * 3L)
  testthat::expect_identical(levels(d$group), c("control", "sjogren", "sle"))

  # prevalences from both sides of the pairs
  testthat::expect_equal(d$prop[d$peptide_id == "pep2" & d$group == "control"], 0.20)
  testthat::expect_equal(d$prop[d$peptide_id == "pep2" & d$group == "sle"], 0.60)
  # peptides missing from prev_tbl at 0
  testthat::expect_true(all(d$prop[d$peptide_id %in% c("pep3", "pep6")] == 0))
  # segments from pos + 1, clipped at the protein end (pep6: 51-80 -> 51-75)
  testthat::expect_equal(unique(d$start[d$peptide_id == "pep2"]), 21)
  testthat::expect_equal(unique(d$end[d$peptide_id == "pep2"]), 50)
  testthat::expect_equal(unique(d$end[d$peptide_id == "pep6"]), 75)

  testthat::expect_identical(
    levels(d$panel),
    c("Species a: protein A (protA)", "Species b: protein B (protB)")
  )

  if (requireNamespace("vdiffr", quietly = TRUE)) {
    vdiffr::expect_doppelganger("profile-plot-basic", p)
  }
})

testthat::test_that("profile_plot orders default groups pair by pair", {
  prev <- data.frame(
    feature = "pep1",
    group1  = c("A", "T1"),
    group2  = c("B", "T2"),
    prop1   = c(0.1, 0.2),
    prop2   = c(0.3, 0.4)
  )
  p <- profile_plot(prev, "protA", peplib = .make_profile_peplib(),
                    show_overlap = FALSE)
  testthat::expect_identical(levels(p$data$group), c("A", "B", "T1", "T2"))
})

testthat::test_that("profile_plot uses names as titles and labels", {
  p <- profile_plot(
    .make_profile_prev(),
    proteins      = c("Protein B" = "protB", "protA"),
    groups        = c(Controls = "control", "sle"),
    peplib        = .make_profile_peplib(),
    custom_colors = c(sle = "#1baf7a", control = "#2a78d6")
  )
  testthat::expect_identical(
    levels(p$data$panel),
    c("Protein B", "Species a: protein A (protA)")
  )
  testthat::expect_identical(levels(p$data$group), c("control", "sle"))

  guide <- ggplot2::get_guide_data(p, "colour")
  testthat::expect_identical(guide$.label, c("Controls", "sle"))
  testthat::expect_identical(guide$colour, c("#2a78d6", "#1baf7a"))
})

testthat::test_that("profile_plot adds DELTA statistics to the titles", {
  delta_tbl <- data.frame(
    rank    = c("protein_seq_id", "protein_seq_id", "species"),
    feature = c("protA", "protA", "protB"),
    group1  = "control",
    group2  = c("sjogren", "sle", "sle"),
    T_obs   = c(3.1, 1.2, 9.9),
    p_perm  = c(0.002, 0.21, 0.001)
  )
  p <- profile_plot(
    .make_profile_prev(),
    proteins  = c("protA", "protB"),
    groups    = c(Controls = "control", "sjogren", "sle"),
    peplib    = .make_profile_peplib(),
    delta_tbl = delta_tbl
  )
  testthat::expect_identical(levels(p$data$panel), c(
    paste0(
      "Species a: protein A (protA)\n",
      "Controls vs sjogren\nDELTA T = 3.10, p = 0.002\n",
      "Controls vs sle\nDELTA T = 1.20, p = 0.21"
    ),
    "Species b: protein B (protB)\nDELTA: not tested"
  ))
})

testthat::test_that("profile_plot fetches the peptides' libraries without peplib", {
  ids_seen <- NULL
  testthat::local_mocked_bindings(
    .ph_library_for_peptides = function(peptide_ids) {
      ids_seen <<- peptide_ids
      .make_profile_peplib()
    }
  )
  p <- profile_plot(.make_profile_prev(), proteins = "protA")
  testthat::expect_setequal(ids_seen, c("pep1", "pep2", "pep4", "pep5"))
  testthat::expect_setequal(unique(p$data$peptide_id), c("pep1", "pep2", "pep3"))
})

testthat::test_that("profile_plot shades the longest run of overlapping hits", {
  p <- profile_plot(
    .make_profile_prev(),
    proteins = c("protA", "protB"),
    peplib   = .make_profile_peplib()
  )
  # sjogren pair: pep1 and pep2 are hits (BH q 0.004), overlap, more common in
  # sjogren; sle pair: only pep2 is a hit, a single peptide is no run
  bands <- p$layers[[1]]$data
  testthat::expect_s3_class(p$layers[[1]]$geom, "GeomRect")
  testthat::expect_equal(nrow(bands), 1L)
  testthat::expect_identical(bands$protein_seq_id, "protA")
  testthat::expect_identical(as.character(bands$group), "sjogren")
  testthat::expect_equal(c(bands$xmin, bands$xmax), c(1, 50))
  testthat::expect_identical(
    p$labels$subtitle,
    "Shaded: runs of overlapping peptides with BH q <= 0.05"
  )

  # BH within each pair: with alpha = 0.003 the sjogren q of 0.004 fails
  p <- profile_plot(.make_profile_prev(), "protA",
                    peplib = .make_profile_peplib(), alpha = 0.003)
  testthat::expect_false(any(vapply(p$layers, function(l) {
    inherits(l$geom, "GeomRect")
  }, logical(1))))
})

testthat::test_that("profile_plot colours a run by the group with higher prevalence", {
  prev <- .make_profile_prev()
  # make pep1 and pep2 more common in control than in sjogren (control's
  # prevalence is set in both pairs, which must agree)
  prev$prop1[c(1, 2, 5, 6)] <- c(0.60, 0.70, 0.60, 0.70)
  p <- profile_plot(prev, "protA", peplib = .make_profile_peplib())
  testthat::expect_identical(as.character(p$layers[[1]]$data$group), "control")
})

testthat::test_that("profile_plot draws bands only for pairs with both groups drawn", {
  p <- profile_plot(.make_profile_prev(), "protA",
                    groups = c("control", "sle"),
                    peplib = .make_profile_peplib())
  testthat::expect_false(any(vapply(p$layers, function(l) {
    inherits(l$geom, "GeomRect")
  }, logical(1))))
  testthat::expect_null(p$labels$subtitle)
})

testthat::test_that("profile_plot without show_overlap needs no p_raw", {
  prev <- .make_profile_prev()
  prev$p_raw <- NULL
  testthat::expect_error(
    profile_plot(prev, "protA", peplib = .make_profile_peplib()),
    "prev_tbl is missing required columns"
  )
  p <- profile_plot(prev, "protA", peplib = .make_profile_peplib(),
                    show_overlap = FALSE)
  testthat::expect_length(p$layers, 1L)
  testthat::expect_s3_class(p$layers[[1]]$geom, "GeomSegment")
})

testthat::test_that("profile_plot validates its inputs", {
  prev <- .make_profile_prev()
  lib <- .make_profile_peplib()

  testthat::expect_error(
    profile_plot(prev[, -5], "protA", peplib = lib),
    "prev_tbl is missing required columns"
  )
  testthat::expect_error(
    profile_plot(prev, "protC", peplib = lib),
    "proteins not found in the peptide library"
  )
  testthat::expect_error(
    profile_plot(prev, "protA", groups = "case", peplib = lib),
    "groups not found in prev_tbl"
  )
  testthat::expect_error(
    profile_plot(prev, "protA", peplib = lib[, -3]),
    "peplib is missing required columns"
  )
  testthat::expect_error(
    profile_plot(prev, "protA", peplib = lib,
                 delta_tbl = data.frame(feature = "protA")),
    "delta_tbl is missing required columns"
  )
  testthat::expect_error(
    profile_plot(prev, "protA", peplib = lib, alpha = 0),
    "alpha must be in"
  )

  # control gets two different prevalences for pep1
  conflict <- prev
  conflict$prop1[5] <- 0.99
  testthat::expect_error(
    profile_plot(conflict, "protA", peplib = lib),
    "Groups have different prevalences"
  )
})

testthat::test_that("profile_plot rejects IDs shared between libraries", {
  # combined (agilent_) and ICAM (icam_) peptides under one protein_seq_id
  lib <- .make_profile_peplib()
  lib$peptide_id <- c("agilent_1", "agilent_2", "icam_3",
                      "agilent_4", "agilent_5", "agilent_6")
  prev <- .make_profile_prev()
  prev$feature <- sub("pep", "agilent_", prev$feature)

  testthat::expect_error(
    profile_plot(prev, "protA", peplib = lib),
    "matches proteins of several peptide libraries"
  )
  testthat::expect_s3_class(profile_plot(prev, "protB", peplib = lib), "ggplot")
})
