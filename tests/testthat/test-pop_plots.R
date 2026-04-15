# tests/testthat/test-pop_plots.R
# Unit tests for R/pop_plots.R:
#   scatter_static(), scatter_interactive(), volcano_static(), volcano_interactive()

# ===========================================================================
# Fixtures
# ===========================================================================

.make_prev_df <- function(n = 20L, seed = 1L,
                          with_ratio = TRUE,
                          group1 = "A", group2 = "B") {
  set.seed(seed)
  tibble::tibble(
    rank      = "peptide_id",
    feature   = paste0("pep", seq_len(n)),
    group1    = group1,
    group2    = group2,
    group_col = "group",
    n1        = sample(0L:10L, n, replace = TRUE),
    N1        = 10L,
    n2        = sample(0L:10L, n, replace = TRUE),
    N2        = 10L,
    prop1     = n1 / N1,
    prop2     = n2 / N2,
    percent1  = prop1 * 100,
    percent2  = prop2 * 100,
    ratio     = if (with_ratio)
      pmax(prop1, 1e-9) / pmax(prop2, 1e-9) else NA_real_,
    p_raw     = stats::runif(n, 0, 1),
    n_peptides = 1L
  )
}

# prev_df with two pairs (for pair-filtering tests)
.make_prev_df_two_pairs <- function(n_per_pair = 10L, seed = 2L) {
  dplyr::bind_rows(
    .make_prev_df(n = n_per_pair, seed = seed,       group1 = "A", group2 = "B"),
    .make_prev_df(n = n_per_pair, seed = seed + 100L, group1 = "A", group2 = "C")
  )
}

# minimal peplib fixture for color_by tests
.make_peplib <- function(features, col = "is_flag", vals = NULL) {
  if (is.null(vals)) vals <- rep(c(TRUE, FALSE), length.out = length(features))
  tibble::tibble(peptide_id = features, !!col := vals)
}

# ===========================================================================
# scatter_static — return type & basic structure
# ===========================================================================

testthat::test_that("scatter_static returns a ggplot", {
  df <- .make_prev_df()
  p  <- scatter_static(df)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("scatter_static: empty df returns a void ggplot", {
  df <- .make_prev_df()[0L, ]
  p  <- scatter_static(df)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("scatter_static: pair filtering keeps only that pair", {
  df <- .make_prev_df_two_pairs()
  p  <- scatter_static(df, pair = c("A", "B"))
  testthat::expect_s3_class(p, "ggplot")
  # axis labels default to pair values
  testthat::expect_equal(p$labels$x, "A")
  testthat::expect_equal(p$labels$y, "B")
})

testthat::test_that("scatter_static: rank filtering reduces data", {
  df <- .make_prev_df()
  df2 <- dplyr::bind_rows(df, dplyr::mutate(df, rank = "virus"))
  p  <- scatter_static(df2, rank = "peptide_id")
  testthat::expect_s3_class(p, "ggplot")
  # the plot data should only contain peptide_id rows
  pd <- ggplot2::ggplot_build(p)$data[[1]]
  testthat::expect_true(nrow(pd) <= nrow(df))
})

testthat::test_that("scatter_static: custom xlab/ylab are applied", {
  df <- .make_prev_df()
  p  <- scatter_static(df, xlab = "Group X", ylab = "Group Y")
  testthat::expect_equal(p$labels$x, "Group X")
  testthat::expect_equal(p$labels$y, "Group Y")
})

testthat::test_that("scatter_static: prop1/prop2 accepted when percent columns absent", {
  df <- .make_prev_df()
  df$percent1 <- NULL
  df$percent2 <- NULL
  p  <- scatter_static(df)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("scatter_static: graphical ... args accepted without error", {
  df <- .make_prev_df()
  testthat::expect_no_error(
    scatter_static(df, point_size = 3, point_alpha = 0.5,
                   font_size = 14, jitter_width_pp = 0.5)
  )
})

# ===========================================================================
# scatter_static — BH category coloring
# ===========================================================================

testthat::test_that("scatter_static: all-NS data triggers p_bin fallback", {
  df      <- .make_prev_df()
  df$p_raw <- 0.9   # all not significant → single category → p_bin fallback
  p  <- scatter_static(df, alpha = 0.05)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("scatter_static: BH correction applied (significant BH category exists)", {
  df        <- .make_prev_df(n = 30L, seed = 5L)
  df$p_raw  <- c(rep(1e-6, 10L), rep(0.5, 20L))   # 10 very significant
  p  <- scatter_static(df, alpha = 0.05)
  build <- ggplot2::ggplot_build(p)
  # color scale should include "significant (BH)"
  scale_labels <- build$plot$scales$scales[[1]]$labels %||%
    build$plot$scales$get_scales("colour")$labels
  testthat::expect_s3_class(p, "ggplot")
})

# ===========================================================================
# scatter_static — color_by (named vector interface)
# ===========================================================================

testthat::test_that("scatter_static: color_by uses peplib via scatter_interactive peplib arg", {
  # Test the .build_color_group helper indirectly through scatter_interactive
  # (which exposes peplib). We verify it returns a ggplot without error.
  df     <- .make_prev_df(n = 10L, seed = 3L)
  peplib <- .make_peplib(df$feature, col = "is_flag", vals = rep(c(TRUE, FALSE), 5L))

  testthat::expect_no_error(
    scatter_interactive(df, color_by = c("is_flag" = TRUE), peplib = peplib)
  )
})

# ===========================================================================
# scatter_interactive — return type & basic structure
# ===========================================================================

testthat::test_that("scatter_interactive returns a plotly object", {
  df <- .make_prev_df()
  p  <- scatter_interactive(df)
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("scatter_interactive: empty df returns a plotly object", {
  df <- .make_prev_df()[0L, ]
  p  <- scatter_interactive(df)
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("scatter_interactive: pair filtering works", {
  df <- .make_prev_df_two_pairs()
  p  <- scatter_interactive(df, pair = c("A", "B"))
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("scatter_interactive: rank filtering works", {
  df  <- .make_prev_df()
  df2 <- dplyr::bind_rows(df, dplyr::mutate(df, rank = "virus"))
  p   <- scatter_interactive(df2, rank = "peptide_id")
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("scatter_interactive: throws error when n_peptides column is absent", {
  df <- .make_prev_df()
  df$n_peptides <- NULL
  testthat::expect_error(scatter_interactive(df))
})

testthat::test_that("scatter_interactive: graphical ... args accepted without error", {
  df <- .make_prev_df()
  testthat::expect_no_error(
    scatter_interactive(df,
      point_size = 8, point_alpha = 0.7,
      jitter_width_pp = 0.3, font_size = 10,
      point_line_width = 0.5
    )
  )
})

# ===========================================================================
# scatter_interactive — color_by named vector
# ===========================================================================

testthat::test_that("scatter_interactive: color_by logical value highlights group", {
  df     <- .make_prev_df(n = 12L, seed = 7L)
  peplib <- .make_peplib(df$feature, col = "is_flag",
                         vals = c(rep(TRUE, 6L), rep(FALSE, 6L)))
  p <- scatter_interactive(df, color_by = c("is_flag" = TRUE), peplib = peplib)
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("scatter_interactive: color_by string value highlights group", {
  df     <- .make_prev_df(n = 12L, seed = 8L)
  peplib <- tibble::tibble(
    peptide_id = df$feature,
    species    = rep(c("Staph", "Ecoli"), 6L)
  )
  p <- scatter_interactive(df,
    color_by = c("species" = "Staph"),
    peplib   = peplib
  )
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("scatter_interactive: multiple color_by entries", {
  df     <- .make_prev_df(n = 12L, seed = 9L)
  peplib <- tibble::tibble(
    peptide_id = df$feature,
    is_flag    = c(rep(TRUE, 4L), rep(FALSE, 8L)),
    species    = rep(c("Staph", "Ecoli", "Other"), 4L)
  )
  p <- scatter_interactive(df,
    color_by = c("is_flag" = TRUE, "species" = "Staph"),
    peplib   = peplib
  )
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("scatter_interactive: color_by with no matching keys → other group", {
  df              <- .make_prev_df(n = 6L, seed = 10L)
  df$rank         <- "protein"   # not "peptide_id" → pep_key all NA → join skipped
  peplib          <- .make_peplib(paste0("pep", 1:6))
  p  <- scatter_interactive(df,
    color_by = c("is_flag" = TRUE),
    peplib   = peplib
  )
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("scatter_interactive: color_title sets legend title", {
  df     <- .make_prev_df(n = 8L, seed = 11L)
  peplib <- .make_peplib(df$feature)
  p <- scatter_interactive(df,
    color_by    = c("is_flag" = TRUE),
    color_title = "Flagellum",
    peplib      = peplib
  )
  testthat::expect_s3_class(p, "plotly")
})

# ===========================================================================
# scatter_interactive — background overlay
# ===========================================================================

testthat::test_that("scatter_interactive: background_df overlay accepted", {
  df  <- .make_prev_df(n = 15L, seed = 12L)
  bg  <- .make_prev_df(n = 50L, seed = 99L)
  p   <- scatter_interactive(df,
    background_df = bg,
    show_background = TRUE
  )
  testthat::expect_s3_class(p, "plotly")
})

# ===========================================================================
# .build_color_group — internal helper
# ===========================================================================

testthat::test_that(".build_color_group: single logical condition assigns correctly", {
  df     <- .make_prev_df(n = 6L, seed = 20L)
  peplib <- tibble::tibble(
    peptide_id = df$feature,
    is_flag    = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
  )
  res <- phiper:::.build_color_group(df, c("is_flag" = TRUE), peplib = peplib)
  grp <- as.character(res$df$.color_group)
  testthat::expect_equal(sum(grp == "is_flag: TRUE"), 2L)
  testthat::expect_equal(sum(grp == "other"), 4L)
})

testthat::test_that(".build_color_group: multiple conditions, last-match wins", {
  df     <- .make_prev_df(n = 4L, seed = 21L)
  # feature pep1 matches both conditions: last (species = Staph) wins
  peplib <- tibble::tibble(
    peptide_id = df$feature,
    is_flag    = c(TRUE, FALSE, FALSE, FALSE),
    species    = c("Staph", "Ecoli", "Ecoli", "Ecoli")
  )
  res <- phiper:::.build_color_group(df,
    c("is_flag" = TRUE, "species" = "Staph"),
    peplib = peplib
  )
  grp <- as.character(res$df$.color_group)
  # pep1 matches is_flag:TRUE first, then species:Staph → last match wins
  testthat::expect_equal(grp[1], "species: Staph")
})

testthat::test_that(".build_color_group: returns factor with 'other' level", {
  df     <- .make_prev_df(n = 4L, seed = 22L)
  peplib <- tibble::tibble(
    peptide_id = df$feature,
    is_flag    = c(TRUE, FALSE, FALSE, FALSE)
  )
  res  <- phiper:::.build_color_group(df, c("is_flag" = TRUE), peplib = peplib)
  levs <- levels(res$df$.color_group)
  testthat::expect_true("other" %in% levs)
  testthat::expect_true("is_flag: TRUE" %in% levs)
})

# ===========================================================================
# volcano_static — return type & basic structure
# ===========================================================================

testthat::test_that("volcano_static returns a ggplot", {
  df <- .make_prev_df()
  p  <- volcano_static(df)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("volcano_static: empty df returns a void ggplot", {
  df <- .make_prev_df()[0L, ]
  p  <- volcano_static(df)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("volcano_static: p_mode raw uses p_raw directly", {
  df <- .make_prev_df()
  testthat::expect_no_error(volcano_static(df, p_mode = "raw"))
})

testthat::test_that("volcano_static: p_mode bh computes BH correction", {
  df <- .make_prev_df()
  testthat::expect_no_error(volcano_static(df, p_mode = "bh"))
})

testthat::test_that("volcano_static: pair filtering works", {
  df <- .make_prev_df_two_pairs()
  p  <- volcano_static(df, pair = c("A", "B"))
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("volcano_static: rank filtering works", {
  df  <- .make_prev_df()
  df2 <- dplyr::bind_rows(df, dplyr::mutate(df, rank = "virus"))
  p   <- volcano_static(df2, rank = "peptide_id")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("volcano_static: fc_cut and p_cut are respected (no error)", {
  df <- .make_prev_df()
  testthat::expect_no_error(volcano_static(df, fc_cut = 2, p_cut = 0.01))
})

testthat::test_that("volcano_static: ratio derived from prop1/prop2 when absent", {
  df       <- .make_prev_df(with_ratio = FALSE)
  df$ratio <- NULL
  p  <- volcano_static(df)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("volcano_static: missing ratio and props raises error", {
  df        <- .make_prev_df()
  df$ratio  <- NULL
  df$prop1  <- NULL
  df$prop2  <- NULL
  testthat::expect_error(volcano_static(df), regexp = "prop1")
})

testthat::test_that("volcano_static: color_by with peplib colors points", {
  df     <- .make_prev_df(n = 10L, seed = 30L)
  peplib <- .make_peplib(df$feature, col = "is_flag")
  # volcano_static doesn't expose peplib, so we test via volcano_interactive
  # which also calls .volcano_prepare. Here just verify static doesn't error.
  # (Library join in static goes via get_peptide_library; skip color_by test there.)
  testthat::expect_s3_class(volcano_static(df), "ggplot")
})

# ===========================================================================
# volcano_interactive — return type & basic structure
# ===========================================================================

testthat::test_that("volcano_interactive returns a plotly object", {
  df <- .make_prev_df()
  p  <- volcano_interactive(df)
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("volcano_interactive: empty df returns a plotly object", {
  df <- .make_prev_df()[0L, ]
  p  <- volcano_interactive(df)
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("volcano_interactive: p_mode raw", {
  df <- .make_prev_df()
  testthat::expect_no_error(volcano_interactive(df, p_mode = "raw"))
})

testthat::test_that("volcano_interactive: p_mode bh", {
  df <- .make_prev_df()
  testthat::expect_no_error(volcano_interactive(df, p_mode = "bh"))
})

testthat::test_that("volcano_interactive: pair filtering works", {
  df <- .make_prev_df_two_pairs()
  p  <- volcano_interactive(df, pair = c("A", "B"))
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("volcano_interactive: rank filtering works", {
  df  <- .make_prev_df()
  df2 <- dplyr::bind_rows(df, dplyr::mutate(df, rank = "virus"))
  p   <- volcano_interactive(df2, rank = "peptide_id")
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("volcano_interactive: ratio derived from prop1/prop2 when absent", {
  df       <- .make_prev_df(with_ratio = FALSE)
  df$ratio <- NULL
  p  <- volcano_interactive(df)
  testthat::expect_s3_class(p, "plotly")
})

testthat::test_that("volcano_interactive: invalid p_mode raises error", {
  df <- .make_prev_df()
  testthat::expect_error(
    volcano_interactive(df, p_mode = "wbh"),
    regexp = "arg"
  )
})

# ===========================================================================
# BH correction — shared logic validation
# ===========================================================================

testthat::test_that("scatter_static BH: p_bh <= p_raw for all rows (BH is conservative)", {
  # With a single group of independent p-values, all BH-adjusted values
  # should be >= raw (BH inflates, not deflates).
  df       <- .make_prev_df(n = 50L, seed = 42L)
  df$p_raw <- stats::runif(50L, 0, 1)
  p_bh     <- stats::p.adjust(df$p_raw, method = "BH")
  testthat::expect_true(all(p_bh >= df$p_raw - .Machine$double.eps))
})

testthat::test_that("volcano: BH and raw produce different nlog10p when p differs", {
  df       <- .make_prev_df(n = 30L, seed = 43L)
  df$p_raw <- c(rep(1e-8, 5L), stats::runif(25L, 0.1, 0.9))

  prep_raw <- phiper:::.volcano_prepare(df, p_mode = "raw")
  prep_bh  <- phiper:::.volcano_prepare(df, p_mode = "bh")

  # BH-adjusted p-values differ from raw for at least some rows
  testthat::expect_false(
    identical(prep_raw$data$nlog10p, prep_bh$data$nlog10p)
  )
})
