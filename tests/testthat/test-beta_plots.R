# tests/testthat/test-beta_plots.R
# Coverage for R/beta_plots.R

# ===========================================================================
# Fixtures
# ===========================================================================

.make_pcoa_res <- function(n = 20L, seed = 1L,
                           with_group = TRUE, with_time = FALSE,
                           with_var_explained = TRUE) {
  set.seed(seed)
  sc <- tibble::tibble(
    sample_id = paste0("s", seq_len(n)),
    PCoA1     = stats::rnorm(n),
    PCoA2     = stats::rnorm(n),
    PCoA3     = stats::rnorm(n)
  )
  if (with_group) sc$group <- rep(c("A", "B"), length.out = n)
  if (with_time)  sc$time  <- rep(c("T1", "T2"), length.out = n)
  ve <- if (with_var_explained) {
    tibble::tibble(`%PCoA1` = 35.2, `%PCoA2` = 18.7, check.names = FALSE)
  } else NULL
  list(sample_coords = sc, var_explained = ve,
       eigenvalues = c(5.0, 2.6, 1.4, 0.9))
}

.make_cap_res <- function(n = 20L, seed = 2L,
                          with_group = TRUE, with_time = FALSE) {
  set.seed(seed)
  sc <- tibble::tibble(
    sample_id = paste0("s", seq_len(n)),
    CAP1      = stats::rnorm(n),
    CAP2      = stats::rnorm(n)
  )
  if (with_group) sc$group <- rep(c("A", "B"), length.out = n)
  if (with_time)  sc$time  <- rep(c("T1", "T2"), length.out = n)
  list(sample_coords = sc, eigenvalues = c(2.5, 1.2, 0.8))
}

.make_disp_res <- function(n_per = 10L, seed = 3L, add_tests = TRUE) {
  set.seed(seed)
  dist_df <- tibble::tibble(
    sample_id = rep(paste0("s", seq_len(n_per)), 2L),
    distance  = c(stats::runif(n_per, 0.2, 0.5),
                  stats::runif(n_per, 0.3, 0.7)),
    level     = rep(c("A", "B"), each = n_per),
    scope     = "group",
    contrast  = "<global>"
  )
  tests_df <- if (add_tests) {
    tibble::tibble(
      scope    = "group",
      contrast = "<global>",
      p_value  = 0.03,
      n_perm   = 999L
    )
  } else tibble::tibble()
  structure(list(distances = dist_df, tests = tests_df),
            class = "beta_dispersion")
}

.make_tsne_df <- function(n = 20L, seed = 5L, dims = 2L) {
  set.seed(seed)
  df <- tibble::tibble(
    tSNE1 = stats::rnorm(n),
    tSNE2 = stats::rnorm(n),
    group = rep(c("A", "B"), length.out = n)
  )
  if (dims >= 3L) df$tSNE3 <- stats::rnorm(n)
  df
}

# ===========================================================================
# plot_pcoa()
# ===========================================================================

testthat::test_that("plot_pcoa: basic output is ggplot", {
  res <- .make_pcoa_res()
  p   <- plot_pcoa(res)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: no group or time — grey points", {
  res <- .make_pcoa_res(with_group = FALSE)
  p   <- plot_pcoa(res, group_col = NULL, time_col = NULL)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: group + time colours + shapes", {
  res <- .make_pcoa_res(with_group = TRUE, with_time = TRUE)
  p   <- plot_pcoa(res, group_col = "group", time_col = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: time only (no group)", {
  res <- .make_pcoa_res(with_group = FALSE, with_time = TRUE)
  p   <- plot_pcoa(res, group_col = NULL, time_col = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: non-default axes (1,3)", {
  res <- .make_pcoa_res()
  p   <- plot_pcoa(res, axes = c(1, 3))
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: show_centroids=FALSE suppresses centroids", {
  res <- .make_pcoa_res()
  p   <- plot_pcoa(res, show_centroids = FALSE)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: centroid_by=group with group column", {
  res <- .make_pcoa_res()
  p   <- plot_pcoa(res, centroid_by = "group")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: centroid_by=time with time column", {
  res <- .make_pcoa_res(with_time = TRUE)
  p   <- plot_pcoa(res, centroid_by = "time", time_col = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: centroid_by=group_time with both columns + connect_centroids=group", {
  res <- .make_pcoa_res(with_group = TRUE, with_time = TRUE)
  p   <- plot_pcoa(res, centroid_by = "group_time", time_col = "time",
                   connect_centroids = "group")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: connect_centroids=time", {
  res <- .make_pcoa_res(with_group = TRUE, with_time = TRUE)
  p   <- plot_pcoa(res, centroid_by = "group_time", time_col = "time",
                   connect_centroids = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: centroid_by=group disabled when no group col", {
  res <- .make_pcoa_res(with_group = FALSE)
  testthat::expect_warning(
    p <- plot_pcoa(res, centroid_by = "group", group_col = "group"),
    regexp = "centroids will be disabled"
  )
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: centroid_by=time disabled when no time col", {
  res <- .make_pcoa_res(with_time = FALSE)
  testthat::expect_warning(
    p <- plot_pcoa(res, centroid_by = "time", time_col = "time"),
    regexp = "centroids will be disabled"
  )
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: centroid_by=group_time disabled when missing a col", {
  res <- .make_pcoa_res(with_time = FALSE)
  testthat::expect_warning(
    p <- plot_pcoa(res, centroid_by = "group_time"),
    regexp = "centroids will be disabled"
  )
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: show_ellipses=TRUE with group ellipses", {
  res <- .make_pcoa_res(n = 30L)
  p   <- plot_pcoa(res, show_ellipses = TRUE, ellipse_by = "group")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: ellipse_by=time with time column", {
  res <- .make_pcoa_res(n = 30L, with_time = TRUE)
  p   <- plot_pcoa(res, show_ellipses = TRUE, ellipse_by = "time",
                   time_col = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: ellipse_by=group_time with both columns", {
  res <- .make_pcoa_res(n = 40L, with_group = TRUE, with_time = TRUE)
  p   <- plot_pcoa(res, show_ellipses = TRUE,
                   ellipse_by = "group_time", time_col = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_pcoa: no var_explained still produces axis labels", {
  res <- .make_pcoa_res(with_var_explained = FALSE)
  p   <- plot_pcoa(res)
  testthat::expect_s3_class(p, "ggplot")
})

# Error paths ----------------------------------------------------------------
testthat::test_that("plot_pcoa: NULL input errors", {
  testthat::expect_error(plot_pcoa(NULL))
})

testthat::test_that("plot_pcoa: missing sample_coords errors", {
  testthat::expect_error(plot_pcoa(list(eigenvalues = 1:5)))
})

testthat::test_that("plot_pcoa: bad axes length errors", {
  res <- .make_pcoa_res()
  testthat::expect_error(plot_pcoa(res, axes = 1L))
  testthat::expect_error(plot_pcoa(res, axes = c(1L, 2L, 3L)))
})

testthat::test_that("plot_pcoa: axis not in data errors", {
  res <- .make_pcoa_res()
  testthat::expect_error(plot_pcoa(res, axes = c(1L, 9L)))
})

# ===========================================================================
# plot_cap()
# ===========================================================================

testthat::test_that("plot_cap: basic output is ggplot", {
  res <- .make_cap_res()
  p   <- plot_cap(res)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: no group or time — grey points", {
  res <- .make_cap_res(with_group = FALSE)
  p   <- plot_cap(res, group_col = NULL, time_col = NULL)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: group + time", {
  res <- .make_cap_res(with_group = TRUE, with_time = TRUE)
  p   <- plot_cap(res, group_col = "group", time_col = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: time only", {
  res <- .make_cap_res(with_group = FALSE, with_time = TRUE)
  p   <- plot_cap(res, group_col = NULL, time_col = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: centroid_by=group", {
  res <- .make_cap_res()
  p   <- plot_cap(res, centroid_by = "group")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: centroid_by=time with time col", {
  res <- .make_cap_res(with_time = TRUE)
  p   <- plot_cap(res, centroid_by = "time", time_col = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: centroid_by=group_time + connect_centroids=group", {
  res <- .make_cap_res(n = 40L, with_group = TRUE, with_time = TRUE)
  p   <- plot_cap(res, centroid_by = "group_time", time_col = "time",
                  connect_centroids = "group")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: centroid_by=group_time + connect_centroids=time", {
  res <- .make_cap_res(n = 40L, with_group = TRUE, with_time = TRUE)
  p   <- plot_cap(res, centroid_by = "group_time", time_col = "time",
                  connect_centroids = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: ellipse_by=group", {
  res <- .make_cap_res(n = 30L)
  p   <- plot_cap(res, show_ellipses = TRUE, ellipse_by = "group")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: ellipse_by=time", {
  res <- .make_cap_res(n = 30L, with_time = TRUE)
  p   <- plot_cap(res, show_ellipses = TRUE, ellipse_by = "time",
                  time_col = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: ellipse_by=group_time", {
  res <- .make_cap_res(n = 40L, with_time = TRUE)
  p   <- plot_cap(res, show_ellipses = TRUE, ellipse_by = "group_time",
                  time_col = "time")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: ellipse_by=none disables ellipses", {
  res <- .make_cap_res()
  p   <- plot_cap(res, show_ellipses = TRUE, ellipse_by = "none")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_cap: axis labels annotated with eigenvalue pct", {
  res <- .make_cap_res()
  p   <- plot_cap(res)
  xlab <- p$labels$x
  testthat::expect_true(grepl("%", xlab))
})

testthat::test_that("plot_cap: no eigenvalues — raw axis labels", {
  res <- .make_cap_res()
  res$eigenvalues <- numeric(0L)
  p <- plot_cap(res)
  testthat::expect_s3_class(p, "ggplot")
})

# Error paths ----------------------------------------------------------------
testthat::test_that("plot_cap: non-list input errors", {
  testthat::expect_error(plot_cap("not a list"))
})

testthat::test_that("plot_cap: missing sample_coords errors", {
  testthat::expect_error(plot_cap(list(eigenvalues = 1:3)))
})

testthat::test_that("plot_cap: bad axes length errors", {
  res <- .make_cap_res()
  testthat::expect_error(plot_cap(res, axes = c(1L, 2L, 3L)))
})

testthat::test_that("plot_cap: axis not in data errors", {
  res <- .make_cap_res()
  testthat::expect_error(plot_cap(res, axes = c(1L, 9L)))
})

# ===========================================================================
# plot_scree()
# ===========================================================================

.scree_res <- list(eigenvalues = c(5.0, 3.2, 2.1, 1.6, 1.1, 0.8, 0.4))

testthat::test_that("plot_scree: bar type returns ggplot", {
  p <- plot_scree(.scree_res, type = "bar")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_scree: line type returns ggplot", {
  p <- plot_scree(.scree_res, type = "line")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_scree: n_axes limits the number of bars", {
  p <- plot_scree(.scree_res, n_axes = 3L)
  built <- ggplot2::ggplot_build(p)
  # geom_col data should have ≤ 3 rows
  n_rows <- nrow(built$data[[1]])
  testthat::expect_lte(n_rows, 3L)
})

testthat::test_that("plot_scree: n_axes=NULL uses min(10, length)", {
  p <- plot_scree(list(eigenvalues = 1:5))
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_scree: negative eigenvalues handled (pos only used for pct)", {
  res <- list(eigenvalues = c(5.0, 3.0, -0.1))
  p   <- plot_scree(res)
  testthat::expect_s3_class(p, "ggplot")
})

# Error paths ----------------------------------------------------------------
testthat::test_that("plot_scree: missing eigenvalues component errors", {
  testthat::expect_error(plot_scree(list(sample_coords = data.frame())))
})

testthat::test_that("plot_scree: empty eigenvalues errors", {
  testthat::expect_error(plot_scree(list(eigenvalues = numeric(0))))
})

# ===========================================================================
# plot_dispersion()
# ===========================================================================

testthat::test_that("plot_dispersion: default output is ggplot", {
  res <- .make_disp_res()
  p   <- plot_dispersion(res)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_dispersion: show_violin=FALSE", {
  res <- .make_disp_res()
  p   <- plot_dispersion(res, show_violin = FALSE)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_dispersion: show_box=FALSE", {
  res <- .make_disp_res()
  p   <- plot_dispersion(res, show_box = FALSE)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_dispersion: show_points=FALSE", {
  res <- .make_disp_res()
  p   <- plot_dispersion(res, show_points = FALSE)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_dispersion: add_pvalue=FALSE skips subtitle", {
  res <- .make_disp_res()
  p   <- plot_dispersion(res, add_pvalue = FALSE)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_dispersion: add_pvalue=TRUE with matching test row adds subtitle", {
  res <- .make_disp_res(add_tests = TRUE)
  p   <- plot_dispersion(res, add_pvalue = TRUE)
  testthat::expect_false(is.null(p$labels$subtitle))
})

testthat::test_that("plot_dispersion: add_pvalue=TRUE with empty tests is silent", {
  res <- .make_disp_res(add_tests = FALSE)
  p   <- plot_dispersion(res, add_pvalue = TRUE)
  testthat::expect_s3_class(p, "ggplot")
})

# Error paths ----------------------------------------------------------------
testthat::test_that("plot_dispersion: non-beta_dispersion input errors", {
  testthat::expect_error(plot_dispersion(list(distances = data.frame())))
})

testthat::test_that("plot_dispersion: missing required columns in distances errors", {
  res <- structure(
    list(distances = data.frame(sample_id = "s1", distance = 0.1),
         tests = tibble::tibble()),
    class = "beta_dispersion"
  )
  testthat::expect_error(plot_dispersion(res))
})

testthat::test_that("plot_dispersion: wrong scope/contrast combination errors", {
  res <- .make_disp_res()
  testthat::expect_error(plot_dispersion(res, scope = "time", contrast = "<global>"))
})

# ===========================================================================
# plot_tsne()
# ===========================================================================

testthat::test_that("plot_tsne: 2D basic returns ggplot", {
  df <- .make_tsne_df()
  p  <- plot_tsne(df, view = "2d")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_tsne: 2D with colour column", {
  df <- .make_tsne_df()
  p  <- plot_tsne(df, view = "2d", colour = "group")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_tsne: 2D with custom palette", {
  df <- .make_tsne_df()
  p  <- plot_tsne(df, view = "2d", colour = "group",
                  palette = c(A = "#FF0000", B = "#0000FF"))
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_tsne: 2D with colour inferred from meta_cols attribute", {
  df <- .make_tsne_df()
  attr(df, "meta_cols") <- "group"
  p  <- plot_tsne(df, view = "2d")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_tsne: 2D with non-existent colour col falls back to no colour", {
  df <- .make_tsne_df()
  p  <- plot_tsne(df, view = "2d", colour = "nonexistent_col")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_tsne: 3D without tSNE3 column errors", {
  df <- .make_tsne_df(dims = 2L)
  testthat::expect_error(plot_tsne(df, view = "3d"))
})

testthat::test_that("plot_tsne: 3D with tSNE3 and plotly returns plotly widget", {
  testthat::skip_if_not_installed("plotly")
  df <- .make_tsne_df(dims = 3L)
  p  <- plot_tsne(df, view = "3d", colour = "group")
  testthat::expect_true(inherits(p, "plotly") || inherits(p, "htmlwidget"))
})

# Error paths ----------------------------------------------------------------
testthat::test_that("plot_tsne: missing tSNE1/tSNE2 errors", {
  df <- data.frame(x = 1:5, y = 1:5)
  testthat::expect_error(plot_tsne(df))
})
