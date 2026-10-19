# tests/testthat/test-plot_utils.R
# Coverage for R/plot_utils.R

# ---------------------------------------------------------------------------
# phip_palette
# ---------------------------------------------------------------------------
testthat::test_that("phip_palette is a character vector of hex colours", {
  testthat::expect_type(phip_palette, "character")
  testthat::expect_true(length(phip_palette) >= 30L)
  testthat::expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", phip_palette)))
})

# ---------------------------------------------------------------------------
# scale_colour_phip / scale_color_phip / scale_fill_phip
# ---------------------------------------------------------------------------
testthat::test_that("scale_colour_phip() returns a ggplot2 Scale object", {
  s <- scale_colour_phip()
  testthat::expect_true(inherits(s, "Scale"))
})

testthat::test_that("scale_color_phip is identical to scale_colour_phip", {
  testthat::expect_identical(scale_color_phip, scale_colour_phip)
})

testthat::test_that("scale_fill_phip() returns a ggplot2 Scale object", {
  s <- scale_fill_phip()
  testthat::expect_true(inherits(s, "Scale"))
})

testthat::test_that("scales apply to a ggplot without error", {
  p <- ggplot2::ggplot(
    data.frame(x = 1:3, g = c("a", "b", "c")),
    ggplot2::aes(x, x, colour = g, fill = g)
  ) +
    ggplot2::geom_point() +
    scale_colour_phip() +
    scale_fill_phip()
  testthat::expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# theme_phip()
# ---------------------------------------------------------------------------
testthat::test_that("theme_phip() returns a ggplot2 theme", {
  th <- theme_phip()
  testthat::expect_s3_class(th, "theme")
})

testthat::test_that("theme_phip() accepts custom base_size and base_family", {
  th <- theme_phip(base_size = 10, base_family = "sans")
  testthat::expect_s3_class(th, "theme")
})

testthat::test_that("theme_phip() applies to a ggplot without error", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    theme_phip()
  testthat::expect_s3_class(p, "ggplot")
})


# ---------------------------------------------------------------------------
# .phip_palette_map()
# ---------------------------------------------------------------------------
testthat::test_that(".phip_palette_map() returns empty named character for 0 levels", {
  out <- .phip_palette_map(character(0))
  testthat::expect_length(out, 0L)
  testthat::expect_named(out, character(0))
})

testthat::test_that(".phip_palette_map() maps a single level", {
  out <- .phip_palette_map("A")
  testthat::expect_length(out, 1L)
  testthat::expect_named(out, "A")
  testthat::expect_identical(out[["A"]], phip_palette[[1L]])
})

testthat::test_that(".phip_palette_map() recycles palette for many levels", {
  lvls <- paste0("G", seq_len(40L))
  out  <- .phip_palette_map(lvls)
  testthat::expect_length(out, 40L)
  testthat::expect_named(out, lvls)
})

# ---------------------------------------------------------------------------
# .ph_overlap_run()
# ---------------------------------------------------------------------------
testthat::test_that(".ph_overlap_run() chains segments overlapping the region so far", {
  # 21-50 overlaps 1-30, 41-70 overlaps 21-50; input order does not matter
  testthat::expect_identical(
    .ph_overlap_run(c(41, 1, 21), c(70, 30, 50)),
    c(n = 3, start = 1, end = 70)
  )
  # 25-30 lies inside 1-60, so 50-80 still joins through the covered region
  testthat::expect_identical(
    .ph_overlap_run(c(1, 25, 50), c(60, 30, 80)),
    c(n = 3, start = 1, end = 80)
  )
})

testthat::test_that(".ph_overlap_run() keeps the longest run, the first on ties", {
  # 1-30 and 31-60 share no residue
  testthat::expect_identical(
    .ph_overlap_run(c(1, 31, 50, 70), c(30, 60, 80, 100)),
    c(n = 3, start = 31, end = 100)
  )
  testthat::expect_identical(
    .ph_overlap_run(c(1, 20, 100, 120), c(30, 50, 130, 150)),
    c(n = 2, start = 1, end = 50)
  )
  testthat::expect_identical(.ph_overlap_run(5, 40), c(n = 1, start = 5, end = 40))
})
