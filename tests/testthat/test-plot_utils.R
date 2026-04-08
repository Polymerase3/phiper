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
# .hex2rgb() / .rgb2hex()
# ---------------------------------------------------------------------------
testthat::test_that(".hex2rgb() converts a known white", {
  r <- .hex2rgb("#FFFFFF")
  testthat::expect_equal(r, c(1, 1, 1), tolerance = 1e-6)
})

testthat::test_that(".hex2rgb() converts a known black", {
  r <- .hex2rgb("#000000")
  testthat::expect_equal(r, c(0, 0, 0), tolerance = 1e-6)
})

testthat::test_that(".hex2rgb() returns NAs for NULL/NA/empty", {
  testthat::expect_true(all(is.na(.hex2rgb(NULL))))
  testthat::expect_true(all(is.na(.hex2rgb(NA_character_))))
  testthat::expect_true(all(is.na(.hex2rgb(""))))
})

testthat::test_that(".rgb2hex() roundtrips a colour through hex2rgb", {
  hex_in <- "#1F78B4"
  hex_out <- .rgb2hex(.hex2rgb(hex_in))
  # round-trip should be close (allow 1 unit rounding in each channel)
  r1 <- grDevices::col2rgb(hex_in)
  r2 <- grDevices::col2rgb(hex_out)
  testthat::expect_true(all(abs(r1 - r2) <= 1L))
})

testthat::test_that(".rgb2hex() returns fallback grey for NA input", {
  testthat::expect_equal(.rgb2hex(c(NA_real_, 0, 0)), "#B3B3B3")
  testthat::expect_equal(.rgb2hex(c(1, 2)), "#B3B3B3")  # too short
})

testthat::test_that(".rgb2hex() clamps values outside [0,1]", {
  out <- .rgb2hex(c(-0.5, 0.5, 1.5))
  testthat::expect_true(grepl("^#[0-9A-Fa-f]{6}$", out))
})

# ---------------------------------------------------------------------------
# .mix_cols() / .tint()
# ---------------------------------------------------------------------------
testthat::test_that(".mix_cols() midpoint between black and white is grey", {
  out <- .mix_cols("#000000", "#FFFFFF", 0.5)
  rgb <- grDevices::col2rgb(out) / 255
  testthat::expect_true(all(abs(rgb - 0.5) < 0.01))
})

testthat::test_that(".mix_cols() with t=0 returns the first colour", {
  out <- .mix_cols("#FF0000", "#FFFFFF", 0)
  r <- grDevices::col2rgb(out)
  ref <- grDevices::col2rgb("#FF0000")
  testthat::expect_true(all(abs(r - ref) <= 1L))
})

testthat::test_that(".mix_cols() returns fallback grey for NA inputs", {
  testthat::expect_equal(.mix_cols(NA, "#FFFFFF", 0.5), "#B3B3B3")
  testthat::expect_equal(.mix_cols("#FF0000", NA, 0.5), "#B3B3B3")
  testthat::expect_equal(.mix_cols("#FF0000", "#FFFFFF", NA), "#B3B3B3")
})

testthat::test_that(".tint() lightens a colour (higher R/G/B than original)", {
  orig <- grDevices::col2rgb("#1F78B4") / 255
  out  <- grDevices::col2rgb(.tint("#1F78B4", 0.5)) / 255
  testthat::expect_true(all(out >= orig - 1e-6))
})

# ---------------------------------------------------------------------------
# .blend_hex()
# ---------------------------------------------------------------------------
testthat::test_that(".blend_hex() w=0 returns the first colour", {
  out <- .blend_hex("#FF0000", "#FFFFFF", 0)
  r   <- grDevices::col2rgb(out)
  ref <- grDevices::col2rgb("#FF0000")
  testthat::expect_true(all(abs(r - ref) <= 1L))
})

testthat::test_that(".blend_hex() w=1 returns the second colour", {
  out <- .blend_hex("#FF0000", "#FFFFFF", 1)
  r   <- grDevices::col2rgb(out)
  ref <- grDevices::col2rgb("#FFFFFF")
  testthat::expect_true(all(abs(r - ref) <= 1L))
})

# ---------------------------------------------------------------------------
# .make_shades()
# ---------------------------------------------------------------------------
testthat::test_that(".make_shades() returns base_hex unchanged for n=1", {
  out <- .make_shades("#1F78B4", 1L)
  testthat::expect_identical(out, "#1F78B4")
})

testthat::test_that(".make_shades() returns n colours for n>1", {
  out <- .make_shades("#1F78B4", 4L)
  testthat::expect_length(out, 4L)
  testthat::expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", out)))
})

# ---------------------------------------------------------------------------
# .build_shaded_map()
# ---------------------------------------------------------------------------
testthat::test_that(".build_shaded_map() returns fill and outline with correct names", {
  out <- .build_shaded_map(c("A", "B"), c("T1", "T2"))
  testthat::expect_named(out, c("fill", "outline"))
  nms <- c("A | T1", "A | T2", "B | T1", "B | T2")
  testthat::expect_true(all(nms %in% names(out$fill)))
  testthat::expect_true(all(nms %in% names(out$outline)))
})

testthat::test_that(".build_shaded_map() accepts custom base_colors", {
  bc <- c(A = "#FF0000", B = "#00FF00")
  out <- .build_shaded_map(c("A", "B"), "T1", base_colors = bc)
  testthat::expect_named(out, c("fill", "outline"))
})

# ---------------------------------------------------------------------------
# .pick_axes()
# ---------------------------------------------------------------------------
testthat::test_that(".pick_axes() finds PCoA columns", {
  df <- data.frame(PCoA1 = 1, PCoA2 = 2)
  out <- .pick_axes(df, c(1, 2))
  testthat::expect_equal(out, c("PCoA1", "PCoA2"))
})

testthat::test_that(".pick_axes() finds CAP columns", {
  df <- data.frame(CAP1 = 1, CAP2 = 2, CAP3 = 3)
  out <- .pick_axes(df, c(1, 2))
  testthat::expect_equal(out, c("CAP1", "CAP2"))
})

testthat::test_that(".pick_axes() finds PC columns", {
  df <- data.frame(PC1 = 1, PC2 = 2)
  out <- .pick_axes(df, c(1, 2))
  testthat::expect_equal(out, c("PC1", "PC2"))
})

testthat::test_that(".pick_axes() aborts when no axis columns found", {
  df <- data.frame(x = 1, y = 2)
  testthat::expect_error(.pick_axes(df, c(1, 2)))
})

# ---------------------------------------------------------------------------
# .axis_labels_with_pct()
# ---------------------------------------------------------------------------
testthat::test_that(".axis_labels_with_pct() returns raw names when var_ex is NULL", {
  out <- .axis_labels_with_pct(NULL, c("PCoA1", "PCoA2"))
  testthat::expect_equal(unname(out), c("PCoA1", "PCoA2"))
})

testthat::test_that(".axis_labels_with_pct() returns raw names when var_ex is empty", {
  out <- .axis_labels_with_pct(data.frame(), c("PCoA1", "PCoA2"))
  testthat::expect_equal(unname(out), c("PCoA1", "PCoA2"))
})

testthat::test_that(".axis_labels_with_pct() formats pct labels correctly", {
  ve <- data.frame(`%PCoA1` = 35.2, `%PCoA2` = 18.7, check.names = FALSE)
  out <- .axis_labels_with_pct(ve, c("PCoA1", "PCoA2"))
  testthat::expect_true(grepl("35.2%", out[["PCoA1"]]))
  testthat::expect_true(grepl("18.7%", out[["PCoA2"]]))
})

testthat::test_that(".axis_labels_with_pct() falls back when pct cols are missing", {
  ve <- data.frame(`%PCoA3` = 10.0, check.names = FALSE)
  out <- .axis_labels_with_pct(ve, c("PCoA1", "PCoA2"))
  testthat::expect_equal(unname(out), c("PCoA1", "PCoA2"))
})

# ---------------------------------------------------------------------------
# .shaded_colors()
# ---------------------------------------------------------------------------
.make_shade_df <- function(n = 10L, seed = 1L, add_time_cont = FALSE,
                           add_time_fac = FALSE) {
  set.seed(seed)
  df <- data.frame(
    group = rep(c("A", "B"), length.out = n),
    stringsAsFactors = FALSE
  )
  if (add_time_cont) df$time_cont <- seq_len(n)
  if (add_time_fac)  df$time_fac  <- rep(c("T1", "T2"), length.out = n)
  df
}

testthat::test_that(".shaded_colors() returns n hex colours for shade='none'", {
  df  <- .make_shade_df(10L)
  out <- .shaded_colors(df, group_col = "group", shade = "none")
  testthat::expect_length(out, 10L)
  testthat::expect_true(all(grepl("^#", out)))
})

testthat::test_that(".shaded_colors() shades by continuous time", {
  df  <- .make_shade_df(10L, add_time_cont = TRUE)
  out <- .shaded_colors(df, group_col = "group", time_cont = "time_cont", shade = "auto")
  testthat::expect_length(out, 10L)
})

testthat::test_that(".shaded_colors() shades by factor time", {
  df  <- .make_shade_df(10L, add_time_fac = TRUE)
  out <- .shaded_colors(df, group_col = "group", time_fac = "time_fac", shade = "auto")
  testthat::expect_length(out, 10L)
})

testthat::test_that(".shaded_colors() shade='group' ignores time columns", {
  df  <- .make_shade_df(6L, add_time_cont = TRUE)
  out <- .shaded_colors(df, group_col = "group", shade = "group")
  testthat::expect_length(out, 6L)
})

# ---------------------------------------------------------------------------
# .make_point_fills()
# ---------------------------------------------------------------------------
testthat::test_that(".make_point_fills() returns base colours in 'group' mode", {
  df  <- .make_shade_df(8L)
  bm  <- c(A = "#1F78B4", B = "#33A02C")
  out <- .make_point_fills(df, bm, group_col = "group", mode = "group")
  testthat::expect_length(out, 8L)
})

testthat::test_that(".make_point_fills() shades with time_cont in 'time' mode", {
  df  <- .make_shade_df(8L, add_time_cont = TRUE)
  bm  <- c(A = "#1F78B4", B = "#33A02C")
  out <- .make_point_fills(df, bm, group_col = "group",
                           time_cont = "time_cont", mode = "time")
  testthat::expect_length(out, 8L)
})

testthat::test_that(".make_point_fills() shades with time_fac in 'time' mode", {
  df  <- .make_shade_df(8L, add_time_fac = TRUE)
  bm  <- c(A = "#1F78B4", B = "#33A02C")
  out <- .make_point_fills(df, bm, group_col = "group",
                           time_fac = "time_fac", mode = "time")
  testthat::expect_length(out, 8L)
})

testthat::test_that(".make_point_fills() auto mode falls back to group when no time", {
  df  <- .make_shade_df(6L)
  bm  <- c(A = "#1F78B4", B = "#33A02C")
  out <- .make_point_fills(df, bm, group_col = "group", mode = "auto")
  testthat::expect_length(out, 6L)
})

# ---------------------------------------------------------------------------
# .first_subview_name()
# ---------------------------------------------------------------------------
testthat::test_that(".first_subview_name() returns the first name", {
  x <- list(alpha = 1, beta = 2)
  testthat::expect_identical(.first_subview_name(x), "alpha")
})

testthat::test_that(".first_subview_name() aborts on unnamed list", {
  testthat::expect_error(.first_subview_name(list(1, 2)))
})
