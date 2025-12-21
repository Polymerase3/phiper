# tests/testthat/test-compute_distance.R

# helper that builds a small real phip_data from the package example data
.get_ps_small_for_distance <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) return(cache)

    phip_path <- phip_example_path()

    ps <- phip_convert(
      data_long_path    = phip_path,
      backend           = "duckdb",
      peptide_library   = TRUE,
      subject_id        = "subject_id",
      peptide_id        = "peptide_id",
      sample_id         = "sample_id",
      exist             = "exist",
      timepoint         = "timepoint_factor",
      fold_change       = "fold_change",
      materialise_table = TRUE,
      auto_expand       = TRUE,
      n_cores           = 2
    )

    cols <- dplyr::tbl_vars(ps$data_long)
    tp_col <- if ("time" %in% cols) {
      "time"
    } else if ("timepoint_factor" %in% cols) {
      "timepoint_factor"
    } else if ("timepoint" %in% cols) {
      "timepoint"
    } else {
      testthat::skip("example data_long has no recognized time column")
    }

    # small subset: 5 peptides at time T1
    keep_pep <- c("16627", "5243", "24799", "16196", "18003")

    ps_filt <- ps %>%
      dplyr::filter(
        peptide_id %in% keep_pep,
        !!rlang::sym(tp_col) == "T1"
      ) %>%
      dplyr::collect()

    cache <<- ps_filt
    cache
  }
})

.pick_existing_value_col <- function(ps, preferred) {
  dat <- if("phip_data" %in% class(ps)) ps$data_long else ps
  cols <- dplyr::tbl_vars(dat)
  if (preferred %in% cols) return(preferred)
}

testthat::test_that("compute_distance works on example phip_data and attaches
                    abundances", {

  ps <- .get_ps_small_for_distance()
  value_col <- .pick_existing_value_col(ps, "fold_change")

  d <- suppressWarnings(compute_distance(
    ps,
    value_col = value_col,
    method_normalization = "relative",
    distance = "euclidean",
    n_threads = 1L
  ))

  testthat::expect_s3_class(d, "dist")
  testthat::expect_true(is.matrix(attr(d, "abundances")))
  testthat::expect_true(is.numeric(as.numeric(d)))
  testthat::expect_true(all(as.numeric(d) >= 0))

  abund <- attr(d, "abundances")
  testthat::expect_equal(attr(d, "Size"), nrow(abund))
  testthat::expect_true(length(rownames(abund)) == nrow(abund))
  testthat::expect_true(all(!is.na(rownames(abund))))

  # relative normalization: nonzero rows should sum to ~1; all-zero rows sum to 0
  rs <- rowSums(abund)
  is_zero_row <- rowSums(abund != 0) == 0
  testthat::expect_true(all(abs(rs[!is_zero_row] - 1) < 1e-10))
  testthat::expect_true(all(rs[is_zero_row] == 0))
})

testthat::test_that("compute_distance auto-detects value_col on example data", {

  ps <- .get_ps_small_for_distance()

  d <- suppressWarnings(compute_distance(
    ps,
    value_col = NULL,
    method_normalization = "none",
    distance = "kulczynski",
    n_threads = 1L
  ))

  testthat::expect_s3_class(d, "dist")
  testthat::expect_true(is.matrix(attr(d, "abundances")))
})

testthat::test_that("compute_distance auto normalization chooses none for binary
                    exist data", {

  ps <- .get_ps_small_for_distance()

  cols <- dplyr::tbl_vars(ps)

  d <- suppressWarnings(compute_distance(
    ps,
    value_col = "exist",
    method_normalization = "auto",
    distance = "manhattan",
    n_threads = 1L
  ))
  a <- attr(d, "abundances")
  vals <- as.vector(a)

  # should remain binary (0/1) if auto picked "none"
  testthat::expect_true(all(vals %in% c(0, 1)))
})

testthat::test_that("compute_distance errors on duplicated (sample_id,
                    peptide_id) pairs", {

  ps <- .get_ps_small_for_distance()
  value_col <- .pick_existing_value_col(ps, "fold_change")

  ps_dup <- ps
  ps_dup <- dplyr::bind_rows(ps, ps[1, , drop = FALSE])

  testthat::expect_error(
    suppressWarnings(compute_distance(
      ps_dup,
      value_col = value_col,
      method_normalization = "none",
      distance = "manhattan",
      n_threads = 1L
    )),
    regexp = "(?i)duplicat|unique|pivot"
  )
})

testthat::test_that("compute_distance validates inputs (chk) on example data", {

  ps <- .get_ps_small_for_distance()
  value_col <- .pick_existing_value_col(ps, "fold_change")

  ps_bad <- ps
  class(ps_bad) <- "not_phip_data"
  testthat::expect_error(
    compute_distance(ps_bad, value_col = value_col),
    regexp = "(?i)phip_data|class"
  )

  ps_nolong <- structure(list(data_long = NULL), class = "phip_data")
  testthat::expect_error(
    compute_distance(ps_nolong, value_col = value_col),
    regexp = "(?i)data_long|null|missing"
  )

  testthat::expect_error(
    compute_distance(ps, value_col = value_col, n_threads = 0L),
    regexp = "(?i)thread|count|gt|positive|chk"
  )

  testthat::expect_error(
    compute_distance(ps, value_col = "definitely_not_a_column"),
    regexp = "(?i)not found|missing|column"
  )
})

testthat::test_that("compute_distance reproducibility: repeated runs match
                    across parameter combinations", {

  ps <- .get_ps_small_for_distance()

  value_fc <- "fold_change"
  value_ex <- "exist"

  combos <- list(
    list(value_col = value_fc, norm = "none",
         dist = "euclidean"),
    list(value_col = value_fc, norm = "relative",
         dist = "manhattan"),
    list(value_col = value_fc, norm = "hellinger",
         dist = "euclidean"),
    list(value_col = value_fc, norm = "log",
         dist = "manhattan")
  )

  combos <- combos[!vapply(combos, function(x) is.null(x$value_col),
                           logical(1))]

  for (cmb in combos) {
    d1 <- withr::with_seed(123, suppressWarnings(compute_distance(
      ps,
      value_col = cmb$value_col,
      method_normalization = cmb$norm,
      distance = cmb$dist,
      n_threads = 1L
    )))
    d2 <- withr::with_seed(123, suppressWarnings(compute_distance(
      ps,
      value_col = cmb$value_col,
      method_normalization = cmb$norm,
      distance = cmb$dist,
      n_threads = 1L
    )))

    testthat::expect_equal(as.numeric(d1), as.numeric(d2), tolerance = 1e-12)
    testthat::expect_equal(attr(d1, "abundances"), attr(d2, "abundances"),
                           tolerance = 1e-12)
    testthat::expect_equal(attr(d1, "Size"), attr(d2, "Size"))
  }
})

testthat::test_that("compute_distance bray agrees with vegan and manual", {

  ps <- .get_ps_small_for_distance()
  value_col <- .pick_existing_value_col(ps, "fold_change")

  d <- suppressWarnings(compute_distance(
    ps,
    value_col = value_col,
    method_normalization = "relative",
    distance = "bray",
    n_threads = 1L
  ))
  testthat::expect_s3_class(d, "dist")

  x <- attr(d, "abundances")


  d_ref <- vegan::vegdist(x, method = "bray")
  testthat::expect_equal(as.numeric(d), as.numeric(d_ref), tolerance = 1e-10)

  # manual bray-curtis in dist upper-triangle order
  n <- nrow(x)
  v <- numeric(n * (n - 1L) / 2L)
  k <- 1L
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      num <- sum(abs(x[i, ] - x[j, ]))
      den <- sum(x[i, ]) + sum(x[j, ])
      v[k] <- if (den == 0) 0 else num / den
      k <- k + 1L
    }
  }
  testthat::expect_equal(as.numeric(d), v, tolerance = 1e-10)
})
