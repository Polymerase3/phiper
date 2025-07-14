skip_if_not_installed("dplyr")
skip_if_not_installed("cli")
skip_if_not_installed("knitr")
skip_if_not_installed("DBI")

# ---------------------------------------------------------------------------
# helper: tiny example data
# ---------------------------------------------------------------------------
counts_tbl <- tibble::tibble(
  peptide      = c("pep1", "pep2", "pep1", "pep2"),
  sample_id    = c("S1",   "S1",   "S2",   "S2"),
  timepoint    = c("T1",   "T1",   "T1",   "T1"),
  fold_change  = c(1.2,    0.8,    1.5,    0.7),
  present      = c(10L,    0L,     8L,     1L)
)

contrasts_df <- tibble::tibble(
  group   = "S1",
  control = "S2"
)

# ---------------------------------------------------------------------------
# constructor + meta flags
# ---------------------------------------------------------------------------
test_that("new_phip_data sets meta flags correctly", {
  pd <- new_phip_data(counts_tbl, contrasts_df, backend = "memory")

  expect_s3_class(pd, "phip_data")
  expect_true(pd$meta$longitudinal)
  expect_true(pd$meta$fold_change)
  expect_equal(get_backend(pd), "memory")
})

# ---------------------------------------------------------------------------
# print method (just make sure it runs and contains certain strings)
# ---------------------------------------------------------------------------
test_that("print.phip_data shows backend and previews", {
  pd <- new_phip_data(counts_tbl, contrasts_df, backend = "memory")

  out <- capture.output(print(pd))
  expect_true(any(grepl("backend: memory", out)))
  expect_true(any(grepl("counts \\(first 5 rows\\):", out)))
  expect_true(any(grepl("contrasts:", out)))
})

# ---------------------------------------------------------------------------
# plain accessors and .check_pd guard
# ---------------------------------------------------------------------------
test_that("accessors work and .check_pd errors on wrong class", {
  pd <- new_phip_data(counts_tbl, contrasts_df, backend = "memory")

  expect_equal(get_counts(pd), counts_tbl)
  expect_equal(get_comparisons(pd), contrasts_df)
  expect_equal(get_meta(pd)$longitudinal, TRUE)

  expect_error(get_counts(list(a = 1)), "`x` must be a <phip_data> object")
})

# ---------------------------------------------------------------------------
# dplyr verb wrappers
# ---------------------------------------------------------------------------
test_that("dplyr wrappers modify data_long lazily", {
  pd <- new_phip_data(counts_tbl, contrasts_df, backend = "memory")

  pd2 <- dplyr::filter(pd, peptide == "pep1")

  expect_s3_class(pd2, "phip_data")
  expect_equal(
    dplyr::collect(pd2$data_long)$peptide,
    c("pep1", "pep1")
  )

  pd3 <- dplyr::select(pd, peptide)
  expect_equal(colnames(pd3$data_long), "peptide")

  pd4 <- dplyr::mutate(pd, new_val = present * 2)
  expect_true("new_val" %in% colnames(pd4$data_long))

  pd5 <- dplyr::arrange(pd, desc(present))
  expect_equal(nrow(pd5$data_long), 4)

  n <- dplyr::summarise(pd, n = dplyr::n())
  expect_equal(n$n, 4)

  collected <- dplyr::collect(pd)
  expect_s3_class(collected, "data.frame")
})

# ---------------------------------------------------------------------------
# disconnect helper (mock DBI)
# ---------------------------------------------------------------------------
test_that("disconnect.phip_data closes duckdb connection if present", {
  skip_if_not_installed("duckdb")

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  pd  <- new_phip_data(counts_tbl, contrasts_df,
                       backend = "duckdb",
                       meta = list(con = con))

  expect_true(DBI::dbIsValid(con))

  disconnect.phip_data(pd)   # should close

  expect_false(DBI::dbIsValid(con))
})
