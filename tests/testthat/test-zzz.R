# tests/testthat/test-zzz.R
# Coverage for R/zzz.R: .onLoad(), re-exports (load_example_data, get_example_path)

# ---------------------------------------------------------------------------
# .onLoad runs without error
# ---------------------------------------------------------------------------
testthat::test_that(".onLoad() runs without error", {
  testthat::expect_no_error(phiper:::.onLoad("", "phiper"))
})

testthat::test_that(".onLoad() is idempotent (safe to call twice)", {
  testthat::expect_no_error({
    phiper:::.onLoad("", "phiper")
    phiper:::.onLoad("", "phiper")
  })
})

# ---------------------------------------------------------------------------
# Re-exported functions are present and callable
# ---------------------------------------------------------------------------
testthat::test_that("load_example_data is exported and callable", {
  testthat::expect_true(is.function(load_example_data))
  obj <- load_example_data("small_mixture")
  testthat::expect_true(!is.null(obj))
})

testthat::test_that("get_example_path is exported and returns a non-empty path", {
  testthat::expect_true(is.function(get_example_path))
  p <- get_example_path("phip_mixture")
  testthat::expect_true(is.character(p))
  testthat::expect_true(nzchar(p))
})
