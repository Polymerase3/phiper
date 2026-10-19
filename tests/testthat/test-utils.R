skip_if_not_installed("chk")

# -------------------------------------------------------------------------
# .ph_check_cond ------------------------------------------------------------
# -------------------------------------------------------------------------
test_that(".ph_check_cond emits errors and warnings correctly", {
  expect_error(
    .ph_check_cond(TRUE,
                   "boom",
                   error = TRUE
    ),
    label = "boom"
  )

  expect_warning(
    .ph_check_cond(TRUE,
                   "soft boom",
                   error = FALSE
    ),
    label = "soft_boom"
  )

  # When condition is FALSE nothing happens
  expect_silent(.ph_check_cond(FALSE, "no-op"))
})

# ------------------------------------------------------------------------------
# .ph_add_quotes ---------------------------------------------------------------
# ------------------------------------------------------------------------------
test_that(".ph_add_quotes behaves", {
  # quotes variations
  expect_equal(.ph_add_quotes("x", quotes = FALSE), "x")
  expect_equal(.ph_add_quotes("x", quotes = TRUE), '"x"')
  expect_equal(.ph_add_quotes("x", quotes = 1), "'x'")
  expect_equal(.ph_add_quotes("x", quotes = 2), '"x"')
})

# ------------------------------------------------------------------------------
# Operators: %||%
# ------------------------------------------------------------------------------
test_that("%||% returns lhs when non-NULL, rhs otherwise", {
  expect_equal("a"  %||% "b", "a")
  expect_equal(NULL %||% "b", "b")
  expect_equal(0L   %||% 99L, 0L)
  expect_null(NULL %||% NULL)
})

# ------------------------------------------------------------------------------
# Logging helpers
# ------------------------------------------------------------------------------
test_that(".ph_opt() reads option with fallback", {
  withr::with_options(list(phiperio.log.verbose = FALSE), {
    expect_false(.ph_opt("verbose", TRUE))
  })
  expect_true(.ph_opt("verbose", TRUE))  # no option set -> default
})

test_that(".ph_now() returns a non-empty character string", {
  out <- .ph_now()
  expect_type(out, "character")
  expect_true(nzchar(out))
})

test_that(".ph_base_prefix() produces a formatted prefix string", {
  out <- .ph_base_prefix("INFO")
  expect_type(out, "character")
  expect_true(grepl("INFO", out))
})

test_that(".ph_wrap() returns character lines containing the text", {
  lines <- .ph_wrap("hello world", "[PREFIX] ")
  expect_type(lines, "character")
  expect_true(any(grepl("hello", lines)))
})

test_that(".ph_compose_lines() includes headline, step, and bullets", {
  out <- .ph_compose_lines("INFO", "headline", "a step", c("bullet one", "bullet two"))
  expect_type(out, "character")
  expect_true(any(grepl("headline", out)))
  expect_true(any(grepl("a step", out)))
  expect_true(any(grepl("bullet one", out)))
})

test_that(".ph_compose_lines() handles NULL step and NULL bullets", {
  out <- .ph_compose_lines("WARN", "just a headline")
  expect_type(out, "character")
  expect_true(any(grepl("headline", out)))
})

test_that(".ph_compose_lines() skips blank/NA bullets", {
  out <- .ph_compose_lines("INFO", "head", bullets = c("good", NA, "", "also good"))
  text <- paste(out, collapse = " ")
  expect_true(grepl("good", text))
  expect_true(grepl("also good", text))
})

test_that(".ph_log_info() returns invisible character when verbose=TRUE", {
  out <- withr::with_options(list(phiperio.log.verbose = TRUE), {
    suppressMessages(.ph_log_info("test headline", verbose = TRUE))
  })
  expect_type(out, "character")
})

test_that(".ph_log_info() returns invisible empty character when verbose=FALSE", {
  out <- .ph_log_info("test", verbose = FALSE)
  expect_length(out, 0L)
})

test_that(".ph_log_ok() returns invisible character when verbose=TRUE", {
  out <- suppressMessages(.ph_log_ok("ok headline", verbose = TRUE))
  expect_type(out, "character")
})

test_that(".ph_log_ok() returns invisible empty character when verbose=FALSE", {
  out <- .ph_log_ok("ok", verbose = FALSE)
  expect_length(out, 0L)
})

test_that(".ph_warn() emits a warning", {
  expect_warning(.ph_warn("watch out"))
})

test_that(".ph_warn() with step and bullets still emits a warning", {
  expect_warning(.ph_warn("warn head", step = "step", bullets = "a bullet"))
})

test_that(".ph_abort() throws an error", {
  expect_error(.ph_abort("something went wrong"))
})

test_that(".ph_abort() with step and bullets throws an error", {
  expect_error(.ph_abort("error head", step = "step", bullets = c("b1", "b2")))
})

test_that(".ph_with_timing() returns the expression result", {
  out <- .ph_with_timing("timing test", expr = 42L, verbose = FALSE)
  expect_equal(out, 42L)
})

test_that(".ph_with_timing() propagates errors from expr", {
  expect_error(
    .ph_with_timing("will fail", expr = stop("inner error"), verbose = FALSE),
    "inner error"
  )
})

# ------------------------------------------------------------------------------
# .ph_add_quotes(): error on invalid quotes argument
# ------------------------------------------------------------------------------
test_that(".ph_add_quotes() with quotes=0 returns original string", {
  expect_equal(.ph_add_quotes("hello", quotes = 0L), "hello")
})

test_that(".ph_add_quotes() with a custom quote char works", {
  expect_equal(.ph_add_quotes("x", quotes = "|"), "|x|")
})

test_that(".ph_add_quotes() with invalid quotes argument throws error", {
  expect_error(.ph_add_quotes("x", quotes = 5L))
})

# -------------------------------------------------------------------------
# .ph_library_for_peptides ------------------------------------------------
# -------------------------------------------------------------------------
test_that(".ph_library_for_peptides fetches the libraries the peptides belong to", {
  lib <- tibble::tibble(peptide_id = c("agilent_1", "icam_2"))
  requested <- NULL
  local_mocked_bindings(
    get_peptide_library = function(library = "combined", force_refresh = FALSE) {
      requested <<- library
      lib
    }
  )

  out <- .ph_library_for_peptides(c("agilent_1", "icam_2", "agilent_3", NA))

  expect_identical(requested, c("combined", "icam"))
  expect_identical(out, lib)
})

test_that(".ph_library_for_peptides returns NULL when no library matches", {
  local_mocked_bindings(
    get_peptide_library = function(...) stop("no library should be fetched")
  )

  expect_null(.ph_library_for_peptides(c("pep1", "pep2")))
})
