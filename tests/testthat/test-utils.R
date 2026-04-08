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

# -------------------------------------------------------------------------
# .ph_check_extension ------------------------------------------------------
# -------------------------------------------------------------------------
test_that(".ph_check_extension passes correct ext and fails otherwise", {
  # good
  expect_silent(
    .ph_check_extension("file.TSV", "file", c("tsv", "csv"))
  )

  # bad
  expect_error(
    .ph_check_extension("file.bad", "file", c("tsv", "csv")),
    class = "chk_error"
  )
})

# -------------------------------------------------------------------------
# .ph_check_null_default ---------------------------------------------------
# -------------------------------------------------------------------------
test_that(".ph_check_null_default returns original or default", {
  expect_equal(
    .ph_check_null_default(5, "x", "m", 10),
    5
  )

  withr::with_options(list(warn = 1), {
    expect_warning(
      out <- .ph_check_null_default(NULL, "x", "m", 10),
      label = "method"
    )
    expect_equal(out, 10)
  })
})

# -------------------------------------------------------------------------
# .ph_check_path -----------------------------------------------------------
# -------------------------------------------------------------------------
test_that(".ph_check_path validates path + ext", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  writeLines("test", tmp)

  # passes
  expect_silent(
    .ph_check_path(tmp, "arg", extension = "csv")
  )

  # fails on extension
  expect_error(
    .ph_check_path(tmp, "arg", extension = "tsv"),
    class = "chk_error"
  )

  # fails on non-existent path
  expect_error(
    .ph_check_path(file.path(tempdir(), "nope.csv"), "arg", extension = "csv"),
    class = "chk_error"
  )
})

test_that(".ph_check_path works for folders", {
  tmp <- withr::local_tempdir()

  # Default is is_dir=FALSE
  expect_error(
    .ph_check_path(tmp, "arg")
  )

  # Explicit *file* check
  expect_error(
    .ph_check_path(tmp, "arg", is_dir=FALSE)
  )

  # Correct folder check
  expect_silent(
    .ph_check_path(tmp, "arg", is_dir=TRUE)
  )

  # Both folder and extension should fail
  expect_error(
    .ph_check_path(tmp, "arg", c("txt", "pqt"), is_dir=TRUE)
  )

  # Non existing directory should fail
  expect_error(
    .ph_check_path("this/path/should/definetly/not/exist", "arg", is_dir=TRUE)
  )
})

# ------------------------------------------------------------------------------
# .ph_word_list / .ph_add_quotes -----------------------------------------------
# ------------------------------------------------------------------------------
test_that(".ph_word_list and .ph_add_quotes behave", {
  expect_equal(
    as.vector(.ph_word_list(c("a", "b", "c"), and_or = "and", quotes = TRUE)),
    '"a", "b", and "c"'
  )

  expect_equal(
    attr(.ph_word_list("a"), "plural"),
    FALSE
  )

  # quotes variations
  expect_equal(.ph_add_quotes("x", quotes = FALSE), "x")
  expect_equal(.ph_add_quotes("x", quotes = TRUE), '"x"')
  expect_equal(.ph_add_quotes("x", quotes = 1), "'x'")
  expect_equal(.ph_add_quotes("x", quotes = 2), '"x"')
})

# ------------------------------------------------------------------------------
# Operators: %nin% and %||%
# ------------------------------------------------------------------------------
test_that("%nin% is the complement of %in%", {
  expect_true(5L %nin% c(1L, 2L, 3L))
  expect_false(2L %nin% c(1L, 2L, 3L))
  expect_equal(c(1, 2, 3) %nin% c(2, 3), c(TRUE, FALSE, FALSE))
  expect_equal(character(0) %nin% "a", logical(0))
})

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
# .ph_check_pd()
# ------------------------------------------------------------------------------
test_that(".ph_check_pd() aborts on a non-phip_data object", {
  expect_error(.ph_check_pd(list(a = 1)))
  expect_error(.ph_check_pd(data.frame()))
})

test_that(".ph_check_pd() passes silently for a phip_data object", {
  pd <- load_example_data("small_mixture")
  if (inherits(pd, "phip_data")) {
    expect_silent(.ph_check_pd(pd))
  } else {
    skip("small_mixture not a phip_data object in this environment")
  }
})

# ------------------------------------------------------------------------------
# .ph_word_list(): additional branches
# ------------------------------------------------------------------------------
test_that(".ph_word_list() handles two-element lists", {
  out <- .ph_word_list(c("a", "b"), and_or = "and")
  expect_equal(as.vector(out), "a and b")
  expect_true(attr(out, "plural"))
})

test_that(".ph_word_list() with and_or=FALSE joins with commas", {
  out <- .ph_word_list(c("x", "y", "z"), and_or = FALSE)
  expect_equal(as.vector(out), "x, y, z")
})

test_that(".ph_word_list() with is_are=TRUE appends is/are", {
  single <- .ph_word_list("x", is_are = TRUE)
  expect_true(grepl("is", single))
  multi <- .ph_word_list(c("x", "y"), is_are = TRUE)
  expect_true(grepl("are", multi))
})

test_that(".ph_word_list() strips NA and empty from word list", {
  out <- .ph_word_list(c("a", NA, "", "b"), and_or = "and")
  expect_equal(as.vector(out), "a and b")
})

test_that(".ph_word_list() returns empty string for NULL input", {
  out <- .ph_word_list(NULL)
  expect_equal(as.vector(out), "")
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

# ------------------------------------------------------------------------------
# .ph_resolve_paths()
# ------------------------------------------------------------------------------
test_that(".ph_resolve_paths() fails without any file source", {
  expect_error(.ph_resolve_paths())
})

test_that(".ph_resolve_paths() fails when input_file provided without hit_file", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  writeLines("test", tmp)
  expect_error(
    .ph_resolve_paths(input_file = tmp, samples_file = tmp)
  )
})

test_that(".ph_resolve_paths() succeeds with exist_file", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  writeLines("test", tmp)
  cfg <- .ph_resolve_paths(exist_file = tmp, samples_file = tmp)
  expect_equal(cfg$exist_file, normalizePath(tmp, winslash = "/"))
  expect_null(cfg$fold_change_file)
})

test_that(".ph_resolve_paths() warns for deprecated output_dir", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  writeLines("test", tmp)
  expect_warning(.ph_resolve_paths(exist_file = tmp, samples_file = tmp, output_dir = "/some/path"))
})

test_that(".ph_resolve_paths() fails when data_long_path given with exist_file", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  writeLines("test", tmp)
  tmpdir <- withr::local_tempdir()
  expect_error(
    .ph_resolve_paths(data_long_path = tmp, exist_file = tmp)
  )
})
