#!/usr/bin/env Rscript
# compute_delta2() validation driver. Not run on CRAN or by R CMD check.
#
#   Rscript run_validation.R [--sections 1,2,3,4,5,tests] [--quick]
#                            [--out results] [--mixture path]
#
# Environment: PHIPER_SOURCE (load phiper from a source checkout with
# pkgload instead of the installed package), VAL_CORES (parallel workers).
# See README.md for the sections, gates and expected run time.

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[i + 1L]
}
quick <- "--quick" %in% args
sections <- strsplit(opt("sections", "1,2,3,4,5,tests"), ",")[[1]]

here <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else getwd()
})
Sys.setenv(VAL_DIR = here)
out_dir <- opt("out", file.path(here, "results"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (f in list.files(file.path(here, "R"), full.names = TRUE)) source(f)
val_load_phiper()

# full run: 1000 datasets x B = 500 (6.1, 6.3, 6.4), 500 (6.2), 200-300 (6.5)
n <- if (quick) {
  list(cal = 50, pow = 30, unb = 50, rho = 20, rho_pow = 30, B = 200L)
} else {
  list(cal = 1000, pow = 500, unb = 1000, rho = 200, rho_pow = 300, B = 500L)
}
mixture <- opt("mixture", Sys.getenv("VAL_MIXTURE", ""))
if (!nzchar(mixture)) {
  src <- Sys.getenv("PHIPER_SOURCE", "")
  cand <- if (nzchar(src)) {
    file.path(src, "inst", "extdata", "phip_mixture.parquet")
  } else {
    system.file("extdata", "phip_mixture.parquet", package = "phiper")
  }
  parquet_ok <- requireNamespace("arrow", quietly = TRUE) ||
    requireNamespace("nanoparquet", quietly = TRUE)
  mixture <- if (file.exists(cand) && parquet_ok) cand else ""
}

timed <- function(label, expr) {
  t0 <- Sys.time()
  res <- force(expr)
  message(sprintf("%s: %.1f min", label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res
}

results <- list()
if ("1" %in% sections) {
  results$s1 <- timed("6.1", val_s1_calibration(n$cal, n$B))
  ggplot2::ggsave(file.path(out_dir, "calibration.png"),
    val_fig_calibration(results$s1$rates), width = 9, height = 4, dpi = 150, bg = val_pal$surface)
}
if ("2" %in% sections) {
  results$s2 <- timed("6.2", val_s2_pooled_power(n$pow, n$B))
  ggplot2::ggsave(file.path(out_dir, "pooled_power.png"),
    val_fig_power(results$s2$power), width = 9, height = 4, dpi = 150, bg = val_pal$surface)
}
if ("3" %in% sections) results$s3 <- timed("6.3", val_s3_unbalanced(n$unb, n$B))
if ("4" %in% sections) results$s4 <- timed("6.4", val_s4_turnover(n$unb, n$B))
if ("5" %in% sections) {
  results$s5 <- timed("6.5", val_s5_rho(n$rho, n$rho_pow, n$B,
    data_path = if (nzchar(mixture)) mixture
  ))
  ggplot2::ggsave(file.path(out_dir, "rho_recovery.png"),
    val_fig_rho(results$s5$recovery), width = 5, height = 5, dpi = 150, bg = val_pal$surface)
}

# 6.6 / 6.7: the brute-force sampler checks and the equivalence suite
if ("tests" %in% sections) {
  src <- Sys.getenv("PHIPER_SOURCE", "")
  test_dir <- if (nzchar(src)) file.path(src, "tests", "testthat") else ""
  files <- c("test-delta2_resample.R", "test-delta2_score.R", "test-delta2_compute.R")
  if (nzchar(test_dir) && all(file.exists(file.path(test_dir, files)))) {
    Sys.setenv(NOT_CRAN = "true")
    tests <- dplyr::bind_rows(lapply(files, function(f) {
      # in phiper's namespace: the tests call internal helpers
      r <- as.data.frame(testthat::test_file(file.path(test_dir, f),
        reporter = "silent", stop_on_failure = FALSE,
        package = "phiper", load_package = "none"
      ))
      tibble::tibble(
        file = f, tests = nrow(r), expectations = sum(r$nb),
        failed = sum(r$failed), errors = sum(r$error)
      )
    }))
    results$tests <- list(table = tests, gates = tibble::tibble(
      section = "6.6/6.7", gate = sprintf("%s passes", tests$file),
      value = tests$failed + tests$errors, lo = 0, hi = 0,
      pass = tests$failed + tests$errors == 0
    ))
  } else {
    message("tests: set PHIPER_SOURCE to a source checkout to run 6.6 / 6.7")
  }
}

for (nm in names(results)) {
  if (!is.null(results[[nm]]$gates)) {
    utils::write.csv(results[[nm]]$gates, file.path(out_dir, sprintf("gates_%s.csv", nm)),
      row.names = FALSE
    )
  }
}
saveRDS(results, file.path(out_dir, "results.rds"))
settings <- sprintf(
  "%s run; B = %d; cores = %d; phiper %s.",
  if (quick) "Quick" else "Full", n$B, val_cores(),
  as.character(utils::packageVersion("phiper"))
)
val_write_report(results, out_dir, settings)
gates <- dplyr::bind_rows(lapply(results, `[[`, "gates"))
message(sprintf("gates: %d of %d pass; report: %s", sum(gates$pass), nrow(gates),
  file.path(out_dir, "report.md")))
