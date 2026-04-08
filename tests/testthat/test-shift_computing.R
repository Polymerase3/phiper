# tests/testthat/test-shift_computing.R
# Exercises uncovered branches in src/shift_computing.cpp via compute_delta():
#   stat_mode: "diff", "score", "mcnemar", "srlr_paired"
#   weight_mode: "se_invvar", "equal"
#   strat_bins: non-zero (stratification path)
#   winsor_z: finite clamping
#   paired design with multiple paired stat modes

# ---------------------------------------------------------------------------
# Guard: restore build_bitset_unpaired to the real C++ wrapper if a prior
# test (test-delta_compute.R) leaked its local_mocked_bindings mock.
# ---------------------------------------------------------------------------
local({
  ns   <- asNamespace("phiper")
  real <- function(hits_by_peptide, N) {
    .Call(`_phiper_build_bitset_unpaired`, hits_by_peptide, N)
  }
  environment(real) <- ns
  tryCatch({
    unlockBinding("build_bitset_unpaired", ns)
    assign("build_bitset_unpaired", real, envir = ns)
    lockBinding("build_bitset_unpaired", ns)
  }, error = function(e) NULL)
})

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
.make_unpaired_df <- function(n_per_grp = 12L, n_pep = 6L, seed = 1L) {
  set.seed(seed)
  # compute_delta requires subject_id; for unpaired each sample is its own subject
  sample_ids <- c(paste0("sA", seq_len(n_per_grp)),
                  paste0("sB", seq_len(n_per_grp)))
  df <- tidyr::crossing(
    sample_id  = sample_ids,
    peptide_id = paste0("pep", seq_len(n_pep))
  ) |>
    dplyr::mutate(
      subject_id = .data$sample_id,
      group      = dplyr::if_else(grepl("^sA", .data$sample_id), "A", "B"),
      exist      = sample(0:1, dplyr::n(), replace = TRUE)
    )
  df
}

.make_paired_df <- function(n_subjects = 8L, n_pep = 5L, seed = 2L) {
  set.seed(seed)
  subs <- paste0("sub", seq_len(n_subjects))
  tidyr::crossing(
    subject_id = subs,
    peptide_id = paste0("pep", seq_len(n_pep)),
    group      = c("T1", "T2")
  ) |>
    dplyr::mutate(
      sample_id = paste0(.data$subject_id, "_", .data$group),
      exist     = sample(0:1, dplyr::n(), replace = TRUE)
    )
}

.mock_peplib <- function(n_pep = 6L) {
  data.frame(
    peptide_id = paste0("pep", seq_len(n_pep)),
    species    = paste0("sp_", seq_len(n_pep)),
    stringsAsFactors = FALSE
  )
}

.run_unpaired <- function(stat_mode, weight_mode = "equal",
                          strat_bins = 0, winsor_z = Inf,
                          n_pep = 6L, seed = 1L) {
  df  <- .make_unpaired_df(n_pep = n_pep, seed = seed)
  pep <- .mock_peplib(n_pep)
  compute_delta(
    x               = df,
    exist_col       = "exist",
    rank_cols       = "species",
    group_cols      = "group",
    peptide_library = pep,
    B_permutations  = 100L,
    stat_mode       = stat_mode,
    weight_mode     = weight_mode,
    strat_bins      = strat_bins,
    winsor_z        = winsor_z,
    log             = FALSE
  )
}

.run_paired <- function(stat_mode, weight_mode = "equal",
                        n_pep = 5L, seed = 2L) {
  df  <- .make_paired_df(n_pep = n_pep, seed = seed)
  pep <- data.frame(
    peptide_id = paste0("pep", seq_len(n_pep)),
    species    = rep("sp1", n_pep),
    stringsAsFactors = FALSE
  )
  compute_delta(
    x               = df,
    exist_col       = "exist",
    rank_cols       = "species",
    group_cols      = "group",
    paired_by       = "subject_id",
    peptide_library = pep,
    B_permutations  = 100L,
    stat_mode       = stat_mode,
    weight_mode     = weight_mode,
    strat_bins      = 0,
    winsor_z        = Inf,
    log             = FALSE
  )
}

# ---------------------------------------------------------------------------
# stat_mode = "diff"  (else branch in combine_T_internal)
# ---------------------------------------------------------------------------
testthat::test_that("stat_mode=diff runs without error and returns valid T_obs", {
  res <- .run_unpaired("diff")
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_false(is.na(res$T_obs[[1L]]))
})

# ---------------------------------------------------------------------------
# stat_mode = "score"
# ---------------------------------------------------------------------------
testthat::test_that("stat_mode=score runs without error", {
  res <- .run_unpaired("score")
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_false(is.na(res$T_obs[[1L]]))
})

# ---------------------------------------------------------------------------
# stat_mode = "srlr"  (default — verify basic result)
# ---------------------------------------------------------------------------
testthat::test_that("stat_mode=srlr runs without error", {
  res <- .run_unpaired("srlr")
  testthat::expect_s3_class(res, "tbl_df")
})

# ---------------------------------------------------------------------------
# stat_mode = "mcnemar" (paired-only stat)
# ---------------------------------------------------------------------------
testthat::test_that("stat_mode=mcnemar requires paired_by — errors without it", {
  testthat::expect_error(
    .run_unpaired("mcnemar"),
    regexp = "paired_by"
  )
})

testthat::test_that("stat_mode=mcnemar runs in paired design", {
  res <- .run_paired("mcnemar")
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_identical(res$design[[1L]], "paired")
  testthat::expect_false(is.na(res$T_obs[[1L]]))
})

# ---------------------------------------------------------------------------
# stat_mode = "srlr_paired" (paired-only stat)
# ---------------------------------------------------------------------------
testthat::test_that("stat_mode=srlr_paired requires paired_by — errors without it", {
  testthat::expect_error(
    .run_unpaired("srlr_paired"),
    regexp = "paired_by"
  )
})

testthat::test_that("stat_mode=srlr_paired runs in paired design", {
  res <- .run_paired("srlr_paired")
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_identical(res$design[[1L]], "paired")
  testthat::expect_false(is.na(res$T_obs[[1L]]))
})

# ---------------------------------------------------------------------------
# weight_mode = "se_invvar" with various stat modes
# ---------------------------------------------------------------------------
testthat::test_that("weight_mode=se_invvar + stat_mode=srlr", {
  res <- .run_unpaired("srlr", weight_mode = "se_invvar")
  testthat::expect_s3_class(res, "tbl_df")
})

testthat::test_that("weight_mode=se_invvar + stat_mode=score", {
  res <- .run_unpaired("score", weight_mode = "se_invvar")
  testthat::expect_s3_class(res, "tbl_df")
})

testthat::test_that("weight_mode=se_invvar + stat_mode=asin", {
  res <- .run_unpaired("asin", weight_mode = "se_invvar")
  testthat::expect_s3_class(res, "tbl_df")
})

testthat::test_that("weight_mode=se_invvar + stat_mode=diff (fallback se_invvar branch)", {
  res <- .run_unpaired("diff", weight_mode = "se_invvar")
  testthat::expect_s3_class(res, "tbl_df")
})

# ---------------------------------------------------------------------------
# weight_mode = "equal"
# ---------------------------------------------------------------------------
testthat::test_that("weight_mode=equal + stat_mode=asin", {
  res <- .run_unpaired("asin", weight_mode = "equal")
  testthat::expect_s3_class(res, "tbl_df")
})

# ---------------------------------------------------------------------------
# strat_bins: non-zero (stratification path in combine_T_internal)
# ---------------------------------------------------------------------------
testthat::test_that("strat_bins=c(0.1,0.5) exercises the stratification path", {
  res <- .run_unpaired("srlr", strat_bins = c(0.1, 0.5))
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_false(is.na(res$T_obs[[1L]]))
})

testthat::test_that("strat_bins with all-boundary values falls back gracefully", {
  # Values outside (0,1) are filtered out -> effectively no stratification
  res <- .run_unpaired("srlr", strat_bins = c(0.0, 1.0))
  testthat::expect_s3_class(res, "tbl_df")
})

# ---------------------------------------------------------------------------
# winsor_z: finite clamping
# ---------------------------------------------------------------------------
testthat::test_that("finite winsor_z clamps extreme z-scores", {
  res <- .run_unpaired("asin", winsor_z = 1.5)
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_false(is.na(res$T_obs[[1L]]))
})

testthat::test_that("very small winsor_z clamps to near zero", {
  res <- .run_unpaired("asin", winsor_z = 0.001)
  testthat::expect_s3_class(res, "tbl_df")
  # T_obs should be very small after aggressive clamping
  testthat::expect_lte(abs(res$T_obs[[1L]]), 0.01 + .Machine$double.eps)
})

# ---------------------------------------------------------------------------
# Paired design with n_eff_sqrt weight + strat_bins
# ---------------------------------------------------------------------------
testthat::test_that("paired + mcnemar + se_invvar runs", {
  res <- .run_paired("mcnemar", weight_mode = "se_invvar")
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_identical(res$design[[1L]], "paired")
})

testthat::test_that("paired + srlr_paired + n_eff_sqrt runs", {
  res <- .run_paired("srlr_paired", weight_mode = "n_eff_sqrt")
  testthat::expect_s3_class(res, "tbl_df")
})
