# tests/testthat/test-compute_alpha_diversity.R

.pick_peplib_rank <- function(ps) {
  lib_cols <- tryCatch(colnames(ps$peptide_library),
                       error = function(e) character(0))
  rank_cols <- setdiff(lib_cols, "peptide_id")
  rank_cols[[1L]]
}

.get_ps_small_for_alpha <- function() {
  cache_env <- getOption("phiper.test.alpha_cache_env")
  if (is.null(cache_env) || !is.environment(cache_env)) {
    cache_env <- new.env(parent = emptyenv())
    withr::local_options(
      list(phiper.test.alpha_cache_env = cache_env),
      .local_envir = parent.frame()
    )
  }
  if (!is.null(cache_env$val)) return(cache_env$val)

  ps <- phip_load_example_data()

  dat_cols <- dplyr::tbl_vars(ps$data_long)
  tp_col <- "timepoint"

  cand_pep <- c("16627", "5243", "24799", "16196", "18003")
  .data <- rlang::.data

  pep_avail <- ps$data_long |>
    dplyr::distinct(.data$peptide_id) |>
    dplyr::collect()

  keep_pep <- intersect(cand_pep, pep_avail$peptide_id)

  ps_small <- ps |>
    dplyr::filter(
      .data$peptide_id %in% keep_pep,
      .data[[tp_col]] == "T1"
    )

  cache_env$val <- list(ps = ps_small, tp_col = tp_col, keep_pep = keep_pep)
  cache_env$val
}

testthat::test_that("compute_alpha_diversity works on example phip_data with interaction", {
  info <- .get_ps_small_for_alpha()
  ps <- info$ps
  tp_col <- info$tp_col

  out <- compute_alpha_diversity(
    ps,
    group_cols = c("group", tp_col),
    ranks = "peptide_id",
    group_interaction = TRUE
  )

  testthat::expect_s3_class(out, "phip_alpha_diversity")
  testthat::expect_true(all(c("group", tp_col) %in% names(out)))

  combo_nm <- paste(c("group", tp_col), collapse = " * ")
  testthat::expect_true(combo_nm %in% names(out))

  req_cols <- c("rank", "sample_id", "richness", "shannon_diversity",
                "simpson_diversity")
  for (nm in names(out)) {
    df <- out[[nm]]
    testthat::expect_s3_class(df, "tbl_df")
    testthat::expect_true(all(req_cols %in% names(df)))
    if (nm == combo_nm) {
      testthat::expect_true("phip_interaction" %in% names(df))
    } else {
      testthat::expect_true(nm %in% names(df))
    }
    testthat::expect_true(all(df$richness >= 0, na.rm = TRUE))
    testthat::expect_true(all(df$shannon_diversity >= 0, na.rm = TRUE))
    testthat::expect_true(all(df$simpson_diversity >= 0, na.rm = TRUE))
    testthat::expect_true(all(df$rank == "peptide_id"))
  }

  testthat::expect_identical(attr(out, "group_cols"), c("group", tp_col))
  testthat::expect_identical(attr(out, "ranks"), "peptide_id")
  testthat::expect_true(isTRUE(attr(out, "interaction")))
  testthat::expect_identical(attr(out, "interaction_sep"), " * ")
})

testthat::test_that("compute_alpha_diversity handles full_cross pruning and
                    peplib ranks", {
  info <- .get_ps_small_for_alpha()
  ps <- info$ps

  ps_full <- ps
  ps_full$meta$full_cross <- TRUE
  ps_full$meta$exist_prop <- 0.25

  out <- compute_alpha_diversity(
    ps_full,
    group_cols = "group",
    ranks = "peptide_id"
  )

  testthat::expect_s3_class(out, "phip_alpha_diversity")

  rank_col <- .pick_peplib_rank(ps_full)

  out2 <- compute_alpha_diversity(
    ps_full,
    group_cols = "group",
    ranks = c("peptide_id", rank_col)
  )

  testthat::expect_true(all(c("peptide_id", rank_col) %in% out2$group$rank))
})

testthat::test_that("compute_alpha_diversity supports interaction_only and
                    fold-change presence", {
  info <- .get_ps_small_for_alpha()
  ps <- info$ps
  tp_col <- info$tp_col

  out <- compute_alpha_diversity(
    ps,
    group_cols = c("group", tp_col),
    ranks = "peptide_id",
    fc_threshold = 1.5,
    shannon_log = "log10",
    group_interaction = TRUE,
    interaction_only = TRUE
  )

  combo_nm <- paste(c("group", tp_col), collapse = " * ")
  testthat::expect_identical(names(out), combo_nm)
  testthat::expect_true(isTRUE(attr(out, "interaction_only")))
  testthat::expect_identical(attr(out, "fc_threshold"), 1.5)
  testthat::expect_identical(attr(out, "shannon_log"), "log10")

  res <- out[[combo_nm]]
  testthat::expect_true(all(c("phip_interaction", "richness") %in% names(res)))
  testthat::expect_true(all(res$rank == "peptide_id"))
})

testthat::test_that("compute_alpha_diversity handles data.frame input without
                    grouping", {
  info <- .get_ps_small_for_alpha()
  df_small <- dplyr::collect(info$ps$data_long)

  out <- compute_alpha_diversity(
    df_small,
    group_cols = NULL,
    ranks = "peptide_id"
  )

  testthat::expect_s3_class(out, "phip_alpha_diversity")
  testthat::expect_identical(names(out), "all_samples")

  res <- out$all_samples
  testthat::expect_true("group" %in% names(res))
  testthat::expect_true(all(res$group == "All samples"))
})

testthat::test_that("compute_alpha_diversity handles empty presence sets", {
  df_empty <- tibble::tibble(
    sample_id = c("s1", "s2"),
    peptide_id = c("p1", "p2"),
    exist = 0L,
    group = "A"
  )

  out <- compute_alpha_diversity(
    df_empty,
    group_cols = NULL,
    ranks = "peptide_id"
  )

  res <- out$all_samples
  testthat::expect_true(all(res$richness == 0L))
  testthat::expect_true(all(res$shannon_diversity == 0))
  testthat::expect_true(all(res$simpson_diversity == 0))
})

testthat::test_that("compute_alpha_diversity carries extra columns on
                    data.frame input", {
  info <- .get_ps_small_for_alpha()
  tp_col <- info$tp_col
  df_small <- dplyr::collect(info$ps$data_long)

  out <- compute_alpha_diversity(
    df_small,
    group_cols = "group",
    ranks = "peptide_id",
    carry_cols = tp_col
  )

  testthat::expect_true("group" %in% names(out))
  testthat::expect_true(tp_col %in% names(out$group))
})

testthat::test_that("compute_alpha_diversity validates inputs and ranks", {
  info <- .get_ps_small_for_alpha()
  tp_col <- info$tp_col
  df_small <- dplyr::collect(info$ps$data_long)

  testthat::expect_error(
    compute_alpha_diversity(
      info$ps,
      group_cols = c("group", tp_col),
      ranks = "peptide_id",
      interaction_only = TRUE
    ),
    regexp = "(?i)interaction"
  )

  testthat::expect_error(
    compute_alpha_diversity(
      df_small,
      group_cols = "not_a_col",
      ranks = "peptide_id"
    ),
    regexp = "(?i)grouping columns"
  )

  testthat::expect_error(
    compute_alpha_diversity(
      info$ps,
      group_cols = "not_a_col",
      ranks = "peptide_id"
    ),
    regexp = "(?i)grouping columns"
  )

  testthat::expect_error(
    compute_alpha_diversity(
      df_small,
      group_cols = NULL,
      ranks = character(0)
    ),
    regexp = "(?i)ranks"
  )

  testthat::expect_warning(
    compute_alpha_diversity(
      df_small,
      group_cols = "group",
      ranks = c("peptide_id", "not_a_rank")
    ),
    regexp = "(?i)rank not found"
  )

  df_no_fc <- dplyr::select(df_small, -fold_change)

  testthat::expect_error(
    compute_alpha_diversity(
      df_no_fc,
      group_cols = NULL,
      ranks = "peptide_id",
      fc_threshold = 1.5
    ),
    regexp = "(?i)missing required columns"
  )

  testthat::expect_error(
    compute_alpha_diversity(
      list(),
      group_cols = NULL,
      ranks = "peptide_id"
    ),
    regexp = "(?i)phip_data|data.frame"
  )
})

testthat::test_that("compute_alpha_diversity handles data.frame rank mapping", {
  info <- .get_ps_small_for_alpha()
  df_small <- dplyr::collect(info$ps$data_long)

  df_small$family <- substr(df_small$peptide_id, 1L, 1L)

  out <- compute_alpha_diversity(
    df_small,
    group_cols = "group",
    ranks = c("peptide_id", "family"),
    shannon_log = "log2"
  )

  testthat::expect_true(all(c("peptide_id", "family") %in% out$group$rank))
})
