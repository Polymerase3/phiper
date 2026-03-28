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

  ps <- load_example_data()

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
    mode = "threshold", threshold = 1.5,
    shannon_base = "log10",
    group_interaction = TRUE,
    interaction_only = TRUE
  )

  combo_nm <- paste(c("group", tp_col), collapse = " * ")
  testthat::expect_identical(names(out), combo_nm)
  testthat::expect_true(isTRUE(attr(out, "interaction_only")))
  testthat::expect_identical(attr(out, "threshold"), 1.5)
  testthat::expect_identical(attr(out, "shannon_base"), "log10")

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
      mode = "threshold", threshold = 1.5
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

testthat::test_that("compute_alpha_diversity `metrics` subset produces only requested columns", {
  info <- .get_ps_small_for_alpha()
  df_small <- dplyr::collect(info$ps$data_long)

  out <- compute_alpha_diversity(
    df_small,
    group_cols = NULL,
    ranks = "peptide_id",
    metrics = c("richness", "shannon")
  )

  res <- out$all_samples
  testthat::expect_true(all(c("richness", "shannon_diversity") %in% names(res)))
  testthat::expect_false("simpson_diversity" %in% names(res))
  testthat::expect_false("pielou_evenness"   %in% names(res))
  testthat::expect_false("berger_parker_dominance" %in% names(res))
  testthat::expect_identical(attr(out, "metrics"), c("richness", "shannon"))
})

testthat::test_that("compute_alpha_diversity single metric works and stores attribute", {
  info <- .get_ps_small_for_alpha()
  df_small <- dplyr::collect(info$ps$data_long)

  out <- compute_alpha_diversity(
    df_small,
    group_cols = "group",
    ranks = "peptide_id",
    metrics = "richness"
  )

  testthat::expect_true("richness" %in% names(out$group))
  testthat::expect_false("shannon_diversity" %in% names(out$group))
  testthat::expect_identical(attr(out, "metrics"), "richness")
})

testthat::test_that("compute_alpha_diversity returns pielou_evenness and berger_parker_dominance", {
  info <- .get_ps_small_for_alpha()
  df_small <- dplyr::collect(info$ps$data_long)

  out <- compute_alpha_diversity(
    df_small,
    group_cols = NULL,
    ranks = "peptide_id"
  )

  res <- out$all_samples
  testthat::expect_true("pielou_evenness"         %in% names(res))
  testthat::expect_true("berger_parker_dominance" %in% names(res))
  # values in [0, 1] for non-empty samples
  present <- res[res$richness > 0, ]
  testthat::expect_true(all(present$pielou_evenness >= 0 & present$pielou_evenness <= 1, na.rm = TRUE))
  testthat::expect_true(all(present$berger_parker_dominance > 0 & present$berger_parker_dominance <= 1, na.rm = TRUE))
})

testthat::test_that("pielou_evenness is NA when richness == 1", {
  df_one <- tibble::tibble(
    sample_id  = "s1",
    peptide_id = "p1",
    exist      = 1L
  )
  out <- compute_alpha_diversity(df_one, group_cols = NULL, ranks = "peptide_id")
  res <- out$all_samples
  testthat::expect_equal(res$richness, 1L)
  testthat::expect_true(is.na(res$pielou_evenness))
  # berger_parker is defined (= 1) when richness = 1
  testthat::expect_equal(res$berger_parker_dominance, 1)
})

testthat::test_that("empty presence set gives NA for pielou_evenness and berger_parker_dominance", {
  df_empty <- tibble::tibble(
    sample_id  = c("s1", "s2"),
    peptide_id = c("p1", "p2"),
    exist      = 0L
  )
  out <- compute_alpha_diversity(df_empty, group_cols = NULL, ranks = "peptide_id")
  res <- out$all_samples
  testthat::expect_true(all(is.na(res$pielou_evenness)))
  testthat::expect_true(all(is.na(res$berger_parker_dominance)))
})

testthat::test_that("compute_alpha_diversity mode = 'abundance' uses raw values", {
  df <- tibble::tibble(
    sample_id  = c("s1", "s1", "s2", "s2"),
    peptide_id = c("p1", "p2", "p1", "p2"),
    exist      = 1L,
    counts_hits = c(100L, 10L, 5L, 5L)
  )

  out <- compute_alpha_diversity(
    df,
    group_cols = NULL,
    ranks = "peptide_id",
    mode = "abundance",
    abundance_col = "counts_hits",
    metrics = c("richness", "berger_parker")
  )

  res <- out$all_samples
  testthat::expect_true("berger_parker_dominance" %in% names(res))
  # s1: max/sum = 100/110; s2: max/sum = 5/10 = 0.5
  s1 <- res[res$sample_id == "s1", ]
  s2 <- res[res$sample_id == "s2", ]
  testthat::expect_equal(s1$berger_parker_dominance, 100 / 110, tolerance = 1e-9)
  testthat::expect_equal(s2$berger_parker_dominance, 0.5,        tolerance = 1e-9)
})

testthat::test_that("compute_alpha_diversity mode = 'abundance' with higher-rank aggregation", {
  df <- tibble::tibble(
    sample_id  = rep("s1", 4),
    peptide_id = c("p1", "p2", "p3", "p4"),
    exist      = 1L,
    counts_hits = c(10L, 20L, 30L, 40L),
    family     = c("A", "A", "B", "B")
  )

  out <- compute_alpha_diversity(
    df,
    group_cols = NULL,
    ranks = "family",
    mode = "abundance",
    abundance_col = "counts_hits",
    abundance_agg = "sum",
    metrics = "richness"
  )

  res <- out$all_samples
  # two families -> richness = 2
  testthat::expect_equal(res$richness, 2L)
})

testthat::test_that("mode validation: threshold requires threshold scalar", {
  df <- tibble::tibble(
    sample_id = "s1", peptide_id = "p1", exist = 1L, fold_change = 2
  )
  testthat::expect_error(
    compute_alpha_diversity(df, ranks = "peptide_id", mode = "threshold"),
    regexp = "(?i)threshold.*scalar"
  )
  testthat::expect_error(
    compute_alpha_diversity(df, ranks = "peptide_id",
                            mode = "threshold", threshold = c(1, 2)),
    regexp = "(?i)threshold.*scalar"
  )
})

testthat::test_that("mode validation: abundance requires abundance_col", {
  df <- tibble::tibble(
    sample_id = "s1", peptide_id = "p1", exist = 1L, counts_hits = 5L
  )
  testthat::expect_error(
    compute_alpha_diversity(df, ranks = "peptide_id", mode = "abundance"),
    regexp = "(?i)abundance_col.*scalar"
  )
})

testthat::test_that("abundance_col not present in data triggers error", {
  df <- tibble::tibble(
    sample_id = "s1", peptide_id = "p1", exist = 1L
  )
  testthat::expect_error(
    compute_alpha_diversity(df, ranks = "peptide_id",
                            mode = "abundance", abundance_col = "no_such_col"),
    regexp = "(?i)not found"
  )
})

testthat::test_that("binary mode warns when abundance_col or threshold supplied", {
  df <- tibble::tibble(
    sample_id = "s1", peptide_id = "p1", exist = 1L, counts_hits = 5L
  )
  testthat::expect_warning(
    compute_alpha_diversity(df, ranks = "peptide_id",
                            mode = "binary", abundance_col = "counts_hits"),
    regexp = "(?i)ignored"
  )
  testthat::expect_warning(
    compute_alpha_diversity(df, ranks = "peptide_id",
                            mode = "binary", threshold = 1),
    regexp = "(?i)ignored"
  )
})

testthat::test_that("abundance mode warns when threshold supplied", {
  df <- tibble::tibble(
    sample_id = "s1", peptide_id = "p1", exist = 1L, counts_hits = 5L
  )
  testthat::expect_warning(
    compute_alpha_diversity(df, ranks = "peptide_id",
                            mode = "abundance", abundance_col = "counts_hits",
                            threshold = 1),
    regexp = "(?i)ignored"
  )
})

testthat::test_that("all-NULL ranks aborts with informative error", {
  df <- tibble::tibble(
    sample_id = "s1", peptide_id = "p1", exist = 1L
  )
  testthat::expect_error(
    compute_alpha_diversity(df, ranks = c("not_col_1", "not_col_2")),
    regexp = "(?i)no valid ranks"
  )
})

testthat::test_that("n_samples attribute is set correctly", {
  info <- .get_ps_small_for_alpha()
  df_small <- dplyr::collect(info$ps$data_long)
  n_expected <- dplyr::n_distinct(df_small$sample_id)

  out <- compute_alpha_diversity(df_small, group_cols = NULL, ranks = "peptide_id")

  testthat::expect_identical(attr(out, "n_samples"), n_expected)
})

testthat::test_that("custom interaction_sep is reflected in list names and column values", {
  info <- .get_ps_small_for_alpha()
  tp_col <- info$tp_col
  df_small <- dplyr::collect(info$ps$data_long)

  out <- compute_alpha_diversity(
    df_small,
    group_cols = c("group", tp_col),
    ranks = "peptide_id",
    group_interaction = TRUE,
    interaction_only = TRUE,
    interaction_sep = ":"
  )

  combo_nm <- paste(c("group", tp_col), collapse = ":")
  testthat::expect_identical(names(out), combo_nm)
  testthat::expect_true(all(grepl(":", out[[combo_nm]]$phip_interaction)))
})

testthat::test_that("compute_alpha_diversity group_cols = NULL on phip_data", {
  info <- .get_ps_small_for_alpha()
  ps <- info$ps

  out <- compute_alpha_diversity(ps, group_cols = NULL, ranks = "peptide_id")

  testthat::expect_s3_class(out, "phip_alpha_diversity")
  testthat::expect_identical(names(out), "all_samples")
  testthat::expect_true("group" %in% names(out$all_samples))
})

testthat::test_that("shannon_base log2 and log10 give different values than ln", {
  df <- tibble::tibble(
    sample_id  = rep("s1", 3),
    peptide_id = c("p1", "p2", "p3"),
    exist      = 1L
  )

  out_ln   <- compute_alpha_diversity(df, group_cols = NULL, ranks = "peptide_id",
                                      shannon_base = "ln",    metrics = "shannon")
  out_log2 <- compute_alpha_diversity(df, group_cols = NULL, ranks = "peptide_id",
                                      shannon_base = "log2",  metrics = "shannon")
  out_log10 <- compute_alpha_diversity(df, group_cols = NULL, ranks = "peptide_id",
                                       shannon_base = "log10", metrics = "shannon")

  h_ln   <- out_ln$all_samples$shannon_diversity
  h_log2 <- out_log2$all_samples$shannon_diversity
  h_log10 <- out_log10$all_samples$shannon_diversity

  testthat::expect_true(h_log2 > h_ln)     # log2 base < e, so H_log2 > H_ln
  testthat::expect_true(h_ln   > h_log10)  # log10 base > e, so H_log10 < H_ln
  # consistency: H_log2 = H_ln / log(2)  ->  H_log2 * log(2) == H_ln
  testthat::expect_equal(h_log2 * log(2), h_ln, tolerance = 1e-9)
})
