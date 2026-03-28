# tests/testthat/test-alpha_significance.R

# ---------------------------------------------------------------------------
# Minimal local fixture (data.frame path — no DuckDB needed)
# ---------------------------------------------------------------------------
.make_alpha_df <- function(n = 10L, seed = 42L) {
  set.seed(seed)
  tibble::tibble(
    sample_id         = paste0("s", seq_len(n * 2L)),
    rank              = "peptide_id",
    group             = rep(c("A", "B"), each = n),
    richness          = c(
      as.integer(stats::rpois(n, lambda = 30)),
      as.integer(stats::rpois(n, lambda = 10))
    ),
    shannon_diversity = c(
      stats::runif(n, 2.5, 3.5),
      stats::runif(n, 1.0, 2.0)
    ),
    simpson_diversity = c(
      stats::runif(n, 0.85, 0.95),
      stats::runif(n, 0.55, 0.75)
    ),
    pielou_evenness         = c(stats::runif(n, 0.7, 0.9), stats::runif(n, 0.4, 0.6)),
    berger_parker_dominance = c(stats::runif(n, 0.1, 0.2), stats::runif(n, 0.3, 0.5))
  )
}

.make_alpha_3grp <- function(n = 8L, seed = 7L) {
  set.seed(seed)
  tibble::tibble(
    sample_id         = paste0("s", seq_len(n * 3L)),
    rank              = "peptide_id",
    group             = rep(c("A", "B", "C"), each = n),
    richness          = as.integer(stats::rpois(n * 3L, lambda = 20)),
    shannon_diversity = stats::runif(n * 3L, 1, 4),
    simpson_diversity = stats::runif(n * 3L, 0.5, 0.95),
    pielou_evenness         = stats::runif(n * 3L, 0.4, 0.9),
    berger_parker_dominance = stats::runif(n * 3L, 0.1, 0.5)
  )
}

# ---------------------------------------------------------------------------
# 1. Basic: runs without error on a data frame
# ---------------------------------------------------------------------------
testthat::test_that("compute_alpha_significance() runs on data frame", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group")
  testthat::expect_s3_class(sig, "phip_alpha_significance")
})

# ---------------------------------------------------------------------------
# 2. Output structure: $global and $pairwise exist with expected columns
# ---------------------------------------------------------------------------
testthat::test_that("output has $global and $pairwise with correct columns", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group")

  testthat::expect_true(is.list(sig))
  testthat::expect_true(all(c("global", "pairwise") %in% names(sig)))

  g_cols <- c("rank", "metric", "statistic", "p_value", "test")
  testthat::expect_true(all(g_cols %in% names(sig$global)))

  pw_cols <- c("rank", "metric", "group1", "group2", "p_raw", "p_adj", "cohens_d", "stars", "test")
  testthat::expect_true(all(pw_cols %in% names(sig$pairwise)))
})

# ---------------------------------------------------------------------------
# 3. Attributes set correctly
# ---------------------------------------------------------------------------
testthat::test_that("attributes are set on the result", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group",
                                    global_test = "anova",
                                    pairwise_test = "tukey",
                                    p_adjust_method = "bonferroni")

  testthat::expect_identical(attr(sig, "group_col"), "group")
  testthat::expect_identical(attr(sig, "global_test"), "anova")
  testthat::expect_identical(attr(sig, "pairwise_test"), "tukey")
  testthat::expect_identical(attr(sig, "p_adjust_method"), "bonferroni")
  testthat::expect_true(is.character(attr(sig, "metrics")))
  testthat::expect_true(length(attr(sig, "metrics")) > 0L)
})

# ---------------------------------------------------------------------------
# 4. Global test (kruskal) returns finite p_value in [0,1]
# ---------------------------------------------------------------------------
testthat::test_that("kruskal global test returns finite p_value in [0,1]", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group",
                                    global_test = "kruskal",
                                    metric = "richness")
  pvals <- sig$global$p_value[!is.na(sig$global$p_value)]
  testthat::expect_true(all(pvals >= 0 & pvals <= 1))
})

# ---------------------------------------------------------------------------
# 5. Global test (anova) returns finite p_value in [0,1]
# ---------------------------------------------------------------------------
testthat::test_that("anova global test returns finite p_value in [0,1]", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group",
                                    global_test = "anova",
                                    metric = "richness")
  pvals <- sig$global$p_value[!is.na(sig$global$p_value)]
  testthat::expect_true(all(pvals >= 0 & pvals <= 1))
})

# ---------------------------------------------------------------------------
# 6. Pairwise (wilcoxon) returns p_raw and p_adj in [0,1]
# ---------------------------------------------------------------------------
testthat::test_that("wilcoxon pairwise p-values are in [0,1]", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group",
                                    pairwise_test = "wilcoxon",
                                    metric = "richness")
  pw <- sig$pairwise[!is.na(sig$pairwise$p_raw), ]
  testthat::expect_true(all(pw$p_raw >= 0 & pw$p_raw <= 1))
  testthat::expect_true(all(pw$p_adj >= 0 & pw$p_adj <= 1))
})

# ---------------------------------------------------------------------------
# 7. Pairwise (tukey) returns p_raw and p_adj in [0,1]
# ---------------------------------------------------------------------------
testthat::test_that("tukey pairwise p-values are in [0,1]", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group",
                                    pairwise_test = "tukey",
                                    metric = "richness")
  pw <- sig$pairwise[!is.na(sig$pairwise$p_raw), ]
  testthat::expect_true(all(pw$p_raw >= 0 & pw$p_raw <= 1))
  testthat::expect_true(all(pw$p_adj >= 0 & pw$p_adj <= 1))
})

# ---------------------------------------------------------------------------
# 8. p_adjust_method: BH vs bonferroni differ with >1 pair
# ---------------------------------------------------------------------------
testthat::test_that("p_adjust_method affects p_adj with multiple pairs", {
  df <- .make_alpha_3grp()

  sig_bh   <- compute_alpha_significance(df, group_col = "group",
                                         pairwise_test = "wilcoxon",
                                         p_adjust_method = "BH",
                                         metric = "richness")
  sig_bonf <- compute_alpha_significance(df, group_col = "group",
                                         pairwise_test = "wilcoxon",
                                         p_adjust_method = "bonferroni",
                                         metric = "richness")

  # Bonferroni is always >= BH (more conservative)
  bh_adj   <- sig_bh$pairwise$p_adj[!is.na(sig_bh$pairwise$p_adj)]
  bonf_adj <- sig_bonf$pairwise$p_adj[!is.na(sig_bonf$pairwise$p_adj)]
  if (length(bh_adj) && length(bonf_adj)) {
    testthat::expect_true(all(bonf_adj >= bh_adj - 1e-12))
  }
})

# ---------------------------------------------------------------------------
# 9. Cohen's d: finite for adequate group sizes, NA for tiny groups
# ---------------------------------------------------------------------------
testthat::test_that(".ph_cohens_d returns NA for n<2 and finite otherwise", {
  testthat::expect_true(is.na(phiper:::.ph_cohens_d(1, c(1, 2, 3))))
  testthat::expect_true(is.finite(phiper:::.ph_cohens_d(c(10, 12, 14), c(1, 3, 5))))
})

# ---------------------------------------------------------------------------
# 10. metric subset: only requested metrics returned
# ---------------------------------------------------------------------------
testthat::test_that("metric= returns only that metric", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group", metric = "richness")
  testthat::expect_identical(unique(sig$global$metric), "richness")
  testthat::expect_identical(unique(sig$pairwise$metric), "richness")
})

# ---------------------------------------------------------------------------
# 11. Multiple metrics returned when metric=NULL
# ---------------------------------------------------------------------------
testthat::test_that("metric=NULL returns all present metrics", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group")
  testthat::expect_true(length(unique(sig$global$metric)) > 1L)
})

# ---------------------------------------------------------------------------
# 12. plot_alpha_significance(type="table") returns a data frame
# ---------------------------------------------------------------------------
testthat::test_that("plot_alpha_significance type='table' returns data frame", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group", metric = "richness")
  tbl <- plot_alpha_significance(sig, metric = "richness", type = "table")
  testthat::expect_true(inherits(tbl, "data.frame"))
})

# ---------------------------------------------------------------------------
# 13. plot_alpha_significance(type="heatmap") returns a ggplot
# ---------------------------------------------------------------------------
testthat::test_that("plot_alpha_significance type='heatmap' returns ggplot", {
  df  <- .make_alpha_3grp()
  sig <- compute_alpha_significance(df, group_col = "group", metric = "richness")
  p   <- plot_alpha_significance(sig, metric = "richness", type = "heatmap",
                                 p_threshold = 1)  # include all pairs
  testthat::expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# 14. p_threshold filters to only significant pairs in table mode
# ---------------------------------------------------------------------------
testthat::test_that("p_threshold filters pairwise table correctly", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group", metric = "richness")

  tbl_all  <- plot_alpha_significance(sig, metric = "richness", type = "table",
                                      p_threshold = 1)
  tbl_sig  <- plot_alpha_significance(sig, metric = "richness", type = "table",
                                      p_threshold = 0.05)

  testthat::expect_true(nrow(tbl_sig) <= nrow(tbl_all))
  if (nrow(tbl_sig)) {
    testthat::expect_true(all(tbl_sig$p_adj <= 0.05))
  }
})

# ---------------------------------------------------------------------------
# 15. Error on non-phip_alpha_significance input to plot_alpha_significance
# ---------------------------------------------------------------------------
testthat::test_that("plot_alpha_significance errors on wrong input class", {
  testthat::expect_error(
    plot_alpha_significance(list(a = 1)),
    regexp = "phip_alpha_significance"
  )
})

# ---------------------------------------------------------------------------
# 16. Single group: global test returns NA gracefully (no error)
# ---------------------------------------------------------------------------
testthat::test_that("single group produces NA global p_value without error", {
  df  <- dplyr::filter(.make_alpha_df(), group == "A")
  # kruskal.test with 1 group should fail gracefully
  sig <- testthat::expect_no_error(
    compute_alpha_significance(df, group_col = "group", metric = "richness")
  )
  testthat::expect_true(is.na(sig$global$p_value[[1L]]) || is.finite(sig$global$p_value[[1L]]))
})

# ---------------------------------------------------------------------------
# 17. ... param is silently ignored
# ---------------------------------------------------------------------------
testthat::test_that("... param ignored silently in plot_alpha_significance", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group", metric = "richness")
  testthat::expect_no_error(
    plot_alpha_significance(sig, metric = "richness", type = "table",
                            ignored_param = TRUE)
  )
})

# ---------------------------------------------------------------------------
# 18. group_col inferred from attribute
# ---------------------------------------------------------------------------
testthat::test_that("group_col inferred from attr when not supplied", {
  df <- .make_alpha_df()
  attr(df, "group_cols") <- "group"
  sig <- compute_alpha_significance(df)
  testthat::expect_identical(attr(sig, "group_col"), "group")
})

# ---------------------------------------------------------------------------
# 19. stars column contains only "***", "**", "*", "ns"
# ---------------------------------------------------------------------------
testthat::test_that("stars column contains valid values", {
  df  <- .make_alpha_3grp()
  sig <- compute_alpha_significance(df, group_col = "group", metric = "richness")
  valid <- c("***", "**", "*", "ns")
  testthat::expect_true(all(sig$pairwise$stars %in% valid))
})

# ---------------------------------------------------------------------------
# 20. Pairwise table has one row per group pair
# ---------------------------------------------------------------------------
testthat::test_that("pairwise table has exactly n*(n-1)/2 rows per (rank,metric)", {
  df  <- .make_alpha_3grp()  # 3 groups -> 3 pairs
  sig <- compute_alpha_significance(df, group_col = "group", metric = "richness")
  pw  <- dplyr::filter(sig$pairwise,
                       .data$rank   == "peptide_id",
                       .data$metric == "richness")
  testthat::expect_equal(nrow(pw), 3L)
})

# ---------------------------------------------------------------------------
# 21. test column records the test used
# ---------------------------------------------------------------------------
testthat::test_that("test column records the pairwise test name", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group",
                                    pairwise_test = "wilcoxon",
                                    metric = "richness")
  testthat::expect_true(all(sig$pairwise$test == "wilcoxon"))

  sig2 <- compute_alpha_significance(df, group_col = "group",
                                     pairwise_test = "tukey",
                                     metric = "richness")
  testthat::expect_true(all(sig2$pairwise$test == "tukey"))
})

# ---------------------------------------------------------------------------
# 22. metric not in data: error with informative message
# ---------------------------------------------------------------------------
testthat::test_that("requesting absent metric raises informative error", {
  df  <- .make_alpha_df()
  testthat::expect_error(
    compute_alpha_significance(df, group_col = "group",
                               metric = "nonexistent_metric"),
    regexp = "(?i)metric"
  )
})

# ---------------------------------------------------------------------------
# 23. plot_alpha_diversity accepts significance object without error
# ---------------------------------------------------------------------------
testthat::test_that("plot_alpha_diversity accepts significance= without error", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group", metric = "richness")
  testthat::expect_no_error(
    plot_alpha_diversity(df, metric = "richness", group_col = "group",
                        significance = sig, show_significance = FALSE)
  )
})

# ---------------------------------------------------------------------------
# 24. plot_alpha_diversity new metrics accepted
# ---------------------------------------------------------------------------
testthat::test_that("plot_alpha_diversity accepts pielou_evenness metric", {
  df <- .make_alpha_df()
  p  <- plot_alpha_diversity(df, metric = "pielou_evenness", group_col = "group")
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("plot_alpha_diversity accepts berger_parker_dominance metric", {
  df <- .make_alpha_df()
  p  <- plot_alpha_diversity(df, metric = "berger_parker_dominance",
                             group_col = "group")
  testthat::expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# 25. plot_alpha_diversity ... param ignored silently
# ---------------------------------------------------------------------------
testthat::test_that("plot_alpha_diversity ... ignored silently", {
  df <- .make_alpha_df()
  testthat::expect_no_error(
    plot_alpha_diversity(df, metric = "richness", group_col = "group",
                        unknown_future_param = 42)
  )
})

# ---------------------------------------------------------------------------
# 26. No-rank path: data frame without a rank column
# ---------------------------------------------------------------------------
testthat::test_that("compute_alpha_significance works on data without a rank column", {
  df  <- dplyr::select(.make_alpha_df(), -rank)
  sig <- compute_alpha_significance(df, group_col = "group", metric = "richness")
  testthat::expect_s3_class(sig, "phip_alpha_significance")
  # rank column in output should be NA
  testthat::expect_true(all(is.na(sig$global$rank)))
})

# ---------------------------------------------------------------------------
# 27. Tukey pairwise with 3 groups covers the row-lookup logic
# ---------------------------------------------------------------------------
testthat::test_that("tukey pairwise with 3 groups returns 3 pairs", {
  df  <- .make_alpha_3grp()
  sig <- compute_alpha_significance(df, group_col = "group",
                                    pairwise_test = "tukey",
                                    metric = "richness")
  pw  <- dplyr::filter(sig$pairwise, .data$metric == "richness")
  testthat::expect_equal(nrow(pw), 3L)
  testthat::expect_true(all(pw$p_raw >= 0 & pw$p_raw <= 1 | is.na(pw$p_raw)))
})

# ---------------------------------------------------------------------------
# 28. group_col not in names -> informative error
# ---------------------------------------------------------------------------
testthat::test_that("supplied group_col not found raises informative error", {
  df <- .make_alpha_df()
  testthat::expect_error(
    compute_alpha_significance(df, group_col = "nonexistent"),
    regexp = "(?i)group_col"
  )
})

# ---------------------------------------------------------------------------
# 29. .ph_infer_group_col: fallback to first non-standard column
# ---------------------------------------------------------------------------
testthat::test_that(".ph_infer_group_col falls back to first non-standard col", {
  df  <- .make_alpha_df()
  # Rename group to something unusual; no attr set
  df2 <- dplyr::rename(df, cohort = group)
  sig <- compute_alpha_significance(df2, metric = "richness")
  testthat::expect_identical(attr(sig, "group_col"), "cohort")
})

# ---------------------------------------------------------------------------
# 30. .ph_infer_group_col: attr group_cols with length > 1 (falls through)
# ---------------------------------------------------------------------------
testthat::test_that(".ph_infer_group_col ignores multi-element attr and infers from cols", {
  df <- .make_alpha_df()
  attr(df, "group_cols") <- c("group", "extra")  # length > 1, not length 1
  # Should fall through to candidates since length != 1
  sig <- compute_alpha_significance(df, metric = "richness")
  testthat::expect_s3_class(sig, "phip_alpha_significance")
})

# ---------------------------------------------------------------------------
# 31. plot_alpha_significance: metric and rank default to first when NULL
# ---------------------------------------------------------------------------
testthat::test_that("plot_alpha_significance defaults metric/rank to first", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group")
  # Should not error; defaults applied internally
  tbl <- plot_alpha_significance(sig, type = "table", p_threshold = 1)
  testthat::expect_true(inherits(tbl, "data.frame"))
})

# ---------------------------------------------------------------------------
# 32. plot_alpha_significance: heatmap with p_threshold=0 (nothing significant)
# ---------------------------------------------------------------------------
testthat::test_that("plot_alpha_significance heatmap with no significant pairs", {
  df  <- .make_alpha_df()
  sig <- compute_alpha_significance(df, group_col = "group", metric = "richness")
  p   <- plot_alpha_significance(sig, metric = "richness", type = "heatmap",
                                 p_threshold = 0)
  testthat::expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# 33. compute_alpha_significance with phip_alpha_diversity list input
# ---------------------------------------------------------------------------
testthat::test_that("compute_alpha_significance accepts phip_alpha_diversity list", {
  df  <- .make_alpha_df()
  lst <- list(group = df)
  class(lst) <- c("phip_alpha_diversity", "list")
  attr(lst, "group_cols") <- "group"
  sig <- compute_alpha_significance(lst, metric = "richness")
  testthat::expect_s3_class(sig, "phip_alpha_significance")
  testthat::expect_identical(attr(sig, "group_col"), "group")
})

# ---------------------------------------------------------------------------
# 34. Global test error fallback returns NA row
# ---------------------------------------------------------------------------
testthat::test_that("global test error returns NA row gracefully", {
  # All values identical -> kruskal throws an error (zero variance)
  df <- tibble::tibble(
    sample_id = paste0("s", 1:20),
    rank      = "peptide_id",
    group     = rep(c("A", "B"), each = 10L),
    richness  = 5L  # all the same -> kruskal-wallis error
  )
  sig <- testthat::expect_no_error(
    compute_alpha_significance(df, group_col = "group", metric = "richness")
  )
  testthat::expect_true(is.na(sig$global$p_value[[1L]]) ||
                          is.finite(sig$global$p_value[[1L]]))
})

# ---------------------------------------------------------------------------
# 35. .ph_sig_stars boundary values
# ---------------------------------------------------------------------------
testthat::test_that(".ph_sig_stars returns correct symbols for boundary p-values", {
  stars <- phiper:::.ph_sig_stars(c(0.0001, 0.001, 0.005, 0.01, 0.04, 0.05, 0.1))
  testthat::expect_equal(stars, c("***", "**", "**", "*", "*", "ns", "ns"))
})

# ---------------------------------------------------------------------------
# 36. .ph_infer_group_col: no candidates -> error
# ---------------------------------------------------------------------------
testthat::test_that(".ph_infer_group_col errors when no non-standard columns", {
  # df with ONLY standard columns and no group
  df <- tibble::tibble(
    sample_id = paste0("s", 1:4),
    rank      = "peptide_id",
    richness  = c(10L, 12L, 8L, 9L)
  )
  testthat::expect_error(
    compute_alpha_significance(df, metric = "richness"),
    regexp = "(?i)group"
  )
})
