# # tests/testthat/test-compute_alpha_diversity.R
# # ==============================================================================
# # unit tests for compute_alpha_diversity using ps_merged_box_bin
# # heavy tests are enabled only when PHIPER_RUN_HEAVY_TESTS=1 and not on cran
# # ==============================================================================
#
# testthat::test_that("compute_alpha_diversity <phip_data> path works on ps_merged_box_bin", {
#   # skip conditions ------------------------------------------------------------
#   testthat::skip_if_not(Sys.getenv("PHIPER_RUN_HEAVY_TESTS", "0") == "1", "heavy tests off")
#   testthat::skip_on_cran()
#   skip_if_no_heavy()
#
#   # sanity on fixture ----------------------------------------------------------
#   pd <- get_heavy_ps()
#   testthat::expect_s3_class(pd, "phip_data")
#   testthat::expect_true(all(c("data_long", "meta") %in% names(pd)))
#   testthat::expect_true(all(c("sample_id", "peptide_id", "exist") %in% colnames(pd$data_long)))
#
#   # call generic two grouping columns default options no interaction -----------
#   out <- compute_alpha_diversity(
#     pd,
#     group_cols = c("big_group", "timepoint_factor"),
#     ranks = "peptide_id"
#   )
#
#   # structure ------------------------------------------------------------------
#   testthat::expect_type(out, "list")
#   testthat::expect_s3_class(out, "phip_alpha_diversity")
#   testthat::expect_setequal(names(out), c("big_group", "timepoint_factor"))
#
#   # each element has the required columns and no negative metrics --------------
#   req_cols <- c("rank", "sample_id", "richness", "shannon_diversity", "simpson_diversity")
#   for (nm in names(out)) {
#     df <- out[[nm]]
#     testthat::expect_s3_class(df, "tbl_df")
#     testthat::expect_true(all(req_cols %in% names(df)))
#     testthat::expect_true(nm %in% names(df)) # grouping column present
#     testthat::expect_true(all(df$richness >= 0, na.rm = TRUE))
#     testthat::expect_true(all(df$shannon_diversity >= 0, na.rm = TRUE))
#     testthat::expect_true(all(df$simpson_diversity >= 0, na.rm = TRUE))
#     testthat::expect_true(all(df$rank == "peptide_id"))
#   }
#
#   # attributes propagated -------------------------------------------------------
#   testthat::expect_identical(attr(out, "group_cols"), c("big_group", "timepoint_factor"))
#   testthat::expect_identical(attr(out, "ranks"), "peptide_id")
# })
#
# testthat::test_that("compute_alpha_diversity data.frame path on filtered ps$data_long", {
#   # skip conditions ------------------------------------------------------------
#   testthat::skip_if_not(Sys.getenv("PHIPER_RUN_HEAVY_TESTS", "0") == "1", "heavy tests off")
#   testthat::skip_on_cran()
#   skip_if_no_heavy()
#
#   pd <- get_heavy_ps()
#   testthat::expect_true("data_long" %in% names(pd))
#
#   # filter to big_group == kid_serum timepoint_factor in T2 T6 -----------------
#   .data <- rlang::.data
#   df_small <- pd$data_long |>
#     dplyr::filter(
#       .data$big_group == "kid_serum",
#       .data$timepoint_factor %in% c("T2", "T6")
#     ) |>
#     dplyr::collect()
#
#   # minimal columns present ----------------------------------------------------
#   testthat::expect_true(all(c("sample_id", "peptide_id", "exist", "timepoint_factor") %in% colnames(df_small)))
#   testthat::expect_true(nrow(df_small) > 0)
#
#   # compute by timepoint_factor via generic -----------------------------------
#   out_df <- compute_alpha_diversity(
#     df_small,
#     group_cols = "timepoint_factor",
#     ranks = "peptide_id"
#   )
#
#   testthat::expect_s3_class(out_df, "phip_alpha_diversity")
#   testthat::expect_identical(names(out_df), "timepoint_factor")
#
#   res <- out_df$timepoint_factor
#   testthat::expect_true(all(c("T2", "T6") %in% unique(res$timepoint_factor)))
#   testthat::expect_true(all(c(
#     "rank", "sample_id",
#     "richness", "shannon_diversity", "simpson_diversity",
#     "timepoint_factor"
#   ) %in% names(res)))
#
#   # metrics are well formed -----------------------------------------------------
#   testthat::expect_true(all(res$richness >= 0, na.rm = TRUE))
#   testthat::expect_true(all(res$shannon_diversity >= 0, na.rm = TRUE))
#   testthat::expect_true(all(res$simpson_diversity >= 0, na.rm = TRUE))
#   testthat::expect_true(all(res$rank == "peptide_id"))
#
#   # there is at least one sample per group -------------------------------------
#   n_by_group <- res |>
#     dplyr::count(.data$timepoint_factor, name = "n")
#   n_map <- stats::setNames(n_by_group$n, n_by_group$timepoint_factor)
#   testthat::expect_true(all(n_map[c("T2", "T6")] > 0))
#
#   # optional sanity not all summaries are na -----------------------------------
#   testthat::expect_false(all(is.na(res$richness)))
#   testthat::expect_false(all(is.na(res$shannon_diversity)))
#   testthat::expect_false(all(is.na(res$simpson_diversity)))
# })
#
# testthat::test_that("direct s3 methods execute for binary_alpha", {
#   # skip conditions ------------------------------------------------------------
#   testthat::skip_if_not(Sys.getenv("PHIPER_RUN_HEAVY_TESTS", "0") == "1", "heavy tests off")
#   testthat::skip_on_cran()
#   skip_if_no_heavy()
#
#   # hit the .phip_data method body directly -----------------------------------
#   ps <- get_heavy_ps()
#   res1 <- phiper:::compute_alpha_diversity(
#     ps,
#     group_cols = c("big_group", "timepoint_factor"),
#     ranks = "peptide_id"
#   )
#   testthat::expect_s3_class(res1, "phip_alpha_diversity")
#   testthat::expect_setequal(names(res1), c("big_group", "timepoint_factor"))
#
#   # build a small df and hit the .data.frame method directly -------------------
#   .data <- rlang::.data
#   df_small <- ps$data_long |>
#     dplyr::filter(
#       .data$big_group == "kid_serum",
#       .data$timepoint_factor %in% c("T2", "T6")
#     ) |>
#     dplyr::collect()
#
#   res2 <- phiper:::compute_alpha_diversity(
#     df_small,
#     group_cols = "timepoint_factor",
#     ranks = "peptide_id"
#   )
#   testthat::expect_s3_class(res2, "phip_alpha_diversity")
#   testthat::expect_identical(names(res2), "timepoint_factor")
# })
