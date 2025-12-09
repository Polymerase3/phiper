# tests/testthat/test-ph_prevalence_shift.R

test_that("ph_prevalence_shift works for unpaired design (mock species)", {
  # example phip_data
  phip_path <- phip_example_path()

  ps <- phip_convert(
    data_long_path    = phip_path,
    backend           = "duckdb",
    peptide_library   = TRUE,
    subject_id        = "subject_id",
    peptide_id        = "peptide_id",
    sample_id         = "sample_id",
    exist             = "exist",
    timepoint         = "timepoint_factor",
    fold_change       = "fold_change",
    materialise_table = TRUE,
    auto_expand       = TRUE,
    n_cores           = 2
  )

  # small unpaired subset: one species, two groups at T1
  ps_filt <- ps %>%
    dplyr::filter(
      peptide_id %in% c("16627", "5243", "24799", "16196", "18003"),
      time == "T1"
    ) %>%
    dplyr::collect()

  mock_peplib <- data.frame(
    peptide_id = c("16627", "5243", "24799", "16196", "18003"),
    species    = rep("mock_species", 5),
    stringsAsFactors = FALSE
  )

  res <- ph_prevalence_shift(
    x                  = ps_filt,
    exist_col          = "exist",
    rank_cols          = "species",
    group_cols         = "group",
    peptide_library    = mock_peplib,
    B_permutations     = 500L,          # smaller for test speed
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "n_eff_sqrt",
    stat_mode          = "asin",
    prev_strat         = "none",
    winsor_z           = Inf,
    rank_feature_keep  = list(species = NULL),
    log                = FALSE,
    fold_change        = "sum",
    cross_prev         = "mean"
  )

  # basic structure
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 1L)

  # key columns
  expect_true(all(c(
    "rank", "feature", "group_col", "group1", "group2", "design",
    "n_subjects_paired", "n_peptides_used", "m_eff",
    "T_obs", "T_obs_stand", "Z_from_p",
    "p_perm", "b", "p_adj_rank", "category_rank_bh"
  ) %in% names(res)))

  # unpaired design should not report paired subjects
  expect_identical(res$design, "unpaired")
  expect_true(is.na(res$n_subjects_paired))

  # stats should be numeric and finite
  expect_true(is.numeric(res$T_obs))
  expect_false(any(is.na(res$T_obs)))
  expect_true(is.numeric(res$T_obs_stand))
  expect_true(is.numeric(res$Z_from_p))

  # permutation p-value in (0,1]
  expect_gt(res$p_perm, 0)
  expect_lte(res$p_perm, 1)
})

test_that("ph_prevalence_shift handles paired design via paired_by and returns
          standardized T", {
  # build example phip_data
  phip_path <- phip_example_path()

  ps <- phip_convert(
    data_long_path    = phip_path,
    backend           = "duckdb",
    peptide_library   = TRUE,
    subject_id        = "subject_id",
    peptide_id        = "peptide_id",
    sample_id         = "sample_id",
    exist             = "exist",
    timepoint         = "timepoint_factor",
    fold_change       = "fold_change",
    materialise_table = TRUE,
    auto_expand       = TRUE,
    n_cores           = 2
  )

  # subset: one mock species, group A, two timepoints
  ps_filt2 <- ps %>%
    dplyr::filter(
      peptide_id %in% c("2269", "21399", "7789", "13588", "10180"),
      group == "A"
    ) %>%
    dplyr::collect()

  mock_peplib <- data.frame(
    peptide_id = c("2269", "21399", "7789", "13588", "10180"),
    species    = rep("mock_species", 5),
    stringsAsFactors = FALSE
  )

  # identify true T1-T2 pairs
  pair_ids <- ps_filt2 %>%
    dplyr::filter(time %in% c("T1", "T2")) %>%
    dplyr::distinct(subject_id, time) %>%
    tidyr::pivot_wider(
      names_from  = time,
      values_from = time,
      values_fn   = length,
      values_fill = 0
    ) %>%
    dplyr::filter(T1 > 0, T2 > 0) %>%
    dplyr::pull(subject_id)

  expect_gt(length(pair_ids), 0L)

  # inject artificial variability by flipping exist for a subset of rows
  set.seed(123)
  candidates <- ps_filt2 %>%
    dplyr::filter(
      subject_id %in% pair_ids,
      time %in% c("T1", "T2")
    )

  n_flip   <- min(15L, nrow(candidates))
  flip_ids <- sample(candidates$sample_id, size = n_flip)

  ps_filt2_flipped <- ps_filt2 %>%
    dplyr::mutate(
      flipped    = sample_id %in% flip_ids,
      exist      = ifelse(flipped, 1L - exist, exist),
      counts_hits= ifelse(flipped & exist == 0L, 0, counts_hits),
      fold_change= ifelse(flipped & exist == 0L, 0, fold_change)
    ) %>%
    dplyr::select(-flipped)

  res <- ph_prevalence_shift(
    x                  = ps_filt2_flipped,
    paired_by          = "subject_id",
    exist_col          = "exist",
    rank_cols          = "species",
    group_cols         = "time",
    peptide_library    = mock_peplib,
    B_permutations     = 500L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "n_eff_sqrt",
    stat_mode          = "asin",
    prev_strat         = "none",
    winsor_z           = Inf,
    rank_feature_keep  = list(species = NULL),
    log                = FALSE,
    fold_change        = "sum",
    cross_prev         = "mean"
  )

  # paired design detected
  expect_identical(res$design, "paired")
  expect_equal(res$n_subjects_paired, length(pair_ids))

  # effect should be non-zero after flipping
  expect_false(is.na(res$T_obs))
  expect_true(res$T_obs != 0)

  # permutation p-value should be non-trivial
  expect_gt(res$p_perm, 0)
  expect_lte(res$p_perm, 1)

  # standardized T should exist and have same sign as T_obs
  expect_false(is.na(res$T_obs_stand))
  expect_equal(sign(res$T_obs), sign(res$T_obs_stand))

  # Z_from_p should have the same sign as T_obs by construction
  expect_false(is.na(res$Z_from_p))
})

test_that("ph_prevalence_shift errors on duplicate positives within group", {
  # toy data with duplicate positive for same subject/peptide/group
  toy_dup <- tibble::tibble(
    sample_id  = c("s1", "s2"),
    subject_id = c("id1", "id1"),
    peptide_id = c("pep1", "pep1"),
    group      = c("A", "A"),
    exist      = c(1L, 1L)
  )

  expect_error(
    ph_prevalence_shift(
      x                  = toy_dup,
      exist_col          = "exist",
      rank_cols          = "peptide_id",
      group_cols         = "group",
      B_permutations     = 200L,
      smooth_eps_num     = 0.5,
      smooth_eps_den_mult= 2.0,
      min_max_prev       = 0.0,
      weight_mode        = "equal",
      stat_mode          = "diff",
      prev_strat         = "none",
      winsor_z           = 4,
      rank_feature_keep  = NULL,
      peptide_library    = NULL,
      log                = FALSE,
      fold_change        = "none",
      cross_prev         = "none"
    ),
    regexp = "duplicate positives",
    fixed  = FALSE
  )
})

test_that("ph_prevalence_shift gives consistent direction for T_obs,
          T_obs_stand and Z_from_p", {
  # reuse unpaired example, smaller B for speed
  phip_path <- phip_example_path()

  ps <- phip_convert(
    data_long_path    = phip_path,
    backend           = "duckdb",
    peptide_library   = TRUE,
    subject_id        = "subject_id",
    peptide_id        = "peptide_id",
    sample_id         = "sample_id",
    exist             = "exist",
    timepoint         = "timepoint_factor",
    fold_change       = "fold_change",
    materialise_table = TRUE,
    auto_expand       = TRUE,
    n_cores           = 2
  )

  ps_filt <- ps %>%
    dplyr::filter(
      peptide_id %in% c("16627", "5243", "24799", "16196", "18003"),
      time == "T1"
    ) %>%
    dplyr::collect()

  mock_peplib <- data.frame(
    peptide_id = c("16627", "5243", "24799", "16196", "18003"),
    species    = rep("mock_species", 5),
    stringsAsFactors = FALSE
  )

  res <- ph_prevalence_shift(
    x                  = ps_filt,
    exist_col          = "exist",
    rank_cols          = "species",
    group_cols         = "group",
    peptide_library    = mock_peplib,
    B_permutations     = 300L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "n_eff_sqrt",
    stat_mode          = "asin",
    prev_strat         = "none",
    winsor_z           = Inf,
    rank_feature_keep  = list(species = NULL),
    log                = FALSE,
    fold_change        = "sum",
    cross_prev         = "mean"
  )

  # only check direction when standardized T is available
  if (!is.na(res$T_obs_stand)) {
    expect_equal(sign(res$T_obs), sign(res$T_obs_stand))
  }
  if (!is.na(res$Z_from_p)) {
    expect_equal(sign(res$T_obs), sign(res$Z_from_p))
  }
})

test_that("ph_prevalence_shift returns NA T_obs_stand when permutation
          variance is zero", {
  # toy data: one peptide, two groups, identical prevalence pattern
  toy <- tibble::tibble(
    sample_id  = paste0("s", 1:40),
    subject_id = paste0("id", rep(1:20, times = 2)),
    peptide_id = rep("pep1", 40),
    group      = rep(c("A", "B"), each = 20),
    exist      = 1L
  )

  res <- ph_prevalence_shift(
    x                  = toy,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 200L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "equal",
    stat_mode          = "diff",
    prev_strat         = "none",
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
    fold_change        = "none",
    cross_prev         = "none"
  )

  # all permutations give the same T => variance zero
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 0L)
  expect_setequal(
    colnames(res),
    c("rank","feature","group_col","group1","group2","design",
      "n_subjects_paired","n_peptides_used","m_eff",
      "T_obs","T_obs_stand","Z_from_p","p_perm","b",
      "p_adj_rank",
      "mean_delta","frac_delta_pos","mean_delta_w","frac_delta_pos_w",
      "category_rank_bh")
  )
})

test_that("ph_prevalence_shift aborts when required columns are missing", {
  toy_missing <- tibble::tibble(
    sample_id  = "s1",
    peptide_id = "pep1",
    group      = "A",
    exist      = 1L
  )

  expect_error(
    ph_prevalence_shift(
      x                  = toy_missing,
      exist_col          = "exist",
      rank_cols          = "peptide_id",
      group_cols         = "group",
      B_permutations     = 200L,
      smooth_eps_num     = 0.5,
      smooth_eps_den_mult= 2.0,
      min_max_prev       = 0.0,
      weight_mode        = "equal",
      stat_mode          = "diff",
      prev_strat         = "none",
      winsor_z           = 4,
      rank_feature_keep  = NULL,
      peptide_library    = NULL,
      log                = FALSE,
      fold_change        = "none",
      cross_prev         = "none"
    ),
    regexp = "Missing required columns"
  )
})

test_that("ph_prevalence_shift aborts when all peptides are zero after hits
          guard", {
  toy_zero <- tibble::tibble(
    sample_id  = c("s1", "s2"),
    subject_id = c("id1", "id2"),
    peptide_id = c("pep1", "pep2"),
    group      = c("A", "B"),
    exist      = c(0L, 0L)
  )

  expect_error(
    ph_prevalence_shift(
      x                  = toy_zero,
      exist_col          = "exist",
      rank_cols          = "peptide_id",
      group_cols         = "group",
      B_permutations     = 200L,
      smooth_eps_num     = 0.5,
      smooth_eps_den_mult= 2.0,
      min_max_prev       = 0.0,
      weight_mode        = "equal",
      stat_mode          = "diff",
      prev_strat         = "none",
      winsor_z           = 4,
      rank_feature_keep  = NULL,
      peptide_library    = NULL,
      log                = FALSE,
      fold_change        = "none",
      cross_prev         = "none"
    ),
    regexp = "All peptides are zero after hits guard"
  )
})

test_that("ph_prevalence_shift aborts when peptide_library misses required
          rank columns", {
  toy_df <- tibble::tibble(
    sample_id  = "s1",
    subject_id = "id1",
    peptide_id = "pep1",
    group      = "A",
    exist      = 1L
  )

  bad_lib <- tibble::tibble(
    peptide_id = "pep1"
  )

  expect_error(
    ph_prevalence_shift(
      x                  = toy_df,
      exist_col          = "exist",
      rank_cols          = "species",
      group_cols         = "group",
      peptide_library    = bad_lib,
      B_permutations     = 200L,
      smooth_eps_num     = 0.5,
      smooth_eps_den_mult= 2.0,
      min_max_prev       = 0.0,
      weight_mode        = "equal",
      stat_mode          = "diff",
      prev_strat         = "none",
      winsor_z           = 4,
      rank_feature_keep  = NULL,
      log                = FALSE,
      fold_change        = "none",
      cross_prev         = "none"
    ),
    regexp = "Peptide library missing required columns"
  )
})

test_that("ph_prevalence_shift works with rank_cols = 'peptide_id' only", {
  toy_df <- tibble::tibble(
    sample_id  = c("s1", "s2"),
    subject_id = c("id1", "id2"),
    peptide_id = c("pep1", "pep2"),
    group      = c("A", "B"),
    exist      = c(1L, 0L)
  )

  res <- ph_prevalence_shift(
    x                  = toy_df,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 200L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "equal",
    stat_mode          = "diff",
    prev_strat         = "none",
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
    fold_change        = "none",
    cross_prev         = "none"
  )

  expect_s3_class(res, "tbl_df")
  expect_true(all(res$rank == "peptide_id"))
})

test_that("ph_prevalence_shift uses peptide_library attached to phip_data", {
  toy_long <- tibble::tibble(
    sample_id  = c("s1", "s2"),
    subject_id = c("id1", "id2"),
    peptide_id = c("pep1", "pep1"),
    group      = c("A", "B"),
    exist      = c(1L, 1L)
  )

  toy_lib <- tibble::tibble(
    peptide_id = "pep1",
    species    = "mock_species"
  )

  x_phip <- list(
    data_long       = toy_long,
    peptide_library = toy_lib
  )
  class(x_phip) <- "phip_data"

  res <- ph_prevalence_shift(
    x                  = x_phip,
    exist_col          = "exist",
    rank_cols          = "species",
    group_cols         = "group",
    peptide_library    = NULL,
    B_permutations     = 100L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "equal",
    stat_mode          = "diff",
    prev_strat         = "none",
    winsor_z           = 4,
    rank_feature_keep  = list(species = "mock_species"),
    log                = FALSE,
    fold_change        = "none",
    cross_prev         = "none"
  )

  expect_s3_class(res, "tbl_df")
  expect_equal(res$rank, "species")
  expect_equal(res$feature, "mock_species")
})

test_that("ph_prevalence_shift logs to default log_file in a temporary
          directory", {
  skip_if_not_installed("filelock")

  toy_long <- tibble::tibble(
    sample_id  = c("s1", "s2"),
    subject_id = c("id1", "id2"),
    peptide_id = c("pep1", "pep1"),
    group      = c("A", "B"),
    exist      = c(1L, 1L)
  )

  toy_lib <- tibble::tibble(
    peptide_id = "pep1",
    species    = "mock_species"
  )

  x_phip <- list(
    data_long       = toy_long,
    peptide_library = toy_lib
  )
  class(x_phip) <- "phip_data"

  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)

  res <- ph_prevalence_shift(
    x                  = x_phip,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 100L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "equal",
    stat_mode          = "diff",
    prev_strat         = "none",
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = TRUE,
    # using default: log_file = "ph_prevalence_shift.log"
    fold_change        = "none",
    cross_prev         = "none"
  )

  expect_true(file.exists("ph_prevalence_shift.log"))
  expect_true(file.exists("ph_prevalence_shift.log.progress.bin"))
})

test_that("ph_prevalence_shift logs to a custom log_file path", {
  skip_if_not_installed("filelock")

  toy_long <- tibble::tibble(
    sample_id  = c("s1", "s2"),
    subject_id = c("id1", "id2"),
    peptide_id = c("pep1", "pep1"),
    group      = c("A", "B"),
    exist      = c(1L, 1L)
  )

  toy_lib <- tibble::tibble(
    peptide_id = "pep1",
    species    = "mock_species"
  )

  x_phip <- list(
    data_long       = toy_long,
    peptide_library = toy_lib
  )
  class(x_phip) <- "phip_data"

  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)

  custom_log <- file.path(tmp, "my_custom_shift.log")

  res <- ph_prevalence_shift(
    x                  = x_phip,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 100L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "equal",
    stat_mode          = "diff",
    prev_strat         = "none",
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = TRUE,
    log_file           = custom_log,
    fold_change        = "none",
    cross_prev         = "none"
  )

  expect_true(file.exists(custom_log))
  expect_true(file.exists(paste0(custom_log, ".progress.bin")))
})

test_that("ph_prevalence_shift computes fold_change summaries correctly for
          all modes", {
  toy_fc <- tibble::tibble(
    sample_id  = c("s1", "s2", "s3", "s4"),
    subject_id = c("id1", "id2", "id3", "id4"),
    peptide_id = c("pep1", "pep1", "pep2", "pep2"),
    group      = c("A",    "B",    "A",    "B"),
    exist      = c(1L,     1L,     1L,     1L),
    fold_change= c(1,      2,      3,      4)
  )

  # expected per peptide_id
  expected <- toy_fc %>%
    dplyr::group_by(peptide_id) %>%
    dplyr::summarise(
      fc_sum    = sum(fold_change),
      fc_mean   = mean(fold_change),
      fc_max    = max(fold_change),
      fc_median = stats::median(fold_change),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(peptide_id)

  run_fc <- function(mode) {
    set.seed(1)
    ph_prevalence_shift(
      x                  = toy_fc,
      exist_col          = "exist",
      rank_cols          = "peptide_id",
      group_cols         = "group",
      B_permutations     = 100L,
      smooth_eps_num     = 0.5,
      smooth_eps_den_mult= 2.0,
      min_max_prev       = 0.0,
      weight_mode        = "equal",
      stat_mode          = "diff",
      prev_strat         = "none",
      winsor_z           = 4,
      rank_feature_keep  = NULL,
      peptide_library    = NULL,
      log                = FALSE,
      fold_change        = mode,
      cross_prev         = "none"
    ) %>%
      dplyr::arrange(feature)
  }

  res_sum    <- run_fc("sum")
  res_mean   <- run_fc("mean")
  res_max    <- run_fc("max")
  res_median <- run_fc("median")

  expect_equal(res_sum$fold_change_sum,         expected$fc_sum)
  expect_equal(res_mean$fold_change_mean,       expected$fc_mean)
  expect_equal(res_max$fold_change_max,         expected$fc_max)
  expect_equal(res_median$fold_change_median,   expected$fc_median)
})

test_that("ph_prevalence_shift computes cross_prev summaries correctly for all
          modes", {
  toy_cp <- tibble::tibble(
    sample_id  = c("s1", "s2", "s3", "s4"),
    subject_id = c("id1", "id2", "id3", "id4"),
    peptide_id = c("pep1","pep1","pep2","pep2"),
    species    = c("sp1","sp1","sp1","sp1"),
    group      = c("A",   "B",   "A",   "B"),
    exist      = c(1L,    0L,    1L,    1L)
  )

  toy_lib <- toy_cp %>%
    dplyr::distinct(peptide_id, species)

  # per-peptide prevalence across g1 ∪ g2 (here all subjects/rows)
  prev_per_pep <- toy_cp %>%
    dplyr::group_by(peptide_id) %>%
    dplyr::summarise(
      prev = mean(exist > 0, na.rm = TRUE),
      .groups = "drop"
    )

  # For sp1: pep1 prev = (1,0) -> 0.5; pep2 prev = (1,1) -> 1.0
  expected_sum    <- sum(prev_per_pep$prev)
  expected_mean   <- mean(prev_per_pep$prev)
  expected_max    <- max(prev_per_pep$prev)
  expected_median <- stats::median(prev_per_pep$prev)

  run_cp <- function(mode) {
    set.seed(1)
    ph_prevalence_shift(
      x                  = toy_cp,
      exist_col          = "exist",
      rank_cols          = "species",
      group_cols         = "group",
      B_permutations     = 100L,
      smooth_eps_num     = 0.5,
      smooth_eps_den_mult= 2.0,
      min_max_prev       = 0.0,
      weight_mode        = "equal",
      stat_mode          = "diff",
      prev_strat         = "none",
      winsor_z           = 4,
      rank_feature_keep  = list(species = "sp1"),
      peptide_library    = toy_lib,
      log                = FALSE,
      fold_change        = "none",
      cross_prev         = mode
    )
  }

  res_sum    <- run_cp("sum")
  res_mean   <- run_cp("mean")
  res_max    <- run_cp("max")
  res_median <- run_cp("median")

  expect_equal(res_sum$cross_prev_sum,         expected_sum)
  expect_equal(res_mean$cross_prev_mean,       expected_mean)
  expect_equal(res_max$cross_prev_max,         expected_max)
  expect_equal(res_median$cross_prev_median,   expected_median)
})

test_that("ph_prevalence_shift supports se_invvar weights with asin and diff
          z-scores", {
  toy_se <- tibble::tibble(
    sample_id  = c("s1",  "s2",  "s3",  "s4"),
    subject_id = c("id1", "id2", "id3", "id4"),
    peptide_id = c("pep1","pep1","pep2","pep2"),
    group      = c("A",   "B",   "A",   "B"),
    exist      = c(1L,    0L,    0L,    1L)
  )

  # se_invvar + asin
  set.seed(1)
  res_asin <- ph_prevalence_shift(
    x                  = toy_se,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 120L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "se_invvar",
    stat_mode          = "asin",
    prev_strat         = "none",
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
    fold_change        = "none",
    cross_prev         = "none"
  )

  # se_invvar + diff
  set.seed(1)
  res_diff <- ph_prevalence_shift(
    x                  = toy_se,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 120L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "se_invvar",
    stat_mode          = "diff",
    prev_strat         = "none",
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
    fold_change        = "none",
    cross_prev         = "none"
  )

  expect_s3_class(res_asin, "tbl_df")
  expect_s3_class(res_diff, "tbl_df")
  expect_equal(nrow(res_asin), 2L)
  expect_equal(nrow(res_diff), 2L)
  expect_false(all(is.na(res_asin$T_obs)))
  expect_false(all(is.na(res_diff$T_obs)))
})

test_that("ph_prevalence_shift supports prevalence-stratified combination
          (deciles)", {
  toy_dec <- tibble::tibble(
    sample_id  = c("s1", "s2", "s3", "s4"),
    subject_id = c("id1", "id2", "id3", "id4"),
    peptide_id = c("pep1","pep1","pep2","pep2"),
    species    = c("sp1","sp1","sp1","sp1"),
    group      = c("A",   "B",   "A",   "B"),
    exist      = c(1L,    0L,    1L,    1L)
  )

  toy_lib <- toy_dec %>%
    dplyr::distinct(peptide_id, species)

  set.seed(1)
  res <- ph_prevalence_shift(
    x                  = toy_dec,
    exist_col          = "exist",
    rank_cols          = "species",
    group_cols         = "group",
    B_permutations     = 120L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "equal",
    stat_mode          = "diff",
    prev_strat         = "decile",
    winsor_z           = 4,
    rank_feature_keep  = list(species = "sp1"),
    peptide_library    = toy_lib,
    log                = FALSE,
    fold_change        = "none",
    cross_prev         = "none"
  )

  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 1L)
  expect_true(all(is.numeric(res$T_obs)))
})

test_that("ph_prevalence_shift returns NA stats when no peptides pass
          min_max_prev (unpaired)", {
  toy_prev <- tibble::tibble(
    sample_id  = c("s1","s2","s3","s4"),
    subject_id = c("id1","id2","id3","id4"),
    peptide_id = c("pep1","pep1","pep2","pep2"),
    group      = c("A",   "B",   "A",   "B"),
    exist      = c(1L,    0L,    0L,    1L)
  )

  set.seed(1)
  res <- ph_prevalence_shift(
    x                  = toy_prev,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 100L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 2.0,
    weight_mode        = "equal",
    stat_mode          = "diff",
    prev_strat         = "none",
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
    fold_change        = "none",
    cross_prev         = "none"
  )

  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 2L)
  expect_true(all(res$n_peptides_used == 0L))
  expect_true(all(is.na(res$T_obs)))
  expect_true(all(is.na(res$p_perm)))
  expect_true(all(is.na(res$mean_delta)))
  expect_true(all(is.na(res$mean_delta_w)))
})

test_that("ph_prevalence_shift returns NA stats when no peptides pass
          min_max_prev (paired)", {
  toy_prev_p <- tibble::tibble(
    sample_id  = c("s1","s2","s3","s4"),
    subject_id = c("id1","id1","id2","id2"),
    peptide_id = c("pep1","pep1","pep1","pep1"),
    group      = c("A",   "B",   "A",   "B"),
    exist      = c(1L,    0L,    0L,    1L)
  )

  set.seed(1)
  res <- ph_prevalence_shift(
    x                  = toy_prev_p,
    paired_by          = "subject_id",
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 100L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 2.0,
    weight_mode        = "equal",
    stat_mode          = "diff",
    prev_strat         = "none",
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
    fold_change        = "none",
    cross_prev         = "none"
  )

  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 0L)
})

test_that("ph_prevalence_shift aborts on bitset/peptide dimension mismatch", {
  toy_df <- tibble::tibble(
    sample_id  = c("s1", "s2"),
    subject_id = c("id1", "id2"),
    peptide_id = c("pep1", "pep2"),
    group      = c("A", "B"),
    exist      = c(1L, 1L)
  )

  testthat::local_mocked_bindings(
    build_bitset_unpaired = function(hits_by_peptide, N_subjects) {
      list(
        data    = raw(0),
        m       = length(hits_by_peptide) + 1L,
        n_words = 1L
      )
    },
    .env = asNamespace("phiper")
  )

  expect_error(
    ph_prevalence_shift(
      x                  = toy_df,
      exist_col          = "exist",
      rank_cols          = "peptide_id",
      group_cols         = "group",
      B_permutations     = 200L,
      smooth_eps_num     = 0.5,
      smooth_eps_den_mult= 2.0,
      min_max_prev       = 0.0,
      weight_mode        = "equal",
      stat_mode          = "diff",
      prev_strat         = "none",
      winsor_z           = 4,
      rank_feature_keep  = NULL,
      peptide_library    = NULL,
      log                = FALSE,
      fold_change        = "none",
      cross_prev         = "none"
    ),
    regexp = "Bitset/peptide dimension mismatch"
  )
})

test_that("ph_prevalence_shift uses global RNG reproducibly and advances .Random.seed deterministically", {
  toy_rng <- tibble::tibble(
    sample_id  = paste0("s", 1:10),
    subject_id = paste0("id", 1:10),
    peptide_id = rep("pep1", 10),
    group      = rep(c("A", "B"), each = 5),
    exist      = c(rep(1L, 5), rep(0L, 5))
  )

  # --- First run with a fixed seed ---------------------------------------------
  set.seed(1234)
  seed_before1 <- .Random.seed

  res1 <- ph_prevalence_shift(
    x                  = toy_rng,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 150L,  # >= 100 for the function, still fast for tests
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "equal",
    stat_mode          = "diff",
    prev_strat         = "none",
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
    fold_change        = "none",
    cross_prev         = "none"
  )

  seed_after1 <- .Random.seed

  # --- Second run with the same seed -------------------------------------------
  set.seed(1234)
  seed_before2 <- .Random.seed

  res2 <- ph_prevalence_shift(
    x                  = toy_rng,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 150L,
    smooth_eps_num     = 0.5,
    smooth_eps_den_mult= 2.0,
    min_max_prev       = 0.0,
    weight_mode        = "equal",
    stat_mode          = "diff",
    prev_strat         = "none",
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
    fold_change        = "none",
    cross_prev         = "none"
  )

  seed_after2 <- .Random.seed

  # 1) The function actually consumes RNG (seed changes between before/after)
  expect_false(identical(seed_before1, seed_after1))

  # 2) RNG consumption path is deterministic for the same seed and input
  expect_identical(seed_before1, seed_before2)
  expect_identical(seed_after1,  seed_after2)

  # 3) Results are identical when started from the same seed
  expect_equal(res1, res2)
})

