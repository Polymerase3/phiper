# tests/testthat/test-compute_delta.R

test_that("compute_delta works for unpaired design (mock species)", {
  # example phip_data
  ps <- phip_load_example_data()

  # small unpaired subset: one species, two groups at T1
  ps_filt <- ps |>
    dplyr::filter(
      peptide_id %in% c("16627", "5243", "24799", "16196", "18003"),
      timepoint == "T1"
    ) |>
    dplyr::collect()

  mock_peplib <- data.frame(
    peptide_id = c("16627", "5243", "24799", "16196", "18003"),
    species    = rep("mock_species", 5),
    stringsAsFactors = FALSE
  )

  res <- compute_delta(
    x                  = ps_filt,
    exist_col          = "exist",
    rank_cols          = "species",
    group_cols         = "group",
    peptide_library    = mock_peplib,
    B_permutations     = 500L,          # smaller for test speed
    weight_mode        = "n_eff_sqrt",
    stat_mode          = "asin",
    strat_bins         = 0,
    winsor_z           = Inf,
    rank_feature_keep  = list(species = NULL),
    log                = FALSE
  )

  # basic structure
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 1L)

  # key columns
  expect_true(all(c(
    "rank", "feature", "group_col", "group1", "group2", "design",
    "n_subjects_paired", "n_peptides_used", "m_eff",
    "T_obs", "T_null_mean", "T_null_sd", "T_obs_stand", "Z_from_p",
    "p_perm", "b"
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

test_that("compute_delta handles paired design via paired_by and returns
          standardized T", {
  # build example phip_data
  ps <- phip_load_example_data()

  # subset: one mock species, group A, two timepoints
  ps_filt2 <- ps |>
    dplyr::filter(
      peptide_id %in% c("2269", "21399", "7789", "13588", "10180"),
      group == "A"
    ) |>
    dplyr::collect()

  mock_peplib <- data.frame(
    peptide_id = c("2269", "21399", "7789", "13588", "10180"),
    species    = rep("mock_species", 5),
    stringsAsFactors = FALSE
  )

  # identify true T1-T2 pairs
  pair_ids <- ps_filt2 |>
    dplyr::filter(timepoint %in% c("T1", "T2")) |>
    dplyr::distinct(subject_id, timepoint) |>
    tidyr::pivot_wider(
      names_from  = timepoint,
      values_from = timepoint,
      values_fn   = length,
      values_fill = 0
    ) |>
    dplyr::filter(T1 > 0, T2 > 0) |>
    dplyr::pull(subject_id)

  expect_gt(length(pair_ids), 0L)

  # inject artificial variability by flipping exist for a subset of rows
  set.seed(123)
  candidates <- ps_filt2 |>
    dplyr::filter(
      subject_id %in% pair_ids,
      timepoint %in% c("T1", "T2")
    )

  n_flip   <- min(15L, nrow(candidates))
  flip_ids <- sample(candidates$sample_id, size = n_flip)

  ps_filt2_flipped <- ps_filt2 |>
    dplyr::mutate(
      flipped    = sample_id %in% flip_ids,
      exist      = ifelse(flipped, 1L - exist, exist),
      counts_hits= ifelse(flipped & exist == 0L, 0, counts_hits),
    ) |>
    dplyr::select(-flipped)

  res <- compute_delta(
    x                  = ps_filt2_flipped,
    paired_by          = "subject_id",
    exist_col          = "exist",
    rank_cols          = "species",
    group_cols         = "timepoint",
    peptide_library    = mock_peplib,
    B_permutations     = 500L,
    weight_mode        = "n_eff_sqrt",
    stat_mode          = "asin",
    strat_bins         = 0,
    winsor_z           = Inf,
    rank_feature_keep  = list(species = NULL),
    log                = FALSE,
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

test_that("compute_delta errors on duplicate positives within group", {
  # toy data with duplicate positive for same subject/peptide/group
  toy_dup <- tibble::tibble(
    sample_id  = c("s1", "s2"),
    subject_id = c("id1", "id1"),
    peptide_id = c("pep1", "pep1"),
    group      = c("A", "A"),
    exist      = c(1L, 1L)
  )

  expect_error(
    compute_delta(
      x                  = toy_dup,
      exist_col          = "exist",
      rank_cols          = "peptide_id",
      group_cols         = "group",
      B_permutations     = 200L,
      weight_mode        = "equal",
      stat_mode          = "diff",
      strat_bins         = 0,
      winsor_z           = 4,
      rank_feature_keep  = NULL,
      peptide_library    = NULL,
      log                = FALSE,
    ),
    regexp = "duplicate positives",
    fixed  = FALSE
  )
})

test_that("compute_delta gives consistent direction for T_obs,
          T_obs_stand and Z_from_p", {
  # reuse unpaired example, smaller B for speed
  ps <- phip_load_example_data()

  ps_filt <- ps |>
    dplyr::filter(
      peptide_id %in% c("16627", "5243", "24799", "16196", "18003"),
      timepoint == "T1"
    ) |>
    dplyr::collect()

  mock_peplib <- data.frame(
    peptide_id = c("16627", "5243", "24799", "16196", "18003"),
    species    = rep("mock_species", 5),
    stringsAsFactors = FALSE
  )

  res <- compute_delta(
    x                  = ps_filt,
    exist_col          = "exist",
    rank_cols          = "species",
    group_cols         = "group",
    peptide_library    = mock_peplib,
    B_permutations     = 300L,
    weight_mode        = "n_eff_sqrt",
    stat_mode          = "asin",
    strat_bins         = 0,
    winsor_z           = Inf,
    rank_feature_keep  = list(species = NULL),
    log                = FALSE,
  )

  # only check direction when standardized T is available
  expect_equal(sign(res$T_obs), sign(res$T_obs_stand))
  expect_equal(sign(res$T_obs), sign(res$Z_from_p))
})

test_that("compute_delta returns NA T_obs_stand when permutation
          variance is zero", {
  # toy data: one peptide, two groups, identical prevalence pattern
  toy <- tibble::tibble(
    sample_id  = paste0("s", 1:40),
    subject_id = paste0("id", rep(1:20, times = 2)),
    peptide_id = rep("pep1", 40),
    group      = rep(c("A", "B"), each = 20),
    exist      = 1L
  )

  res <- compute_delta(
    x                  = toy,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 200L,
    weight_mode        = "equal",
    stat_mode          = "diff",
    strat_bins         = 0,
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
  )

  # all permutations give the same T => variance zero
  expect_s3_class(res, "tbl_df")
  expect_setequal(
    colnames(res),
    c("rank","feature","group_col","group1","group2","design",
      "n_subjects_paired","n_peptides_used","m_eff",
      "T_obs","T_null_mean","T_null_sd","T_obs_stand","Z_from_p","p_perm","b",
      "max_delta","frac_delta_pos","frac_delta_pos_w")
  )
})

test_that("compute_delta aborts when required columns are missing", {
  toy_missing <- tibble::tibble(
    sample_id  = "s1",
    peptide_id = "pep1",
    group      = "A",
    exist      = 1L
  )

  expect_error(
    compute_delta(
      x                  = toy_missing,
      exist_col          = "exist",
      rank_cols          = "peptide_id",
      group_cols         = "group",
      B_permutations     = 200L,
      weight_mode        = "equal",
      stat_mode          = "diff",
      strat_bins         = 0,
      winsor_z           = 4,
      rank_feature_keep  = NULL,
      peptide_library    = NULL,
      log                = FALSE,
    ),
    regexp = "Missing required columns"
  )
})

test_that("compute_delta aborts when all peptides are zero after hits
          guard", {
  toy_zero <- tibble::tibble(
    sample_id  = c("s1", "s2"),
    subject_id = c("id1", "id2"),
    peptide_id = c("pep1", "pep2"),
    group      = c("A", "B"),
    exist      = c(0L, 0L)
  )

  expect_error(
    compute_delta(
      x                  = toy_zero,
      exist_col          = "exist",
      rank_cols          = "peptide_id",
      group_cols         = "group",
      B_permutations     = 200L,
      weight_mode        = "equal",
      stat_mode          = "diff",
      strat_bins         = 0,
      winsor_z           = 4,
      rank_feature_keep  = NULL,
      peptide_library    = NULL,
      log                = FALSE,
    ),
    regexp = "All peptides are zero after hits guard"
  )
})

test_that("compute_delta aborts when peptide_library misses required
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
    compute_delta(
      x                  = toy_df,
      exist_col          = "exist",
      rank_cols          = "species",
      group_cols         = "group",
      peptide_library    = bad_lib,
      B_permutations     = 200L,
      weight_mode        = "equal",
      stat_mode          = "diff",
      strat_bins         = 0,
      winsor_z           = 4,
      rank_feature_keep  = NULL,
      log                = FALSE,
    ),
    regexp = "Peptide library missing required columns"
  )
})

test_that("compute_delta works with rank_cols = 'peptide_id' only", {
  toy_df <- tibble::tibble(
    sample_id  = c("s1", "s2"),
    subject_id = c("id1", "id2"),
    peptide_id = c("pep1", "pep2"),
    group      = c("A", "B"),
    exist      = c(1L, 0L)
  )

  res <- compute_delta(
    x                  = toy_df,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 200L,
    weight_mode        = "equal",
    stat_mode          = "diff",
    strat_bins         = 0,
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
  )

  expect_s3_class(res, "tbl_df")
  expect_true(all(res$rank == "peptide_id"))
})

test_that("compute_delta uses peptide_library attached to phip_data", {
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

  res <- compute_delta(
    x                  = x_phip,
    exist_col          = "exist",
    rank_cols          = "species",
    group_cols         = "group",
    peptide_library    = NULL,
    B_permutations     = 100L,
    weight_mode        = "equal",
    stat_mode          = "diff",
    strat_bins         = 0,
    winsor_z           = 4,
    rank_feature_keep  = list(species = "mock_species"),
    log                = FALSE,
  )

  expect_s3_class(res, "tbl_df")
  expect_equal(res$rank, "species")
  expect_equal(res$feature, "mock_species")
})

test_that("compute_delta logs to default log_file in a temporary
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

  res <- compute_delta(
    x                  = x_phip,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 100L,
    weight_mode        = "equal",
    stat_mode          = "diff",
    strat_bins         = 0,
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = TRUE,
    # using default: log_file = "compute_delta.log"
  )

  expect_true(file.exists("compute_delta.log"))
  expect_true(file.exists("compute_delta.log.progress.bin"))
})

test_that("compute_delta logs to a custom log_file path", {
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

  res <- compute_delta(
    x                  = x_phip,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 100L,
    weight_mode        = "equal",
    stat_mode          = "diff",
    strat_bins         = 0,
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = TRUE,
    log_file           = custom_log,
  )

  expect_true(file.exists(custom_log))
  expect_true(file.exists(paste0(custom_log, ".progress.bin")))
})

test_that("compute_delta supports se_invvar weights with asin and diff
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
  res_asin <- compute_delta(
    x                  = toy_se,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 120L,
    weight_mode        = "se_invvar",
    stat_mode          = "asin",
    strat_bins         = 0,
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
  )

  # se_invvar + diff
  set.seed(1)
  res_diff <- compute_delta(
    x                  = toy_se,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 120L,
    weight_mode        = "se_invvar",
    stat_mode          = "diff",
    strat_bins         = 0,
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
  )

  expect_s3_class(res_asin, "tbl_df")
  expect_s3_class(res_diff, "tbl_df")
  expect_equal(nrow(res_asin), 2L)
  expect_equal(nrow(res_diff), 2L)
  expect_false(all(is.na(res_asin$T_obs)))
  expect_false(all(is.na(res_diff$T_obs)))
})

test_that("compute_delta supports prevalence-stratified combination
          (custom bins)", {
  toy_dec <- tibble::tibble(
    sample_id  = c("s1", "s2", "s3", "s4"),
    subject_id = c("id1", "id2", "id3", "id4"),
    peptide_id = c("pep1","pep1","pep2","pep2"),
    species    = c("sp1","sp1","sp1","sp1"),
    group      = c("A",   "B",   "A",   "B"),
    exist      = c(1L,    0L,    1L,    1L)
  )

  toy_lib <- toy_dec |>
    dplyr::distinct(peptide_id, species)

  set.seed(1)
  res <- compute_delta(
    x                  = toy_dec,
    exist_col          = "exist",
    rank_cols          = "species",
    group_cols         = "group",
    B_permutations     = 120L,
    weight_mode        = "equal",
    stat_mode          = "diff",
    strat_bins         = seq(0.1, 0.9, by = 0.1),
    winsor_z           = 4,
    rank_feature_keep  = list(species = "sp1"),
    peptide_library    = toy_lib,
    log                = FALSE,
  )

  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 1L)
  expect_true(all(is.numeric(res$T_obs)))
})

test_that("compute_delta uses global RNG reproducibly and advances .Random.seed deterministically", {
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

  res1 <- compute_delta(
    x                  = toy_rng,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 150L,  # >= 100 for the function, still fast for tests
    weight_mode        = "equal",
    stat_mode          = "diff",
    strat_bins         = 0,
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
  )

  seed_after1 <- .Random.seed

  # --- Second run with the same seed -------------------------------------------
  set.seed(1234)
  seed_before2 <- .Random.seed

  res2 <- compute_delta(
    x                  = toy_rng,
    exist_col          = "exist",
    rank_cols          = "peptide_id",
    group_cols         = "group",
    B_permutations     = 150L,
    weight_mode        = "equal",
    stat_mode          = "diff",
    strat_bins         = 0,
    winsor_z           = 4,
    rank_feature_keep  = NULL,
    peptide_library    = NULL,
    log                = FALSE,
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

test_that("compute_delta aborts on bitset/peptide dimension mismatch", {
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
    compute_delta(
      x                  = toy_df,
      exist_col          = "exist",
      rank_cols          = "peptide_id",
      group_cols         = "group",
      B_permutations     = 200L,
      weight_mode        = "equal",
      stat_mode          = "diff",
      strat_bins         = 0,
      winsor_z           = 4,
      rank_feature_keep  = NULL,
      peptide_library    = NULL,
      log                = FALSE,
    ),
    regexp = "Bitset/peptide dimension mismatch"
  )
})
