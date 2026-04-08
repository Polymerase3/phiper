# tests/testthat/test-pop_compute.R

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.make_pop_df <- function(n_per_group = 5L, n_peptides = 3L, seed = 1L) {
  set.seed(seed)
  n <- n_per_group * 2L * n_peptides
  tibble::tibble(
    sample_id  = rep(paste0("s", seq_len(n_per_group * 2L)), each = n_peptides),
    peptide_id = rep(paste0("pep", seq_len(n_peptides)), times = n_per_group * 2L),
    group      = rep(c("A", "B"), each = n_per_group * n_peptides),
    exist      = sample(0L:1L, n, replace = TRUE)
  )
}

# ---------------------------------------------------------------------------
# 1. Basic unpaired – data.frame, rank_cols = "peptide_id"
# ---------------------------------------------------------------------------

test_that("compute_pop returns correct columns for unpaired peptide-level", {
  df <- .make_pop_df()

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist"
  )

  expect_s3_class(res, "data.frame")
  expected_cols <- c(
    "rank", "feature", "group_col",
    "group1", "n1", "N1", "prop1", "percent1",
    "group2", "n2", "N2", "prop2", "percent2",
    "ratio", "delta_ratio", "p_raw"
  )
  expect_true(all(expected_cols %in% names(res)))
})

test_that("compute_pop: unpaired has one row per (rank, feature)", {
  df <- .make_pop_df(n_per_group = 4L, n_peptides = 4L)

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist"
  )

  expect_true(nrow(res) >= 1L)
  expect_true(all(res$rank == "peptide_id"))
  expect_equal(res$group_col, rep("group", nrow(res)))
})

test_that("compute_pop: groups are labelled correctly (group1 <= group2)", {
  df <- .make_pop_df()

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist"
  )

  expect_true(all(res$group1 <= res$group2))
  expect_equal(unique(res$group1), "A")
  expect_equal(unique(res$group2), "B")
})

# ---------------------------------------------------------------------------
# 2. Fisher p-value sanity
# ---------------------------------------------------------------------------

test_that("compute_pop Fisher p-value is in (0,1]", {
  df <- .make_pop_df(n_per_group = 10L, n_peptides = 2L, seed = 42L)

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist"
  )

  finite_p <- res$p_raw[!is.na(res$p_raw)]
  expect_true(all(finite_p > 0))
  expect_true(all(finite_p <= 1))
})

test_that("compute_pop detects strong signal (p << 0.05 when groups perfectly differ)", {
  n <- 20L
  df <- tibble::tibble(
    sample_id  = paste0("s", seq_len(2L * n)),
    peptide_id = rep("pep1", 2L * n),
    group      = rep(c("A", "B"), each = n),
    exist      = c(rep(1L, n), rep(0L, n))
  )

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist"
  )

  expect_equal(nrow(res), 1L)
  expect_lt(res$p_raw, 0.001)
  expect_equal(res$n1, n)
  expect_equal(res$n2, 0L)
})

test_that("compute_pop returns p = 1 when prevalence is identical in both groups", {
  n <- 10L
  df <- tibble::tibble(
    sample_id  = paste0("s", seq_len(2L * n)),
    peptide_id = rep("pep1", 2L * n),
    group      = rep(c("A", "B"), each = n),
    exist      = 1L
  )

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist"
  )

  expect_equal(nrow(res), 1L)
  expect_equal(res$p_raw, 1)
})

# ---------------------------------------------------------------------------
# 3. Non-peptide rank via peptide_library
# ---------------------------------------------------------------------------

test_that("compute_pop works with species rank via peptide_library", {
  df <- tibble::tibble(
    sample_id  = c("s1","s1","s2","s2","s3","s3","s4","s4"),
    peptide_id = rep(c("pep1","pep2"), 4L),
    group      = rep(c("A","B"), each = 4L),
    exist      = c(1L,0L, 1L,1L, 0L,1L, 0L,0L)
  )
  lib <- tibble::tibble(
    peptide_id = c("pep1","pep2"),
    species    = c("sp_alpha","sp_alpha")
  )

  res <- compute_pop(
    x               = df,
    rank_cols       = "species",
    group_cols      = "group",
    exist_col       = "exist",
    peptide_library = lib
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_equal(res$rank, "species")
  expect_equal(res$feature, "sp_alpha")
})

test_that("compute_pop errors when non-peptide rank requested but no peptide_library", {
  df <- tibble::tibble(
    sample_id  = c("s1","s2"),
    peptide_id = c("pep1","pep1"),
    group      = c("A","B"),
    exist      = c(1L,0L)
  )

  expect_error(
    compute_pop(
      x          = df,
      rank_cols  = "species",
      group_cols = "group",
      exist_col  = "exist"
    ),
    regexp = "peptide_library"
  )
})

test_that("compute_pop errors when peptide_library is missing required rank column", {
  df <- tibble::tibble(
    sample_id  = c("s1","s2"),
    peptide_id = c("pep1","pep1"),
    group      = c("A","B"),
    exist      = c(1L,0L)
  )
  bad_lib <- tibble::tibble(peptide_id = "pep1")

  expect_error(
    compute_pop(
      x               = df,
      rank_cols       = "species",
      group_cols      = "group",
      exist_col       = "exist",
      peptide_library = bad_lib
    ),
    regexp = "species"
  )
})

# ---------------------------------------------------------------------------
# 4. pop_k_min threshold
# ---------------------------------------------------------------------------

test_that("compute_pop: pop_k_min = 2 requires at least 2 positive peptides", {
  # 3 peptides per sample; only samples with >= 2 positives are "present"
  df <- tibble::tibble(
    sample_id  = c(rep("s1", 3L), rep("s2", 3L), rep("s3", 3L), rep("s4", 3L)),
    peptide_id = rep(c("pep1","pep2","pep3"), 4L),
    group      = rep(c("A","B"), each = 6L),
    exist      = c(1L,1L,0L,  # s1: 2 positives -> present
                   1L,0L,0L,  # s2: 1 positive  -> absent with k=2
                   1L,1L,1L,  # s3: 3 positives -> present
                   0L,0L,0L)  # s4: 0 positives -> absent
  )

  res_k1 <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist",
    pop_k_min  = 1L
  )

  res_k2 <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist",
    pop_k_min  = 2L
  )

  # k=2 should have fewer or equal present counts per feature
  expect_true(all(res_k2$n1 + res_k2$n2 <= res_k1$n1 + res_k1$n2))
})

# ---------------------------------------------------------------------------
# 5. Multiple rank_cols and multiple group_cols
# ---------------------------------------------------------------------------

test_that("compute_pop handles multiple rank_cols", {
  df <- tibble::tibble(
    sample_id  = c("s1","s1","s2","s2","s3","s3","s4","s4"),
    peptide_id = rep(c("pep1","pep2"), 4L),
    group      = rep(c("A","B"), each = 4L),
    exist      = c(1L,0L, 0L,1L, 1L,1L, 0L,0L)
  )
  lib <- tibble::tibble(
    peptide_id = c("pep1","pep2"),
    species    = c("sp1","sp2")
  )

  res <- compute_pop(
    x               = df,
    rank_cols       = c("peptide_id", "species"),
    group_cols      = "group",
    exist_col       = "exist",
    peptide_library = lib
  )

  expect_true(all(c("peptide_id", "species") %in% res$rank))
})

test_that("compute_pop handles multiple group_cols", {
  df <- tibble::tibble(
    sample_id  = paste0("s", 1:8),
    peptide_id = rep("pep1", 8L),
    grp1       = rep(c("A","B"), each = 4L),
    grp2       = rep(c("X","Y"), 4L),
    exist      = c(1L,1L,0L,0L, 0L,0L,1L,1L)
  )

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = c("grp1","grp2"),
    exist_col  = "exist"
  )

  expect_true(all(c("grp1","grp2") %in% res$group_col))
})

# ---------------------------------------------------------------------------
# 6. Argument validation errors
# ---------------------------------------------------------------------------

test_that("compute_pop errors when group has more than 2 levels", {
  df <- tibble::tibble(
    sample_id  = paste0("s", 1:3),
    peptide_id = rep("pep1", 3L),
    group      = c("A","B","C"),
    exist      = 1L
  )

  expect_error(
    compute_pop(
      x          = df,
      rank_cols  = "peptide_id",
      group_cols = "group",
      exist_col  = "exist"
    ),
    regexp = "binary"
  )
})

test_that("compute_pop errors when group has only 1 level", {
  df <- tibble::tibble(
    sample_id  = paste0("s", 1:3),
    peptide_id = rep("pep1", 3L),
    group      = rep("A", 3L),
    exist      = 1L
  )

  expect_error(
    compute_pop(
      x          = df,
      rank_cols  = "peptide_id",
      group_cols = "group",
      exist_col  = "exist"
    ),
    regexp = "binary"
  )
})

test_that("compute_pop errors on invalid rank_cols argument type", {
  df <- .make_pop_df()
  expect_error(
    compute_pop(x = df, rank_cols = 123, group_cols = "group"),
    regexp = NULL
  )
})

test_that("compute_pop errors on invalid pop_k_min < 1", {
  df <- .make_pop_df()
  expect_error(
    compute_pop(x = df, rank_cols = "peptide_id", group_cols = "group", pop_k_min = 0L),
    regexp = NULL
  )
})

# ---------------------------------------------------------------------------
# 7. Paired design – McNemar test
# ---------------------------------------------------------------------------

test_that("compute_pop paired design returns p_raw without ratio/delta_ratio columns", {
  df <- tibble::tibble(
    sample_id  = c("s1a","s1b","s2a","s2b","s3a","s3b","s4a","s4b"),
    subject_id = rep(paste0("sub", 1:4), each = 2L),
    peptide_id = rep("pep1", 8L),
    timepoint  = rep(c("T1","T2"), 4L),
    exist      = c(1L,0L, 0L,1L, 1L,1L, 0L,0L)
  )

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "timepoint",
    exist_col  = "exist",
    paired     = "subject_id"
  )

  expect_s3_class(res, "data.frame")
  expect_true("p_raw" %in% names(res))
  # paired path does NOT return ratio or delta_ratio
  expect_false("ratio" %in% names(res))
  expect_false("delta_ratio" %in% names(res))
})

test_that("compute_pop paired: p_raw is in (0,1]", {
  df <- tibble::tibble(
    sample_id  = c("s1a","s1b","s2a","s2b","s3a","s3b","s4a","s4b",
                   "s5a","s5b","s6a","s6b"),
    subject_id = rep(paste0("sub", 1:6), each = 2L),
    peptide_id = rep("pep1", 12L),
    timepoint  = rep(c("T1","T2"), 6L),
    exist      = c(1L,0L, 1L,0L, 1L,0L, 0L,1L, 0L,1L, 1L,1L)
  )

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "timepoint",
    exist_col  = "exist",
    paired     = "subject_id"
  )

  expect_true(all(res$p_raw[!is.na(res$p_raw)] > 0))
  expect_true(all(res$p_raw[!is.na(res$p_raw)] <= 1))
})

test_that("compute_pop paired: N1 == N2 (same paired count)", {
  df <- tibble::tibble(
    sample_id  = c("s1a","s1b","s2a","s2b","s3a","s3b","s4a","s4b"),
    subject_id = rep(paste0("sub", 1:4), each = 2L),
    peptide_id = rep("pep1", 8L),
    timepoint  = rep(c("T1","T2"), 4L),
    exist      = c(1L,0L, 0L,1L, 1L,1L, 0L,0L)
  )

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "timepoint",
    exist_col  = "exist",
    paired     = "subject_id"
  )

  expect_equal(res$N1, res$N2)
})

test_that("compute_pop errors when paired column is missing from data", {
  df <- tibble::tibble(
    sample_id  = c("s1","s2"),
    peptide_id = rep("pep1", 2L),
    group      = c("A","B"),
    exist      = c(1L,0L)
  )

  expect_error(
    compute_pop(
      x          = df,
      rank_cols  = "peptide_id",
      group_cols = "group",
      exist_col  = "exist",
      paired     = "subject_id"
    ),
    regexp = "not found"
  )
})

test_that("compute_pop errors when paired is not FALSE or a string", {
  df <- .make_pop_df()
  expect_error(
    compute_pop(
      x          = df,
      rank_cols  = "peptide_id",
      group_cols = "group",
      paired     = 42L
    ),
    regexp = NULL
  )
})

# ---------------------------------------------------------------------------
# 8. phip_data input
# ---------------------------------------------------------------------------

test_that("compute_pop works with phip_data-like list input", {
  long_df <- tibble::tibble(
    sample_id  = c("s1","s1","s2","s2","s3","s3","s4","s4"),
    peptide_id = rep(c("pep1","pep2"), 4L),
    group      = rep(c("A","B"), each = 4L),
    exist      = c(1L,0L, 0L,1L, 1L,0L, 0L,1L)
  )
  lib <- tibble::tibble(
    peptide_id = c("pep1","pep2"),
    species    = c("sp1","sp2")
  )

  x_phip <- list(data_long = long_df, peptide_library = lib)
  class(x_phip) <- "phip_data"

  res <- compute_pop(
    x          = x_phip,
    rank_cols  = "species",
    group_cols = "group",
    exist_col  = "exist"
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(res$rank == "species"))
})

test_that("compute_pop: phip_data peptide_library used when peptide_library arg is NULL", {
  long_df <- tibble::tibble(
    sample_id  = c("s1","s2","s3","s4"),
    peptide_id = rep("pep1", 4L),
    group      = c("A","A","B","B"),
    exist      = c(1L,0L,0L,1L)
  )
  lib <- tibble::tibble(peptide_id = "pep1", genus = "Genus_X")

  x_phip <- list(data_long = long_df, peptide_library = lib)
  class(x_phip) <- "phip_data"

  res <- compute_pop(
    x               = x_phip,
    rank_cols       = "genus",
    group_cols      = "group",
    exist_col       = "exist",
    peptide_library = NULL
  )

  expect_equal(res$rank[1], "genus")
  expect_equal(res$feature[1], "Genus_X")
})

# ---------------------------------------------------------------------------
# 9. Output proportions and percent are consistent
# ---------------------------------------------------------------------------

test_that("compute_pop: percent = 100 * prop and prop = n / N", {
  df <- tibble::tibble(
    sample_id  = paste0("s", 1:8),
    peptide_id = rep("pep1", 8L),
    group      = rep(c("A","B"), each = 4L),
    exist      = c(1L,1L,0L,0L, 1L,0L,0L,0L)
  )

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist"
  )

  expect_equal(nrow(res), 1L)
  expect_equal(res$prop1, res$n1 / res$N1)
  expect_equal(res$percent1, 100 * res$prop1)
  expect_equal(res$prop2, res$n2 / res$N2)
  expect_equal(res$percent2, 100 * res$prop2)
})

# ---------------------------------------------------------------------------
# 10. Cohort sizes (N1, N2)
# ---------------------------------------------------------------------------

test_that("compute_pop: N1 and N2 reflect correct group sizes", {
  n_a <- 6L; n_b <- 4L
  df <- tibble::tibble(
    sample_id  = paste0("s", seq_len(n_a + n_b)),
    peptide_id = rep("pep1", n_a + n_b),
    group      = c(rep("A", n_a), rep("B", n_b)),
    exist      = 1L
  )

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist"
  )

  expect_equal(res$N1, n_a)
  expect_equal(res$N2, n_b)
})

# ---------------------------------------------------------------------------
# 11. NA group values are ignored
# ---------------------------------------------------------------------------

test_that("compute_pop ignores NA group values", {
  df <- tibble::tibble(
    sample_id  = paste0("s", 1:6),
    peptide_id = rep("pep1", 6L),
    group      = c("A","A","B","B", NA_character_, NA_character_),
    exist      = 1L
  )

  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist"
  )

  expect_equal(res$N1 + res$N2, 4L)
})

# ---------------------------------------------------------------------------
# 12. exist_col coercion of NA to 0
# ---------------------------------------------------------------------------

test_that("compute_pop coerces NA exist to 0", {
  df <- tibble::tibble(
    sample_id  = c("s1","s2","s3","s4"),
    peptide_id = rep("pep1", 4L),
    group      = c("A","A","B","B"),
    exist      = c(NA_integer_, 1L, 1L, 0L)
  )

  # should not error; NA treated as absent
  res <- compute_pop(
    x          = df,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist"
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  # group A: NA->0, 1 => n_present = 1
  expect_equal(res$n1, 1L)
})

# ---------------------------------------------------------------------------
# 13. view column prepended for phip_data with view attribute
# ---------------------------------------------------------------------------

test_that("compute_pop prepends view column when phip_data has view attribute", {
  long_df <- tibble::tibble(
    sample_id  = c("s1","s2","s3","s4"),
    peptide_id = rep("pep1", 4L),
    group      = c("A","A","B","B"),
    exist      = c(1L,0L,0L,1L)
  )

  x_phip <- list(
    data_long = long_df,
    meta      = list(view = "my_view")
  )
  class(x_phip) <- "phip_data"

  res <- compute_pop(
    x          = x_phip,
    rank_cols  = "peptide_id",
    group_cols = "group",
    exist_col  = "exist"
  )

  expect_true("view" %in% names(res))
  expect_equal(names(res)[1], "view")
  expect_true(all(res$view == "my_view"))
})
