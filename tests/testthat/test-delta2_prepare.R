# tests/testthat/test-delta2_prepare.R

# bitset column j as a logical vector over the first N rows
.bitset_col <- function(prep, j, N) {
  bytes_per_col <- prep$bitset_words * 8L
  bytes <- prep$bitset_raw[(j - 1L) * bytes_per_col + seq_len(bytes_per_col)]
  as.logical(rawToBits(bytes))[seq_len(N)]
}

# subject-level bitset exactly as compute_delta() builds it
.v1_bitset <- function(df, exist_col = "exist") {
  subjects_order <- sort(unique(df$subject_id))
  peptides_order <- sort(unique(df$peptide_id))
  pos <- unique(df[df[[exist_col]] > 0L, c("subject_id", "peptide_id")])
  hits_split <- split(
    match(pos$subject_id, subjects_order),
    match(pos$peptide_id, peptides_order)
  )
  peptides_order <- peptides_order[as.integer(names(hits_split))]
  bs <- build_bitset_unpaired(
    lapply(unname(hits_split), as.integer), length(subjects_order)
  )
  list(
    subjects_order = subjects_order, peptides_order = peptides_order,
    bitset_raw = bs$data, bitset_m = bs$m, bitset_words = bs$n_words
  )
}

# 6 subjects, 11 samples, between-subject group, two timepoints
.toy_delta2 <- function() {
  samples <- data.frame(
    sample_id = c(
      "s1_T1", "s1_T2", "s2_T1", "s2_T1b", "s2_T2", "s3_T1",
      "s4_T1", "s4_T2", "s5_T2", "s6_T1", "s6_T2"
    ),
    subject_id = c(
      "s1", "s1", "s2", "s2", "s2", "s3", "s4", "s4", "s5", "s6", "s6"
    ),
    timepoint = c(
      "T1", "T2", "T1", "T1", "T2", "T1", "T1", "T2", "T2", "T1", "T2"
    ),
    group = c("A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  pos <- list(
    p1 = c("s1_T1", "s1_T2", "s2_T1", "s2_T1b", "s4_T2", "s6_T1"),
    p2 = c("s3_T1", "s5_T2", "s6_T2"),
    p3 = "s2_T2",
    p0 = character()
  )
  df <- merge(samples, data.frame(peptide_id = names(pos)), by = NULL)
  df$exist <- as.integer(mapply(
    function(s, p) s %in% pos[[p]], df$sample_id, df$peptide_id
  ))
  df
}

test_that("one sample per subject reproduces the v1 bitset bit-for-bit", {
  fx <- readRDS(test_path("fixtures", "delta_regression.rds"))
  df <- fx$data[fx$data$timepoint == "T1", ]

  prep <- .ph_delta2_prepare(df, group_cols = "group", time_col = "timepoint")
  v1 <- .v1_bitset(df)

  expect_identical(prep$bitset_raw, v1$bitset_raw)
  expect_identical(prep$bitset_m, v1$bitset_m)
  expect_identical(prep$bitset_words, v1$bitset_words)
  expect_identical(prep$peptides_order, v1$peptides_order)
  expect_identical(prep$subjects_order, v1$subjects_order)
  expect_identical(prep$sample_subject, seq_along(v1$subjects_order))
  expect_true(all(prep$subject_n == 1L))
})

test_that("several samples per subject give hand-computed counts", {
  prep <- .ph_delta2_prepare(
    .toy_delta2(),
    group_cols = "group", time_col = "timepoint"
  )
  N <- length(prep$samples_order)

  expect_identical(prep$samples_order, c(
    "s1_T1", "s1_T2", "s2_T1", "s2_T1b", "s2_T2", "s3_T1",
    "s4_T1", "s4_T2", "s5_T2", "s6_T1", "s6_T2"
  ))
  expect_identical(
    prep$sample_subject,
    c(1L, 1L, 2L, 2L, 2L, 3L, 4L, 4L, 5L, 6L, 6L)
  )
  expect_identical(
    prep$sample_time,
    c(1L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 2L)
  )
  expect_identical(prep$subject_n, c(2L, 3L, 1L, 2L, 1L, 2L))
  expect_identical(prep$time_levels, c("T1", "T2"))
  expect_identical(prep$subject_group$group, c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_identical(prep$sample_group$group, rep(1:2, c(6L, 5L)))

  # p0 has no positives and is dropped, as in compute_delta()
  expect_identical(prep$peptides_order, c("p1", "p2", "p3"))
  expect_identical(prep$bitset_m, 3L)

  # positive samples per peptide (rows) and timepoint (cols)
  counts <- t(vapply(seq_len(prep$bitset_m), function(j) {
    bits <- .bitset_col(prep, j, N)
    as.integer(tabulate(prep$sample_time[bits], nbins = 2L))
  }, integer(2)))
  expect_identical(counts, rbind(c(4L, 2L), c(1L, 2L), c(0L, 1L)))
})

test_that("a subject with a single sample is kept in the bitset", {
  prep <- .ph_delta2_prepare(
    .toy_delta2(),
    group_cols = "group", time_col = "timepoint"
  )
  N <- length(prep$samples_order)

  s3 <- match("s3", prep$subjects_order)
  expect_identical(prep$subject_n[s3], 1L)
  row <- which(prep$sample_subject == s3)
  expect_identical(prep$samples_order[row], "s3_T1")
  expect_true(.bitset_col(prep, match("p2", prep$peptides_order), N)[row])
})

test_that("within-subject group columns are NA at subject level", {
  prep <- .ph_delta2_prepare(
    .toy_delta2(),
    group_cols = c("group", "timepoint")
  )

  expect_identical(prep$time_levels, 1L)
  expect_true(all(prep$sample_time == 1L))
  # s3 and s5 have a single timepoint, the others vary within subject
  expect_identical(
    prep$subject_group$timepoint,
    c(NA, NA, 1L, NA, 2L, NA)
  )
})

test_that("a sample_id mapped to two subjects is rejected", {
  df <- .toy_delta2()
  df$subject_id[df$sample_id == "s1_T2" & df$peptide_id == "p1"] <- "s2"
  expect_error(
    .ph_delta2_prepare(df, group_cols = "group", time_col = "timepoint"),
    "single subject_id"
  )
})
