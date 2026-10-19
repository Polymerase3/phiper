# ==============================================================================
# compute_delta2(): sample-level presence structure (phase 1)
# ==============================================================================
# compute_delta() collapses a subject's samples into one bitset row
# (distinct(subject_id, peptide_id)), so positivity is an OR over all of a
# subject's samples and subjects with more samples get more chances to look
# positive. compute_delta2() indexes samples instead. This file only builds the
# data structure; it shares build_bitset_unpaired() with compute_delta() and
# leaves compute_delta() itself untouched.

#' @title Internal helper: .ph_delta2_prepare
#' @description Build the sample-level presence bitset and the index vectors
#'   that map samples to subjects, timepoints and group levels.
#'
#' @param x A `phip_data` object or a data frame in long format with columns
#'   `sample_id`, `subject_id`, `peptide_id`, `exist_col`, `group_cols` and,
#'   if given, `time_col`.
#' @param group_cols Character vector of grouping columns.
#' @param exist_col Name of the 0/1 presence column. Default `"exist"`.
#' @param time_col Optional name of the timepoint column. If `NULL`, every
#'   sample gets the single timepoint level 1.
#'
#' @return A list with
#'   - `samples_order`: sample IDs ordered by `(subject_id, time, sample_id)`;
#'     the bitset rows.
#'   - `sample_subject`: 1-based index into `subjects_order`, per sample.
#'   - `sample_time`: 1-based index into `time_levels`, per sample.
#'   - `sample_group`: named list (one entry per `group_cols` element) of
#'     1-based indices into `group_levels[[gc]]`, per sample.
#'   - `subject_group`: named list as `sample_group`, per subject; `NA` where
#'     the group value varies within the subject (e.g. a timepoint column).
#'   - `subject_n`: number of samples per subject.
#'   - `subjects_order`, `time_levels`, `group_levels`: the level vectors the
#'     indices above point into.
#'   - `peptides_order`: peptide IDs with at least one positive sample; the
#'     bitset columns.
#'   - `bitset_raw`, `bitset_m`, `bitset_words`: packed 64-bit bitset of
#'     `bitset_m` peptide columns with `bitset_words` words each, as returned
#'     by [build_bitset_unpaired()].
#' @keywords internal
.ph_delta2_prepare <- function(x, group_cols, exist_col = "exist",
                               time_col = NULL) {
  chk::chk_character(group_cols)
  chk::chk_true(length(group_cols) >= 1)
  chk::chk_string(exist_col)
  if (!is.null(time_col)) chk::chk_string(time_col)

  sample_cols <- unique(c("sample_id", "subject_id", time_col, group_cols))
  need_cols <- unique(c(sample_cols, "peptide_id", exist_col))

  if (inherits(x, "phip_data")) {
    df_long <- x$data_long |>
      dplyr::select(dplyr::any_of(need_cols))
  } else {
    chk::chk_data(x)
    miss <- setdiff(need_cols, colnames(x))
    if (length(miss)) {
      .ph_abort("Missing required columns",
        bullets = paste("-", miss)
      )
    }
    df_long <- tibble::as_tibble(x) |>
      dplyr::select(dplyr::any_of(need_cols))
  }

  # --- SAMPLE TABLE: one row per sample ---------------------------------------
  # No hits guard as in compute_delta(): rows are samples, so several positive
  # samples of one subject are kept apart instead of rejected.
  samples <- df_long |>
    dplyr::distinct(!!!rlang::syms(sample_cols)) |>
    dplyr::collect()

  if (anyNA(samples$sample_id) || anyNA(samples$subject_id)) {
    .ph_abort("sample_id and subject_id must not contain NA.")
  }
  dup_samples <- unique(samples$sample_id[duplicated(samples$sample_id)])
  if (length(dup_samples)) {
    .ph_abort(
      sprintf(
        "Each sample_id must have a single subject_id, %s value.",
        paste(c(time_col, group_cols), collapse = ", ")
      ),
      bullets = paste("- sample_id =", utils::head(dup_samples, 10))
    )
  }

  # --- FREEZE ORDERS ----------------------------------------------------------
  # Same arrange() as compute_delta(), so subjects (and, with one sample per
  # subject, samples) come out in the same order.
  subjects_order <- samples |>
    dplyr::distinct(subject_id) |>
    dplyr::arrange(subject_id) |>
    dplyr::pull(subject_id)

  if (is.null(time_col)) {
    time_levels <- 1L
    samples$.time <- 1L
  } else {
    time_levels <- samples |>
      dplyr::distinct(!!rlang::sym(time_col)) |>
      dplyr::arrange(!!rlang::sym(time_col)) |>
      dplyr::pull(!!rlang::sym(time_col))
    samples$.time <- match(samples[[time_col]], time_levels)
  }

  samples <- samples |>
    dplyr::arrange(subject_id, .data$.time, sample_id)

  samples_order <- samples$sample_id
  sample_subject <- match(samples$subject_id, subjects_order)
  sample_time <- as.integer(samples$.time)
  N_samples <- length(samples_order)
  N_subjects <- length(subjects_order)
  subject_n <- tabulate(sample_subject, nbins = N_subjects)

  group_levels <- list()
  sample_group <- list()
  subject_group <- list()
  for (gc in group_cols) {
    lv <- samples |>
      dplyr::distinct(!!rlang::sym(gc)) |>
      dplyr::arrange(!!rlang::sym(gc)) |>
      dplyr::pull(!!rlang::sym(gc))
    idx <- match(samples[[gc]], lv)
    # one level per subject, NA if the column varies within the subject
    lo <- vapply(split(idx, sample_subject), min, integer(1))
    hi <- vapply(split(idx, sample_subject), max, integer(1))
    group_levels[[gc]] <- lv
    sample_group[[gc]] <- idx
    subject_group[[gc]] <- unname(dplyr::if_else(lo == hi, lo, NA_integer_))
  }

  peptides_order <- df_long |>
    dplyr::distinct(peptide_id) |>
    dplyr::arrange(peptide_id) |>
    dplyr::collect() |>
    dplyr::pull(peptide_id)

  # --- BUILD HITS FROM DISTINCT POSITIVE SAMPLES ------------------------------
  pos_pairs <- df_long |>
    dplyr::filter(!!rlang::sym(exist_col) > 0L) |>
    dplyr::distinct(sample_id, peptide_id) |>
    dplyr::collect()

  samp_index <- match(pos_pairs$sample_id, samples_order)
  pep_index <- match(pos_pairs$peptide_id, peptides_order)

  hits_split <- split(samp_index, pep_index)

  nonempty_pep_ids <- as.integer(names(hits_split))
  if (length(nonempty_pep_ids) == 0L) {
    .ph_abort("All peptides are zero (no positives).")
  }
  peptides_order <- peptides_order[nonempty_pep_ids]
  hits_by_peptide <- lapply(unname(hits_split), function(i) sort(as.integer(i)))

  # --- BUILD 64-bit BITSET ONCE -----------------------------------------------
  bs <- build_bitset_unpaired(hits_by_peptide, N_samples)
  if (bs$m != length(peptides_order)) {
    .ph_abort("Bitset/peptide dimension mismatch.")
  }

  list(
    samples_order  = samples_order,
    sample_subject = sample_subject,
    sample_time    = sample_time,
    sample_group   = sample_group,
    subject_group  = subject_group,
    subject_n      = subject_n,
    subjects_order = subjects_order,
    time_levels    = time_levels,
    group_levels   = group_levels,
    peptides_order = peptides_order,
    bitset_raw     = bs$data,
    bitset_m       = bs$m,
    bitset_words   = bs$n_words
  )
}
