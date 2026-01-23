#' @title Global Shift in Peptide-level Prevalence via Subject-level Permutation
#'
#' @description
#' Test for a **global** (antigen-wide) shift in peptide-level
#' prevalence between two groups within each `(rank, feature, group_col)`
#' stratum by aggregating per-peptide effects into a single Stouffer-type
#' statistic and assessing significance via **label permutation**.
#'
#' @details
#' For each stratum `(rank, feature, group_col)` that defines an ordered pair of
#' groups `(g1, g2)`, the procedure is:
#'
#' 1. **Smoothed peptide-level prevalences.**
#'    For each peptide, smoothed prevalences in `g1` and `g2` are
#'    \deqn{
#'      \hat p_k = \frac{x_k + \varepsilon}{n_k + \lambda \varepsilon},\quad k \in \{g1, g2\},
#'    }
#'    where `x_k` is the number of samples with presence (`exist > 0`), `n_k`
#'    is the number of distinct samples (or paired subjects), `\varepsilon` is
#'    `smooth_eps_num`, and `\lambda` is `smooth_eps_den_mult`. Peptides with
#'    \eqn{\max(\hat p_{g1}, \hat p_{g2}) <} `min_max_prev` are discarded.
#'
#' 2. **Per-peptide effect and z-score (`stat_mode`).**
#'
#'    - `"diff"`: effect \eqn{\Delta = \hat p_{g2} - \hat p_{g1}} with
#'      \deqn{
#'        SE = \sqrt{\hat p_{g1}(1 - \hat p_{g1}) / n_{g1} +
#'                   \hat p_{g2}(1 - \hat p_{g2}) / n_{g2}} ,
#'      }
#'      and z-score \eqn{z = \Delta / SE} (with a small floor on `SE`).
#'
#'    - `"asin"`: variance-stabilized difference of angles,
#'      \deqn{
#'        z = \frac{\arcsin \sqrt{\hat p_{g2}} - \arcsin \sqrt{\hat p_{g1}}}
#'                 {\sqrt{1 / (4 n_{g1}) + 1 / (4 n_{g2})}} .
#'      }
#'
#'    Each \eqn{z_i} is then **winsorized** to \eqn{\pm} `winsor_z`.
#'
#' 3. **Weights (`weight_mode`).**
#'
#'    - `"equal"`: all peptides get weight 1.
#'    - `"se_invvar"`: weights proportional to the inverse standard error
#'      (or the analogous denominator for `"asin"`).
#'    - `"n_eff_sqrt"`:
#'      \deqn{
#'        w_i = \sqrt{ n_{g1} \hat p_{g1} + n_{g2} \hat p_{g2} },
#'      }
#'      i.e. the square root of the expected number of positives across both
#'      groups.
#'
#' 4. **Combine into a single test statistic (Stouffer).**
#'    Let \eqn{w_i} be the weights and \eqn{z_i} the peptide-level z-scores.
#'
#'    - If `prev_strat = "none"`:
#'      \deqn{
#'        T_{\mathrm{obs}} =
#'          \frac{\sum_i w_i z_i}{\sqrt{\sum_i w_i^2}} .
#'      }
#'
#'    - If `prev_strat = "decile"`: peptides are binned into deciles of pooled
#'      prevalence
#'      \deqn{
#'        \frac{n_{g1}\hat p_{g1} + n_{g2}\hat p_{g2}}{n_{g1} + n_{g2}},
#'      }
#'      a Stouffer statistic \eqn{T_b} is computed within each bin, and the
#'      **mean** of the bin-level z-scores is reported as
#'      \eqn{T_{\mathrm{obs}}}.
#'
#' 5. **Permutation scheme and p-value.**
#'    A **two-sided** permutation p-value for the global shift is computed:
#'
#'    - **Paired design**: if both groups have measurements for the same
#'      `subject_id`, each subject’s labels (\code{g1} \eqn{\leftrightarrow}
#'      \code{g2}) are independently
#'      flipped with probability 1/2 and steps (1–4) are recomputed to obtain
#'      \eqn{T_b}.
#'
#'    - **Unpaired design**: sample labels are shuffled while preserving group
#'      sizes (random split into `n1` / `n2`), and steps (1–4) are recomputed
#'      to obtain \eqn{T_b}.
#'
#'    With \eqn{B} permutations, the p-value is
#'    \deqn{
#'      p =
#'      \frac{1 + \sum_{b=1}^{B} \mathbf{1}\{|T_b| \ge |T_{\mathrm{obs}}|\}}
#'           {1 + B},
#'    }
#'    the standard add-one estimator that remains non-zero under the global
#'    null.
#'
#' 6. **Multiplicity.**
#'    No multiple-testing correction is performed; `p_perm` is returned as
#'    computed.
#'
#' **Input requirements.**
#'
#' - `exist_col` is treated as 0/1 presence.
#' - There must be **at most one positive** per
#'   (`subject_id`, `peptide_id`, `group_col`, `group_value`); paired designs
#'   can have up to two positives across the two group levels. Violations
#'   trigger an error. Example (group levels A/B): for a single subject and
#'   peptide, you may have A=1 and B=0 (or A=0 and B=1, or A=1 and B=1), but you
#'   cannot have two rows both with A=1 (or two rows both with B=1).
#' - Non-peptide ranks specified in `rank_cols` must be resolvable from a
#'   peptide library (see `peptide_library` below).
#'
#' **Peptide library resolution.**
#'
#' When `rank_cols` includes non-`"peptide_id"` ranks, a peptide library is
#' required. It is resolved in the following order:
#'
#' 1. The explicit `peptide_library` argument.
#' 2. `x$peptide_library` if `x` is a `phip_data` with an attached library.
#' 3. `phiper::get_peptide_meta()` if available (i.e. `phiper` is installed and
#'    exports this function).
#'
#' If none of these are available, the function aborts with an informative
#' error.
#'
#' @section Parallelization:
#' The permutation contrasts are evaluated either sequentially or in parallel
#' using the **current** `future` backend:
#'
#' - The function does **not** change the global `future::plan()`.
#' - If `future` is installed, the number of workers is inferred from
#'   `future::nbrOfWorkers()`. To run in parallel, set up a `future` plan
#'   (e.g. `future::plan(multisession, workers = 8)`) **outside** the function.
#' - Internal C++ code is constrained to a single thread per R worker
#'   (`RcppParallel::setThreadOptions(numThreads = 1)` and BLAS/OpenMP limits),
#'   to avoid oversubscription when using multiple workers.
#'
#' Reproducibility of permutations is controlled via the global RNG:
#' call `set.seed()` before `compute_delta()` to obtain reproducible
#' results; each contrast draws its own seed from the global RNG state.
#'
#' @param x A `phip_data` object (recommended to attach a `peptide_library`
#'   when `rank_cols` includes non-`"peptide_id"` and not providing
#' `peptide_library` argument) or a data frame with the necessary columns.
#' @param rank_cols Character vector of rank columns (e.g. `"peptide_id"`,
#'   `"species"`, ...). For non-peptide ranks, annotations are joined from a
#'   peptide library (see Details).
#' @param group_cols Character vector of grouping columns that define
#'   universes; pairwise group contrasts are constructed within each
#'   `(rank, feature, group_col)` stratum.
#' @param exist_col Name of the 0/1 presence column. Default `"exist"`.
#' @param paired_by Optional single column name identifying paired samples
#'   (e.g. subject ID). When provided, paired tests are used where applicable.
#' @param interaction Logical; reserved for future use (currently ignored). If
#'   `TRUE`, an interaction of the first two `group_cols` may be used as an
#'   additional `group_col` in future versions.
#' @param combine_cols Optional length-2 character vector; reserved for future
#'   use. If non-`NULL`, only this interaction would be used as `group_col`
#'   in future versions.
#' @param interaction_sep Separator for interaction values. Default `"::"`.
#' @param B_permutations Number of permutations `B` used for the permutation
#'   p-value. Default `2000L`. Must be at least 100.
#' @param smooth_eps_num Laplace numerator epsilon \eqn{\varepsilon} used to
#'   smooth prevalence estimates. Default `0.5`.
#' @param smooth_eps_den_mult Multiplicative factor \eqn{\lambda} applied to
#'   the denominator epsilon in the smoothing formula. Default `2.0`.
#' @param min_max_prev Minimum required value of `max(p1, p2)` for a peptide to
#'   be retained in the analysis. Default `0.0` (no prevalence filter).
#' @param weight_mode One of `c("equal", "se_invvar", "n_eff_sqrt")`,
#'   controlling how peptide-level z-scores are weighted in the Stouffer
#'   combination (see Details).
#' @param stat_mode One of `c("diff", "asin")`, determining whether z-scores
#'   are based on prevalence differences or the arcsin–sqrt transform.
#' @param prev_strat One of `c("none", "decile")`. If `"decile"`, the Stouffer
#'   statistic is computed within prevalence deciles and the mean of
#'   decile-level z-scores is used as the global test statistic.
#' @param winsor_z Winsorization threshold applied to peptide-level z-scores.
#'   Values beyond `±winsor_z` are truncated. Default `4.0`.
#' @param rank_feature_keep Optional **named list** mapping `rank` to a vector
#'   of `feature` values to keep. Only rank–feature strata in this list are
#'   tested; others are dropped after the peptide-level pivot.
#' @param peptide_library Optional data frame providing peptide annotations for
#'   non-peptide ranks. Must at least contain `peptide_id` and all requested
#'   `rank_cols` besides `"peptide_id"`. If `NULL`, the function falls back to
#'   `x$peptide_library` or `phiper::get_peptide_meta()` if available.
#' @param log Logical; if `TRUE`, write progress messages (per contrast and
#'   overall) using the package's logging helpers.
#' @param log_file Path to a log file used by the logging helpers if `log` is
#'   `TRUE`. Default `"compute_delta.log"`.
#'
#' @return
#' A tibble with one row per tested stratum:
#'
#' - `rank`, `feature`: rank and feature identifiers.
#' - `group_col`, `group1`, `group2`: grouping column and the ordered pair of
#'   groups compared.
#' - `design`: `"paired"` or `"unpaired"`.
#' - `n_subjects_paired`: number of paired subjects used in a paired design
#'   (or `NA` for unpaired).
#' - `n_peptides_used`: number of peptides contributing to the test.
#' - `m_eff`: effective number of peptides after prevalence filtering and
#'   weighting (as returned by the C++ helper).
#' - `T_obs`: observed Stouffer-type test statistic.
#' - `p_perm`: two-sided permutation p-value.
#' - `b`: number of permutations actually used (may be `< B_permutations` if
#'   early stopping is implemented in the C++ helper).
#' - `max_delta`, `frac_delta_pos`, `frac_delta_pos_w`:
#'   maximum absolute peptide-level prevalence difference
#'   and unweighted/weighted fractions of positive peptide-level deltas.
#'
#' @examples
#' \donttest{
#' # Load example PhIP-Seq data shipped with the package
#' pd <- phip_load_example_data()
#'
#' # Small unpaired subset with a mock peptide library
#' pd_filt <- pd |>
#'   dplyr::filter(
#'     peptide_id %in% c("16627", "5243", "24799", "16196", "18003"),
#'     timepoint == "T1"
#'   ) |>
#'   dplyr::collect()
#'
#' mock_peplib <- data.frame(
#'   peptide_id = c("16627", "5243", "24799", "16196", "18003"),
#'   species    = rep("mock_species", 5),
#'   stringsAsFactors = FALSE
#' )
#'
#' res <- compute_delta(
#'   x                  = pd_filt,
#'   exist_col          = "exist",
#'   rank_cols          = "species",
#'   group_cols         = "group",
#'   peptide_library    = mock_peplib,
#'   B_permutations     = 500L,  # smaller for speed
#'   weight_mode        = "n_eff_sqrt",
#'   stat_mode          = "asin",
#'   prev_strat         = "none",
#'   winsor_z           = Inf,
#'   rank_feature_keep  = list(species = NULL),
#'   log                = FALSE
#' )
#' res
#' }
#'
#' @export
compute_delta <- function(
  x, rank_cols, group_cols,
  exist_col = "exist",
  paired_by = NULL,
  interaction = FALSE,
  combine_cols = NULL,
  interaction_sep = "::",
  B_permutations = 2000L,
  smooth_eps_num = 0.5,
  smooth_eps_den_mult = 2.0,
  min_max_prev = 0.0,
  weight_mode = c("equal", "se_invvar", "n_eff_sqrt"),
  stat_mode = c("diff", "asin"),
  prev_strat = c("none", "decile"),
  winsor_z = 4.0,
  rank_feature_keep = NULL,
  peptide_library = NULL,
  log = FALSE,
  log_file = "compute_delta.log"
) {
  # --- 0) Argument validation -------------------------------------------------
  chk::chk_character(rank_cols)
  chk::chk_true(length(rank_cols) >= 1)
  chk::chk_character(group_cols)
  chk::chk_true(length(group_cols) >= 1)
  chk::chk_string(exist_col)
  if (!is.null(paired_by)) chk::chk_string(paired_by)
  chk::chk_number(B_permutations)
  chk::chk_true(B_permutations >= 100)

  # --- 1) Prepare data once ---------------------------------------------------
  # Required columns from `x`
  need_cols <- c("sample_id", "subject_id", "peptide_id", exist_col, group_cols)
  if (!is.null(paired_by)) {
    need_cols <- c(need_cols, paired_by)
  }
  need_cols <- unique(need_cols)

  if (inherits(x, "phip_data")) {
    df_long <- x$data_long |>
      dplyr::select(tidyselect::any_of(need_cols))
  } else {
    chk::chk_data(x)
    miss <- setdiff(need_cols, colnames(x))
    if (length(miss)) {
      .ph_abort("Missing required columns",
        bullets = paste("-", miss)
      )
    }
    df_long <- tibble::as_tibble(x) |>
      dplyr::select(tidyselect::any_of(need_cols))
  }

  # --- STRICT HITS GUARD: at most one positive per (subject_id, peptide_id,
  # group value) ---
  dup_pos <- df_long |>
    dplyr::filter(!!rlang::sym(exist_col) > 0L) |>
    tidyr::pivot_longer(tidyselect::all_of(group_cols),
      names_to = "group_col", values_to = "group_value"
    ) |>
    dplyr::count(subject_id, peptide_id, group_col, group_value,
      name = "n_pos"
    ) |>
    dplyr::filter(n_pos > 1L) |>
    dplyr::collect()

  if (nrow(dup_pos) > 0L) {
    eg <- dup_pos |>
      dplyr::slice_head(n = 10) |>
      dplyr::mutate(example = paste0(
        "subject=", subject_id,
        ", peptide=", peptide_id,
        ", group_col=", group_col,
        ", group_value=", group_value,
        ", n_pos=", n_pos
      )) |>
      dplyr::pull(example)
    .ph_abort(
      "Invalid input: duplicate positives within the SAME group for some
      (subject_id, peptide_id). One positive per group is allowed (paired
      designs can have up to 2 across groups).",
      bullets = c(eg, if (nrow(dup_pos) > 10) {
        sprintf(
          "... and %d more.",
          nrow(dup_pos) - 10
        )
      })
    )
  }


  # --- SUBJECT METADATA (1 row per subject; lazy-safe) ------------------------
  subjects_meta <- df_long |>
    dplyr::group_by(subject_id) |>
    dplyr::summarise(
      dplyr::across(
        tidyselect::all_of(group_cols),
        dplyr::first
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(subject_id) |>
    dplyr::collect()

  # --- FREEZE ORDERS (lazy-safe) ----------------------------------------------
  subjects_order <- df_long |>
    dplyr::distinct(subject_id) |>
    dplyr::arrange(subject_id) |>
    dplyr::collect() |>
    dplyr::pull(subject_id)

  peptides_order <- df_long |>
    dplyr::distinct(peptide_id) |>
    dplyr::arrange(peptide_id) |>
    dplyr::collect() |>
    dplyr::pull(peptide_id)

  # --- BUILD HITS DIRECTLY FROM DISTINCT POSITIVES ----------------------------
  pos_pairs <- df_long |>
    dplyr::filter(!!rlang::sym(exist_col) > 0L) |>
    dplyr::distinct(subject_id, peptide_id) |>
    dplyr::collect()

  subj_index <- match(pos_pairs$subject_id, subjects_order)
  pep_index <- match(pos_pairs$peptide_id, peptides_order)

  hits_split <- split(subj_index, pep_index)

  nonempty_pep_ids <- as.integer(names(hits_split))
  if (length(nonempty_pep_ids) == 0L) {
    .ph_abort("All peptides are zero after hits guard (no positives).")
  }
  peptides_order <- peptides_order[nonempty_pep_ids]

  m <- length(peptides_order)
  hits_by_peptide <- vector("list", m)
  hits_by_peptide[seq_along(nonempty_pep_ids)] <- lapply(hits_split, as.integer)

  # --- BUILD 64-bit BITSET ONCE -----------------------------------------------
  N_subjects <- length(subjects_order)
  bs <- build_bitset_unpaired(hits_by_peptide, N_subjects)

  bitset_raw <- bs$data
  bitset_m <- bs$m
  bitset_words <- bs$n_words
  if (bitset_m != length(peptides_order)) {
    .ph_abort("Bitset/peptide dimension mismatch.")
  }

  # Ready for downstream steps:
  # subjects_order, peptides_order, subjects_meta, bitset_raw, bitset_m,
  # bitset_words


  # Objects ready for downstream use:
  # - subjects_order: character/ID vector (row map)
  # - peptides_order: character/ID vector (col map; all-zero peptides removed)
  # - subjects_meta : per-subject metadata (for contrasts)
  # - bitset_raw, bitset_m, bitset_words : packed 64-bit bitset (col-major)

  # --- 2) Peptide library attach (simple, in R) -------------------------------
  ranks_need <- setdiff(rank_cols, "peptide_id")

  get_lib_tbl <- function() {
    if (length(ranks_need) == 0L) {
      # Only peptide_id requested --> trivial long map
      return(
        tibble::tibble(
          peptide_id = peptides_order, rank = "peptide_id",
          feature = peptides_order
        )
      )
    }

    lib_src <- NULL
    if (!is.null(peptide_library)) {
      lib_src <- peptide_library
    } else if (inherits(x, "phip_data") && !is.null(x$peptide_library)) {
      lib_src <- x$peptide_library
    } else if (rlang::is_installed("phiper") &&
      "get_peptide_meta" %in% getNamespaceExports("phiper")) {
      # auto-fetch from phiper if available
      lib_src <- phiper::get_peptide_meta()
    } else {
      .ph_abort(
        "Peptide library required for non-peptide ranks.",
        bullets = c(
          "- Provide `peptide_library` with the needed columns,",
          "- Or attach a peptide_library to `x` (phip_data),",
          "- Or ensure phiper::get_peptide_meta() is available."
        )
      )
    }

    # Select needed columns and collect to R (handles DuckDB/lazy)
    lib_needed <- c("peptide_id", ranks_need)
    miss_tax <- setdiff(ranks_need, colnames(lib_src))
    if (length(miss_tax)) {
      .ph_abort("Peptide library missing required columns:",
        bullets = paste("-", miss_tax)
      )
    }

    lib_small <- lib_src |>
      dplyr::select(tidyselect::any_of(lib_needed)) |>
      dplyr::distinct()

    if (inherits(lib_small, "tbl_sql") || inherits(lib_small, "tbl_lazy")) {
      lib_small <- lib_small |> dplyr::collect()
    }

    # Keep only peptides we actually retained (peptides_order), then pivot to
    # long
    lib_small <- lib_small |>
      dplyr::filter(.data$peptide_id %in% peptides_order)

    rank_map_long <- lib_small |>
      tidyr::pivot_longer(
        cols      = tidyselect::all_of(ranks_need),
        names_to  = "rank",
        values_to = "feature"
      ) |>
      dplyr::mutate(feature = as.character(.data$feature)) |>
      dplyr::filter(!is.na(.data$feature)) |>
      dplyr::select(peptide_id, rank, feature) |>
      dplyr::distinct()

    # Also add the direct peptide_id rank if requested
    if ("peptide_id" %in% rank_cols) {
      rank_map_long <- dplyr::bind_rows(
        rank_map_long,
        tibble::tibble(
          peptide_id = peptides_order, rank = "peptide_id",
          feature = peptides_order
        )
      )
    }

    rank_map_long
  }

  rank_map_long <- get_lib_tbl()
  # rank_map_long columns: peptide_id, rank, feature  (in R memory)

  # --- 3) Initialize progress logging -----------------------------------------
  progress_path <- .progress_file(log_file)
  if (isTRUE(log)) .progress_init(progress_path)

  # --- 4) Construct grouping universes (contrasts) ----------------------------

  # (a) Filter rank_map_long by rank_feature_keep (if provided)
  if (!is.null(rank_feature_keep) && length(rank_feature_keep)) {
    for (rk in names(rank_feature_keep)) {
      vals <- rank_feature_keep[[rk]]
      if (is.null(vals)) {
        next
      } else {
        keep_vals <- as.character(vals)
        rank_map_long <- rank_map_long |>
          dplyr::filter(!(rank == rk) | feature %in% keep_vals)
      }
    }
  }

  # (b) rank-feature counts after peptide filtering
  rf_counts <- rank_map_long |>
    dplyr::count(rank, feature, name = "n_peptides")

  # (c) levels per group_col across all rows (not per-subject)
  # Using subjects_meta can drop within-subject levels (e.g., timepoints),
  # which yields zero contrasts and an empty result.
  lvl <- df_long |>
    tidyr::pivot_longer(tidyselect::all_of(group_cols),
      names_to = "group_col", values_to = "group_value"
    ) |>
    dplyr::distinct(group_col, group_value) |>
    dplyr::collect()

  # (d) all ordered pairs per group_col (no many-to-many warning)
  all_pairs <- lvl |>
    dplyr::group_by(group_col) |>
    dplyr::reframe({
      vals <- sort(unique(group_value))
      if (length(vals) >= 2L) {
        cmb <- utils::combn(vals, 2)
        tibble::tibble(group1 = cmb[1, ], group2 = cmb[2, ])
      } else {
        tibble::tibble(group1 = character(), group2 = character())
      }
    })

  # (e) paired/unpaired by subject overlap
  if (!is.null(paired_by)) {
    pair_sym <- rlang::sym(paired_by)

    subj_groups <- df_long |>
      tidyr::pivot_longer(
        tidyselect::all_of(group_cols),
        names_to  = "group_col",
        values_to = "group_value"
      ) |>
      dplyr::distinct(!!pair_sym, group_col, group_value) |>
      dplyr::collect() |>
      dplyr::rename(pair_id = !!pair_sym)

    pairs_shared <- subj_groups |>
      dplyr::rename(group1 = group_value) |>
      dplyr::inner_join(
        subj_groups |> dplyr::rename(group2 = group_value),
        by = c("pair_id", "group_col")
      ) |>
      dplyr::filter(.data$group1 < .data$group2) |>
      dplyr::count(group_col, group1, group2, name = "n_shared")

    contrasts_base <- all_pairs |>
      dplyr::left_join(pairs_shared, by = c("group_col", "group1", "group2")) |>
      dplyr::mutate(
        design = dplyr::if_else(
          !is.na(.data$n_shared) & .data$n_shared > 0L,
          "paired",
          "unpaired"
        )
      ) |>
      dplyr::select(group_col, group1, group2, design)
  } else {
    contrasts_base <- all_pairs |>
      dplyr::mutate(design = "unpaired") |>
      dplyr::select(group_col, group1, group2, design)
  }

  # (f) master plan = rank-feature strata × contrasts (only strata with
  # peptides)
  master_plan <- tidyr::crossing(
    rf_counts |> dplyr::select(rank, feature, n_peptides),
    contrasts_base
  ) |>
    dplyr::filter(.data$n_peptides > 0L)

  # --- 5) Work weighting (by usable peptides) ---------------------------------
  master_plan <- master_plan |>
    dplyr::mutate(work_weight = .data$n_peptides)

  # --- 7) Header line (CPP-only) ----------------------------------------------
  n_contrasts <- nrow(master_plan)
  B_num <- as.numeric(B_permutations)
  pep_sum <- sum(as.numeric(master_plan$n_peptides), na.rm = TRUE)
  work_total <- pep_sum * B_num
  peps_total <- length(peptides_order)
  hdr_fast <- "CPP"
  log_to_file <- isTRUE(log) && (!is.null(log_file) && nzchar(log_file))

  if (isTRUE(log)) {
    hdr_bullets <- c(
      sprintf("total contrasts: %d", n_contrasts),
      sprintf("B per contrast: %d", B_permutations),
      sprintf("peptides after filtering: %d", peps_total),
      sprintf("weighted work total (sum B*peptides): %.0f", work_total),
      sprintf(
        "mode: %s",
        if (rlang::is_installed("future") &&
          isTRUE(tryCatch(future::nbrOfWorkers() > 1L,
            error = function(...) FALSE
          ))) {
          "parallel (future backend)"
        } else {
          "sequential"
        }
      ),
      sprintf("engine: %s", hdr_fast)
    )
    if (log_to_file) {
      .ph_log_info_file(log_file, "Starting permutation computations",
        bullets = hdr_bullets
      )
    } else {
      .ph_log_info("Starting permutation computations", bullets = hdr_bullets)
    }
  }

  # --- Prep indices for CPP ---------------------------------------------------
  # subject row index map (1..N)
  subj_row_map <- tibble::tibble(
    subject_id = subjects_order,
    row = seq_along(subjects_order)
  )

  # group membership with subject row indices (lazy-safe collect already done
  # earlier)
  subj_groups_idx <- df_long |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(group_cols),
      names_to = "group_col", values_to = "group_value"
    ) |>
    dplyr::distinct(subject_id, group_col, group_value) |>
    dplyr::collect() |>
    dplyr::left_join(subj_row_map, by = "subject_id")

  # peptide column index map (within bitset columns 1..m)
  pep_col_map <- stats::setNames(seq_along(peptides_order), peptides_order)

  get_pep_cols <- function(rk, ft) {
    pid <- rank_map_long |>
      dplyr::filter(.data$rank == rk, .data$feature == ft) |>
      dplyr::pull(peptide_id)
    as.integer(pep_col_map[pid])
  }

  get_group_rows <- function(gc, gval) {
    subj_groups_idx |>
      dplyr::filter(.data$group_col == gc, .data$group_value == gval) |>
      dplyr::pull(row) |>
      as.integer()
  }

  get_paired_hits <- function(gc, g1, g2) {
    pair_sym <- rlang::sym(paired_by)

    s1 <- df_long |>
      dplyr::filter(.data[[gc]] == g1) |>
      dplyr::pull(!!pair_sym)
    s2 <- df_long |>
      dplyr::filter(.data[[gc]] == g2) |>
      dplyr::pull(!!pair_sym)

    pair_ids <- sort(intersect(unique(s1), unique(s2)))
    if (!length(pair_ids)) {
      return(NULL)
    }

    id_map <- stats::setNames(seq_along(pair_ids), pair_ids)

    # G1: block-level presence
    g1_hits <- df_long |>
      dplyr::filter(
        !!rlang::sym(exist_col) > 0L,
        .data[[gc]] == g1,
        .data[[paired_by]] %in% pair_ids
      ) |>
      dplyr::transmute(
        peptide_id,
        idx = id_map[.data[[paired_by]]]
      ) |>
      dplyr::distinct(peptide_id, idx) |>
      dplyr::group_by(peptide_id) |>
      dplyr::summarise(idx_hits = list(as.integer(sort(idx))), .groups = "drop")

    # G2:
    g2_hits <- df_long |>
      dplyr::filter(
        !!rlang::sym(exist_col) > 0L,
        .data[[gc]] == g2,
        .data[[paired_by]] %in% pair_ids
      ) |>
      dplyr::transmute(
        peptide_id,
        idx = id_map[.data[[paired_by]]]
      ) |>
      dplyr::distinct(peptide_id, idx) |>
      dplyr::group_by(peptide_id) |>
      dplyr::summarise(idx_hits = list(as.integer(sort(idx))), .groups = "drop")

    hh <- dplyr::full_join(g1_hits, g2_hits,
      by = "peptide_id",
      suffix = c("_g1", "_g2")
    ) |>
      dplyr::mutate(
        idx_hits_g1 = lapply(idx_hits_g1, function(x) x %||% integer(0)),
        idx_hits_g2 = lapply(idx_hits_g2, function(x) x %||% integer(0))
      ) |>
      dplyr::arrange(peptide_id)

    list(
      pair_ids = pair_ids,
      hits_g1  = hh$idx_hits_g1,
      hits_g2  = hh$idx_hits_g2,
      pep_ids  = hh$peptide_id
    )
  }

  # --- core runner for one contrast (calls CPP helper + logs) ----------
  .run_one <- function(i) {
    st <- master_plan[i, , drop = FALSE]
    pep_cols <- get_pep_cols(st$rank, st$feature)
    if (length(pep_cols) == 0L) {
      return(NULL)
    }

    g1_rows <- get_group_rows(st$group_col, st$group1)
    g2_rows <- get_group_rows(st$group_col, st$group2)
    if (length(g1_rows) == 0L || length(g2_rows) == 0L) {
      return(NULL)
    }

    # draw per-contrast seed from global RNG; reproducible via set.seed()
    # outside
    seed_i <- as.integer(sample.int(.Machine$integer.max, 1L))

    # Defaults for UNPAIRED call
    hits_g1 <- list()
    hits_g2 <- list()
    P <- 0L

    if (identical(st$design, "paired")) {
      paired <- get_paired_hits(st$group_col, st$group1, st$group2)
      if (is.null(paired)) {
        return(NULL)
      }

      # restrict to peptides in this (rank, feature)
      pep_cols_all <- get_pep_cols(st$rank, st$feature)
      pep_keep_ids <- names(pep_col_map)[pep_cols_all]
      keep_mask <- paired$pep_ids %in% pep_keep_ids

      if (!any(keep_mask)) {
        return(NULL)
      }

      hits_g1 <- paired$hits_g1[keep_mask]
      hits_g2 <- paired$hits_g2[keep_mask]
      P <- length(paired$pair_ids)
    }

    res <- cpp_shift_contrast(
      bitset_raw       = bitset_raw,
      n_words          = bitset_words,
      pep_cols         = pep_cols,
      g1_rows          = g1_rows,
      g2_rows          = g2_rows,
      hits_g1_paired   = hits_g1,
      hits_g2_paired   = hits_g2,
      P                = as.integer(P),
      B                = as.integer(B_permutations),
      seed             = seed_i,
      smooth_eps_num   = smooth_eps_num,
      smooth_eps_den   = smooth_eps_den_mult,
      min_max_prev     = min_max_prev,
      weight_mode      = weight_mode,
      stat_mode        = stat_mode,
      prev_strat       = prev_strat,
      winsor_z         = winsor_z,
      design           = st$design
    )
    if (is.null(res)) {
      return(NULL)
    }

    row <- dplyr::bind_cols(
      st[, c("rank", "feature", "group_col", "group1", "group2", "design")],
      tibble::tibble(
        n_subjects_paired = if (identical(st$design, "paired")) {
          as.integer(P)
        } else {
          NA_integer_
        },
        n_peptides_used = as.integer(res$n_peptides_used),
        m_eff = as.numeric(res$m_eff),
        T_obs = as.numeric(res$T_obs),
        T_null_mean = as.numeric(res$T_null_mean),
        T_null_sd = as.numeric(res$T_null_sd),
        T_obs_stand = {
          mean_null <- as.numeric(res$T_null_mean)
          sd_null <- as.numeric(res$T_null_sd)
          ifelse(
            !is.na(mean_null) & !is.na(sd_null) & sd_null > 0,
            (as.numeric(res$T_obs) - mean_null) / sd_null,
            NA_real_
          )
        },
        Z_from_p = {
          pval <- as.numeric(res$p_perm)
          T_val <- as.numeric(res$T_obs)
          ifelse(
            is.na(pval),
            NA_real_,
            sign(T_val) * stats::qnorm(pval / 2, lower.tail = FALSE)
          )
        },
        b = as.integer(res$b),
        p_perm = as.numeric(res$p_perm),
        max_delta       = as.numeric(res$max_delta),
        frac_delta_pos  = as.numeric(res$frac_delta_pos),
        frac_delta_pos_w = as.numeric(res$frac_delta_pos_w)
      )
    )

    if (isTRUE(log)) {
      msg <- sprintf(
        "%s / %s | n_pep=%d b=%d p=%g",
        row$rank, row$feature, row$n_peptides_used, row$b, row$p_perm
      )
      if (log_to_file) .ph_log_info_file(log_file, msg) else .ph_log_info(msg)

      inc <- as.numeric(row$n_peptides_used) * B_num
      done_w <- .progress_add(progress_path, inc)

      prog <- if (is.finite(done_w) && is.finite(work_total) &&
        work_total > 0) {
        pmin(100, 100 * (done_w / work_total))
      } else {
        NA_real_
      }

      if (log_to_file) {
        .ph_log_info_file(log_file, sprintf(
          "progress=%s [C]",
          if (is.na(prog)) "n/a" else sprintf("%.1f%%", prog)
        ))
      } else {
        .ph_log_info(sprintf(
          "progress=%s [C]",
          if (is.na(prog)) "n/a" else sprintf("%.1f%%", prog)
        ))
      }
    }

    row
  }

  # --- 7/8) Computation path: sequential vs parallel --------------------------
  result_rows <- list()

  # decide workers
  # decide workers based on current future backend
  n_workers <- 1L
  if (rlang::is_installed("future")) {
    n_workers <- tryCatch(future::nbrOfWorkers(), error = function(...) 1L)
  }
  n_workers <- max(1L, min(as.integer(n_workers), nrow(master_plan)))

  # ---- dynamic scheduling for future.apply ----
  op_old <- options(
    future.scheduling   = Inf,
    future.chunk.size   = 1
  )
  on.exit(options(op_old), add = TRUE)

  # ---- avoid inner oversubscription (TBB inside RcppParallel) ----
  if (rlang::is_installed("RcppParallel")) {
    RcppParallel::setThreadOptions(numThreads = 1) # fairness > peak single-task
  } #            speed
  env_threads <- c(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1"
  )
  result_rows <- withr::with_envvar(env_threads, {
    if (n_workers == 1L) {
      # --- sequential
      for (i in seq_len(nrow(master_plan))) {
        row <- .run_one(i)
        if (!is.null(row)) result_rows[[length(result_rows) + 1L]] <- row
      }
      result_rows
    } else {
      # --- parallel via future.apply (respect existing future plan)
      idx <- order(master_plan$work_weight, decreasing = TRUE)
      future.apply::future_lapply(
        X = idx, FUN = .run_one,
        future.seed = TRUE,
        future.scheduling = Inf,
        future.chunk.size = 1,
        future.packages = c("dplyr", "tibble", "tidyr")
      )
    }
  })

  # ---- Final log -------------------------------------------------------------
  if (isTRUE(log)) {
    # figure out ranks for footer without failing on empty res
    ranks_footer <- tryCatch(
      paste(sort(unique(master_plan$rank)), collapse = ", "),
      error = function(...) NA_character_
    )
    headline <- "Done (permutation shift)"
    bullets <- c(
      sprintf("ranks: %s", ranks_footer),
      sprintf("B: %d", B_permutations),
      sprintf("workers: %d", ifelse(n_workers > 1L, n_workers, 1L))
    )
    if (log_to_file) {
      .ph_log_ok_file(log_file, headline, bullets = bullets)
    } else {
      .ph_log_ok(headline, bullets = bullets)
    }
  }

  # ---- Collect results -------------------------------------------------------
  rr <- Filter(Negate(is.null), result_rows)
  res <- if (length(rr)) dplyr::bind_rows(rr) else tibble::tibble()

  # Return early if empty
  if (nrow(res) == 0L) {
    return(tibble::tibble(
      rank = character(), feature = character(), group_col = character(),
      group1 = character(), group2 = character(), design = character(),
      n_subjects_paired = integer(), n_peptides_used = integer(),
      m_eff = numeric(),
      T_obs = numeric(), T_null_mean = numeric(), T_null_sd = numeric(),
      T_obs_stand = numeric(), Z_from_p = numeric(),
      p_perm = numeric(), b = integer(),
      max_delta = numeric(), frac_delta_pos = numeric(),
      frac_delta_pos_w = numeric()
    ))
  }

  # Ensure numeric types we rely on
  res <- res |>
    dplyr::mutate(
      p_perm = as.numeric(p_perm),
      m_eff  = as.numeric(m_eff)
    )

  # ---- final select/arrange ---------------------------------------
  res <- res |>
    dplyr::mutate(
      n_subjects_paired = dplyr::coalesce(.data$n_subjects_paired, NA_integer_)
    ) |>
    dplyr::select(
      rank, feature, group_col, group1, group2, design,
      n_subjects_paired, n_peptides_used, m_eff,
      T_obs, T_null_mean, T_null_sd, T_obs_stand, Z_from_p, p_perm, b,
      max_delta, frac_delta_pos, frac_delta_pos_w
    ) |>
    dplyr::arrange(rank, feature, group_col, group1, group2)

  return(res)
}

# ====== LOG-TO-FILE HELPERS (phiper-style, same layout) =======================
.progress_file <- function(log_file = NULL, stream_path = NULL) {
  base <- if (!is.null(log_file) && nzchar(log_file)) {
    log_file
  } else {
    stream_path %||% "ph_prev_shift"
  }
  paste0(base, ".progress.bin")
}

.progress_init <- function(path) {
  if (file.exists(path)) unlink(path)
  con <- file(path, open = "wb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  serialize(0, con, xdr = FALSE) # start at 0
  invisible(TRUE)
}

.progress_lock_path <- function(p) paste0(p, ".lock")

.progress_lock <- function(lock_path, timeout_ms = 60000, poll_sec = 0.05) {
  start_time <- Sys.time()
  repeat {
    if (file.create(lock_path)) {
      return(lock_path)
    }
    if (as.numeric(difftime(Sys.time(), start_time, units = "secs")) * 1000 >=
        timeout_ms) {
      .ph_abort("Timed out waiting for progress lock. Disable logging or try
                again.")
    }
    Sys.sleep(poll_sec)
  }
}

.progress_add <- function(path, delta) {
  lock_path <- .progress_lock(.progress_lock_path(path), timeout_ms = 60000)
  on.exit(try(unlink(lock_path), silent = TRUE), add = TRUE)

  cur <- 0
  if (file.exists(path)) {
    conr <- file(path, open = "rb")
    on.exit(try(close(conr), silent = TRUE), add = TRUE)
    cur <- tryCatch(unserialize(conr), error = function(e) 0)
  }
  newv <- cur + as.numeric(delta)
  conw <- file(path, open = "wb")
  on.exit(try(close(conw), silent = TRUE), add = TRUE)
  serialize(newv, conw, xdr = FALSE)
  newv
}

.ph_log_to_file <- function(path, level = "INFO", headline,
                            step = NULL, bullets = NULL) {
  lines <- .ph_compose_lines(level, headline, step, bullets)
  cat(paste0(lines, collapse = "\n"), "\n",
    sep = "", file = path,
    append = TRUE
  )
  invisible(lines)
}

.ph_log_ok_file <- function(path, headline, step = NULL, bullets = NULL) {
  .ph_log_to_file(path, "OK", headline, step, bullets)
}

.ph_log_info_file <- function(path, headline, step = NULL, bullets = NULL) {
  .ph_log_to_file(path, "INFO", headline, step, bullets)
}

.ph_log_warn_file <- function(path, headline, step = NULL, bullets = NULL) {
  .ph_log_to_file(path, "WARN", headline, step, bullets)
}
