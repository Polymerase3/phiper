# ==============================================================================
#' Global shift in peptide-level prevalence via subject-level permutation
#'
#' @importFrom foreach %dopar% foreach registerDoSEQ
#' @importFrom doParallel registerDoParallel
#'
#' @description
#' Tests for a **global** (genome/antigen-wide) shift in peptide-level prevalence
#' between two groups within each `(rank, feature, group_col)` stratum by
#' aggregating per-peptide effects into a single Stouffer-type statistic and
#' assessing significance via **label permutation** (paired or unpaired, chosen
#' automatically from the data).
#'
#' @details
#' **Overview.** For each stratum `(rank, feature, group_col)` that yields a
#' unique ordered pair of groups `(g1, g2)`, the method:
#'
#' 1) **Peptide-level prevalence with smoothing.**
#'    For each peptide, compute smoothed prevalences in `g1` and `g2`:
#'    \deqn{ \hat p_k = \frac{x_k + \varepsilon}{n_k + \lambda \varepsilon}, \; k \in \{g1,g2\} , }
#'    where `x_k` is the number of samples with presence (`exist > 0`), `n_k`
#'    is the number of (distinct) samples (or paired subjects), `\varepsilon`
#'    is `smooth_eps_num`, and `\lambda` is `smooth_eps_den_mult`. Peptides with
#'    `max( \hat p_{g1}, \hat p_{g2} ) < min_max_prev` are discarded.
#'
#' 2) **Per-peptide effect and z-score (`stat_mode`).**
#'    - `"diff"`: effect \eqn{\Delta = \hat p_{g2} - \hat p_{g1}} with
#'      \eqn{ SE = \sqrt{\hat p_{g1}(1-\hat p_{g1})/n_{g1} + \hat p_{g2}(1-\hat p_{g2})/n_{g2}} } and
#'      \eqn{ z = \Delta / SE } (guarded by a tiny floor).
#'    - `"asin"`: variance-stabilized difference of angles:
#'      \eqn{ z = \frac{\arcsin\sqrt{\hat p_{g2}} - \arcsin\sqrt{\hat p_{g1}}}{ \sqrt{ 1/(4n_{g1}) + 1/(4n_{g2}) } } .}
#'    Each \eqn{z} is **winsorized** to \eqn{\pm} `winsor_z`.
#'
#' 3) **Weights (`weight_mode`).**
#'    - `"equal"`: all peptides weight 1.
#'    - `"se_invvar"`: inverse standard error (or the analogous denominator for
#'      `"asin"`).
#'    - `"n_eff_sqrt"`: \eqn{ \sqrt{ n_{g1}\hat p_{g1} + n_{g2}\hat p_{g2} } } (proxy for expected positives).
#'
#' 4) **Combine into one test statistic (Stouffer).**
#'    Let \eqn{w_i} be weights and \eqn{z_i} peptide z-scores.
#'    - If `prev_strat = "none"`:
#'      \deqn{ T_{\mathrm{obs}} = \frac{ \sum_i w_i z_i }{ \sqrt{ \sum_i w_i^2 } } . }
#'    - If `prev_strat = "decile"`: bin peptides by deciles of pooled prevalence
#'      \eqn{ (\;n_{g1}\hat p_{g1} + n_{g2}\hat p_{g2}\;)/(\;n_{g1}+n_{g2}\;) } and combine a
#'      Stouffer \eqn{T_b} within each bin; report the **mean** of bin-level z’s
#'      as \eqn{T_{\mathrm{obs}}} (stabilizes across prevalence regimes).
#'
#' 5) **Permutation scheme and p-value.**
#'    We compute a **two-sided** permutation p-value for the global shift:
#'
#'    - **Paired design:** when both groups have measurements for the same
#'      `subject_id`, we independently flip each subject’s labels (`g1`↔`g2`)
#'      with 50/50 probability, recompute all steps (1–4) to get \eqn{T_b}.
#'    - **Unpaired design:** we shuffle sample labels while **preserving group
#'      sizes** (draw a random split of samples into `n1`/`n2`), recompute (1–4)
#'      to get \eqn{T_b}.
#'
#'    With \eqn{B} permutations, the p-value is
#'    \deqn{ p = \frac{1 + \sum_{b=1}^{B} \mathbf{1}\{ |T_b| \ge |T_{\mathrm{obs}}| \} }{1 + B} , }
#'    which is the standard *add-one* estimate ensuring non-zero p under the
#'    global null.
#'
#' 6) **Multiplicity.** For each `rank`, apply BH to `p_perm` → `p_adj_rank`.
#'
#' **Notes.**
#' - If a stratum doesn’t contain exactly two group levels, it’s skipped.
#' - Progress can be printed every `progress_perm_every` permutations.
#'
#' @param x A `phip_data` (recommended; must carry `peptide_library` when
#'   `rank_cols` include non-`peptide_id`) or a data.frame with the necessary
#'   columns.
#' @param rank_cols Character vector of rank columns (e.g. `"peptide_id"`,
#'   `"species"`, ...). Non-peptide ranks are joined from `x$peptide_library`.
#' @param group_cols Character vector of grouping columns that define universes;
#'   pairwise contrasts are constructed within each universe.
#' @param exist_col Name of the 0/1 presence column. Default `"exist"`.
#' @param interaction Logical; if `TRUE`, also create the interaction of the
#'   first two `group_cols`. Ignored if `combine_cols` is set.
#' @param combine_cols Optional length-2 vector; if provided, only that
#'   interaction is used as `group_col`.
#' @param interaction_sep Separator for interaction values. Default `"::"`.
#' @param B_permutations Number of permutations `B`. Default `2000L`.
#' @param seed RNG seed. Default `1L`.
#' @param smooth_eps_num Laplace numerator epsilon. Default `0.5`.
#' @param smooth_eps_den_mult Multiplier for denominator epsilon. Default `2.0`.
#' @param min_max_prev Keep peptides with `max(p1, p2) >=` this. Default `0.0`.
#' @param weight_mode One of `c("equal","se_invvar","n_eff_sqrt")`.
#' @param stat_mode One of `c("diff","asin")`.
#' @param prev_strat One of `c("none","decile")` for prevalence-stratified
#'   combining.
#' @param winsor_z Winsorization threshold for per-peptide z. Default `4.0`.
#' @param parallel Reserved (not used internally). Default `NULL`.
#' @param collect Keep `TRUE`; materialization is not required here.
#' @param register_name Reserved; no-op here.
#' @param progress_perm_every Print a small log every k permutations. Default `1L`.
#' @param rank_feature_keep Optional **named list** `rank -> values to keep`;
#'   filters `(rank, feature)` after pivot.
#'
#' @return A tibble with one row per tested stratum:
#'   `rank, feature, group_col, group1, group2, design, n_subjects_paired,
#'   n_peptides_used, T_obs, p_perm, p_adj_rank, mean_delta, frac_delta_pos,
#'   mean_delta_w, frac_delta_pos_w, category_rank_bh`.
#'
#' @export
ph_prevalence_shift2 <- function(
    x, rank_cols, group_cols,
    exist_col            = "exist",
    interaction          = FALSE,
    combine_cols         = NULL,
    interaction_sep      = "::",
    # permutation options
    B_permutations       = 2000L,
    seed                 = 1L,
    smooth_eps_num       = 0.5,
    smooth_eps_den_mult  = 2.0,
    min_max_prev         = 0.0,
    weight_mode          = c("equal","se_invvar","n_eff_sqrt"),
    stat_mode            = c("diff","asin"),
    prev_strat           = c("none","decile"),
    winsor_z             = 4.0,
    # parallel
    parallel             = NULL,   # integer cores; NULL => auto (1 by default)
    rank_feature_keep    = NULL,
    peptide_library      = NULL,
    auto_fetch_library   = FALSE,
    # logging
    log                  = FALSE,
    log_file             = "ph_prevalence_shift.log",
    # STREAM/LOG — new arguments
    stream_path          = NULL,          # path to RDS multi-object stream (serialize() appends)
    return_results       = TRUE,          # read back the stream and return tibble
    append_stream        = FALSE,          # if TRUE, don't re-init header in stream
    fold_change          = c("none","sum","mean","max","median"),
    cross_prev           = c("none","sum","mean","max","median")
) {
  # --- 0) Argument validation ----------------------------------------------------
  chk::chk_character(rank_cols)
  chk::chk_true(length(rank_cols) >= 1)
  chk::chk_character(group_cols)
  chk::chk_true(length(group_cols) >= 1)
  chk::chk_string(exist_col)
  chk::chk_number(B_permutations)
  chk::chk_true(B_permutations >= 100)
  fold_change <- match.arg(fold_change)
  cross_prev  <- match.arg(cross_prev)

  # --- 1) Prepare data once ------------------------------------------------------
  # Required columns from `x`
  need_cols <- c("sample_id","subject_id","peptide_id", exist_col, group_cols)
  if (!identical(fold_change, "none")) need_cols <- c(need_cols, "fold_change")
  if (inherits(x, "phip_data")) {
    df_long <- x$data_long |>
      dplyr::select(tidyselect::any_of(need_cols))
  } else {
    chk::chk_data(x)
    miss <- setdiff(need_cols, colnames(x))
    if (length(miss)) .ph_abort("Missing required columns", bullets = paste("-", miss))
    df_long <- tibble::as_tibble(x) |>
      dplyr::select(tidyselect::any_of(need_cols))
  }

  # --- STRICT HITS GUARD: at most one positive per (subject_id, peptide_id, group value) ---
  dup_pos <- df_long |>
    dplyr::filter(!!rlang::sym(exist_col) > 0L) |>
    tidyr::pivot_longer(tidyselect::all_of(group_cols),
      names_to = "group_col", values_to = "group_value"
    ) |>
    dplyr::count(subject_id, peptide_id, group_col, group_value, name = "n_pos") |>
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
      "Invalid input: duplicate positives within the SAME group for some (subject_id, peptide_id). One positive per group is allowed (paired designs can have up to 2 across groups).",
      bullets = c(eg, if (nrow(dup_pos) > 10) sprintf("... and %d more.", nrow(dup_pos) - 10))
    )
  }


  # --- SUBJECT METADATA (1 row per subject; lazy-safe) --------------------------
  subjects_meta <- df_long |>
    dplyr::group_by(subject_id) |>
    dplyr::summarise(dplyr::across(tidyselect::all_of(group_cols), dplyr::first),
      .groups = "drop"
    ) |>
    dplyr::arrange(subject_id) |>
    dplyr::collect()

  # --- FREEZE ORDERS (lazy-safe) ------------------------------------------------
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

  # --- BUILD HITS DIRECTLY FROM DISTINCT POSITIVES ------------------------------
  pos_pairs <- df_long |>
    dplyr::filter(!!rlang::sym(exist_col) > 0L) |>
    dplyr::distinct(subject_id, peptide_id) |>
    dplyr::collect()

  subj_index <- match(pos_pairs$subject_id, subjects_order)
  pep_index <- match(pos_pairs$peptide_id, peptides_order)

  hits_dt <- data.table::data.table(pep = pep_index, subj = subj_index)
  hits_split <- split(hits_dt$subj, hits_dt$pep)

  nonempty_pep_ids <- as.integer(names(hits_split))
  if (length(nonempty_pep_ids) == 0L) .ph_abort("All peptides are zero after hits guard (no positives).")
  peptides_order <- peptides_order[nonempty_pep_ids]

  m <- length(peptides_order)
  hits_by_peptide <- vector("list", m)
  hits_by_peptide[seq_along(nonempty_pep_ids)] <- lapply(hits_split, as.integer)

  # --- BUILD 64-bit BITSET ONCE --------------------------------------------------
  N_subjects <- length(subjects_order)
  bs <- build_bitset_unpaired(hits_by_peptide, N_subjects)

  bitset_raw <- bs$data
  bitset_m <- bs$m
  bitset_words <- bs$n_words
  if (bitset_m != length(peptides_order)) .ph_abort("Bitset/peptide dimension mismatch.")

  # Ready for downstream steps:
  # subjects_order, peptides_order, subjects_meta, bitset_raw, bitset_m, bitset_words


  # Objects ready for downstream use:
  # - subjects_order: character/ID vector (row map)
  # - peptides_order: character/ID vector (col map; all-zero peptides removed)
  # - subjects_meta : per-subject metadata (for contrasts)
  # - bitset_raw, bitset_m, bitset_words : packed 64-bit bitset (col-major)

  # --- 2) Peptide library attach (simple, in R) ---------------------------------
  ranks_need <- setdiff(rank_cols, "peptide_id")

  get_lib_tbl <- function() {
    if (length(ranks_need) == 0L) {
      # Only peptide_id requested → trivial long map
      return(
        tibble::tibble(peptide_id = peptides_order, rank = "peptide_id", feature = peptides_order)
      )
    }

    lib_src <- NULL
    if (inherits(x, "phip_data") && !is.null(x$peptide_library)) {
      lib_src <- x$peptide_library
    } else if (!is.null(peptide_library)) {
      lib_src <- peptide_library
    } else if (isTRUE(auto_fetch_library)) {
      if (!rlang::is_installed("phiper") || !("get_peptide_meta" %in% getNamespaceExports("phiper"))) {
        .ph_abort("auto_fetch_library=TRUE but phiper::get_peptide_meta() is not available.")
      }
      lib_src <- phiper::get_peptide_meta()
    } else {
      .ph_abort(
        "Peptide library required for non-peptide ranks.",
        bullets = c(
          "- Provide `peptide_library` with the needed columns,",
          "- Or set `auto_fetch_library = TRUE`,",
          "- Or include those rank columns in `x`."
        )
      )
    }

    # Select needed columns and collect to R (handles DuckDB/lazy)
    lib_needed <- c("peptide_id", ranks_need)
    miss_tax <- setdiff(ranks_need, colnames(lib_src))
    if (length(miss_tax)) .ph_abort("Peptide library missing required columns:", bullets = paste("-", miss_tax))

    lib_small <- lib_src |>
      dplyr::select(tidyselect::any_of(lib_needed)) |>
      dplyr::distinct()

    if (inherits(lib_small, "tbl_sql") || inherits(lib_small, "tbl_lazy")) {
      lib_small <- lib_small |> dplyr::collect()
    }

    # Keep only peptides we actually retained (peptides_order), then pivot to long
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
        tibble::tibble(peptide_id = peptides_order, rank = "peptide_id", feature = peptides_order)
      )
    }

    rank_map_long
  }

  rank_map_long <- get_lib_tbl()
  # rank_map_long columns: peptide_id, rank, feature  (in R memory)

  # --- 3) Initialize streamer ----------------------------------------------------
  if (!is.null(stream_path) && !isTRUE(append_stream)) {
    header <- list(
      type = "ph_prevalence_shift_stream_header",
      version = 1L,
      created = as.character(Sys.time()),
      engine = "CPP",
      args = list(
        weight_mode     = weight_mode,
        stat_mode       = stat_mode,
        prev_strat      = prev_strat,
        B_permutations  = B_permutations,
        min_max_prev    = min_max_prev,
        winsor_z        = winsor_z
      ),
      n_subjects = length(subjects_order),
      n_peptides_total = length(peptides_order),
      n_contrasts = NA_integer_ # will be filled after we construct universes
    )
    .stream_init(stream_path, header)
  }

  progress_path <- .progress_file(log_file, stream_path)
  if (isTRUE(log)) .progress_init(progress_path)

  # --- 4) Construct grouping universes (contrasts) -------------------------------

  # (a) Filter rank_map_long by rank_feature_keep (if provided)
  if (!is.null(rank_feature_keep) && length(rank_feature_keep)) {
    for (rk in names(rank_feature_keep)) {
      vals <- rank_feature_keep[[rk]]
      if (is.null(vals)) {
        # NULL => brak filtrowania dla tej rangi
        next
      } else {
        keep_vals <- as.character(vals)
        rank_map_long <- rank_map_long %>%
          dplyr::filter(!(rank == rk) | feature %in% keep_vals)
      }
    }
  }

  # (b) rank-feature counts after peptide filtering
  rf_counts <- rank_map_long |>
    dplyr::count(rank, feature, name = "n_peptides")

  # (c) levels per group_col from subjects_meta
  lvl <- subjects_meta |>
    tidyr::pivot_longer(tidyselect::all_of(group_cols),
      names_to = "group_col", values_to = "group_value"
    ) |>
    dplyr::distinct(group_col, group_value)

  # (d) all ordered pairs per group_col (no many-to-many warning)
  all_pairs <- lvl %>%
    dplyr::group_by(group_col) %>%
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
  subj_groups <- df_long |>
    tidyr::pivot_longer(tidyselect::all_of(group_cols),
      names_to = "group_col", values_to = "group_value"
    ) |>
    dplyr::distinct(subject_id, group_col, group_value) |>
    dplyr::collect()

  pairs_shared <- subj_groups |>
    dplyr::rename(group1 = group_value) |>
    dplyr::inner_join(subj_groups |> dplyr::rename(group2 = group_value),
      by = c("subject_id", "group_col")
    ) |>
    dplyr::filter(.data$group1 < .data$group2) |>
    dplyr::count(group_col, group1, group2, name = "n_shared")

  contrasts_base <- all_pairs |>
    dplyr::left_join(pairs_shared, by = c("group_col", "group1", "group2")) |>
    dplyr::mutate(design = dplyr::if_else(!is.na(.data$n_shared) & .data$n_shared > 0L, "paired", "unpaired")) |>
    dplyr::select(group_col, group1, group2, design)

  # (f) master plan = rank-feature strata × contrasts (only strata with peptides)
  master_plan <- tidyr::crossing(
    rf_counts |> dplyr::select(rank, feature, n_peptides),
    contrasts_base
  ) |>
    dplyr::filter(.data$n_peptides > 0L)

  # --- 5) Work weighting (by usable peptides) -----------------------------------
  master_plan <- master_plan |>
    dplyr::mutate(work_weight = .data$n_peptides)

  # --- 7) Header line (CPP-only) -------------------------------------------------
  n_contrasts <- nrow(master_plan)
  B_num <- as.numeric(B_permutations) # <- double
  pep_sum <- sum(as.numeric(master_plan$n_peptides), na.rm = TRUE) # <- double
  work_total <- pep_sum * B_num # <- double (no as.integer!)
  peps_total <- length(peptides_order)
  hdr_fast <- "CPP"
  log_to_file <- isTRUE(log) && (!is.null(log_file) && nzchar(log_file))

  if (isTRUE(log)) {
    hdr_bullets <- c(
      sprintf("total contrasts: %d", n_contrasts),
      sprintf("B per contrast: %d", B_permutations),
      sprintf("peptides after filtering: %d", peps_total),
      sprintf("weighted work total (sum B*peptides): %.0f", work_total), # <- %.0f
      sprintf("mode: %s", if (!is.null(parallel) && parallel > 1L) "parallel (requested)" else "auto/sequential"),
      sprintf("engine: %s", hdr_fast)
    )
    if (log_to_file) {
      .ph_log_info_file(log_file, "Starting permutation computations", bullets = hdr_bullets)
    } else {
      .ph_log_info("Starting permutation computations", bullets = hdr_bullets)
    }
  }


  # Prepare result container only if not streaming
  result_rows <- if (is.null(stream_path)) list() else NULL

  # --- Prep indices for CPP ------------------------------------------------------
  # subject row index map (1..N)
  subj_row_map <- tibble::tibble(subject_id = subjects_order, row = seq_along(subjects_order))

  # group membership with subject row indices (lazy-safe collect already done earlier)
  subj_groups_idx <- df_long |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(group_cols),
      names_to = "group_col", values_to = "group_value"
    ) |>
    dplyr::distinct(subject_id, group_col, group_value) |>
    dplyr::collect() |>
    dplyr::left_join(subj_row_map, by = "subject_id")

  # peptide column index map (within bitset columns 1..m)
  pep_col_map <- setNames(seq_along(peptides_order), peptides_order)

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
    s1 <- subj_groups_idx |>
      dplyr::filter(group_col == gc, group_value == g1) |>
      dplyr::pull(subject_id)
    s2 <- subj_groups_idx |>
      dplyr::filter(group_col == gc, group_value == g2) |>
      dplyr::pull(subject_id)
    subj_pair <- sort(intersect(s1, s2))
    if (!length(subj_pair)) {
      return(NULL)
    }
    id_map <- setNames(seq_along(subj_pair), subj_pair)

    # build per-peptide integer indices (1..P) for g1 and g2 in the SAME subject order
    g1_hits <- df_long |>
      dplyr::filter(!!rlang::sym(exist_col) > 0L, .data$subject_id %in% subj_pair) |>
      dplyr::select(subject_id, peptide_id, dplyr::all_of(gc)) |>
      dplyr::filter(.data[[gc]] == g1) |>
      dplyr::mutate(idx = id_map[.data$subject_id]) |>
      dplyr::distinct(peptide_id, idx) |>
      dplyr::group_by(peptide_id) |>
      dplyr::summarise(idx_hits = list(as.integer(sort(idx))), .groups = "drop")

    g2_hits <- df_long |>
      dplyr::filter(!!rlang::sym(exist_col) > 0L, .data$subject_id %in% subj_pair) |>
      dplyr::select(subject_id, peptide_id, dplyr::all_of(gc)) |>
      dplyr::filter(.data[[gc]] == g2) |>
      dplyr::mutate(idx = id_map[.data$subject_id]) |>
      dplyr::distinct(peptide_id, idx) |>
      dplyr::group_by(peptide_id) |>
      dplyr::summarise(idx_hits = list(as.integer(sort(idx))), .groups = "drop")

    # align peptide sets and keep only peptides present in either group
    hh <- dplyr::full_join(g1_hits, g2_hits, by = "peptide_id", suffix = c("_g1", "_g2")) |>
      dplyr::mutate(
        idx_hits_g1 = purrr::map(idx_hits_g1, ~ .x %||% integer(0)),
        idx_hits_g2 = purrr::map(idx_hits_g2, ~ .x %||% integer(0))
      ) |>
      dplyr::arrange(peptide_id)

    list(
      subj_pair = subj_pair,
      hits_g1 = hh$idx_hits_g1,
      hits_g2 = hh$idx_hits_g2,
      pep_ids = hh$peptide_id
    )
  }

  # --- core runner for one contrast (calls CPP helper; streams + logs) ----------
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

    seed_i <- as.integer(seed + i - 1L)

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
      P <- length(paired$subj_pair)
    }

    res <- cpp_shift_contrast( # <-- use exported symbol name
      bitset_raw       = bitset_raw,
      n_words          = bitset_words,
      pep_cols         = pep_cols,
      g1_rows          = g1_rows,
      g2_rows          = g2_rows,
      hits_g1_paired   = hits_g1, # empty list for unpaired
      hits_g2_paired   = hits_g2, # empty list for unpaired
      P                = as.integer(P), # 0 for unpaired
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

    # ---- optional fold_change summary (lazy-safe) ----
    fc_val <- NA_real_
    if (!identical(fold_change, "none")) {
      # peptide IDs in this stratum:
      pep_ids_here <- names(pep_col_map)[pep_cols]

      # subjects in this contrast (both groups):
      subj_rows_both <- c(g1_rows, g2_rows)
      subj_ids_here <- subjects_order[unique(subj_rows_both)]

      # lazy filter; if fold_change column is absent, this select will drop it silently
      fc_tbl <- df_long |>
        dplyr::filter(
          .data$peptide_id %in% pep_ids_here,
          .data$subject_id %in% subj_ids_here
        ) |>
        dplyr::select(tidyselect::any_of(c("fold_change"))) # no schema probing

      # summarize only if fold_change actually flowed through
      if ("fold_change" %in% names(fc_tbl)) {
        # Prefer DB-side summarize; for 'median' fall back to R if backend lacks it
        if (fold_change == "sum") fc_tbl <- dplyr::summarise(fc_tbl, v = sum(.data$fold_change, na.rm = TRUE))
        if (fold_change == "mean") fc_tbl <- dplyr::summarise(fc_tbl, v = mean(.data$fold_change, na.rm = TRUE))
        if (fold_change == "max") fc_tbl <- dplyr::summarise(fc_tbl, v = max(.data$fold_change, na.rm = TRUE))
        if (fold_change == "median") {
          # dbplyr->DuckDB supports median; if not, collect minimal vector
          fc_try <- try(dplyr::summarise(fc_tbl, v = median(.data$fold_change, na.rm = TRUE)), silent = TRUE)
          fc_tbl <- if (inherits(fc_try, "try-error")) {
            tibble::tibble(v = stats::median(dplyr::collect(fc_tbl)$fold_change, na.rm = TRUE))
          } else {
            fc_try
          }
        }
        fc_val <- dplyr::pull(dplyr::collect(fc_tbl), v)[1] %||% NA_real_
      }
    }

    # ---- optional cross-contrast prevalence summary (lazy-safe) ----
    cp_val <- NA_real_
    if (!identical(cross_prev, "none")) {
      # peptides in this stratum and subjects in the contrast (both groups)
      pep_ids_here <- names(pep_col_map)[pep_cols]
      subj_rows_both <- c(g1_rows, g2_rows)
      subj_ids_here  <- subjects_order[unique(subj_rows_both)]

      # build per-peptide prevalence pooled over g1 ∪ g2
      cp_tbl <- df_long |>
        dplyr::filter(.data$peptide_id %in% pep_ids_here,
                      .data$subject_id %in% subj_ids_here) |>
        dplyr::mutate(.exist = !!rlang::sym(exist_col) > 0L) |>
        dplyr::group_by(.data$peptide_id) |>
        dplyr::summarise(prev = mean(.exist, na.rm = TRUE), .groups = "drop")

      # summarize prevalence vector by requested reducer
      if (cross_prev == "sum")   cp_tbl <- dplyr::summarise(cp_tbl, v = sum(.data$prev, na.rm = TRUE))
      if (cross_prev == "mean")  cp_tbl <- dplyr::summarise(cp_tbl, v = mean(.data$prev, na.rm = TRUE))
      if (cross_prev == "max")   cp_tbl <- dplyr::summarise(cp_tbl, v = max(.data$prev, na.rm = TRUE))
      if (cross_prev == "median") {
        cp_try <- try(dplyr::summarise(cp_tbl, v = median(.data$prev, na.rm = TRUE)), silent = TRUE)
        cp_tbl <- if (inherits(cp_try, "try-error")) {
          tibble::tibble(v = stats::median(dplyr::collect(cp_tbl)$prev, na.rm = TRUE))
        } else cp_try
      }
      cp_val <- dplyr::pull(dplyr::collect(cp_tbl), v)[1] %||% NA_real_
    }

    row <- dplyr::bind_cols(
      st[, c("rank", "feature", "group_col", "group1", "group2", "design")],
      tibble::tibble(
        n_subjects_paired = if (identical(st$design, "paired")) as.integer(P) else NA_integer_,
        n_peptides_used   = as.integer(res$n_peptides_used),
        m_eff             = as.numeric(res$m_eff),
        T_obs             = as.numeric(res$T_obs),
        b                 = as.integer(res$b),
        p_perm            = as.numeric(res$p_perm),
        mean_delta        = as.numeric(res$mean_delta),
        frac_delta_pos    = as.numeric(res$frac_delta_pos),
        mean_delta_w      = as.numeric(res$mean_delta_w),
        frac_delta_pos_w  = as.numeric(res$frac_delta_pos_w),
        !!paste0("fold_change_", fold_change) := fc_val,
        !!paste0("cross_prev_",  cross_prev)  := cp_val
      )
    )

    if (!is.null(stream_path)) .stream_write_object(stream_path, row)

    if (isTRUE(log)) {
      msg <- sprintf(
        "%s / %s | n_pep=%d b=%d p=%g",
        row$rank, row$feature, row$n_peptides_used, row$b, row$p_perm
      )
      if (log_to_file) .ph_log_info_file(log_file, msg) else .ph_log_info(msg)

      inc <- as.numeric(row$n_peptides_used) * B_num
      done_w <- .progress_add(progress_path, inc)

      prog <- if (is.finite(done_w) && is.finite(work_total) && work_total > 0) {
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

  # --- 7/8) Computation path: sequential vs parallel ----------------------------
  result_rows <- if (is.null(stream_path)) list() else NULL

  # decide workers
  n_workers <- 1L
  if (!is.null(parallel) && !is.na(parallel) && parallel > 1L) {
    n_workers <- min(as.integer(parallel), nrow(master_plan))
  } else if (rlang::is_installed("future") && isTRUE(tryCatch(future::nbrOfWorkers() > 1L, error = function(...) FALSE))) {
    n_workers <- min(future::nbrOfWorkers(), nrow(master_plan))
  }

  # ---- dynamic scheduling for future.apply ----
  op_old <- options(
    future.scheduling   = Inf, # każdy element osobno, dynamicznie
    future.chunk.size   = 1 # bez bundlowania kilku kontrastów na raz
  )
  on.exit(options(op_old), add = TRUE)

  # ---- avoid inner oversubscription (TBB inside RcppParallel) ----
  if (rlang::is_installed("RcppParallel")) {
    RcppParallel::setThreadOptions(numThreads = 1) # fairness > peak single-task speed
  }
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")


  if (n_workers == 1L) {
    # --- sequential
    for (i in seq_len(nrow(master_plan))) {
      row <- .run_one(i)
      if (is.null(stream_path) && !is.null(row)) result_rows[[length(result_rows) + 1L]] <- row
    }
  } else {
    # --- parallel via future.apply
    original_plan <- NULL
    plan_changed <- FALSE
    if (rlang::is_installed("future")) {
      original_plan <- future::plan()
      if (future::nbrOfWorkers() != n_workers) {
        if (.Platform$OS.type == "windows") {
          future::plan(future::multisession, workers = n_workers)
        } else {
          future::plan(future::multicore, workers = n_workers)
        }
        plan_changed <- TRUE
      }
    }
    idx <- order(master_plan$work_weight, decreasing = TRUE)
    res_list <- future.apply::future_lapply(
      X = idx, FUN = .run_one,
      future.seed = TRUE,
      future.scheduling = Inf,
      future.chunk.size = 1,
      future.packages = c("dplyr", "tibble", "tidyr")
    )
    if (is.null(stream_path)) result_rows <- res_list
    if (plan_changed && !is.null(original_plan)) future::plan(original_plan)
  }

  # ---- Final log --------------------------------------------------------------
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
    if (log_to_file) .ph_log_ok_file(log_file, headline, bullets = bullets) else .ph_log_ok(headline, bullets = bullets)
  }

  # ---- Restore original future plan -------------------------------------------
  if (exists("original_plan", inherits = FALSE) && isTRUE(plan_changed) && !is.null(original_plan)) {
    try(future::plan(original_plan), silent = TRUE)
  }

  # ---- Collect results ---------------------------------------------------------
  res <- if (is.null(stream_path)) {
    # in-memory: compact & bind; tolerate all-NULL
    rr <- purrr::compact(result_rows)
    if (length(rr)) dplyr::bind_rows(rr) else tibble::tibble()
  } else if (isTRUE(return_results)) {
    .stream_read_all(stream_path) %>%
      dplyr::filter(.data$type != "ph_prevalence_shift_stream_header")
  } else {
    return(invisible(NULL))
  }

  # Return early if empty
  if (nrow(res) == 0L) {
    return(tibble::tibble(
      rank = character(), feature = character(), group_col = character(),
      group1 = character(), group2 = character(), design = character(),
      n_subjects_paired = integer(), n_peptides_used = integer(), m_eff = numeric(),
      T_obs = numeric(), p_perm = numeric(), b = integer(),
      p_adj_rank = numeric(),
      mean_delta = numeric(), frac_delta_pos = numeric(),
      mean_delta_w = numeric(), frac_delta_pos_w = numeric(),
      category_rank_bh = character()
    ))
  }

  # Ensure numeric types we rely on
  res <- res %>%
    dplyr::mutate(
      p_perm = as.numeric(p_perm),
      m_eff  = as.numeric(m_eff)
    )

  # ---- p.adjust + final select/arrange -----------------------------------------
  res <- res %>%
    dplyr::group_by(rank) %>%
    dplyr::mutate(p_adj_rank = p.adjust(p_perm, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      category_rank_bh = dplyr::case_when(
        !is.na(p_adj_rank) & p_adj_rank < 0.05 ~ "significant (BH, per rank)",
        !is.na(p_perm) & p_perm < 0.05 ~ "nominal only",
        TRUE ~ "not significant"
      ),
      # if you didn't populate it in .run_one(), keep NA for unpaired
      n_subjects_paired = dplyr::coalesce(.data$n_subjects_paired, NA_integer_)
    ) %>%
    dplyr::select(
      rank, feature, group_col, group1, group2, design,
      n_subjects_paired, n_peptides_used, m_eff,
      T_obs, p_perm, b,
      p_adj_rank,
      mean_delta, frac_delta_pos, mean_delta_w, frac_delta_pos_w,
      dplyr::any_of(paste0("fold_change_", fold_change)),
      dplyr::any_of(paste0("cross_prev_",  cross_prev)),
      category_rank_bh
    ) %>%
    dplyr::arrange(rank, feature, group_col, group1, group2)

  return(res)
}
