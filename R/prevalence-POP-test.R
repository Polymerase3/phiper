# ==============================================================================
# Helpers: constructor + printer for class `ph_prev_result`
# ==============================================================================

#' build a `ph_prev_result` object
#'
#' @param data tibble/data.frame with results
#' @param meta named list with accounting info (m_by_rank, pools, pairs, etc.)
#' @return `data` with extra attributes and class `ph_prev_result`
.as_prev_result <- function(data, meta) {
  attr(data, "prev_meta") <- meta
  class(data) <- unique(c("ph_prev_result", class(data)))
  data
}

#' @export
print.ph_prev_result <- function(x, ...) {
  meta <- attr(x, "prev_meta")
  cat("<ph_prev_result>\n")
  if (!is.null(meta$paired)) {
    cat("  design      :", if (isTRUE(meta$paired)) "paired (mcnemar)" else "unpaired (fisher)", "\n")
  }
  if (!is.null(meta$weight_mode)) cat("  weights     :", meta$weight_mode, "\n")
  if (!is.null(meta$pop_k_min)) cat("  pop k-min   :", meta$pop_k_min, "\n")
  if (!is.null(meta$fdr_scope)) cat("  fdr scope   :", meta$fdr_scope, "\n")
  if (!is.null(meta$m_by_rank)) {
    cat(
      "  m per rank  : ",
      paste(sprintf("%s=%s", names(meta$m_by_rank), format(meta$m_by_rank, big.mark = ",")), collapse = "; "),
      "\n"
    )
  }
  NextMethod()
}

# ==============================================================================
# summary() for `ph_prev_result`
# ==============================================================================

#' Summarize FDR accounting and design for `ph_prev_result`
#'
#' @description
#' Provides a compact overview of the testing design and FDR accounting:
#' - per-rank pool sizes (`POOL_r`), per-universe level counts and pair counts,
#' - the number of pairwise comparisons aggregated across universes (`PAIRS`),
#' - the per-rank family sizes `m_r = POOL_r * PAIRS`,
#' - design flags (paired/unpaired), weight mode, `pop_k_min`, `fdr_scope`, and `view`.
#'
#' @param object A `ph_prev_result` object returned by `ph_prevalence_compare()`.
#' @param ... Unused, kept for S3 compatibility.
#'
#' @return An object of class `summary.ph_prev_result` containing:
#' \itemize{
#'   \item \code{overview}: tibble with one row per rank (`rank`, `POOL`, `PAIRS`, `m_r`).
#'   \item \code{pairs_by_universe}: tibble with `group_col`, `k_levels`, `n_pairs`.
#'   \item \code{pool_by_rank}: tibble with `rank`, `POOL`.
#'   \item \code{totals}: named list with `PAIRS_total` and `m_total`.
#'   \item \code{design}: named list with `paired`, `weight_mode`, `pop_k_min`, `fdr_scope`, `view`.
#' }
#' The object has a custom `print()` for nice console output.
#'
#' @examples
#' \dontrun{
#' s <- summary(res)   # where res is a ph_prev_result
#' s
#' }
#'
#' @export
#' @method summary ph_prev_result
summary.ph_prev_result <- function(object, ...) {
  if (!inherits(object, "ph_prev_result")) {
    stop("summary.ph_prev_result: object must be a 'ph_prev_result'.")
  }

  meta <- attr(object, "prev_meta") %||% list()

  pool_by_rank <- meta$pool_by_rank %||% tibble::tibble(rank = character(), POOL = integer())
  pairs_by_universe <- meta$pairs_by_universe %||% tibble::tibble(group_col = character(), k_levels = integer(), n_pairs = integer())
  m_by_rank <- meta$m_by_rank %||% setNames(integer(0), character(0))

  # compute pairs_total (sum across universes) and overview per rank
  pairs_total <- if (nrow(pairs_by_universe)) sum(pairs_by_universe$n_pairs, na.rm = TRUE) else 0L
  if (nrow(pool_by_rank)) {
    ov <- pool_by_rank |>
      dplyr::mutate(
        PAIRS = pairs_total,
        m_r = POOL * pairs_total
      ) |>
      dplyr::relocate(rank, POOL, PAIRS, m_r)
  } else {
    ov <- tibble::tibble(rank = character(), POOL = integer(), PAIRS = integer(), m_r = integer())
  }

  # totals
  m_total <- if (length(m_by_rank)) sum(as.integer(m_by_rank), na.rm = TRUE) else sum(ov$m_r, na.rm = TRUE)

  out <- list(
    overview = tibble::as_tibble(ov),
    pairs_by_universe = tibble::as_tibble(pairs_by_universe),
    pool_by_rank = tibble::as_tibble(pool_by_rank),
    totals = list(PAIRS_total = pairs_total, m_total = m_total),
    design = list(
      paired = isTRUE(meta$paired),
      weight_mode = meta$weight_mode %||% NA_character_,
      pop_k_min = meta$pop_k_min %||% NA_integer_,
      fdr_scope = meta$fdr_scope %||% "per-rank",
      view = meta$view %||% NA_character_
    )
  )
  class(out) <- "summary.ph_prev_result"
  out
}

#' @export
print.summary.ph_prev_result <- function(x, ...) {
  # header
  cat("<summary.ph_prev_result>\n")
  des <- x$design %||% list()
  cat("  design      :", if (isTRUE(des$paired)) "paired (mcnemar)" else "unpaired (fisher)", "\n")
  cat("  weights     :", des$weight_mode, "\n")
  cat("  pop k-min   :", des$pop_k_min, "\n")
  cat("  fdr scope   :", des$fdr_scope, "\n")
  if (!is.null(des$view) && !is.na(des$view)) cat("  view        :", des$view, "\n")

  # totals
  tots <- x$totals %||% list()
  cat("  pairs total :", format(tots$PAIRS_total %||% 0L, big.mark = ","), "\n")
  cat("  m total     :", format(tots$m_total %||% 0L, big.mark = ","), "\n")

  # overview table
  if (nrow(x$overview)) {
    cat("\n  per-rank overview:\n")
    print(x$overview, n = min(10L, nrow(x$overview)))
  } else {
    cat("\n  per-rank overview: <empty>\n")
  }

  # universes table (short)
  if (nrow(x$pairs_by_universe)) {
    cat("\n  universes (levels & pairs):\n")
    print(x$pairs_by_universe, n = min(10L, nrow(x$pairs_by_universe)))
  } else {
    cat("\n  universes (levels & pairs): <empty>\n")
  }

  invisible(x)
}

# ==============================================================================
# filter pairs helper for `ph_prev_result`
# ==============================================================================

#' Filter pairwise results by groups/ranks/features with optional q-value gates
#'
#' @description
#' Order-agnostic filter for rows comparing two group levels (e.g. "B" vs "M12"),
#' z opcjonalnym zawężeniem do ranków/feature'ów/uniwersów oraz progów istotności.
#' Działa zarówno na obiekcie klasy `ph_prev_result`, jak i zwykłym `data.frame`.
#'
#' @details
#' **pary grup (symetrycznie).** Dla każdej zadanej pary \{gA, gB\} dobierane są
#' wiersze, w których `(min(group1,group2), max(group1,group2))` zgadza się z
#' `sort(c(gA, gB))`, niezależnie od kolejności w danych.
#'
#' **filtrowanie.**
#' - `ranks`: wektor rang do zatrzymania.
#' - `features`: wektor lub regex; jeśli `features_regex=TRUE`, wszystkie wzorce
#'   są łączone logicznym OR (dopasowanie nieczułe na wielkość liter).
#' - `group_universe`: dokładne wartości `group_col` albo regex (gdy
#'   `universe_regex=TRUE`).
#' - progi istotności: `p_raw_max`, `q_bh_max`, `q_wbh_max` oraz `passed_only`
#'   (korzysta z kolumn `passed_rank_bh`/`passed_rank_wbh`, jeśli są dostępne).
#'
#' **klasa.** Jeżeli wejście ma klasę `ph_prev_result`, wyjście zachowuje ją oraz
#' atrybut `prev_meta` (dodaje `subsetted = TRUE`, oraz ewentualną notkę o filtrze).
#'
#' @param df A `ph_prev_result` or `data.frame` with columns `group1`, `group2`
#'   and optionally `group_col`, `rank`, `feature`, p/q columns.
#' @param gA,gB character vector or length-1 character (pojedyncza para). Można
#'   też przekazać listę par (np. `list(c("B","M12"), c("T0","T2"))`) albo
#'   dwu-kolumnowy `data.frame`/`tibble` z nazwami kolumn `gA`, `gB`.
#' @param ranks Optional character vector of ranks to keep.
#' @param features Optional character vector or regex pattern(s).
#' @param features_regex Logical; treat `features` as regex pattern(s) (OR).
#' @param group_universe Optional character vector or regex pattern(s) for `group_col`.
#' @param universe_regex Logical; treat `group_universe` as regex pattern(s) (OR).
#' @param col_rank,col_feature,col_g1,col_g2,col_groupcol Column names in `df`.
#' @param p_raw_max,q_bh_max,q_wbh_max Optional numeric thresholds (e.g. 0.05).
#' @param passed_only Logical; if TRUE, keep rows with any `passed_* == TRUE` found.
#' @param drop_na Logical; drop rows with NA in group1/group2 before pair match.
#' @param keep_cols Optional character vector of columns to retain (NULL = all).
#'
#' @return Filtered object of the same base type as `df`; if input is a
#'   `ph_prev_result`, output keeps the class and augmented metadata.
#'
#' @examples
#' \dontrun{
#' # single pair, regex on features:
#' # prev_filter_pairs(res, "B", "M12", ranks = "species",
#' #                   features = "flagellin|fliC", features_regex = TRUE,
#' #                   q_bh_max = 0.1)
#'
#' # multiple pairs:
#' # pairs <- list(c("B","M12"), c("T0","T2"))
#' # prev_filter_pairs(res, pairs[[1]][1], pairs[[1]][2], ranks="peptide_id")
#' }
#'
#' @export
#' Filter pairwise results by groups/ranks/features with optional q-value gates
#' (robust row subsetting; works whether df is data.table or data.frame)
#' @export
prev_filter_pairs <- function(
  df,
  gA, gB,
  ranks = NULL,
  features = NULL,
  features_regex = FALSE,
  group_universe = NULL,
  universe_regex = FALSE,
  col_rank = "rank",
  col_feature = "feature",
  col_g1 = "group1",
  col_g2 = "group2",
  col_groupcol = "group_col",
  p_raw_max = NULL,
  q_bh_max = NULL,
  q_wbh_max = NULL,
  passed_only = FALSE,
  drop_na = TRUE,
  keep_cols = NULL
) {
  # -- normalize pairs input (allow list/data.frame of pairs) -----------------
  normalize_pairs <- function(gA, gB) {
    if (is.character(gA) && is.character(gB) && length(gA) == 1L && length(gB) == 1L) {
      return(list(c(gA, gB)))
    }
    if (is.list(gA) && missing(gB)) {
      ok <- vapply(gA, function(x) is.character(x) && length(x) == 2L, logical(1))
      if (!all(ok)) stop("prev_filter_pairs: all list elements must be length-2 character vectors.")
      return(gA)
    }
    if (is.data.frame(gA) && missing(gB)) {
      if (!all(c("gA", "gB") %in% names(gA))) stop("prev_filter_pairs: data.frame must have columns 'gA' and 'gB'.")
      return(split(as.matrix(gA[, c("gA", "gB"), drop = FALSE]), seq_len(nrow(gA))))
    }
    stop("prev_filter_pairs: provide either scalars gA,gB; a list of pairs; or a data.frame with columns gA,gB.")
  }
  pairs_list <- normalize_pairs(gA, gB)

  # -- basic column checks ----------------------------------------------------
  needed <- c(col_g1, col_g2)
  if (!all(needed %in% names(df))) stop("prev_filter_pairs: df is missing required columns: ", paste(setdiff(needed, names(df)), collapse = ", "))

  # no mutation of original
  dt <- data.table::as.data.table(df)

  # --- optional rank pre-filter (robust to scoping) --------------------------
  if (!is.null(ranks)) {
    if (!col_rank %in% names(dt)) stop("prev_filter_pairs: column '", col_rank, "' not found.")
    ranks <- as.character(ranks)
    col_vec <- as.character(dt[[col_rank]])
    mask <- col_vec %in% ranks
    dt <- dt[mask, , drop = FALSE]
  }

  # --- optional feature pre-filter (safe for factors/NA) ---------------------
  if (!is.null(features)) {
    if (!col_feature %in% names(dt)) stop("prev_filter_pairs: column '", col_feature, "' not found.")
    fvec <- as.character(dt[[col_feature]])
    fvec[is.na(fvec)] <- ""
    if (isTRUE(features_regex)) {
      pats <- as.character(features)
      keep <- rep(FALSE, nrow(dt))
      for (p in pats) keep <- keep | grepl(p, fvec, ignore.case = TRUE, perl = TRUE)
      dt <- dt[keep, , drop = FALSE]
    } else {
      dt <- dt[tolower(fvec) %in% tolower(as.character(features)), , drop = FALSE]
    }
  }

  # --- optional universe filter (group_col exact or regex or) ----------------
  if (!is.null(group_universe) && col_groupcol %in% names(dt)) {
    gvec <- as.character(dt[[col_groupcol]])
    gvec[is.na(gvec)] <- ""
    if (isTRUE(universe_regex)) {
      pats <- as.character(group_universe)
      keep <- rep(FALSE, nrow(dt))
      for (p in pats) keep <- keep | grepl(p, gvec, ignore.case = TRUE, perl = TRUE)
      dt <- dt[keep, , drop = FALSE]
    } else {
      dt <- dt[tolower(gvec) %in% tolower(as.character(group_universe)), , drop = FALSE]
    }
  }

  # --- order-agnostic pair match (two comparisons per row) -------------------
  c1 <- dt[[col_g1]]
  c2 <- dt[[col_g2]]
  ok <- if (drop_na) (!is.na(c1) & !is.na(c2)) else rep(TRUE, nrow(dt))
  lo <- data.table::fifelse(c1 <= c2, c1, c2)
  hi <- data.table::fifelse(c1 <= c2, c2, c1)
  keep_any <- rep(FALSE, nrow(dt))
  for (pp in pairs_list) {
    tgt <- sort(as.character(pp))
    keep_any <- keep_any | (ok & lo == tgt[1] & hi == tgt[2])
  }
  dt <- dt[keep_any, , drop = FALSE]

  # --- optional significance filters -----------------------------------------
  if (!is.null(p_raw_max) && "p_raw" %in% names(dt)) {
    mask <- is.na(dt[["p_raw"]]) | dt[["p_raw"]] <= p_raw_max
    dt <- dt[mask, , drop = FALSE]
  }
  if (!is.null(q_bh_max) && "p_adj_rank" %in% names(dt)) {
    mask <- is.na(dt[["p_adj_rank"]]) | dt[["p_adj_rank"]] <= q_bh_max
    dt <- dt[mask, , drop = FALSE]
  }
  if (!is.null(q_wbh_max) && "p_adj_rank_wbh" %in% names(dt)) {
    mask <- is.na(dt[["p_adj_rank_wbh"]]) | dt[["p_adj_rank_wbh"]] <= q_wbh_max
    dt <- dt[mask, , drop = FALSE]
  }
  if (isTRUE(passed_only)) {
    pass_cols <- intersect(c("passed_rank_bh", "passed_rank_wbh"), names(dt))
    if (length(pass_cols)) {
      pass_any <- Reduce(`|`, lapply(pass_cols, function(cc) data.table::fcoalesce(dt[[cc]], FALSE)))
      dt <- dt[pass_any, , drop = FALSE]
    }
  }

  # --- select columns if requested -------------------------------------------
  if (!is.null(keep_cols)) {
    keep_cols <- intersect(keep_cols, names(dt))
    dt <- dt[, ..keep_cols]
  }

  res <- data.table::setDF(dt)

  # --- preserve class/metadata for ph_prev_result ----------------------------
  if (inherits(df, "ph_prev_result")) {
    meta <- attr(df, "prev_meta")
    if (is.null(meta)) meta <- list()
    meta$subsetted <- TRUE
    meta$subset_n <- nrow(res)
    attr(res, "prev_meta") <- meta
    class(res) <- class(df)
  }

  res
}

# ==============================================================================
#' Prevalence by group with pairwise tests (POP; per-rank FDR; BH & weighted BH)
#'
#' @description
#' Computes prevalence (counts, proportions, percentages) for features defined by
#' the requested `rank_cols` across one or more grouping "universes" (`group_cols`,
#' optional interaction), and performs **pairwise** statistical tests between all
#' group levels. Presence is computed with a k-of-n POP rule within each sample
#' (`pop_k_min`, default 1). Two testing modes are supported:
#' - **Unpaired:** Fisher’s exact test (2×2) per (rank, feature, group pair).
#' - **Paired:** McNemar’s exact (binomial) per subject (`paired = TRUE` requires
#'   `subject_id`).
#'
#' P-values are adjusted **per rank** (single FDR family for each rank across all
#' requested universes/pairs/features) using BH and optional weighted BH
#' (`weight_mode = "peptide_count"`).
#'
#' @details
#' **universe construction.** For each `group_col` we create a universe of levels
#' present in the data (non-missing). If `interaction = TRUE` (or `combine_cols` is
#' provided), we also build a combined universe where the group value is
#' `<col1>::<col2>`. Denominators `N` are distinct sample counts per
#' `(group_col, group_value)`.
#'
#' **presence (pop).** For each sample and each `(rank, feature)`, we count the
#' number of positive peptides contributing to that feature; a feature is marked
#' present if `k >= pop_k_min` (default 1). This yields `n_present` per
#' `(group_col, group_value, rank, feature)` and prevalence `prop = n_present / N`.
#'
#' **pairwise tests.** For each universe with `K` levels we form all unordered
#' pairs `K*(K-1)/2`. For each `(rank, feature, pair)`:
#' - unpaired: build the 2×2 table from (`n1`, `N1 - n1`, `n2`, `N2 - n2`) and run
#'   Fisher’s exact test (two-sided).
#' - paired: compute discordant counts `(n01, n10)` per subject and run
#'   McNemar’s exact binomial test (`binom.test(n01, n01+n10)`).
#'
#' **fdr families and the number of comparisons.**
#' Let `POOL_r` be the number of unique features for a given **rank** `r`
#' **across all requested universes** (after presence aggregation and any missing
#' value filtering). Let `PAIRS = sum_u K_u*(K_u-1)/2` be the total number of
#' unordered level-pairs summed over all universes `u` (each `u` is one
#' `group_col` or an interaction universe). Then the size of the **single FDR
#' family** for rank `r` is:
#'
#' \deqn{ m_r \;=\; POOL_r \times PAIRS. }
#'
#' All p-values produced for rank `r` (covering **all** universes and **all**
#' level pairs) are adjusted together as one family of size `m_r`. This makes the
#' per-rank FDR **stricter** when (i) there are many features for that rank, and/or
#' (ii) many universes or levels generate many pairwise comparisons.
#'
#' **bh (benjamini–hochberg).** Within each rank `r` (and `view` if present),
#' BH q-values are computed via `p.adjust(method="BH")` on the vector of p-values
#' of length `m_r` (excluding NAs). Reported columns:
#' `p_adj_rank`, `passed_rank_bh`, `category_rank_bh`.
#'
#' **weighted bh.** If `weight_mode="peptide_count"`, each `(rank,feature)` gets a
#' base weight equal to the number of distinct peptides mapping to that feature.
#' Within each rank `r` (and `view` if present), **weights are scaled to sum to**
#' `m_r`:
#'
#' \deqn{ w_i^\* \;=\; w_i \cdot \frac{m_r}{\sum_j w_j}. }
#'
#' We adjust using the standard weighted step-up rule on `p_i / w_i^\*`. The
#' resulting q-values are reported in `p_adj_rank_wbh`, with flags
#' `passed_rank_wbh` and labels `category_rank_wbh`. If `weight_mode="none"`,
#' all weights are 1 and wBH reduces to BH.
#'
#' **logging of comparisons.** The function logs:
#' - `POOL_r` per rank,
#' - number of levels `K_u` and pair counts per universe,
#' - `PAIRS = sum_u K_u*(K_u-1)/2`,
#' - `m_r = POOL_r * PAIRS` for each rank.
#' These values are also returned in the object metadata (see Value).
#'
#' @param x A `phip_data` object with DuckDB backend **or** a data.frame/tibble
#'   with at least: `sample_id`, `peptide_id`, `exist_col`, all `group_cols`, and
#'   if `paired=TRUE` also `subject_id`. If `rank_cols` include non-peptide taxa,
#'   `x` must provide `peptide_library` with those columns.
#' @param rank_cols Character vector of rank columns, e.g. `c("peptide_id","species")`.
#' @param group_cols Character vector of grouping columns defining universes.
#' @param exist_col Name of the binary presence column (default `"exist"`).
#' @param weight_mode `"peptide_count"` (default) or `"none"`.
#' @param parallel Logical; compute Fisher p-values in parallel if possible (default `NULL` = auto).
#' @param compute_ratios_db Logical; compute simple ratios/delta ratios in SQL (unpaired only).
#' @param interaction Logical; also create an interaction universe of the first two `group_cols`.
#' @param combine_cols Optional length-2 character vector to build only that interaction universe.
#' @param interaction_sep Separator for interaction labels (default `"::"`).
#' @param collect Logical; if `TRUE` (default) return a collected tibble; otherwise a lazy table.
#' @param register_name Optional DuckDB table name for materialization (unpaired path).
#' @param pop_k_min Integer ≥1; k-of-n POP threshold per sample (default 1).
#' @param paired Logical; use paired design (McNemar exact) with `subject_id` (default `FALSE`).
#'        NOTE: can also be a character scalar naming the column that links related samples
#'        (e.g. "subject_id" or "dyade"). If so, only samples present in both groups
#'        for that identifier will be used for paired McNemar tests.
#'
#' @return An object of class `ph_prev_result`, i.e., a tibble (or lazy table if
#'   `collect = FALSE` on the unpaired path) with attributes:
#'   - `prev_meta$m_by_rank`: named integer vector of `m_r` per rank,
#'   - `prev_meta$pairs_by_universe`: tibble with `group_col`, `k_levels`, `n_pairs`,
#'   - `prev_meta$pool_by_rank`: tibble with `rank`, `POOL`,
#'   - other bookkeeping: `paired`, `weight_mode`, `pop_k_min`, `fdr_scope = "per-rank"`,
#'   - `register_name` (unpaired path) and `view` (if available).
#'
#' Columns include (subset may differ between paths): `view`, `rank`, `feature`,
#' `n_peptides`, `group_col`, `group1`, `group2`, `n1`, `N1`, `prop1`, `percent1`,
#' `n2`, `N2`, `prop2`, `percent2`, optional `ratio`, `delta_ratio`, `p_raw`,
#' `p_adj_rank`, `passed_rank_bh`, `category_rank_bh`, `p_adj_rank_wbh`,
#' `passed_rank_wbh`, `category_rank_wbh`.
#'
#' @examples
#' \dontrun{
#' res <- ph_prevalence_compare(pd, rank_cols=c("species"), group_cols=c("big_group"))
#' print(res)
#' }
#'
#' @export
# ==============================================================================
ph_prevalence_compare <- function(x,
                                  rank_cols,
                                  group_cols,
                                  exist_col = "exist",
                                  weight_mode = c("peptide_count", "none"),
                                  parallel = NULL,
                                  compute_ratios_db = TRUE,
                                  interaction       = FALSE,
                                  combine_cols      = NULL,
                                  interaction_sep   = "::",
                                  collect           = TRUE,
                                  register_name     = NULL,
                                  pop_k_min         = 1L,
                                  paired            = FALSE,
                                  peptide_library   = NULL) {

  weight_mode <- match.arg(weight_mode)

  .ph_with_timing(
    headline = "prevalence_compare (per-rank fdr)",
    step = NULL,
    bullets = NULL,
    expr = {
      .q <- function(con, nm) as.character(DBI::dbQuoteIdentifier(con, nm))
      .sym <- rlang::sym
      chunk_n <- getOption("phiper.prev.chunk", 1e6L)

      # --- choose library handle: explicit arg > x$peptide_library > x$meta$peptide_library
      lib_handle <- peptide_library %||%
        tryCatch(x$peptide_library, error = function(...) NULL) %||%
        tryCatch(x$meta$peptide_library, error = function(...) NULL)

      # ---- basic checks / logging -------------------------------------------
      tryCatch(
        {
          chk::chk_character(rank_cols)
          chk::chk_true(length(rank_cols) >= 1)
          chk::chk_character(group_cols)
          chk::chk_true(length(group_cols) >= 1)
          chk::chk_string(exist_col)
          chk::chk_number(pop_k_min)
          chk::chk_true(pop_k_min >= 1)
        },
        error = function(e) .ph_abort("invalid arguments", bullets = e$message)
      )

      # paired: FALSE or single column name
      paired_col <- NULL
      if (!identical(paired, FALSE)) {
        if (!chk::vld_string(paired)) .ph_abort("paired must be FALSE or a single column name (string).")
        paired_col <- as.character(paired)
      }

      .ph_log_info(
        "preparing input data",
        bullets = c(
          paste0("ranks: ", paste(rank_cols, collapse = ", ")),
          paste0("group_cols: ", paste(group_cols, collapse = ", ")),
          if (!is.null(combine_cols)) {
            paste0("combine_cols: ", paste(combine_cols, collapse = " + "))
          } else if (interaction) "interaction: true (additive to per-column)",
          paste0("exist_col: ", exist_col),
          paste0("weight_mode: ", weight_mode),
          paste0("collect: ", collect),
          paste0("pop_k_min: ", pop_k_min),
          paste0("paired: ", ifelse(is.null(paired_col), "FALSE", paired_col))
        )
      )

      # ---- normalize to long df ---------------------------------------------
      df_long <- try(
        {
          if (inherits(x, "phip_data")) {
            x$data_long |>
              dplyr::select(tidyselect::any_of(c(
                "sample_id", "subject_id", "peptide_id",
                exist_col, group_cols, paired_col
              )))
          } else {
            chk::chk_data(x)
            need <- c("sample_id", "peptide_id", exist_col, group_cols)
            tibble::as_tibble(x) |>
              dplyr::select(tidyselect::any_of(c(need, paired_col)))
          }
        },
        silent = TRUE
      )
      if (inherits(df_long, "try-error")) .ph_abort("could not prepare input data.")

      df_long <- df_long |>
        dplyr::mutate(!!.sym(exist_col) := dplyr::coalesce(as.integer(!!.sym(exist_col)), 0L))

      # ---- paired column detection ------------------------------------------
      found_paired_col <- NULL
      if (!is.null(paired_col)) {
        cols_here <- tryCatch(colnames(df_long), error = function(...) character(0))
        if (paired_col %in% cols_here) {
          found_paired_col <- paired_col
        } else {
          cand <- paste0(paired_col, "_recoded")
          if (cand %in% cols_here) {
            found_paired_col <- cand
            .ph_warn(headline = sprintf("paired column '%s' not found; using '%s' instead.", paired_col, cand))
          } else {
            approx <- grep(paired_col, cols_here, value = TRUE)
            if (length(approx) == 1L) {
              found_paired_col <- approx
              .ph_warn(headline = sprintf("paired column '%s' not found; using similar column '%s'.", paired_col, approx))
            } else if (length(approx) > 1L) {
              .ph_abort(
                headline = sprintf("paired specified as '%s' but multiple matching columns found.", paired_col),
                bullets = approx
              )
            } else {
              .ph_abort(headline = sprintf(
                "paired specified as '%s' but this column not found in the input after filtering.",
                paired_col
              ))
            }
          }
        }
      }

      # ---- connection (data.frame-friendly) ---------------------------------
      con <- NULL
      if (inherits(x, "phip_data")) con <- tryCatch(x$meta$con, error = function(...) NULL)
      if (is.null(con) && inherits(df_long, "tbl_sql")) con <- dbplyr::remote_con(df_long)

      if (is.null(con)) {
        if (!requireNamespace("duckdb", quietly = TRUE)) {
          .ph_abort("no duckdb connection found and package {duckdb} is not installed.")
        }
        con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
        on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

        tmp_name <- paste0("ph_tmp_long_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        DBI::dbWriteTable(con, tmp_name, tibble::as_tibble(df_long), temporary = TRUE)
        df_long <- dplyr::tbl(con, tmp_name)
      }

      view_const <- if (inherits(x, "phip_data")) attr(x, "view") %||% (x$meta$view %||% NA_character_) else NA_character_

      # ---- ranks & library ---------------------------------------------------
      ranks_needing_lib <- setdiff(rank_cols, "peptide_id")

      # validate lib (only if needed)
      if (length(ranks_needing_lib)) {
        if (is.null(lib_handle)) .ph_abort("peptide library required for non-peptide ranks (not found).")

        # get column names safely for both df and tbl_sql
        lib_cols <- if (inherits(lib_handle, "tbl_sql")) {
          colnames(dplyr::collect(dplyr::slice_head(lib_handle, n = 0)))
        } else {
          colnames(lib_handle)
        }

        miss_tax <- setdiff(c("peptide_id", ranks_needing_lib), lib_cols)
        if (length(miss_tax)) .ph_abort("requested taxonomy columns not in peptide_library.", bullets = paste("-", miss_tax))
      }

      # ensure library is a tbl on the same connection if we need it
      lib_tbl_for_join <- NULL
      if (length(ranks_needing_lib)) {
        if (inherits(lib_handle, "tbl_sql")) {
          lib_tbl_for_join <- lib_handle
        } else {
          lib_tbl_for_join <- tibble::as_tibble(lib_handle)
          lib_name <- paste0("ph_tmp_lib_", format(Sys.time(), "%Y%m%d_%H%M%S"))
          DBI::dbWriteTable(con, lib_name, lib_tbl_for_join, temporary = TRUE)
          lib_tbl_for_join <- dplyr::tbl(con, lib_name)
        }
      }

      df_ranked <- df_long
      if (length(ranks_needing_lib)) {
        lib_cols <- c("peptide_id", ranks_needing_lib)
        lib_min <- lib_tbl_for_join |> dplyr::select(tidyselect::all_of(lib_cols)) |> dplyr::distinct()
        df_ranked <- df_ranked |> dplyr::left_join(lib_min, by = "peptide_id", copy = TRUE)
      }

      available_ranks <- intersect(rank_cols, colnames(df_ranked))
      if (!length(available_ranks)) .ph_abort("none of the requested rank_cols are available.")
      .ph_log_info("ranks resolved", bullets = paste("- available:", paste(available_ranks, collapse = ", ")))

      df_ranked_long <- df_ranked %>%
        tidyr::pivot_longer(
          cols = tidyselect::all_of(available_ranks),
          names_to = "rank", values_to = "feature"
        ) %>%
        dplyr::mutate(feature = as.character(feature))

      # ---- grouping universes -----------------------------------------------
      make_interaction <- function(tbl, c1, c2, sep) {
        comb_name <- paste(c1, c2, sep = " + ")
        tbl |>
          dplyr::filter(!is.na(.data[[c1]]) & !is.na(.data[[c2]])) |>
          dplyr::mutate(
            group_col = comb_name,
            group_value = paste0(.data[[c1]], sep, .data[[c2]])
          )
      }

      per_column <- df_ranked_long |>
        tidyr::pivot_longer(
          cols = tidyselect::all_of(group_cols),
          names_to = "group_col", values_to = "group_value"
        ) |>
        dplyr::filter(!is.na(group_value))

      if (!is.null(combine_cols)) {
        gs_view <- make_interaction(df_ranked_long, combine_cols[1], combine_cols[2], interaction_sep)
        .ph_log_info("grouping universes",
          bullets = c(paste0("- only interaction of: ", paste(combine_cols, collapse = " + ")))
        )
      } else if (isTRUE(interaction)) {
        if (length(group_cols) < 2L) .ph_abort("interaction=TRUE needs at least two group_cols.")
        inter_view <- make_interaction(df_ranked_long, group_cols[1], group_cols[2], interaction_sep)
        gs_view <- dplyr::bind_rows(per_column, inter_view)
        .ph_log_info("grouping universes",
          bullets = c(
            paste0("- per-column: ", paste(group_cols, collapse = ", ")),
            paste0("- plus interaction: ", paste(group_cols[1:2], collapse = " + "))
          )
        )
      } else {
        gs_view <- per_column
        .ph_log_info("grouping universes",
          bullets = c(paste0("- per-column only: ", paste(group_cols, collapse = ", ")))
        )
      }

      # ---- cohort sizes (n) per universe ------------------------------------
      .ph_log_info("computing cohort sizes (n) per universe")
      group_sizes <- gs_view |>
        dplyr::distinct(group_col, group_value, sample_id) |>
        dplyr::count(group_col, group_value, name = "N")
      gs_n <- group_sizes |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::collect()
      if (gs_n$n == 0L) .ph_abort("no non-missing group values; cannot compute denominators.")

      # ---- k-of-n presence per (sample, rank, feature) ----------------------
      .ph_log_info("computing presence per sample via k-of-n rule")
      group_by_cols <- c("group_col", "group_value", "sample_id")
      if (!is.null(found_paired_col)) group_by_cols <- c(group_by_cols, found_paired_col)
      group_by_cols <- c(group_by_cols, "rank", "feature")
      k_tbl <- gs_view %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_by_cols))) %>%
        dplyr::summarise(k = sum(!!.sym(exist_col) > 0L), .groups = "drop") %>%
        dplyr::mutate(present = k >= !!pop_k_min)

      # ============================= non-paired path ==========================
      if (is.null(found_paired_col)) {
        .ph_log_info("counting present samples per feature (pop, non-paired)")
        present_counts <- k_tbl %>%
          dplyr::filter(present) %>%
          dplyr::distinct(group_col, group_value, rank, feature, sample_id) %>%
          dplyr::count(group_col, group_value, rank, feature, name = "n_present")
        features_per_rank <- present_counts |> dplyr::distinct(rank, feature)

        # ---- FDR accounting ---------------------------------------------------
        pool_tbl <- features_per_rank |>
          dplyr::count(rank, name = "POOL") |>
          dplyr::arrange(rank) |>
          dplyr::collect()
        lev_tbl <- group_sizes |>
          dplyr::count(group_col, name = "k_levels") |>
          dplyr::mutate(n_pairs = dplyr::if_else(k_levels >= 2, (k_levels * (k_levels - 1)) / 2, 0)) |>
          dplyr::collect()
        sum_pairs <- sum(lev_tbl$n_pairs, na.rm = TRUE)
        m_by_rank <- stats::setNames(pool_tbl$POOL * sum_pairs, pool_tbl$rank)
        .ph_log_info("fdr accounting",
                     bullets = c(
                       paste0("pool per rank: ",
                              paste(paste0(pool_tbl$rank, "=", pool_tbl$POOL), collapse = "; ")),
                       paste0("universes: ",
                              paste(paste0(lev_tbl$group_col, " (k=", lev_tbl$k_levels,
                                           ", pairs=", round(lev_tbl$n_pairs, 2), ")"), collapse = "; ")),
                       paste0("pairs across universes (sum): ", sum_pairs),
                       paste0("total tests m per rank = pool * pairs: ",
                              paste(paste0(names(m_by_rank), "=", m_by_rank), collapse = "; "))
                     ))

        # ---- base grid & stats -----------------------------------------------
        base_grid <- group_sizes |>
          dplyr::select(group_col, group_value) |>
          dplyr::distinct() |>
          dplyr::mutate(.dummy = 1L) |>
          dplyr::inner_join(features_per_rank |> dplyr::mutate(.dummy = 1L), by = ".dummy") |>
          dplyr::select(group_col, group_value, rank, feature)

        stats_long <- base_grid |>
          dplyr::left_join(present_counts, by = c("group_col", "group_value", "rank", "feature")) |>
          dplyr::left_join(group_sizes, by = c("group_col", "group_value")) |>
          dplyr::mutate(
            n_present = dplyr::coalesce(n_present, 0L),
            prop = dplyr::if_else(N > 0, n_present / N, NA_real_),
            percent = 100 * prop
          )
        if (!is.na(view_const)) {
          stats_long <- stats_long |> dplyr::mutate(view = !!view_const, .before = 1)
        }

        # ---- pairwise comparisons --------------------------------------------
        .ph_log_info("building pairwise comparisons")
        by_cols <- intersect(c("view", "rank", "feature", "group_col"), colnames(stats_long))
        has_view <- "view" %in% colnames(stats_long)
        pairs_joined <- stats_long |>
          dplyr::inner_join(stats_long, by = by_cols, suffix = c("_1", "_2")) |>
          dplyr::filter(group_value_1 != group_value_2) |>
          dplyr::mutate(
            group1   = dplyr::if_else(group_value_1 <= group_value_2, group_value_1, group_value_2),
            group2   = dplyr::if_else(group_value_1 <= group_value_2, group_value_2, group_value_1),
            n1_val   = dplyr::if_else(group_value_1 <= group_value_2, n_present_1, n_present_2),
            n1_tot   = dplyr::if_else(group_value_1 <= group_value_2, N_1, N_2),
            p1_val   = dplyr::if_else(group_value_1 <= group_value_2, prop_1, prop_2),
            pct1val  = dplyr::if_else(group_value_1 <= group_value_2, percent_1, percent_2),
            n2_val   = dplyr::if_else(group_value_1 <= group_value_2, n_present_2, n_present_1),
            n2_tot   = dplyr::if_else(group_value_1 <= group_value_2, N_2, N_1),
            p2_val   = dplyr::if_else(group_value_1 <= group_value_2, prop_2, prop_1),
            pct2val  = dplyr::if_else(group_value_1 <= group_value_2, percent_2, percent_1)
          ) |>
          dplyr::filter(group_value_1 == group1)

        pairs_lazy <- (if (has_view) {
          pairs_joined |>
            dplyr::transmute(
              view, rank, feature, group_col,
              group1 = group1, n1 = n1_val, n1_tot = n1_tot, prop1 = p1_val, percent1 = pct1val,
              group2 = group2, n2 = n2_val, n2_tot = n2_tot, prop2 = p2_val, percent2 = pct2val
            )
        } else {
          pairs_joined |>
            dplyr::transmute(
              rank, feature, group_col,
              group1 = group1, n1 = n1_val, n1_tot = n1_tot, prop1 = p1_val, percent1 = pct1val,
              group2 = group2, n2 = n2_val, n2_tot = n2_tot, prop2 = p2_val, percent2 = pct2val
            )
        }) |>
          dbplyr::window_order(group_col, rank, feature, group1, group2) |>
          dplyr::mutate(row_id = dplyr::row_number()) |>
          dplyr::filter(n1_tot > 0, n2_tot > 0)

        if (isTRUE(compute_ratios_db)) {
          pairs_lazy <- pairs_lazy |>
            dplyr::mutate(
              prop1_eps = dplyr::if_else(n1 == 0L, (n1 + 0.5) / n1_tot, prop1),
              prop2_eps = dplyr::if_else(n2 == 0L, (n2 + 0.5) / n2_tot, prop2),
              ratio = prop1_eps / prop2_eps,
              d1 = dplyr::if_else(n1 == 0L, (n1 + 1.0) / n1_tot, prop1),
              d2 = dplyr::if_else(n2 == 0L, (n2 + 1.0) / n2_tot, prop2),
              delta_ratio = dplyr::if_else(d1 >= d2, d1 / d2 - 1, -(d2 / d1 - 1))
            )
        }

        # ---- materialize + compute p-values ----------------------------------
        if (is.null(register_name)) register_name <- paste0("ph_prev_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        tbl_q <- .q(con, register_name)
        sql_body <- dbplyr::sql_render(pairs_lazy, con)
        DBI::dbExecute(con, paste0("drop table if exists ", tbl_q))
        DBI::dbExecute(con, paste0("create table ", tbl_q, " as ", sql_body))
        .ph_log_ok("materialized duckdb table",
          bullets = c(
            paste0("name: ", register_name),
            "computing p-values (fisher-only); then fdr per rank (bh / wbh)"
          )
        )

        p_core_fisher <- function(n1, N1, n2, N2) {
          n1 <- as.double(n1)
          N1 <- as.double(N1)
          n2 <- as.double(n2)
          N2 <- as.double(N2)
          if (any(!is.finite(c(n1, N1, n2, N2)))) {
            return(NA_real_)
          }
          if (N1 <= 0 || N2 <= 0) {
            return(NA_real_)
          }
          a <- max(0, min(n1, N1))
          b <- max(0, min(n2, N2))
          c1 <- N1 - a
          c2 <- N2 - b
          out <- try(stats::fisher.test(matrix(c(a, c1, b, c2), 2, byrow = TRUE))$p.value, silent = TRUE)
          if (inherits(out, "try-error")) NA_real_ else as.numeric(out)
        }

        for (col in c(
          "p_raw", "p_adj_rank", "p_adj_rank_wbh",
          "passed_rank_bh", "passed_rank_wbh",
          "category_rank_bh", "category_rank_wbh"
        )) {
          DBI::dbExecute(
            con,
            paste0(
              "alter table ", tbl_q, " add column if not exists ", col, " ",
              if (grepl("^p_", col)) "double" else if (grepl("^passed", col)) "boolean" else "varchar"
            )
          )
        }

        tmp_q <- .q(con, paste0(register_name, "_p_tmp_", format(Sys.time(), "%H%M%S")))
        DBI::dbExecute(con, paste0("create temp table ", tmp_q, " (row_id bigint, p_raw double)"))

        rs <- DBI::dbSendQuery(con, paste0(
          "select row_id, n1, n1_tot, n2, n2_tot from ", tbl_q, " order by row_id"
        ))
        on.exit(try(DBI::dbClearResult(rs), silent = TRUE), add = TRUE)

        want_parallel <- if (is.null(parallel)) {
          requireNamespace("future", quietly = TRUE) &&
            requireNamespace("future.apply", quietly = TRUE) &&
            tryCatch(future::nbrOfWorkers() > 1L, error = function(...) FALSE)
        } else {
          isTRUE(parallel)
        }
        if (want_parallel && (!requireNamespace("future", quietly = TRUE) ||
          !requireNamespace("future.apply", quietly = TRUE))) {
          .ph_abort("parallel requested but {future}/{future.apply} not available.")
        }

        repeat {
          chunk <- DBI::dbFetch(rs, n = chunk_n)
          if (!nrow(chunk)) break
          colnames(chunk) <- tolower(colnames(chunk))
          req <- c("row_id", "n1", "n1_tot", "n2", "n2_tot")
          miss <- setdiff(req, colnames(chunk))
          if (length(miss)) .ph_abort("fetched chunk is missing required columns.", bullets = paste("-", miss))
          chunk$row_id <- as.numeric(chunk$row_id)
          n <- nrow(chunk)

          if (want_parallel && tryCatch(future::nbrOfWorkers() > 1L, error = function(...) FALSE)) {
            splits <- split(seq_len(n), cut(seq_len(n), breaks = future::nbrOfWorkers(), labels = FALSE))
            parts <- future.apply::future_lapply(
              splits,
              FUN = function(ii, dat) {
                vapply(ii, function(i) p_core_fisher(dat$n1[i], dat$n1_tot[i], dat$n2[i], dat$n2_tot[i]), numeric(1))
              },
              dat = chunk, future.seed = TRUE, future.scheduling = 1, future.globals = FALSE
            )
            pvals <- unlist(parts, use.names = FALSE)
          } else {
            pvals <- vapply(
              seq_len(n),
              function(i) p_core_fisher(chunk$n1[i], chunk$n1_tot[i], chunk$n2[i], chunk$n2_tot[i]),
              numeric(1)
            )
          }

          DBI::dbAppendTable(
            con,
            name = gsub('^\"|\"$', "", tmp_q),
            value = data.frame(row_id = chunk$row_id, p_raw = as.numeric(pvals))
          )
        }

        DBI::dbExecute(con, paste0(
          "update ", tbl_q, " t set p_raw = s.p_raw from ", tmp_q, " s where t.row_id = s.row_id"
        ))
        DBI::dbExecute(con, paste0("drop table ", tmp_q))

        # ---- weights per (rank,feature) --------------------------------------
        weight_tbl_name <- paste0(register_name, "_weights")
        wq <- .q(con, weight_tbl_name)
        DBI::dbExecute(con, paste0("drop table if exists ", wq))

        if (weight_mode == "peptide_count" && length(ranks_needing_lib)) {
          if (is.null(lib_handle)) .ph_abort("peptide library required for peptide_count weights (not found).")
          lib_cols_needed <- c("peptide_id", ranks_needing_lib)
          lib_src <- if (inherits(lib_handle, "tbl_sql")) lib_handle else tibble::as_tibble(lib_handle)
          tmp_lib <- .q(con, paste0(register_name, "_libtmp"))
          DBI::dbExecute(con, paste0("drop table if exists ", tmp_lib))
          if (inherits(lib_src, "tbl_sql")) {
            lib_df <- lib_src %>% dplyr::select(tidyselect::all_of(lib_cols_needed)) %>% dplyr::distinct() %>% dplyr::collect()
            DBI::dbWriteTable(con, name = gsub('^\"|\"$', '', tmp_lib),
                              value = tibble::as_tibble(lib_df), temporary = TRUE)
          } else {
            DBI::dbWriteTable(con, name = gsub('^\"|\"$', '', tmp_lib),
                              value = tibble::as_tibble(lib_src[, lib_cols_needed, drop = FALSE] |> dplyr::distinct()),
                              temporary = TRUE)
          }
          bad <- ranks_needing_lib[!grepl("^[A-Za-z][A-Za-z0-9_]*$", ranks_needing_lib)]
          if (length(bad)) {
            .ph_abort("taxa column names must be alphanumeric/underscore (start with a letter).",
              bullets = paste("-", bad)
            )
          }
          selects <- vapply(
            ranks_needing_lib,
            function(rc) {
              rc_q <- .q(con, rc)
              paste0(
                "select '", rc, "' as rank, ", rc_q, " as feature, peptide_id from ", tmp_lib,
                " where ", rc_q, " is not null"
              )
            }, character(1)
          )
          sql_union <- paste(selects, collapse = " union all ")
          DBI::dbExecute(con, paste0(
            "create table ", wq, " as ",
            "select rank, feature, count(distinct peptide_id) as n_peptides ",
            "from (", sql_union, ") ",
            "group by rank, feature"
          ))

          if ("peptide_id" %in% available_ranks) {
            DBI::dbExecute(con, paste0(
              "insert into ", wq, " (rank, feature, n_peptides)
               select 'peptide_id' as rank, feature, 1 as n_peptides
               from (select distinct feature from ", tbl_q, " where rank = 'peptide_id')"
            ))
          }

          DBI::dbExecute(con, paste0("drop table if exists ", tmp_lib))
        } else {
          DBI::dbExecute(con, paste0(
            " create table ", wq, " as
              select distinct rank, feature, 1::integer as n_peptides from ", tbl_q
          ))
        }

        # ---- collect or return lazy ------------------------------------------
        if (isTRUE(collect)) {
          wq_name <- gsub('^\"|\"$', "", wq)
          out_df <- dplyr::tbl(con, register_name) |>
            dplyr::left_join(dplyr::tbl(con, wq_name), by = c("rank", "feature")) |>
            dplyr::arrange(rank, feature, group_col, group1, group2) |>
            dplyr::collect()

          # BH / wBH per (view?, rank)
          do_bh <- function(df) {
            df |>
              dplyr::mutate(
                p_adj_rank = {
                  ok <- !is.na(p_raw)
                  out <- rep(NA_real_, length(p_raw))
                  if (any(ok)) out[ok] <- p.adjust(p_raw[ok], method = "BH")
                  out
                },
                passed_rank_bh = !is.na(p_adj_rank) & p_adj_rank < 0.05,
                category_rank_bh = dplyr::case_when(
                  !is.na(p_adj_rank) & p_adj_rank < 0.05 ~ "significant (BH, per rank)",
                  !is.na(p_raw) & p_raw < 0.05 ~ "nominal only",
                  TRUE ~ "not significant"
                )
              )
          }

          do_wbh <- function(df) {
            df2 <- df |> dplyr::mutate(n_peptides = dplyr::coalesce(n_peptides, 1.0))
            split_vars <- intersect(c("view", "rank"), names(df2))
            pieces <- split(df2, df2[split_vars], drop = TRUE)
            pieces_adj <- lapply(pieces, function(dd) {
              idx <- which(!is.na(dd$p_raw))
              if (!length(idx)) {
                dd$p_adj_rank_wbh <- NA_real_
                dd$passed_rank_wbh <- NA
                dd$category_rank_wbh <- "not significant"
                return(dd)
              }
              m <- length(idx)
              w_base <- dd$n_peptides[idx]
              w_scaled <- w_base * (m / sum(w_base))
              p_over_w <- dd$p_raw[idx] / w_scaled
              ord <- order(p_over_w, na.last = NA)
              ranks <- seq_along(ord)
              raw <- m * p_over_w[ord] / ranks
              adj <- cummin(rev(raw))
              adj <- rev(adj)
              q <- rep(NA_real_, nrow(dd))
              q[idx[ord]] <- pmin(1.0, adj)
              dd$p_adj_rank_wbh <- q
              dd$passed_rank_wbh <- !is.na(dd$p_adj_rank_wbh) & dd$p_adj_rank_wbh < 0.05
              dd$category_rank_wbh <- dplyr::case_when(
                !is.na(dd$p_adj_rank_wbh) & dd$p_adj_rank_wbh < 0.05 ~ "significant (wBH, per rank)",
                !is.na(dd$p_raw)          & dd$p_raw          < 0.05 ~ "nominal only",
                TRUE ~ "not significant"
              )
              dd
            })
            dplyr::bind_rows(pieces_adj)
          }

          split_vars <- intersect(c("view", "rank"), names(out_df))
          out_df_bh <- do.call(rbind, lapply(split(out_df, out_df[split_vars], drop = TRUE), do_bh))
          out_df_wbh <- do_wbh(out_df_bh)

          out_df <- out_df_wbh |>
            dplyr::select(
              tidyselect::any_of("view"), rank, feature, n_peptides,
              group_col, group1, group2,
              n1,
              N1 = n1_tot, prop1, percent1,
              n2, N2 = n2_tot, prop2, percent2,
              tidyselect::any_of(c("ratio", "delta_ratio")),
              p_raw, p_adj_rank, passed_rank_bh, category_rank_bh,
              p_adj_rank_wbh, passed_rank_wbh, category_rank_wbh
            )

          meta <- list(
            fdr_scope = "per-rank",
            pool_by_rank = tibble::as_tibble(pool_tbl),
            pairs_by_universe = tibble::as_tibble(lev_tbl),
            m_by_rank = stats::setNames(pool_tbl$POOL * sum(lev_tbl$n_pairs, na.rm = TRUE), pool_tbl$rank),
            paired = FALSE,
            weight_mode = weight_mode,
            pop_k_min = pop_k_min,
            register_name = register_name,
            view = view_const
          )

          # store whatever was provided as library (optional)
          meta$peptide_library      <- lib_handle
          meta$peptide_library_cols <- tryCatch({
            if (inherits(lib_handle, "tbl_sql")) colnames(dplyr::collect(dplyr::slice_head(lib_handle, n = 0)))
            else colnames(lib_handle)
          }, error = function(...) NULL)

          out <- .as_prev_result(out_df, meta)
          return(out)
        } else {
          .ph_log_ok("materialization complete (TABLE; p computed). returning lazy table without BH/wBH (collect to compute).",
            bullets = register_name
          )

          meta <- list(
            fdr_scope = "per-rank",
            pool_by_rank = tibble::as_tibble(pool_tbl),
            pairs_by_universe = tibble::as_tibble(lev_tbl),
            m_by_rank = stats::setNames(pool_tbl$POOL * sum(lev_tbl$n_pairs, na.rm = TRUE), pool_tbl$rank),
            paired = FALSE,
            weight_mode = weight_mode,
            pop_k_min = pop_k_min,
            register_name = register_name,
            view = view_const
          )
          meta$peptide_library      <- lib_handle
          meta$peptide_library_cols <- tryCatch({
            if (inherits(lib_handle, "tbl_sql")) colnames(dplyr::collect(dplyr::slice_head(lib_handle, n = 0)))
            else colnames(lib_handle)
          }, error = function(...) NULL)
          lazy_tbl <- dplyr::tbl(con, register_name)
          out <- .as_prev_result(lazy_tbl, meta)
          return(out)
        }
      }

      # ============================= paired path ==============================
      .ph_log_info("paired design detected: running mcnemar exact (binomial)")

      present_counts <- k_tbl %>%
        dplyr::filter(present) %>%
        dplyr::distinct(group_col, group_value, rank, feature, sample_id) %>%
        dplyr::count(group_col, group_value, rank, feature, name = "n_present")
      features_per_rank <- present_counts |> dplyr::distinct(rank, feature)

      pool_tbl <- features_per_rank |>
        dplyr::count(rank, name = "POOL") |>
        dplyr::arrange(rank) |>
        dplyr::collect()
      lev_tbl <- group_sizes |>
        dplyr::count(group_col, name = "k_levels") |>
        dplyr::mutate(n_pairs = dplyr::if_else(k_levels >= 2, (k_levels * (k_levels - 1)) / 2, 0)) |>
        dplyr::collect()
      sum_pairs <- sum(lev_tbl$n_pairs, na.rm = TRUE)
      m_by_rank <- stats::setNames(pool_tbl$POOL * sum_pairs, pool_tbl$rank)
      .ph_log_info("fdr accounting (paired)",
        bullets = c(
          paste0("pool per rank: ", paste(paste0(pool_tbl$rank, "=", pool_tbl$POOL), collapse = "; ")),
          paste0(
            "universes: ",
            paste(paste0(lev_tbl$group_col, " (k=", lev_tbl$k_levels, ", pairs=", lev_tbl$n_pairs, ")"),
              collapse = "; "
            )
          ),
          paste0("pairs across universes (sum): ", sum_pairs),
          paste0(
            "total tests m per rank = pool * pairs: ",
            paste(paste0(names(m_by_rank), "=", m_by_rank), collapse = "; ")
          )
        )
      )

      base_grid <- group_sizes |>
        dplyr::select(group_col, group_value) |>
        dplyr::distinct() |>
        dplyr::mutate(.dummy = 1L) |>
        dplyr::inner_join(features_per_rank |> dplyr::mutate(.dummy = 1L), by = ".dummy") |>
        dplyr::select(group_col, group_value, rank, feature)

      stats_long <- base_grid |>
        dplyr::left_join(present_counts, by = c("group_col", "group_value", "rank", "feature")) |>
        dplyr::left_join(group_sizes, by = c("group_col", "group_value")) |>
        dplyr::mutate(
          n_present = dplyr::coalesce(n_present, 0L),
          prop = dplyr::if_else(N > 0, n_present / N, NA_real_),
          percent = 100 * prop
        )
      if (!is.na(view_const)) {
        stats_long <- stats_long |> dplyr::mutate(view = !!view_const, .before = 1)
      }

      by_cols <- intersect(c("view", "rank", "feature", "group_col"), colnames(stats_long))
      has_view <- "view" %in% colnames(stats_long)

      pairs_joined <- stats_long %>%
        dplyr::inner_join(stats_long, by = by_cols, suffix = c("_1","_2")) %>%
        dplyr::filter(group_value_1 != group_value_2) %>%
        dplyr::mutate(
          group1 = dplyr::if_else(group_value_1 < group_value_2, group_value_1, group_value_2),
          group2 = dplyr::if_else(group_value_1 < group_value_2, group_value_2, group_value_1),
          n1 = dplyr::if_else(group_value_1 < group_value_2, n_present_1, n_present_2),
          n2 = dplyr::if_else(group_value_1 < group_value_2, n_present_2, n_present_1),
          prop1 = dplyr::if_else(group_value_1 < group_value_2, prop_1, prop_2),
          prop2 = dplyr::if_else(group_value_1 < group_value_2, prop_2, prop_1),
          percent1 = dplyr::if_else(group_value_1 < group_value_2, percent_1, percent_2),
          percent2 = dplyr::if_else(group_value_1 < group_value_2, percent_2, percent_1)
        ) %>%
        dplyr::filter(group_value_1 == group1) %>%
        { if (has_view) dplyr::transmute(., view, rank, feature, group_col, group1, group2, n1, n2, prop1, percent1, prop2, percent2)
          else dplyr::transmute(., rank, feature, group_col, group1, group2, n1, n2, prop1, percent1, prop2, percent2) } %>%
        dplyr::distinct()

      sdat <- k_tbl %>%
        dplyr::select(group_col, group_value, rank, feature, dplyr::all_of(found_paired_col), present) %>%
        dplyr::collect() %>%
        dplyr::rename(subject_id = !!rlang::sym(found_paired_col))

      disc <- sdat %>%
        dplyr::inner_join(., ., by = c("group_col", "rank", "feature", "subject_id"), suffix = c("_1", "_2")) %>%
        dplyr::filter(group_value_1 < group_value_2) %>%
        dplyr::group_by(group_col, rank, feature, group1 = group_value_1, group2 = group_value_2) %>%
        dplyr::summarise(
          n01 = sum((!present_1) & (present_2), na.rm = TRUE),
          n10 = sum((present_1) & (!present_2), na.rm = TRUE),
          .groups = "drop"
        ) %>%
        {
          n_vec <- .$n01 + .$n10
          p_vals <- vapply(seq_along(n_vec), function(i) {
            n <- n_vec[i]
            x <- .$n01[i]
            if (is.na(n) || n <= 0) {
              return(NA_real_)
            }
            stats::binom.test(x, n, alternative = "two.sided")$p.value
          }, numeric(1))
          dplyr::mutate(., p_raw = p_vals)
        }

      paired_summary <- sdat %>%
        dplyr::inner_join(., ., by = c("group_col", "rank", "feature", "subject_id"), suffix = c("_1", "_2")) %>%
        dplyr::filter(group_value_1 < group_value_2) %>%
        dplyr::group_by(group_col, rank, feature, group1 = group_value_1, group2 = group_value_2) %>%
        dplyr::summarise(
          N_paired  = dplyr::n_distinct(subject_id),
          n1_paired = sum(present_1, na.rm = TRUE),
          n2_paired = sum(present_2, na.rm = TRUE),
          .groups   = "drop"
        ) %>%
        dplyr::mutate(
          prop1_paired = dplyr::if_else(N_paired > 0, n1_paired / N_paired, NA_real_),
          prop2_paired = dplyr::if_else(N_paired > 0, n2_paired / N_paired, NA_real_)
        )

      pairs_joined_df <- pairs_joined %>% dplyr::collect()
      res <- pairs_joined_df %>%
        dplyr::left_join(disc, by = c("group_col", "rank", "feature", "group1", "group2")) %>%
        dplyr::left_join(paired_summary, by = c("group_col", "rank", "feature", "group1", "group2")) %>%
        dplyr::mutate(
          N1 = N_paired, N2 = N_paired,
          n1 = dplyr::coalesce(n1_paired, n1),
          n2 = dplyr::coalesce(n2_paired, n2),
          prop1 = dplyr::coalesce(prop1_paired, prop1),
          prop2 = dplyr::coalesce(prop2_paired, prop2),
          percent1 = 100 * prop1,
          percent2 = 100 * prop2
        )

      w_tbl <- if (weight_mode == "peptide_count" && length(ranks_needing_lib)) {
        if (is.null(lib_handle)) .ph_abort("peptide library required for peptide_count weights (not found).")

        lib_min <- lib_handle %>%
          dplyr::select(tidyselect::all_of(c("peptide_id", ranks_needing_lib))) %>%
          dplyr::distinct()

        if (inherits(lib_min, "tbl_sql") || inherits(lib_min, "tbl_dbi") || inherits(lib_min, "tbl_lazy")) {
          .ph_log_info("collecting peptide_library into memory for paired wBH weights",
            bullets = paste0(
              "rows (pre-collect): ",
              tryCatch(dplyr::count(lib_min) %>% dplyr::pull(n),
                error = function(e) "unknown"
              )
            )
          )
          lib_min <- lib_min %>% dplyr::collect()
          .ph_log_info("collection complete", bullets = paste0("rows (post-collect): ", nrow(lib_min)))
        }

        base_tbl <- purrr::map_dfr(ranks_needing_lib, function(rc) {
          lib_min %>%
            dplyr::filter(!is.na(.data[[rc]])) %>%
            dplyr::distinct(.data[[rc]], peptide_id) %>%
            dplyr::count(rank = rc, feature = .data[[rc]], name = "n_peptides")
        })

        if ("peptide_id" %in% available_ranks) {
          pid_vals <- unique(res$feature[res$rank == "peptide_id"])
          base_tbl <- dplyr::bind_rows(
            base_tbl,
            tibble::tibble(rank = "peptide_id", feature = pid_vals, n_peptides = 1L)
          )
        }
        base_tbl
      } else {
        res %>%
          dplyr::distinct(rank, feature) %>%
          dplyr::mutate(n_peptides = 1L)
      }

      do_bh <- function(df) {
        df %>%
          dplyr::mutate(
            p_adj_rank = {
              ok <- !is.na(p_raw)
              out <- rep(NA_real_, length(p_raw))
              if (any(ok)) out[ok] <- p.adjust(p_raw[ok], method = "BH")
              out
            },
            passed_rank_bh = !is.na(p_adj_rank) & p_adj_rank < 0.05,
            category_rank_bh = dplyr::case_when(
              !is.na(p_adj_rank) & p_adj_rank < 0.05 ~ "significant (BH, per rank)",
              !is.na(p_raw) & p_raw < 0.05 ~ "nominal only",
              TRUE ~ "not significant"
            )
          )
      }

      do_wbh <- function(df, w_tbl) {
        df2 <- df %>%
          dplyr::left_join(w_tbl, by = c("rank", "feature")) %>%
          dplyr::mutate(n_peptides = dplyr::coalesce(n_peptides, 1.0))
        split_vars <- intersect(c("view", "rank"), names(df2))
        pieces <- split(df2, df2[split_vars], drop = TRUE)
        pieces_adj <- lapply(pieces, function(dd) {
          idx <- which(!is.na(dd$p_raw))
          if (!length(idx)) {
            dd$p_adj_rank_wbh <- NA_real_
            dd$passed_rank_wbh <- NA
            dd$category_rank_wbh <- "not significant"
            return(dd)
          }
          m <- length(idx)
          w_base <- dd$n_peptides[idx]
          w_scaled <- w_base * (m / sum(w_base))
          p_over_w <- dd$p_raw[idx] / w_scaled
          ord <- order(p_over_w, na.last = NA)
          ranks <- seq_along(ord)
          raw <- m * p_over_w[ord] / ranks
          adj <- cummin(rev(raw))
          adj <- rev(adj)
          q <- rep(NA_real_, nrow(dd))
          q[idx[ord]] <- pmin(1.0, adj)
          dd$p_adj_rank_wbh <- q
          dd$passed_rank_wbh <- !is.na(dd$p_adj_rank_wbh) & dd$p_adj_rank_wbh < 0.05
          dd$category_rank_wbh <- dplyr::case_when(
            !is.na(dd$p_adj_rank_wbh) & dd$p_adj_rank_wbh < 0.05 ~ "significant (wBH, per rank)",
            !is.na(dd$p_raw & dd$p_raw < 0.05) ~ "nominal only",
            TRUE ~ "not significant"
          )
          dd
        })
        dplyr::bind_rows(pieces_adj)
      }

      split_vars <- intersect(c("view", "rank"), names(res))
      res_bh <- do.call(rbind, lapply(split(res, res[split_vars], drop = TRUE), do_bh))
      res_wbh <- do_wbh(res_bh, w_tbl)

      out_df <- res_wbh %>%
        dplyr::arrange(rank, feature, group_col, group1, group2) %>%
        dplyr::select(
          tidyselect::any_of("view"), rank, feature,
          n_peptides = tidyselect::any_of("n_peptides"),
          group_col, group1, group2,
          n1, N1, prop1, percent1,
          n2, N2, prop2, percent2,
          tidyselect::any_of(c("ratio", "delta_ratio")),
          p_raw, p_adj_rank, passed_rank_bh, category_rank_bh,
          p_adj_rank_wbh, passed_rank_wbh, category_rank_wbh
        )

      .ph_log_ok("done (paired mcnemar per-rank fdr)",
        bullets = c(
          paste0("rows: ", nrow(out_df)),
          paste0("ranks: ", paste(unique(out_df$rank), collapse = ", ")),
          paste0("pop_k_min: ", pop_k_min)
        )
      )

      meta <- list(
        fdr_scope = "per-rank",
        pool_by_rank = tibble::as_tibble(pool_tbl),
        pairs_by_universe = tibble::as_tibble(lev_tbl),
        m_by_rank = m_by_rank,
        paired = TRUE,
        paired_col = found_paired_col,
        weight_mode = weight_mode,
        pop_k_min = pop_k_min,
        register_name = register_name,
        view = view_const
      )
      meta$peptide_library      <- lib_handle
      meta$peptide_library_cols <- tryCatch({
        if (inherits(lib_handle, "tbl_sql")) colnames(dplyr::collect(dplyr::slice_head(lib_handle, n = 0)))
        else colnames(lib_handle)
      }, error = function(...) NULL)

      out <- .as_prev_result(out_df, meta)
      return(out)
    }
  )
}


# phiper-style writers for results (generic + ph_prev_result)
# - Generic dispatcher: write_result()
# - Method: write_result.ph_prev_result()
#
# Uses phiper logging/check utilities (.ph_log_info/.ph_log_ok/.ph_abort/.ph_warn,
# .chk_extension, .chk_null_default, .ph_opt) which are expected to be present in
# the environment (phiper package namespace).
#
# Behavior summary:
# - Supports formats inferred from path or explicitly via `format`: "xlsx", "csv", "parquet"
# - filter: "none" | "only_significant" (p_raw < 0.05) | "only_significant_fdr" (passed_rank_wbh or fallback passed_rank_bh)
# - sheet_by_rank: for xlsx -> one sheet per rank (when TRUE); for csv/parquet -> one file per rank
# - sorts each output by p_adj_rank_wbh (NA last), then p_adj_rank, then p_raw
# - uses phiper logging for progress and errors; does NOT install packages

#' Generic writer for phiper result objects
#'
#' @export
write_result <- function(x, path, ...) {
  UseMethod("write_result")
}

#' Write method for ph_prev_result (ph_prevalence_compare output)
#'
#' @param x
#'   A `ph_prev_result` object (tibble/data.frame with attribute `prev_meta`).
#' @param path File path. Extension determines format (.xlsx/.csv/.parquet) unless format is provided.
#' @param filter One of c("none","only_significant","only_significant_fdr").
#' @param format Optional override for format: "xlsx","csv","parquet".
#' @param sheet_by_rank Logical; if TRUE create one sheet/file per rank when multiple ranks present.
#' @param overwrite Logical; allow overwriting existing files (default FALSE).
#' @param ... Reserved for future use.
#' @export
write_result.ph_prev_result <- function(x,
                                        path,
                                        filter = c("none", "only_significant", "only_significant_fdr"),
                                        format = NULL,
                                        sheet_by_rank = TRUE,
                                        overwrite = FALSE,
                                        ...) {
  filter <- match.arg(filter)

  .ph_with_timing(headline = "write_result (ph_prev_result)", expr = {
    # ---- basic checks ----------------------------------------------------
    if (!inherits(x, "ph_prev_result")) {
      .ph_abort("x must be a ph_prev_result object")
    }

    if (missing(path) || !nzchar(path)) {
      .ph_abort("please provide a valid 'path' argument")
    }

    # infer/validate format
    ext <- tolower(tools::file_ext(path))
    if (is.null(format)) {
      format <- switch(ext,
        "xlsx" = "xlsx",
        "csv" = "csv",
        "parquet" = "parquet",
        ""
      )
      if (format == "") {
        .ph_abort("cannot infer format from path; supply extension .xlsx/.csv/.parquet or set format argument")
      }
    } else {
      format <- match.arg(tolower(format), c("xlsx", "csv", "parquet"))
    }

    .ph_log_info(
      headline = "writer configuration",
      bullets = c(
        sprintf("path: %s", path),
        sprintf("format: %s", format),
        sprintf("filter: %s", filter),
        sprintf("sheet_by_rank: %s", ifelse(isTRUE(sheet_by_rank), "TRUE", "FALSE")),
        sprintf("overwrite: %s", ifelse(isTRUE(overwrite), "TRUE", "FALSE"))
      )
    )

    # ---- required packages ------------------------------------------------
    if (format == "xlsx" && !requireNamespace("openxlsx", quietly = TRUE)) {
      .ph_abort("package 'openxlsx' required for .xlsx output (install it with install.packages('openxlsx'))")
    }
    if (format == "parquet" && !requireNamespace("arrow", quietly = TRUE)) {
      .ph_abort("package 'arrow' required for .parquet output (install it with install.packages('arrow'))")
    }
    if (format == "csv" && !requireNamespace("readr", quietly = TRUE)) {
      .ph_warn("package 'readr' not available; falling back to utils::write.csv")
    }

    # ---- coerce to data.frame & ensure columns exist ----------------------
    df <- tryCatch(as.data.frame(x), error = function(e) {
      .ph_abort("failed to coerce result to data.frame", bullets = e$message)
    })

    # ensure canonical columns exist (create NA defaults to avoid errors)
    for (col in c("p_raw", "p_adj_rank_wbh", "p_adj_rank", "passed_rank_wbh", "passed_rank_bh", "rank")) {
      if (!col %in% colnames(df)) df[[col]] <- NA
    }

    # ---- filtering --------------------------------------------------------
    .ph_log_info("filtering rows", bullets = sprintf("mode: %s", filter))

    if (filter == "only_significant") {
      keep_idx <- which(!is.na(df$p_raw) & df$p_raw < 0.05)
      .ph_log_info("only_significant: selecting rows with p_raw < 0.05",
        bullets = sprintf("rows kept: %d", length(keep_idx))
      )
      df <- df[keep_idx, , drop = FALSE]
    } else if (filter == "only_significant_fdr") {
      # prefer wbh flag, fallback to bh
      keep_flag <- ifelse(!is.na(df$passed_rank_wbh), df$passed_rank_wbh,
        ifelse(!is.na(df$passed_rank_bh), df$passed_rank_bh, FALSE)
      )
      keep_idx <- which(as.logical(keep_flag))
      .ph_log_info("only_significant_fdr: selecting rows passing FDR",
        bullets = sprintf("rows kept: %d", length(keep_idx))
      )
      df <- df[keep_idx, , drop = FALSE]
    } else {
      .ph_log_info("no filtering applied (filter = 'none')")
    }

    if (nrow(df) == 0L) {
      .ph_warn("result contains zero rows after filtering; will still write empty file(s)")
    }

    # ---- sorting helper --------------------------------------------------
    sort_key <- function(d) {
      # NAs pushed to end -> convert to Inf
      a <- ifelse(is.na(d$p_adj_rank_wbh), Inf, d$p_adj_rank_wbh)
      b <- ifelse(is.na(d$p_adj_rank), Inf, d$p_adj_rank)
      c <- ifelse(is.na(d$p_raw), Inf, d$p_raw)
      order(a, b, c, na.last = TRUE)
    }

    # ---- rank splitting --------------------------------------------------
    ranks <- unique(as.character(df$rank))
    if (length(ranks) == 0) ranks <- character(0)

    sanitize_sheet_name <- function(n) {
      # replace characters invalid in Excel sheet names: : \ / ? * [ ]
      n2 <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", n)
      n2 <- substr(n2, 1, 31)
      if (nchar(n2) == 0) n2 <- "sheet"
      n2
    }


    base_path_sans_ext <- tools::file_path_sans_ext(path)

    # ---- write depending on format ---------------------------------------
    if (format == "xlsx") {
      wb <- openxlsx::createWorkbook()

      if (length(ranks) <= 1 || !isTRUE(sheet_by_rank)) {
        out_df <- df
        if (nrow(out_df)) out_df <- out_df[sort_key(out_df), , drop = FALSE]
        openxlsx::addWorksheet(wb, "results")
        openxlsx::writeData(wb, sheet = "results", out_df)
      } else {
        for (r in ranks) {
          df_r <- df[df$rank == r, , drop = FALSE]
          if (nrow(df_r)) df_r <- df_r[sort_key(df_r), , drop = FALSE]
          sheet_name <- sanitize_sheet_name(r)
          # avoid duplicate sheet names
          i <- 1
          orig <- sheet_name
          while (sheet_name %in% openxlsx::sheets(wb)) {
            i <- i + 1
            sheet_name <- substr(paste0(orig, "_", i), 1, 31)
          }
          openxlsx::addWorksheet(wb, sheet_name)
          openxlsx::writeData(wb, sheet = sheet_name, df_r)
        }
      }

      if (file.exists(path) && !isTRUE(overwrite)) {
        .ph_abort("path exists; set overwrite = TRUE to replace it", bullets = path)
      }
      openxlsx::saveWorkbook(wb, file = path, overwrite = TRUE)
      .ph_log_ok("wrote xlsx file", bullets = path)
      return(invisible(path))
    }

    if (format == "csv") {
      if (length(ranks) <= 1 || !isTRUE(sheet_by_rank)) {
        out_df <- df
        if (nrow(out_df)) out_df <- out_df[sort_key(out_df), , drop = FALSE]
        if (file.exists(path) && !isTRUE(overwrite)) {
          .ph_abort("path exists; set overwrite = TRUE to replace it", bullets = path)
        }
        if (requireNamespace("readr", quietly = TRUE)) {
          readr::write_csv(out_df, path)
        } else {
          utils::write.csv(out_df, path, row.names = FALSE)
        }
        .ph_log_ok("wrote csv file", bullets = path)
        return(invisible(path))
      } else {
        out_files <- character(length(ranks))
        for (i in seq_along(ranks)) {
          r <- ranks[i]
          df_r <- df[df$rank == r, , drop = FALSE]
          if (nrow(df_r)) df_r <- df_r[sort_key(df_r), , drop = FALSE]
          safe <- gsub("\\s+", "_", r)
          fname <- paste0(base_path_sans_ext, "_", safe, ".csv")
          if (file.exists(fname) && !isTRUE(overwrite)) {
            .ph_abort("path exists; set overwrite = TRUE to replace it", bullets = fname)
          }
          if (requireNamespace("readr", quietly = TRUE)) {
            readr::write_csv(df_r, fname)
          } else {
            utils::write.csv(df_r, fname, row.names = FALSE)
          }
          out_files[i] <- fname
          .ph_log_info("wrote csv for rank", bullets = sprintf("%s -> %s", r, fname))
        }
        .ph_log_ok("wrote csv files", bullets = paste(out_files, collapse = ", "))
        return(invisible(out_files))
      }
    }

    if (format == "parquet") {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        .ph_abort("arrow required for parquet writing")
      }
      if (length(ranks) <= 1 || !isTRUE(sheet_by_rank)) {
        out_df <- df
        if (nrow(out_df)) out_df <- out_df[sort_key(out_df), , drop = FALSE]
        if (file.exists(path) && !isTRUE(overwrite)) {
          .ph_abort("path exists; set overwrite = TRUE to replace it", bullets = path)
        }
        arrow::write_parquet(out_df, path)
        .ph_log_ok("wrote parquet file", bullets = path)
        return(invisible(path))
      } else {
        out_files <- character(length(ranks))
        for (i in seq_along(ranks)) {
          r <- ranks[i]
          df_r <- df[df$rank == r, , drop = FALSE]
          if (nrow(df_r)) df_r <- df_r[sort_key(df_r), , drop = FALSE]
          safe <- gsub("\\s+", "_", r)
          fname <- paste0(base_path_sans_ext, "_", safe, ".parquet")
          if (file.exists(fname) && !isTRUE(overwrite)) {
            .ph_abort("path exists; set overwrite = TRUE to replace it", bullets = fname)
          }
          arrow::write_parquet(df_r, fname)
          out_files[i] <- fname
          .ph_log_info("wrote parquet for rank", bullets = sprintf("%s -> %s", r, fname))
        }
        .ph_log_ok("wrote parquet files", bullets = paste(out_files, collapse = ", "))
        return(invisible(out_files))
      }
    }

    .ph_abort("unhandled format internally")
  }) # ph_with_timing
}
