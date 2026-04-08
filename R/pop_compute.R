# ==============================================================================
# Helper: constructor for class `ph_prev_result`
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
#' \deqn{ w_i^{*} \;=\; w_i \cdot \frac{m_r}{\sum_j w_j}. }
#'
#' We adjust using the standard weighted step-up rule on \eqn{p_i / w_i^{*}}. The
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
#' @param pop_k_min Integer \eqn{\ge} 1; k-of-n POP threshold per sample (default 1).
#' @param paired Logical; use paired design (McNemar exact) with `subject_id` (default `FALSE`).
#'        NOTE: can also be a character scalar naming the column that links related samples
#'        (e.g. "subject_id" or "dyade"). If so, only samples present in both groups
#'        for that identifier will be used for paired McNemar tests.
#' @param peptide_library Optional peptide metadata table used to map
#'   non-peptide rank columns. If `NULL`, the function will use
#'   `x$peptide_library`.
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
          colnames(dplyr::collect(head(lib_handle, 0)))
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
                  if (any(ok)) out[ok] <- stats::p.adjust(p_raw[ok], method = "BH")
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

          do_wbh_simple <- function(df) {
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
          out_df_wbh <- do_wbh_simple(out_df_bh)

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
            if (inherits(lib_handle, "tbl_sql")) colnames(dplyr::collect(head(lib_handle, 0)))
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
            if (inherits(lib_handle, "tbl_sql")) colnames(dplyr::collect(head(lib_handle, 0)))
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

        base_tbl <- dplyr::bind_rows(lapply(ranks_needing_lib, function(rc) {
          lib_min %>%
            dplyr::filter(!is.na(.data[[rc]])) %>%
            dplyr::distinct(.data[[rc]], peptide_id) %>%
            dplyr::count(rank = rc, feature = .data[[rc]], name = "n_peptides")
        }))

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
              if (any(ok)) out[ok] <- stats::p.adjust(p_raw[ok], method = "BH")
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

      do_wbh_weighted <- function(df, w_tbl) {
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
      res_wbh <- do_wbh_weighted(res_bh, w_tbl)

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
        if (inherits(lib_handle, "tbl_sql")) colnames(dplyr::collect(head(lib_handle, 0)))
        else colnames(lib_handle)
      }, error = function(...) NULL)

      out <- .as_prev_result(out_df, meta)
      return(out)
    }
  )
}

