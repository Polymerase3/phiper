#' Prevalence comparison by group (POP framework)
#'
#' @description
#' Computes feature-level prevalence across binary group columns and performs
#' pairwise statistical tests:
#' - **Unpaired:** Fisher's exact test (2x2) per (rank, feature, group pair).
#' - **Paired:** McNemar's exact binomial test per subject.
#'
#' Presence per sample is determined by a k-of-n rule: a feature is considered
#' present in a sample if at least `pop_k_min` of its contributing peptides are
#' positive. Each `group_col` must have exactly 2 non-missing levels in the data.
#'
#' @param x A `phip_data` object or a data.frame/tibble with at least columns
#'   `sample_id`, `peptide_id`, `exist_col`, and all `group_cols`. If
#'   `paired` is set, the paired linking column must also be present. If
#'   `rank_cols` include non-peptide taxa, `x` must provide `peptide_library`.
#' @param rank_cols Character vector of rank columns, e.g. `c("peptide_id", "species")`.
#' @param group_cols Character vector of binary grouping columns. Each column must
#'   have exactly 2 non-missing levels in the data.
#' @param exist_col Name of the binary presence column (default `"exist"`).
#' @param pop_k_min Integer >= 1; k-of-n POP threshold per sample (default 1).
#' @param paired `FALSE` (default) or a single string naming the column that
#'   links related samples (e.g. `"subject_id"`). When set, paired McNemar
#'   tests are used instead of Fisher.
#' @param peptide_library Optional peptide metadata table with `peptide_id` and
#'   any non-peptide rank columns. If `NULL`, taken from `x$peptide_library`.
#'
#' @return A data.frame with one row per (rank, feature, group_col) comparison:
#'   `rank`, `feature`, `group_col`, `group1`, `group2`,
#'   `n1`, `N1`, `prop1`, `percent1`,
#'   `n2`, `N2`, `prop2`, `percent2`,
#'   `ratio`, `delta_ratio` (unpaired only), `p_raw`.
#'   A `view` column is prepended when the input carries a view attribute.
#'
#' @examples
#' \dontrun{
#' res <- compute_pop(pd, rank_cols = "species", group_cols = "group")
#' res
#' }
#'
#' @export
compute_pop <- function(x,
                        rank_cols,
                        group_cols,
                        exist_col       = "exist",
                        pop_k_min       = 1L,
                        paired          = FALSE,
                        peptide_library = NULL) {

  .ph_with_timing(
    headline = "compute_pop",
    step     = NULL,
    bullets  = NULL,
    expr     = {
      .q      <- function(con, nm) as.character(DBI::dbQuoteIdentifier(con, nm))
      .sym    <- rlang::sym
      chunk_n <- getOption("phiper.prev.chunk", 1e6L)

      # --- library handle ---------------------------------------------------
      lib_handle <- peptide_library %||%
        tryCatch(x$peptide_library,      error = function(...) NULL) %||%
        tryCatch(x$meta$peptide_library, error = function(...) NULL)

      # --- argument checks --------------------------------------------------
      tryCatch(
        {
          chk::chk_character(rank_cols);  chk::chk_true(length(rank_cols)  >= 1)
          chk::chk_character(group_cols); chk::chk_true(length(group_cols) >= 1)
          chk::chk_string(exist_col)
          chk::chk_number(pop_k_min); chk::chk_true(pop_k_min >= 1)
        },
        error = function(e) .ph_abort("invalid arguments", bullets = e$message)
      )

      paired_col <- NULL
      if (!identical(paired, FALSE)) {
        if (!chk::vld_string(paired)) .ph_abort("paired must be FALSE or a single column name.")
        paired_col <- as.character(paired)
      }

      .ph_log_info(
        "compute_pop",
        bullets = c(
          paste0("ranks     : ", paste(rank_cols,  collapse = ", ")),
          paste0("group_cols: ", paste(group_cols, collapse = ", ")),
          paste0("exist_col : ", exist_col),
          paste0("pop_k_min : ", pop_k_min),
          paste0("paired    : ", ifelse(is.null(paired_col), "FALSE", paired_col))
        )
      )

      # --- normalize input to long df ---------------------------------------
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
            tibble::as_tibble(x) |>
              dplyr::select(tidyselect::any_of(
                c("sample_id", "peptide_id", exist_col, group_cols, paired_col)
              ))
          }
        },
        silent = TRUE
      )
      if (inherits(df_long, "try-error")) .ph_abort("could not prepare input data.")

      df_long <- df_long |>
        dplyr::mutate(!!.sym(exist_col) := dplyr::coalesce(as.integer(!!.sym(exist_col)), 0L))

      # --- paired column detection ------------------------------------------
      found_paired_col <- NULL
      if (!is.null(paired_col)) {
        cols_here <- tryCatch(colnames(df_long), error = function(...) character(0))
        if (paired_col %in% cols_here) {
          found_paired_col <- paired_col
        } else {
          cand <- paste0(paired_col, "_recoded")
          if (cand %in% cols_here) {
            found_paired_col <- cand
            .ph_warn(headline = sprintf("paired column '%s' not found; using '%s'.", paired_col, cand))
          } else {
            approx <- grep(paired_col, cols_here, value = TRUE)
            if (length(approx) == 1L) {
              found_paired_col <- approx
              .ph_warn(headline = sprintf("paired column '%s' not found; using '%s'.", paired_col, approx))
            } else if (length(approx) > 1L) {
              .ph_abort(sprintf("multiple columns matching '%s'.", paired_col), bullets = approx)
            } else {
              .ph_abort(sprintf("paired column '%s' not found in input.", paired_col))
            }
          }
        }
      }

      # --- DuckDB connection ------------------------------------------------
      con <- NULL
      if (inherits(x, "phip_data")) con <- tryCatch(x$meta$con, error = function(...) NULL)
      if (is.null(con) && inherits(df_long, "tbl_sql")) con <- dbplyr::remote_con(df_long)

      if (is.null(con)) {
        if (!requireNamespace("duckdb", quietly = TRUE))
          .ph_abort("no DuckDB connection found and package {duckdb} is not installed.")
        con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
        on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
        tmp_name <- paste0("ph_tmp_long_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        DBI::dbWriteTable(con, tmp_name, tibble::as_tibble(df_long), temporary = TRUE)
        df_long <- dplyr::tbl(con, tmp_name)
      }

      view_const <- if (inherits(x, "phip_data")) {
        attr(x, "view") %||% (x$meta$view %||% NA_character_)
      } else {
        NA_character_
      }

      # --- library / rank columns -------------------------------------------
      ranks_needing_lib <- setdiff(rank_cols, "peptide_id")

      if (length(ranks_needing_lib)) {
        if (is.null(lib_handle)) .ph_abort("peptide_library required for non-peptide ranks (not found).")
        lib_cols <- if (inherits(lib_handle, "tbl_sql")) {
          colnames(dplyr::collect(head(lib_handle, 0)))
        } else {
          colnames(lib_handle)
        }
        miss_tax <- setdiff(c("peptide_id", ranks_needing_lib), lib_cols)
        if (length(miss_tax))
          .ph_abort("requested rank columns missing from peptide_library.", bullets = paste("-", miss_tax))
      }

      lib_tbl_for_join <- NULL
      if (length(ranks_needing_lib)) {
        if (inherits(lib_handle, "tbl_sql")) {
          lib_tbl_for_join <- lib_handle
        } else {
          lib_name <- paste0("ph_tmp_lib_", format(Sys.time(), "%Y%m%d_%H%M%S"))
          DBI::dbWriteTable(con, lib_name, tibble::as_tibble(lib_handle), temporary = TRUE)
          lib_tbl_for_join <- dplyr::tbl(con, lib_name)
        }
      }

      df_ranked <- df_long
      if (length(ranks_needing_lib)) {
        lib_min <- lib_tbl_for_join |>
          dplyr::select(tidyselect::all_of(c("peptide_id", ranks_needing_lib))) |>
          dplyr::distinct()
        df_ranked <- df_ranked |> dplyr::left_join(lib_min, by = "peptide_id", copy = TRUE)
      }

      available_ranks <- intersect(rank_cols, colnames(df_ranked))
      if (!length(available_ranks)) .ph_abort("none of the requested rank_cols are available.")
      .ph_log_info("ranks resolved",
        bullets = paste("available:", paste(available_ranks, collapse = ", "))
      )

      df_ranked_long <- df_ranked |>
        tidyr::pivot_longer(
          cols      = tidyselect::all_of(available_ranks),
          names_to  = "rank",
          values_to = "feature"
        ) |>
        dplyr::mutate(feature = as.character(feature))

      # --- grouping universes (per-column only) -----------------------------
      gs_view <- df_ranked_long |>
        tidyr::pivot_longer(
          cols      = tidyselect::all_of(group_cols),
          names_to  = "group_col",
          values_to = "group_value"
        ) |>
        dplyr::filter(!is.na(group_value))

      # --- cohort sizes + binary validation ---------------------------------
      .ph_log_info("computing cohort sizes and validating binary group_cols")
      group_sizes <- gs_view |>
        dplyr::distinct(group_col, group_value, sample_id) |>
        dplyr::count(group_col, group_value, name = "N")

      gs_n <- group_sizes |> dplyr::summarise(n = dplyr::n()) |> dplyr::collect()
      if (gs_n$n == 0L) .ph_abort("no non-missing group values; cannot compute denominators.")

      lev_check <- group_sizes |>
        dplyr::count(group_col, name = "k_levels") |>
        dplyr::collect()
      non_binary <- lev_check$group_col[lev_check$k_levels != 2L]
      if (length(non_binary)) {
        .ph_abort(
          "each group_col must have exactly 2 non-missing levels (binary).",
          bullets = sprintf(
            "%s (%d levels)",
            non_binary,
            lev_check$k_levels[match(non_binary, lev_check$group_col)]
          )
        )
      }

      # --- k-of-n presence per (sample, rank, feature) ----------------------
      .ph_log_info("computing presence per sample via k-of-n rule")
      group_by_cols <- c("group_col", "group_value", "sample_id")
      if (!is.null(found_paired_col)) group_by_cols <- c(group_by_cols, found_paired_col)
      group_by_cols <- c(group_by_cols, "rank", "feature")

      k_tbl <- gs_view |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group_by_cols))) |>
        dplyr::summarise(k = sum(!!.sym(exist_col) > 0L), .groups = "drop") |>
        dplyr::mutate(present = k >= !!pop_k_min)

      # ========================= non-paired path ============================
      if (is.null(found_paired_col)) {
        .ph_log_info("counting present samples per feature (pop, unpaired)")
        present_counts <- k_tbl |>
          dplyr::filter(present) |>
          dplyr::distinct(group_col, group_value, rank, feature, sample_id) |>
          dplyr::count(group_col, group_value, rank, feature, name = "n_present")
        features_per_rank <- present_counts |> dplyr::distinct(rank, feature)

        base_grid <- group_sizes |>
          dplyr::select(group_col, group_value) |>
          dplyr::distinct() |>
          dplyr::mutate(.dummy = 1L) |>
          dplyr::inner_join(features_per_rank |> dplyr::mutate(.dummy = 1L), by = ".dummy") |>
          dplyr::select(group_col, group_value, rank, feature)

        stats_long <- base_grid |>
          dplyr::left_join(present_counts, by = c("group_col", "group_value", "rank", "feature")) |>
          dplyr::left_join(group_sizes,    by = c("group_col", "group_value")) |>
          dplyr::mutate(
            n_present = dplyr::coalesce(n_present, 0L),
            prop      = dplyr::if_else(N > 0, n_present / N, NA_real_),
            percent   = 100 * prop
          )

        if (!is.na(view_const))
          stats_long <- stats_long |> dplyr::mutate(view = !!view_const, .before = 1)

        .ph_log_info("building pairwise comparisons")
        by_cols  <- intersect(c("view", "rank", "feature", "group_col"), colnames(stats_long))
        has_view <- "view" %in% colnames(stats_long)

        pairs_joined <- stats_long |>
          dplyr::inner_join(stats_long, by = by_cols, suffix = c("_1", "_2")) |>
          dplyr::filter(group_value_1 != group_value_2) |>
          dplyr::mutate(
            group1  = dplyr::if_else(group_value_1 <= group_value_2, group_value_1, group_value_2),
            group2  = dplyr::if_else(group_value_1 <= group_value_2, group_value_2, group_value_1),
            n1_val  = dplyr::if_else(group_value_1 <= group_value_2, n_present_1, n_present_2),
            n1_tot  = dplyr::if_else(group_value_1 <= group_value_2, N_1,         N_2),
            p1_val  = dplyr::if_else(group_value_1 <= group_value_2, prop_1,      prop_2),
            pct1val = dplyr::if_else(group_value_1 <= group_value_2, percent_1,   percent_2),
            n2_val  = dplyr::if_else(group_value_1 <= group_value_2, n_present_2, n_present_1),
            n2_tot  = dplyr::if_else(group_value_1 <= group_value_2, N_2,         N_1),
            p2_val  = dplyr::if_else(group_value_1 <= group_value_2, prop_2,      prop_1),
            pct2val = dplyr::if_else(group_value_1 <= group_value_2, percent_2,   percent_1)
          ) |>
          dplyr::filter(group_value_1 == group1)

        pairs_lazy <- (if (has_view) {
          pairs_joined |>
            dplyr::transmute(
              view, rank, feature, group_col,
              group1, n1 = n1_val, n1_tot, prop1 = p1_val, percent1 = pct1val,
              group2, n2 = n2_val, n2_tot, prop2 = p2_val, percent2 = pct2val
            )
        } else {
          pairs_joined |>
            dplyr::transmute(
              rank, feature, group_col,
              group1, n1 = n1_val, n1_tot, prop1 = p1_val, percent1 = pct1val,
              group2, n2 = n2_val, n2_tot, prop2 = p2_val, percent2 = pct2val
            )
        }) |>
          dbplyr::window_order(group_col, rank, feature, group1, group2) |>
          dplyr::mutate(row_id = dplyr::row_number()) |>
          dplyr::filter(n1_tot > 0, n2_tot > 0) |>
          dplyr::mutate(
            prop1_eps   = dplyr::if_else(n1 == 0L, (n1 + 0.5) / n1_tot, prop1),
            prop2_eps   = dplyr::if_else(n2 == 0L, (n2 + 0.5) / n2_tot, prop2),
            ratio       = prop1_eps / prop2_eps,
            d1          = dplyr::if_else(n1 == 0L, (n1 + 1.0) / n1_tot, prop1),
            d2          = dplyr::if_else(n2 == 0L, (n2 + 1.0) / n2_tot, prop2),
            delta_ratio = dplyr::if_else(d1 >= d2, d1 / d2 - 1, -(d2 / d1 - 1))
          )

        # materialize for chunked Fisher p-value computation
        tbl_name <- paste0("ph_pop_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        tbl_q    <- .q(con, tbl_name)
        sql_body <- dbplyr::sql_render(pairs_lazy, con)
        DBI::dbExecute(con, paste0("drop table if exists ", tbl_q))
        DBI::dbExecute(con, paste0("create table ", tbl_q, " as ", sql_body))
        DBI::dbExecute(con, paste0("alter table ", tbl_q, " add column if not exists p_raw double"))

        .ph_log_ok("materialized; computing Fisher p-values",
          bullets = sprintf("table: %s", tbl_name)
        )

        p_core_fisher <- function(n1, N1, n2, N2) {
          n1 <- as.double(n1); N1 <- as.double(N1)
          n2 <- as.double(n2); N2 <- as.double(N2)
          if (any(!is.finite(c(n1, N1, n2, N2))) || N1 <= 0 || N2 <= 0) return(NA_real_)
          a   <- max(0, min(n1, N1)); b <- max(0, min(n2, N2))
          out <- try(stats::fisher.test(matrix(c(a, N1 - a, b, N2 - b), 2, byrow = TRUE))$p.value,
                     silent = TRUE)
          if (inherits(out, "try-error")) NA_real_ else as.numeric(out)
        }

        use_parallel <- requireNamespace("future",       quietly = TRUE) &&
                        requireNamespace("future.apply", quietly = TRUE) &&
                        tryCatch(future::nbrOfWorkers() > 1L, error = function(...) FALSE)

        tmp_q <- .q(con, paste0(tbl_name, "_ptmp"))
        DBI::dbExecute(con, paste0("create temp table ", tmp_q, " (row_id bigint, p_raw double)"))

        rs <- DBI::dbSendQuery(con, paste0(
          "select row_id, n1, n1_tot, n2, n2_tot from ", tbl_q, " order by row_id"
        ))
        on.exit(try(DBI::dbClearResult(rs), silent = TRUE), add = TRUE)

        repeat {
          chunk <- DBI::dbFetch(rs, n = chunk_n)
          if (!nrow(chunk)) break
          colnames(chunk) <- tolower(colnames(chunk))
          chunk$row_id <- as.numeric(chunk$row_id)
          n <- nrow(chunk)

          if (use_parallel) {
            splits <- split(seq_len(n),
              cut(seq_len(n), breaks = future::nbrOfWorkers(), labels = FALSE))
            parts <- future.apply::future_lapply(
              splits,
              FUN = function(ii, dat) {
                vapply(ii, function(i)
                  p_core_fisher(dat$n1[i], dat$n1_tot[i], dat$n2[i], dat$n2_tot[i]),
                  numeric(1))
              },
              dat = chunk,
              future.seed = TRUE, future.scheduling = 1, future.globals = FALSE
            )
            pvals <- unlist(parts, use.names = FALSE)
          } else {
            pvals <- vapply(seq_len(n), function(i)
              p_core_fisher(chunk$n1[i], chunk$n1_tot[i], chunk$n2[i], chunk$n2_tot[i]),
              numeric(1))
          }

          DBI::dbAppendTable(con,
            name  = gsub('^\"|\"$', "", tmp_q),
            value = data.frame(row_id = chunk$row_id, p_raw = as.numeric(pvals))
          )
        }

        DBI::dbExecute(con, paste0(
          "update ", tbl_q, " t set p_raw = s.p_raw from ", tmp_q, " s where t.row_id = s.row_id"
        ))
        DBI::dbExecute(con, paste0("drop table ", tmp_q))

        out_df <- dplyr::tbl(con, tbl_name) |>
          dplyr::arrange(rank, feature, group_col, group1, group2) |>
          dplyr::select(
            tidyselect::any_of("view"), rank, feature, group_col,
            group1, n1, N1 = n1_tot, prop1, percent1,
            group2, n2, N2 = n2_tot, prop2, percent2,
            ratio, delta_ratio, p_raw
          ) |>
          dplyr::collect() |>
          as.data.frame()

        .ph_log_ok("done (compute_pop, unpaired)",
          bullets = c(
            sprintf("rows  : %d", nrow(out_df)),
            sprintf("ranks : %s", paste(unique(out_df$rank), collapse = ", ")),
            sprintf("k_min : %d", pop_k_min)
          )
        )
        return(out_df)
      }

      # ========================= paired path ================================
      .ph_log_info("paired design: running McNemar exact (binomial)")

      present_counts <- k_tbl |>
        dplyr::filter(present) |>
        dplyr::distinct(group_col, group_value, rank, feature, sample_id) |>
        dplyr::count(group_col, group_value, rank, feature, name = "n_present")
      features_per_rank <- present_counts |> dplyr::distinct(rank, feature)

      base_grid <- group_sizes |>
        dplyr::select(group_col, group_value) |>
        dplyr::distinct() |>
        dplyr::mutate(.dummy = 1L) |>
        dplyr::inner_join(features_per_rank |> dplyr::mutate(.dummy = 1L), by = ".dummy") |>
        dplyr::select(group_col, group_value, rank, feature)

      stats_long <- base_grid |>
        dplyr::left_join(present_counts, by = c("group_col", "group_value", "rank", "feature")) |>
        dplyr::left_join(group_sizes,    by = c("group_col", "group_value")) |>
        dplyr::mutate(
          n_present = dplyr::coalesce(n_present, 0L),
          prop      = dplyr::if_else(N > 0, n_present / N, NA_real_),
          percent   = 100 * prop
        )

      if (!is.na(view_const))
        stats_long <- stats_long |> dplyr::mutate(view = !!view_const, .before = 1)

      by_cols  <- intersect(c("view", "rank", "feature", "group_col"), colnames(stats_long))
      has_view <- "view" %in% colnames(stats_long)

      pairs_joined <- stats_long |>
        dplyr::inner_join(stats_long, by = by_cols, suffix = c("_1", "_2")) |>
        dplyr::filter(group_value_1 != group_value_2) |>
        dplyr::mutate(
          group1   = dplyr::if_else(group_value_1 < group_value_2, group_value_1, group_value_2),
          group2   = dplyr::if_else(group_value_1 < group_value_2, group_value_2, group_value_1),
          n1       = dplyr::if_else(group_value_1 < group_value_2, n_present_1, n_present_2),
          n2       = dplyr::if_else(group_value_1 < group_value_2, n_present_2, n_present_1),
          prop1    = dplyr::if_else(group_value_1 < group_value_2, prop_1,    prop_2),
          prop2    = dplyr::if_else(group_value_1 < group_value_2, prop_2,    prop_1),
          percent1 = dplyr::if_else(group_value_1 < group_value_2, percent_1, percent_2),
          percent2 = dplyr::if_else(group_value_1 < group_value_2, percent_2, percent_1)
        ) |>
        dplyr::filter(group_value_1 == group1) |>
        (\(d) if (has_view) {
          dplyr::transmute(d, view, rank, feature, group_col, group1, group2,
                           n1, n2, prop1, percent1, prop2, percent2)
        } else {
          dplyr::transmute(d, rank, feature, group_col, group1, group2,
                           n1, n2, prop1, percent1, prop2, percent2)
        })() |>
        dplyr::distinct()

      sdat <- k_tbl |>
        dplyr::select(group_col, group_value, rank, feature,
                      dplyr::all_of(found_paired_col), present) |>
        dplyr::collect() |>
        dplyr::rename(subject_id = !!rlang::sym(found_paired_col))

      disc <- sdat |>
        dplyr::inner_join(sdat,
          by     = c("group_col", "rank", "feature", "subject_id"),
          suffix = c("_1", "_2")
        ) |>
        dplyr::filter(group_value_1 < group_value_2) |>
        dplyr::group_by(group_col, rank, feature,
                        group1 = group_value_1, group2 = group_value_2) |>
        dplyr::summarise(
          n01 = sum((!present_1) & (present_2), na.rm = TRUE),
          n10 = sum(( present_1) & (!present_2), na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          p_raw = {
            n_vec <- n01 + n10
            vapply(seq_along(n_vec), function(i) {
              n <- n_vec[i]
              if (is.na(n) || n <= 0) NA_real_
              else stats::binom.test(n01[i], n, alternative = "two.sided")$p.value
            }, numeric(1))
          }
        )

      paired_summary <- sdat |>
        dplyr::inner_join(sdat,
          by     = c("group_col", "rank", "feature", "subject_id"),
          suffix = c("_1", "_2")
        ) |>
        dplyr::filter(group_value_1 < group_value_2) |>
        dplyr::group_by(group_col, rank, feature,
                        group1 = group_value_1, group2 = group_value_2) |>
        dplyr::summarise(
          N_paired  = dplyr::n_distinct(subject_id),
          n1_paired = sum(present_1, na.rm = TRUE),
          n2_paired = sum(present_2, na.rm = TRUE),
          .groups   = "drop"
        ) |>
        dplyr::mutate(
          prop1_paired = dplyr::if_else(N_paired > 0, n1_paired / N_paired, NA_real_),
          prop2_paired = dplyr::if_else(N_paired > 0, n2_paired / N_paired, NA_real_)
        )

      out_df <- pairs_joined |>
        dplyr::collect() |>
        dplyr::left_join(disc,           by = c("group_col", "rank", "feature", "group1", "group2")) |>
        dplyr::left_join(paired_summary, by = c("group_col", "rank", "feature", "group1", "group2")) |>
        dplyr::mutate(
          N1       = N_paired, N2 = N_paired,
          n1       = dplyr::coalesce(n1_paired, n1),
          n2       = dplyr::coalesce(n2_paired, n2),
          prop1    = dplyr::coalesce(prop1_paired, prop1),
          prop2    = dplyr::coalesce(prop2_paired, prop2),
          percent1 = 100 * prop1,
          percent2 = 100 * prop2
        ) |>
        dplyr::arrange(rank, feature, group_col, group1, group2) |>
        dplyr::select(
          tidyselect::any_of("view"), rank, feature, group_col,
          group1, n1, N1, prop1, percent1,
          group2, n2, N2, prop2, percent2,
          p_raw
        ) |>
        as.data.frame()

      .ph_log_ok("done (compute_pop, paired)",
        bullets = c(
          sprintf("rows  : %d", nrow(out_df)),
          sprintf("ranks : %s", paste(unique(out_df$rank), collapse = ", ")),
          sprintf("k_min : %d", pop_k_min)
        )
      )
      out_df
    }
  )
}
