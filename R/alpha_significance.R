# ==============================================================================
# compute_alpha_significance() + plot_alpha_significance()
# Statistical comparison of alpha diversity between groups
# ==============================================================================

# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

#' Cohen's d (pooled SD) for two numeric vectors
#' @noRd
.ph_cohens_d <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  nx <- length(x); ny <- length(y)
  if (nx < 2L || ny < 2L) return(NA_real_)
  pooled_sd <- sqrt(((nx - 1L) * stats::var(x) + (ny - 1L) * stats::var(y)) / (nx + ny - 2L))
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  (mean(x) - mean(y)) / pooled_sd
}

#' Convert p-values to significance stars
#' @noRd
.ph_sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

#' Infer the group column from a phip_alpha_diversity object or data frame
#'
#' Uses (in order): supplied value, `group_col` attribute, first column that is
#' not one of the standard alpha-diversity column names.
#' @noRd
.ph_infer_group_col <- function(alpha_df, supplied = NULL) {
  if (!is.null(supplied)) return(supplied)
  # Try attribute first
  attr_gc <- attr(alpha_df, "group_cols", exact = TRUE)
  if (!is.null(attr_gc) && length(attr_gc) == 1L && attr_gc %in% names(alpha_df)) {
    return(attr_gc)
  }
  # Fallback: first non-standard column
  standard_cols <- c(
    "sample_id", "rank",
    "richness", "shannon_diversity", "simpson_diversity",
    "pielou_evenness", "berger_parker_dominance"
  )
  candidates <- setdiff(names(alpha_df), standard_cols)
  if (!length(candidates)) {
    .ph_abort(
      "Cannot infer group column from alpha data.",
      step    = "group_col inference",
      bullets = "Supply group_col explicitly."
    )
  }
  candidates[[1L]]
}

# ------------------------------------------------------------------------------
# compute_alpha_significance()
# ------------------------------------------------------------------------------

#' Compute statistical significance of alpha diversity between groups
#'
#' @description
#' Runs global (Kruskal-Wallis or one-way ANOVA) and pairwise (Wilcoxon or
#' Tukey HSD) hypothesis tests on each `(rank, metric)` combination present in
#' the output of [compute_alpha()]. Pairwise p-values are adjusted with
#' [p.adjust()] and Cohen's d effect sizes are appended.
#'
#' @param x A `"phip_alpha_diversity"` list (output of [compute_alpha()]) or a
#'   single alpha-diversity data frame.
#' @param group_col Character scalar; name of the grouping column. Inferred
#'   from `attr(x, "group_cols")` when `NULL` (default).
#' @param metric Character vector; subset of metrics to test. `NULL` (default)
#'   uses all numeric metric columns present in the data.
#' @param global_test One of `"kruskal"` (Kruskal-Wallis, default) or `"anova"`
#'   (one-way ANOVA).
#' @param pairwise_test One of `"wilcoxon"` (Wilcoxon rank-sum, default) or
#'   `"tukey"` (Tukey HSD from `aov()`).
#' @param p_adjust_method Method passed to [p.adjust()]. Default `"BH"`.
#'
#' @return A named list of class `"phip_alpha_significance"` with two tibbles:
#'   \describe{
#'     \item{`$global`}{One row per `(rank, metric)`: `rank`, `metric`,
#'       `statistic`, `p_value`, `test`.}
#'     \item{`$pairwise`}{One row per `(rank, metric, pair)`: `rank`, `metric`,
#'       `group1`, `group2`, `p_raw`, `p_adj`, `cohens_d`, `stars`, `test`.}
#'   }
#'   Attributes: `group_col`, `global_test`, `pairwise_test`,
#'   `p_adjust_method`, `metrics`, `ranks`.
#'
#' @examples
#' \dontrun{
#' alpha <- compute_alpha(phip_obj, group_cols = "group")
#' sig   <- compute_alpha_significance(alpha)
#' sig$global
#' sig$pairwise
#' }
#' @export
compute_alpha_significance <- function(
    x,
    group_col       = NULL,
    metric          = NULL,
    global_test     = c("kruskal", "anova"),
    pairwise_test   = c("wilcoxon", "tukey"),
    p_adjust_method = "BH"
) {
  global_test   <- match.arg(global_test)
  pairwise_test <- match.arg(pairwise_test)

  # ---- normalise x to a single data frame ------------------------------------
  gc_attr <- NULL
  if (inherits(x, "phip_alpha_diversity") || (is.list(x) && !is.data.frame(x))) {
    gc_attr  <- attr(x, "group_cols", exact = TRUE)
    inner    <- x[vapply(x, inherits, logical(1L), "data.frame")]
    alpha_df <- suppressWarnings(dplyr::bind_rows(inner))
    if (!is.null(gc_attr) && length(gc_attr) >= 1L) {
      attr(alpha_df, "group_cols") <- gc_attr
    }
  } else {
    alpha_df <- tibble::as_tibble(x)
  }

  # ---- infer group column ----------------------------------------------------
  group_col <- .ph_infer_group_col(alpha_df, supplied = group_col)

  .ph_check_cond(
    !group_col %in% names(alpha_df),
    error_message = sprintf("group_col %s not found in alpha data.", .ph_add_quotes(group_col, 1L)),
    step   = "input validation",
    error  = TRUE
  )

  # ---- determine metrics to test ---------------------------------------------
  known_metrics <- c(
    "richness", "shannon_diversity", "simpson_diversity",
    "pielou_evenness", "berger_parker_dominance"
  )
  present_metrics <- intersect(known_metrics, names(alpha_df))
  if (!is.null(metric)) {
    .ph_check_cond(
      !all(metric %in% present_metrics),
      error_message = "Some requested metrics are not present in the alpha data.",
      step   = "metric selection",
      bullets = sprintf("missing: %s",
        paste(.ph_add_quotes(setdiff(metric, present_metrics), 1L), collapse = ", ")),
      error  = TRUE
    )
    present_metrics <- intersect(metric, present_metrics)
  }

  .ph_check_cond(
    length(present_metrics) == 0L,
    error_message = "No recognised metric columns found in alpha data.",
    step  = "metric selection",
    error = TRUE
  )

  # ---- determine ranks -------------------------------------------------------
  ranks_present <- if ("rank" %in% names(alpha_df)) unique(alpha_df$rank) else NA_character_

  # ---- run tests per (rank, metric) ------------------------------------------
  global_rows   <- list()
  pairwise_rows <- list()

  for (rk in ranks_present) {
    if (!is.na(rk)) {
      sub_df <- dplyr::filter(alpha_df, .data$rank == !!rk)
    } else {
      sub_df <- alpha_df
    }

    for (mt in present_metrics) {
      vals  <- sub_df[[mt]]
      grps  <- as.character(sub_df[[group_col]])

      # drop NA pairs
      ok    <- !is.na(vals) & !is.na(grps)
      vals  <- vals[ok]
      grps  <- grps[ok]

      grp_levels <- unique(grps)

      # --- global test --------------------------------------------------------
      g_row <- tryCatch({
        if (global_test == "kruskal") {
          kt <- stats::kruskal.test(vals ~ factor(grps))
          tibble::tibble(
            rank      = rk,
            metric    = mt,
            statistic = kt$statistic[[1L]],
            p_value   = kt$p.value,
            test      = "kruskal-wallis"
          )
        } else {
          at <- stats::aov(vals ~ factor(grps))
          sm <- summary(at)[[1L]]
          tibble::tibble(
            rank      = rk,
            metric    = mt,
            statistic = sm[["F value"]][[1L]],
            p_value   = sm[["Pr(>F)"]][[1L]],
            test      = "anova"
          )
        }
      }, error = function(e) {
        tibble::tibble(rank = rk, metric = mt, statistic = NA_real_,
                       p_value = NA_real_, test = global_test)
      })
      global_rows <- c(global_rows, list(g_row))

      # --- pairwise test ------------------------------------------------------
      if (length(grp_levels) < 2L) next
      pairs <- utils::combn(grp_levels, 2L, simplify = FALSE)
      if (!length(pairs)) next

      pw_list <- lapply(pairs, function(pr) {
        v1 <- vals[grps == pr[[1L]]]
        v2 <- vals[grps == pr[[2L]]]

        p_raw_val <- tryCatch({
          if (pairwise_test == "wilcoxon") {
            wt <- stats::wilcox.test(v1, v2, exact = FALSE)
            wt$p.value
          } else {
            at2 <- stats::aov(vals ~ factor(grps))
            tk  <- stats::TukeyHSD(at2)[[1L]]
            nm  <- paste(sort(c(pr[[2L]], pr[[1L]])), collapse = "-")
            nm2 <- paste(c(pr[[2L]], pr[[1L]]), collapse = "-")
            row <- tk[rownames(tk) == nm | rownames(tk) == nm2, , drop = FALSE]
            if (nrow(row) > 0L) row[1L, "p adj"] else NA_real_
          }
        }, error = function(e) NA_real_)

        tibble::tibble(
          rank    = rk,
          metric  = mt,
          group1  = pr[[1L]],
          group2  = pr[[2L]],
          p_raw   = p_raw_val,
          cohens_d = .ph_cohens_d(v1, v2),
          test    = pairwise_test
        )
      })
      pw_df <- dplyr::bind_rows(pw_list)

      # p-adjust within (rank, metric)
      pw_df$p_adj <- stats::p.adjust(pw_df$p_raw, method = p_adjust_method)
      pw_df$stars <- .ph_sig_stars(pw_df$p_adj)

      pairwise_rows <- c(pairwise_rows, list(pw_df))
    }
  }

  global_tbl   <- dplyr::bind_rows(global_rows)
  pairwise_tbl <- dplyr::bind_rows(pairwise_rows)
  if (nrow(pairwise_tbl)) {
    pairwise_tbl <- dplyr::select(
      pairwise_tbl,
      rank, metric, group1, group2, p_raw, p_adj, cohens_d, stars, test
    )
  }

  out <- structure(
    list(global = global_tbl, pairwise = pairwise_tbl),
    class           = "phip_alpha_significance",
    group_col       = group_col,
    global_test     = global_test,
    pairwise_test   = pairwise_test,
    p_adjust_method = p_adjust_method,
    metrics         = present_metrics,
    ranks           = ranks_present
  )
  out
}

# ------------------------------------------------------------------------------
# plot_alpha_significance()
# ------------------------------------------------------------------------------

#' Visualise alpha diversity significance results
#'
#' @description
#' Summarises the pairwise significance table from [compute_alpha_significance()]
#' either as a filtered tibble (`type = "table"`) or a Cohen's d heatmap with
#' significance annotations (`type = "heatmap"`).
#'
#' @param x A `"phip_alpha_significance"` object from
#'   [compute_alpha_significance()].
#' @param metric Character scalar; metric to display. Defaults to the first
#'   metric present in `x`.
#' @param rank Character scalar; rank to display. Defaults to the first rank
#'   present in `x`. Ignored when the data has no `rank` column.
#' @param type One of `"table"` (default) or `"heatmap"`.
#' @param p_threshold Numeric; pairs with `p_adj > p_threshold` are shown in
#'   grey in the heatmap and excluded from the table. Default `0.05`.
#' @param ... Reserved; ignored.
#'
#' @return
#' - `type = "table"`: a [`tibble`] of significant pairwise results.
#' - `type = "heatmap"`: a `ggplot` heatmap (Cohen's d fill, stars overlaid).
#'
#' @examples
#' \dontrun{
#' sig <- compute_alpha_significance(alpha_list)
#' plot_alpha_significance(sig, type = "table")
#' plot_alpha_significance(sig, type = "heatmap")
#' }
#' @export
plot_alpha_significance <- function(
    x,
    metric      = NULL,
    rank        = NULL,
    type        = c("table", "heatmap"),
    p_threshold = 0.05,
    ...
) {
  type <- match.arg(type)

  .ph_check_cond(
    !inherits(x, "phip_alpha_significance"),
    error_message = "`x` must be a 'phip_alpha_significance' object (output of compute_alpha_significance()).",
    step  = "input validation",
    error = TRUE
  )

  # ---- defaults --------------------------------------------------------------
  if (is.null(metric)) metric <- attr(x, "metrics")[[1L]]
  if (is.null(rank))   rank   <- attr(x, "ranks")[[1L]]

  # ---- filter pairwise table -------------------------------------------------
  pw <- x$pairwise

  if ("metric" %in% names(pw)) {
    pw <- dplyr::filter(pw, .data$metric == !!metric)
  }
  if ("rank" %in% names(pw) && !is.na(rank)) {
    pw <- dplyr::filter(pw, .data$rank == !!rank)
  }

  # ---- table mode ------------------------------------------------------------
  if (type == "table") {
    return(dplyr::filter(pw, .data$p_adj <= !!p_threshold))
  }

  # ---- heatmap mode ----------------------------------------------------------
  # Build symmetric matrix of pairs (upper + lower triangle)
  all_groups <- unique(c(pw$group1, pw$group2))

  # Create symmetric long form
  pw_sym <- dplyr::bind_rows(
    pw,
    dplyr::rename(pw, group1 = group2, group2 = group1)
  )

  # Fill in reverse so the tile is always populated
  pw_full <- tidyr::crossing(group1 = all_groups, group2 = all_groups) |>
    dplyr::left_join(
      dplyr::select(pw_sym, group1, group2, p_adj, cohens_d, stars),
      by = c("group1", "group2")
    ) |>
    dplyr::mutate(
      sig       = !is.na(.data$p_adj) & .data$p_adj <= !!p_threshold,
      label     = dplyr::if_else(.data$sig, .data$stars, ""),
      fill_d    = dplyr::if_else(.data$sig, abs(.data$cohens_d), NA_real_),
      diagonal  = .data$group1 == .data$group2
    )

  p <- ggplot2::ggplot(
    pw_full,
    ggplot2::aes(x = .data$group2, y = .data$group1, fill = .data$fill_d)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      size   = 4,
      colour = "black"
    ) +
    ggplot2::scale_fill_gradient(
      low      = "#d6eaf8",
      high     = "#1a5276",
      na.value = "grey90",
      name     = "|Cohen's d|",
      limits   = c(0, NA)
    ) +
    ggplot2::labs(
      x     = attr(x, "group_col"),
      y     = attr(x, "group_col"),
      title = sprintf("Alpha significance: %s (rank: %s)", metric, rank)
    ) +
    theme_phip() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  p
}
