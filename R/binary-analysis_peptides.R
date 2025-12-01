#' Calculate peptide prevalence per group (phiper-ready; DuckDB/Arrow friendly)
#'
#' @param x phip_data object
#' @param group_cols character vector of grouping columns (must be present in x$data_long)
#' @param prevalence_threshold numeric; keep peptides reaching this % in ≥1 group
#' @param join_meta logical; join peptide metadata from x$peptide_tbl if available
#' @param exist_col optional character; presence/absence column name. If NULL, auto-detects "exist" or "present".
#' @param return "wide" or "long"
#'
#' @return named list of tibbles, one per group_col
#' @export
phip_calc_prevalence <- function(x,
                                 group_cols,
                                 prevalence_threshold = 0,
                                 join_meta = TRUE,
                                 exist_col = NULL,
                                 return = c("wide", "long")) {
  stopifnot(inherits(x, "phip_data"))
  return <- match.arg(return)

  dl <- x$data_long
  vars <- tryCatch(dplyr::tbl_vars(dl), error = function(e) colnames(dl))

  # presence column autodetect
  if (is.null(exist_col)) {
    if ("exist" %in% vars) {
      exist_col <- "exist"
    } else if ("present" %in% vars) {
      exist_col <- "present"
    } else {
      stop("Couldn't find a presence column. Tried 'exist' and 'present'. Available: ", paste(vars, collapse = ", "))
    }
  } else if (!exist_col %in% vars) {
    stop("exist_col '", exist_col, "' not found. Available: ", paste(vars, collapse = ", "))
  }
  evar <- rlang::sym(exist_col)

  if (!all(group_cols %in% vars)) {
    missing <- setdiff(group_cols, vars)
    stop(
      "These group_cols are not present on x$data_long: ",
      paste(missing, collapse = ", "), "\nAvailable: ", paste(vars, collapse = ", ")
    )
  }

  pep_meta <- if (isTRUE(join_meta) && !is.null(x$peptide_tbl)) x$peptide_tbl else NULL
  meta_cols <- if (!is.null(pep_meta)) intersect(c("peptide_id", "Description", "class", "order", "family", "genus", "species"), names(pep_meta)) else character(0)

  .pivot_wider_safe <- function(tbl, ...) {
    tryCatch(tidyr::pivot_wider(tbl, ...),
      error = function(e) tidyr::pivot_wider(dplyr::collect(tbl), ...)
    )
  }

  res_list <- lapply(group_cols, function(gcol) {
    # 1) counts, N, percent per (peptide, group) — compute Percent directly (dbplyr-safe)
    df_long <- dl %>%
      dplyr::filter(!is.na(.data[[gcol]])) %>%
      dplyr::group_by(.data$peptide_id, Group = .data[[gcol]]) %>%
      dplyr::summarise(
        count   = sum(!!evar, na.rm = TRUE),
        N       = dplyr::n_distinct(.data$sample_id),
        Percent = 100 * sum(!!evar, na.rm = TRUE) / dplyr::n_distinct(.data$sample_id),
        .groups = "drop"
      )

    # 2) hard prevalence filter
    if (!is.null(prevalence_threshold) && prevalence_threshold > 0) {
      df_long <- df_long %>%
        dplyr::group_by(.data$peptide_id) %>%
        dplyr::filter(max(.data$Percent, na.rm = TRUE) >= prevalence_threshold) %>%
        dplyr::ungroup()
    }

    # 3) optional metadata
    if (length(meta_cols)) {
      df_long <- df_long %>%
        dplyr::left_join(dplyr::select(pep_meta, dplyr::all_of(meta_cols)), by = "peptide_id")
    }

    if (return == "long") {
      df_long %>%
        dplyr::relocate(.data$peptide_id, dplyr::all_of(meta_cols), .before = .data$Group) %>%
        dplyr::arrange(.data$peptide_id, .data$Group)
    } else {
      wide_core <- df_long %>%
        dplyr::select(.data$peptide_id, .data$Group, .data$count, .data$Percent) %>%
        .pivot_wider_safe(
          names_from  = .data$Group,
          values_from = c(.data$count, .data$Percent),
          names_glue  = "{Group}_{.value}"
        ) %>%
        dplyr::rename_with(~ sub("_Percent$", "", .x), tidyselect::ends_with("_Percent"))

      if (length(meta_cols)) {
        wide_core <- df_long %>%
          dplyr::distinct(.data$peptide_id, dplyr::across(dplyr::all_of(meta_cols))) %>%
          dplyr::right_join(wide_core, by = "peptide_id")
      }

      wide_core %>%
        dplyr::relocate(.data$peptide_id, dplyr::all_of(meta_cols)) %>%
        dplyr::arrange(.data$peptide_id)
    }
  })

  names(res_list) <- group_cols
  res_list
}

#' Pairwise prevalence tests (Fisher + ratios) for each grouping variable
#'
#' @param x           phip_data object (used to compute per-group N = n_distinct(sample_id))
#' @param prev_list   named list from phip_calc_prevalence(..., return = "wide")
#' @param group_cols  character vector; defaults to names(prev_list)
#' @param prevalence_threshold numeric; optional pairwise filter (keep peptide if >= threshold in ≥1 of the two groups)
#' @param p_adjust    method for p-value adjustment ("BH" default)
#' @param epsilon_prop  numeric; pseudocount for proportion ratio (applied only to zeros), default 0.5
#' @param epsilon_delta numeric; pseudocount for delta-ratio (applied only to zeros), default 1
#' @param alternative  fisher.test alternative ("two.sided", "greater", "less")
#'
#' @return nested list: \code{results[[group_col]][["G1_vs_G2"]]} = list(comparison_df=..., N=c(N1,N2), groups=c(group1,g2), group_col=gcol)
#' @export
phip_test_prevalence <- function(
  x,
  prev_list,
  group_cols = names(prev_list),
  prevalence_threshold = NULL,
  p_adjust = "BH",
  epsilon_prop = 0.5,
  epsilon_delta = 1,
  alternative = "two.sided"
) {
  stopifnot(inherits(x, "phip_data"))
  dl <- x$data_long

  results <- list()

  for (gcol in group_cols) {
    prev_wide <- prev_list[[gcol]]
    if (is.null(prev_wide)) next

    # collect wide table if it's lazy
    prev_wide <- tryCatch(dplyr::collect(prev_wide), error = function(e) prev_wide)

    # group levels present in data
    levs <- dl %>%
      dplyr::filter(!is.na(.data[[gcol]])) %>%
      dplyr::distinct(Group = .data[[gcol]]) %>%
      dplyr::arrange(Group) %>%
      dplyr::pull(Group)
    if (length(levs) < 2L) next

    # per-group N
    N_tbl <- dl %>%
      dplyr::filter(!is.na(.data[[gcol]])) %>%
      dplyr::group_by(Group = .data[[gcol]]) %>%
      dplyr::summarise(N = dplyr::n_distinct(.data$sample_id), .groups = "drop") %>%
      dplyr::collect()
    N_map <- stats::setNames(N_tbl$N, N_tbl$Group)

    cols_wide <- colnames(prev_wide)
    count_cols <- grep("_count$", cols_wide, value = TRUE)
    levs_present <- intersect(levs, sub("_count$", "", count_cols))
    if (length(levs_present) < 2L) next

    pairs <- utils::combn(levs_present, 2, simplify = FALSE)
    pair_res <- list()

    for (pair in pairs) {
      g1 <- pair[1]
      g2 <- pair[2]
      cnt1 <- paste0(g1, "_count")
      cnt2 <- paste0(g2, "_count")
      if (!all(c(cnt1, cnt2, g1, g2) %in% names(prev_wide))) next

      N1 <- as.integer(N_map[g1])
      N2 <- as.integer(N_map[g2])
      if (is.na(N1) || is.na(N2) || N1 <= 0L || N2 <= 0L) next

      # pairwise subset (apply optional pairwise prevalence filter)
      df_pair <- prev_wide
      if (!is.null(prevalence_threshold) && prevalence_threshold > 0) {
        df_pair <- df_pair %>%
          dplyr::filter((.data[[g1]] >= prevalence_threshold) | (.data[[g2]] >= prevalence_threshold))
      }
      df_pair <- tryCatch(dplyr::collect(df_pair), error = function(e) df_pair)
      if (!nrow(df_pair)) next

      # counts (coalesce NA -> 0)
      v1 <- df_pair[[cnt1]]
      v1[is.na(v1)] <- 0L
      v2 <- df_pair[[cnt2]]
      v2[is.na(v2)] <- 0L

      # Fisher per peptide
      pvals <- vapply(seq_len(nrow(df_pair)), function(i) {
        a <- as.integer(v1[i])
        b <- as.integer(v2[i])
        mat <- matrix(
          c(
            a + 1, N1 - a + 1,
            b + 1, N2 - b + 1
          ),
          nrow = 2, byrow = TRUE
        )
        suppressWarnings(stats::fisher.test(mat, alternative = alternative)$p.value)
      }, numeric(1))

      # Ratios (pseudocounts only when zero)
      v1_prop <- (v1 + ifelse(v1 == 0, epsilon_prop, 0)) / N1
      v2_prop <- (v2 + ifelse(v2 == 0, epsilon_prop, 0)) / N2
      ratio <- v1_prop / v2_prop
      ratio[!is.finite(ratio)] <- NA_real_

      v1d <- (v1 + ifelse(v1 == 0, epsilon_delta, 0)) / N1
      v2d <- (v2 + ifelse(v2 == 0, epsilon_delta, 0)) / N2
      dr <- ifelse(v1d >= v2d, v1d / v2d - 1, -(v2d / v1d - 1))
      dr[!is.finite(dr)] <- 0

      padj <- stats::p.adjust(pvals, method = p_adjust)
      sig_raw <- pvals < 0.05
      sig_bh <- padj < 0.05

      categories <- dplyr::case_when(
        sig_bh ~ "significant post FDR correction",
        sig_raw & !sig_bh ~ "significant prior correction",
        TRUE ~ "not significant"
      )

      # build output **from df_pair** (not prev_wide)
      keep_meta <- intersect(
        c("Description", "class", "order", "family", "genus", "species"),
        names(df_pair)
      )
      comparison_df <- dplyr::tibble(Peptide = df_pair$peptide_id)
      if (length(keep_meta)) {
        comparison_df <- dplyr::bind_cols(
          comparison_df,
          df_pair[, keep_meta, drop = FALSE]
        )
      }
      comparison_df <- dplyr::bind_cols(
        comparison_df,
        df_pair[, c(g1, g2, cnt1, cnt2), drop = FALSE]
      ) %>%
        dplyr::mutate(
          Delta_ratio = dr,
          ratio = ratio,
          pvals_not_adj = pvals,
          passed_not_adj = sig_raw,
          pvals_bh = padj,
          passed_bh = sig_bh,
          categories = categories
        ) %>%
        dplyr::arrange(
          Delta_ratio,
          dplyr::desc(passed_bh),
          pvals_bh,
          pvals_not_adj
        )

      pair_res[[paste0(g1, "_vs_", g2)]] <- list(
        comparison_df = comparison_df,
        N = c(N1 = N1, N2 = N2),
        groups = c(group1 = g1, group2 = g2),
        group_col = gcol
      )
    }

    results[[gcol]] <- pair_res
  }

  results
}

#' Plot prevalence comparison scatter (standalone)
#'
#' @param comparison One element from \code{phip_test_prevalence()[[gcol]][[pair]]}
#' @param highlight_cols   Optional character vector of boolean columns in comparison_df to highlight
#' @param highlight_colors Optional named vector (names=highlight_cols)
#' @param default_color    Color for non-highlighted points
#' @param significant_colors Named colors used when no highlight_cols are given
#' @param interactive      TRUE -> plotly, FALSE -> ggplot2
#'
#' @return plotly object (interactive=TRUE) or ggplot
#' @export
phip_plot_prevalence_scatter <- function(
  comparison,
  highlight_cols = NULL,
  highlight_colors = NULL,
  default_color = "gray70",
  significant_colors = c(
    "not significant"                 = "dodgerblue",
    "significant prior correction"    = "forestgreen",
    "significant post FDR correction" = "firebrick"
  ),
  interactive = TRUE
) {
  # ---- basic checks ----------------------------------------------------------
  stopifnot(is.list(comparison), !is.null(comparison$comparison_df))
  df <- comparison$comparison_df
  g1 <- unname(comparison$groups[[1]])
  g2 <- unname(comparison$groups[[2]])
  N <- unname(comparison$N)

  if (!all(c(g1, g2) %in% names(df))) {
    stop("Axis columns '", g1, "' or '", g2, "' not found in comparison_df.")
  }
  cnt1 <- paste0(g1, "_count")
  cnt2 <- paste0(g2, "_count")
  if (!all(c(cnt1, cnt2) %in% names(df))) {
    stop("Count columns '", cnt1, "' or '", cnt2, "' not found in comparison_df.")
  }

  has_desc <- "Description" %in% names(df)
  has_species <- "species" %in% names(df)

  # Ensure finite plotting values
  df_plot <- df %>%
    dplyr::mutate(
      log2ratio = suppressWarnings(log2(.data$ratio))
    ) %>%
    dplyr::filter(
      is.finite(.data[[g1]]),
      is.finite(.data[[g2]])
    )

  # ---- highlighting logic ----------------------------------------------------
  use_highlight <- !is.null(highlight_cols) && length(highlight_cols) > 0
  if (use_highlight) {
    missing <- setdiff(highlight_cols, names(df_plot))
    if (length(missing)) {
      stop("These highlight_cols are not in the data: ", paste(missing, collapse = ", "))
    }

    df_plot <- df_plot %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        .trues    = list(highlight_cols[as.logical(c_across(dplyr::all_of(highlight_cols)))]),
        highlight = if (length(.trues) == 0) "none" else .trues[[1]]
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.trues)

    df_plot$highlight <- forcats::fct_explicit_na(df_plot$highlight, na_level = "none")
    df_plot$highlight <- factor(df_plot$highlight, levels = c("none", highlight_cols))

    # legend labels via Wilcoxon tests (BH-adjusted) comparing flagged vs non-flagged log2ratio
    pvals <- vapply(highlight_cols, function(flag) {
      x <- df_plot %>%
        dplyr::filter(.data[[flag]] %in% TRUE) %>%
        dplyr::pull(.data$log2ratio)
      y <- df_plot %>%
        dplyr::filter(!(.data[[flag]] %in% TRUE)) %>%
        dplyr::pull(.data$log2ratio)
      if (length(x) == 0 || length(y) == 0) {
        return(NA_real_)
      }
      stats::wilcox.test(x, y)$p.value
    }, numeric(1))
    pvals_adj <- stats::p.adjust(pvals, method = "BH")
    fmt_p <- ifelse(is.na(pvals_adj), "NA", format.pval(pvals_adj, digits = 1, eps = 0.001))
    legend_labels <- setNames(paste0(highlight_cols, " (P=", fmt_p, ")"), highlight_cols)

    # colors
    if (is.null(highlight_colors)) {
      pal <- grDevices::hcl.colors(max(length(highlight_cols), 8), palette = "Set2")
      highlight_colors <- stats::setNames(pal[seq_along(highlight_cols)], highlight_cols)
    }
    manual_vals <- c(none = default_color, highlight_colors)

    # tooltip (only for highlighted points; 'none' -> NA)
    df_plot <- df_plot %>%
      dplyr::mutate(
        tooltip_txt = dplyr::if_else(
          .data$highlight == "none",
          NA_character_,
          paste0(
            "Peptide: ", .data$Peptide, "<br>",
            if (has_desc) paste0("Desc: ", .data$Description, "<br>") else "",
            if (has_species) paste0("Species: ", .data$species, "<br>") else "",
            g1, ": ", .data[[g1]], " / ",
            g2, ": ", .data[[g2]], "<br>",
            "Highlight: ", as.character(.data$highlight)
          )
        )
      ) %>%
      dplyr::arrange(.data$highlight) # draw greys first

    color_aes <- ggplot2::aes(color = .data$highlight, text = .data$tooltip_txt)
    color_scale <- ggplot2::scale_color_manual(
      name   = NULL,
      values = manual_vals,
      breaks = highlight_cols,
      labels = unname(legend_labels)
    )
    show_legend <- TRUE
    legend_theme <- ggplot2::theme(
      legend.position      = c(0, 1),
      legend.justification = c(0, 1),
      legend.background    = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.8), color = "gray80"),
      legend.key.size      = grid::unit(10, "pt"),
      legend.text          = ggplot2::element_text(size = 9),
      legend.title         = ggplot2::element_text(size = 9, face = "bold")
    )
  } else {
    # fallback to significance categories
    if (!"categories" %in% names(df_plot)) {
      stop("Column 'categories' not found; needed when highlight_cols is missing.")
    }
    df_plot <- df_plot %>%
      dplyr::mutate(
        tooltip_txt = dplyr::if_else(
          .data$categories == "not significant",
          NA_character_,
          paste0(
            "Peptide: ", .data$Peptide, "<br>",
            if (has_desc) paste0("Desc: ", .data$Description, "<br>") else "",
            if (has_species) paste0("Species: ", .data$species, "<br>") else "",
            g1, ": ", .data[[g1]], " / ",
            g2, ": ", .data[[g2]]
          )
        )
      )
    color_aes <- ggplot2::aes(color = .data$categories, text = .data$tooltip_txt)
    color_scale <- ggplot2::scale_color_manual(
      values = significant_colors, name = NULL,
      labels = c("ns", "significant", "significant FDR")
    )
    show_legend <- FALSE
    legend_theme <- ggplot2::theme(legend.position = "none")
  }

  # ---- build plot ------------------------------------------------------------
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data[[g1]], y = .data[[g2]])) +
    ggplot2::geom_point(color_aes, alpha = 0.65) +
    color_scale +
    ggplot2::labs(
      x = paste0("% ", g1, " in whom\na peptide is significantly bound\n(n = ", N[1], ")"),
      y = paste0("% ", g2, " in whom\na peptide is significantly bound\n(n = ", N[2], ")")
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(colour = "black", fill = NA),
      plot.margin      = ggplot2::margin(t = 10, r = 15, b = 15, l = 10, unit = "pt"),
      axis.text.y.left = ggplot2::element_text(size = 10, face = "italic"),
      axis.title.x     = ggplot2::element_text(face = "bold"),
      axis.title.y     = ggplot2::element_text(face = "bold")
    ) +
    legend_theme

  # ---- interactive wrapper (optional) ----------------------------------------
  if (isTRUE(interactive) && requireNamespace("plotly", quietly = TRUE)) {
    pl <- plotly::ggplotly(p, tooltip = "text", width = 550, height = 550)

    if (use_highlight) {
      # hide 'none' in legend and relabel traces to BH-adjusted labels
      for (i in seq_along(pl$x$data)) {
        tr <- pl$x$data[[i]]
        if (!is.null(tr$name)) {
          if (tr$name %in% c("none")) {
            tr$showlegend <- FALSE
          } else if (tr$name %in% names(legend_labels)) {
            tr$name <- unname(legend_labels[[tr$name]])
          }
        }
        pl$x$data[[i]] <- tr
      }
    }

    pl <- plotly::layout(
      pl,
      showlegend = show_legend,
      legend = list(
        x = 0, xanchor = "left",
        y = 1, yanchor = "top",
        font = list(size = nine <- 9) # tiny helper to avoid magic numbers
      ),
      margin = list(l = 80, r = 80, b = 80, t = 80, pad = 0),
      hoverlabel = list(font = list(size = 10)),
      xaxis = list(scaleratio = 1, scaleanchor = "y"),
      yaxis = list(scaleratio = 1, scaleanchor = "x")
    )
    return(pl)
  }

  p
}
