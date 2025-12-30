#' @title Δ prevalence vs pooled prevalence (per peptide)
#'
#' @description
#' Build a static ggplot showing the per-peptide shift in prevalence
#' (\eqn{\Delta = group2 - group1}) as a function of pooled prevalence
#' (\eqn{(group1 + group2)/2}). The input should be a tibble/data frame produced
#' by \code{ph_compute_prevalence()} or equivalent with columns
#' \code{group1}, \code{group2}, \code{prop1}, and \code{prop2}.
#'
#' @param prev_tbl Data frame with columns \code{group1}, \code{group2},
#'   \code{prop1}, \code{prop2}. Optional \code{feature} is used for row
#'   identity only.
#' @param group_pair_values Optional length-2 character vector
#'   \code{c(group1, group2)}. Use this when \code{prev_tbl} contains multiple
#'   group pairs.
#' @param group_labels Optional length-2 character vector of display labels
#'   \code{c(label_group1, label_group2)}. Defaults to \code{group1}/\code{group2}.
#' @param point_jitter_width,point_jitter_height Jitter amounts for the points.
#'   Defaults 0.005.
#' @param point_alpha Point transparency. Default 0.25.
#' @param point_size Point size. Default 0.6.
#' @param add_smooth Add a GAM smooth curve (\code{mgcv}). Default \code{TRUE}.
#' @param smooth_k Basis dimension \code{k} for the smooth. Default 5.
#' @param arrow_color Color for the directional arrows and labels. Default
#'   \code{"red"}.
#' @param arrow_head_length_mm Arrow head length in mm. Default 4.
#' @param arrow_x_frac Arrow X position as a fraction of the max pooled
#'   prevalence. Default 0.97 (near the right edge).
#' @param plot_title,plot_subtitle Optional plot labels for the title/subtitle.
#' @param x_label,y_label Optional axis labels. Defaults are generated from the
#'   group labels.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' prev_tbl <- make_prev_tbl(x_kid1, rank="is_flagellum", g1="T2", g2="T8", feature_value=TRUE)
#' p <- deltaplot_prevalence(prev_tbl, group_pair_values = c("B","M12"))
#' }
#'
#' @export
deltaplot_prevalence <- function(
    prev_tbl,
    group_pair_values = NULL,
    group_labels     = NULL,
    # Styling
    point_jitter_width  = 0.005,
    point_jitter_height = 0.005,
    point_alpha    = 0.25,
    point_size     = 0.6,
    add_smooth     = TRUE,
    smooth_k       = 5,
    arrow_color    = "red",
    arrow_head_length_mm = 4,
    arrow_x_frac   = 0.97,
    plot_title     = NULL,
    plot_subtitle  = NULL,
    x_label        = NULL,
    y_label        = NULL
) {
  # ---- Basic checks -----------------------------------------------------------
  d <- prev_tbl
  need <- c("group1","group2","prop1","prop2")
  miss <- setdiff(need, names(d))
  if (length(miss)) {
    stop("deltaplot_prevalence(): missing columns from ph_compute_prevalence(): ",
         paste(miss, collapse = ", "))
  }

  # ---- Select exactly one (group1, group2) pair ------------------------------
  pairs <- unique(d[, c("group1","group2")])
  if (!is.null(group_pair_values)) {
    if (length(group_pair_values) != 2L) {
      stop("group_pair_values must be length-2 vector: c(g1, g2).")
    }
    d <- d[d$group1 == group_pair_values[1] &
             d$group2 == group_pair_values[2], , drop = FALSE]
    if (!nrow(d)) stop("No rows left after filtering to group_pair_values = c('",
                       group_pair_values[1], "', '", group_pair_values[2], "').")
    g1_raw <- group_pair_values[1]
    g2_raw <- group_pair_values[2]
  } else {
    if (nrow(pairs) != 1L) {
      stop("Multiple (group1, group2) pairs detected. ",
           "Filter your tibble to one pair or pass group_pair_values = c(g1, g2). ",
           "Found pairs: ",
           paste0(utils::capture.output(print(pairs)), collapse = " "))
    }
    g1_raw <- pairs$group1[1]
    g2_raw <- pairs$group2[1]
  }

  # ---- Labels for display -----------------------------------------------------
  if (is.null(group_labels)) {
    g1_lab <- as.character(g1_raw)
    g2_lab <- as.character(g2_raw)
  } else {
    if (length(group_labels) != 2L) {
      stop("group_labels must be length-2 vector: c(label_g1, label_g2).")
    }
    g1_lab <- group_labels[1]
    g2_lab <- group_labels[2]
  }

  # ---- Compute pooled and delta ----------------------------------------------
  w <- d %>%
    dplyr::transmute(
      id     = if ("feature" %in% names(.)) .data$feature else dplyr::row_number(),
      pooled = (prop1 + prop2) / 2,
      delta  =  prop2 - prop1
    ) %>%
    dplyr::filter(is.finite(pooled), is.finite(delta)) %>%
    dplyr::mutate(
      pooled_clip = pmin(pmax(as.numeric(pooled), 1e-6), 1 - 1e-6)
    )

  if (!nrow(w)) stop("No finite rows to plot after pooled/delta computation.")

  # Arrow near right edge of observed range
  arrow_x <- max(w$pooled_clip, na.rm = TRUE) * arrow_x_frac

  # ---- Plot -------------------------------------------------------------------
  p <- ggplot2::ggplot(w, ggplot2::aes(pooled_clip, delta)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_jitter(alpha = point_alpha,
                         height = point_jitter_height,
                         width  = point_jitter_width,
                         size   = point_size)

  if (isTRUE(add_smooth)) {
    p <- p + ggplot2::geom_smooth(
      method  = "gam",
      formula = y ~ s(x, k = smooth_k),
      se      = FALSE
    )
  }

  p <- p +
    ggplot2::annotate(
      "segment",
      x = arrow_x, xend = arrow_x, y = 0, yend =  0.06,
      colour = arrow_color, linewidth = 0.6,
      arrow = grid::arrow(length = grid::unit(arrow_head_length_mm, "mm"))
    ) +
    ggplot2::annotate(
      "text",
      x = arrow_x, y = 0.065,
      label = paste0("More in ", g2_lab),
      colour = arrow_color, fontface = "bold", vjust = 0
    ) +
    ggplot2::annotate(
      "segment",
      x = arrow_x, xend = arrow_x, y = 0, yend = -0.06,
      colour = arrow_color, linewidth = 0.6,
      arrow = grid::arrow(length = grid::unit(arrow_head_length_mm, "mm"))
    ) +
    ggplot2::annotate(
      "text",
      x = arrow_x, y = -0.065,
      label = paste0("More in ", g1_lab),
      colour = arrow_color, fontface = "bold", vjust = 1
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::labs(
      title    = plot_title %||% paste0("Per-peptide shift vs pooled prevalence (", g2_lab, " - ", g1_lab, ")"),
      subtitle = plot_subtitle,
      x        = x_label %||% paste0("Pooled prevalence (", g1_lab, " & ", g2_lab, ")"),
      y        = y_label %||% paste0("\u0394 prevalence (", g2_lab, " - ", g1_lab, ")")
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 20, 10, 10),
      text = ggplot2::element_text(size = 20)
    )

  p
}

# ==============================================================================
# Interactive Δ prevalence vs pooled prevalence (per peptide)
# Input: tibble from ph_compute_prevalence() with:
#   group1, group2, prop1, prop2, (n1, N1, percent1), (n2, N2, percent2),
#   optional: feature, p_adj_rank_wbh
# ==============================================================================
#' @title Interactive Δ prevalence vs pooled prevalence (per peptide)
#' @description Interactive plotly version of \code{deltaplot_prevalence()}.
#' @param prev_tbl A data frame with columns \code{group1}, \code{group2},
#'   \code{prop1}, \code{prop2} (and optional \code{feature}).
#' @param group_pair_values Optional length-2 vector \code{c(g1, g2)} selecting a pair.
#' @param group_labels Optional length-2 vector of display labels.
#' @param point_alpha Point transparency. Default 0.6.
#' @param point_size Point size. Default 6.
#' @param add_smooth Add a GAM smooth curve. Default \code{TRUE}.
#' @param smooth_k Basis dimension \code{k} for the smooth. Default 5.
#' @param arrow_color Color for arrows and labels. Default \code{"red"}.
#' @param arrow_x_frac Arrow X position as a fraction of the x-range. Default 0.97.
#' @param arrow_length_frac Arrow length as a fraction of the y-range. Default 0.30.
#' @param label_x_gap_frac Horizontal label offset as a fraction of the x-range.
#' @param label_y_gap_frac Vertical label offset as a fraction of the y-range.
#' @param plot_title,plot_subtitle Optional plot labels.
#' @param add_jitter Logical; if \code{TRUE}, apply display-only jitter.
#' @param point_jitter_width,point_jitter_height Jitter amounts. Defaults 0.005.
#' @param jitter_seed Optional seed for jitter reproducibility.
#' @return A plotly object.
#' @export
deltaplot_prevalence_interactive <- function(
    prev_tbl,
    group_pair_values = NULL,
    group_labels     = NULL,
    # Styling
    point_alpha     = 0.6,
    point_size      = 6,
    add_smooth      = TRUE,
    smooth_k        = 5,
    arrow_color     = "red",
    # stała geometria strzałek/etykiet
    arrow_x_frac    = 0.97,  # 0..1 po osi X
    arrow_length_frac = 0.30,  # 0..1 wysokości osi Y
    label_x_gap_frac  = 0.03,  # ile osi X w lewo od strzałki
    label_y_gap_frac  = 0.02,  # pionowy margines etykiet
    plot_title      = NULL,
    plot_subtitle   = NULL,
    # Jitter (display-only)
    add_jitter          = TRUE,
    point_jitter_width  = 0.005,
    point_jitter_height = 0.005,
    jitter_seed         = NULL
){
  need <- c("group1","group2","prop1","prop2")
  miss <- setdiff(need, names(prev_tbl))
  if (length(miss)) stop("deltaplot_prevalence_interactive(): missing: ", paste(miss, collapse=", "))
  d <- prev_tbl

  # wybór pary
  if (!is.null(group_pair_values)) {
    stopifnot(length(group_pair_values) == 2L)
    d <- d[d$group1 == group_pair_values[1] &
             d$group2 == group_pair_values[2], , drop = FALSE]
    if (!nrow(d)) stop("No rows after filtering to group_pair_values.")
    g1_raw <- group_pair_values[1]; g2_raw <- group_pair_values[2]
  } else {
    pairs <- unique(d[, c("group1","group2")])
    if (nrow(pairs) != 1L) stop("Multiple (group1,group2) pairs – pass group_pair_values.")
    g1_raw <- pairs$group1[1]; g2_raw <- pairs$group2[1]
  }
  if (is.null(group_labels)) { g1_lab <- as.character(g1_raw); g2_lab <- as.character(g2_raw)
  } else { stopifnot(length(group_labels) == 2L); g1_lab <- group_labels[1]; g2_lab <- group_labels[2] }

  # kolumny pomocnicze
  if (!("feature" %in% names(d))) {
    d$feature <- if ("peptide_id" %in% names(d)) as.character(d$peptide_id) else as.character(seq_len(nrow(d)))
  }
  if (!("percent1" %in% names(d))) d$percent1 <- d$prop1 * 100
  if (!("percent2" %in% names(d))) d$percent2 <- d$prop2 * 100
  pick_or <- function(df, nm, default = NA_real_) if (nm %in% names(df)) df[[nm]] else default

  w <- d %>%
    dplyr::transmute(
      feature   = .data$feature,
      pooled    = (prop1 + prop2) / 2,
      delta     =  prop2 - prop1,
      n1 = pick_or(d,"n1"), N1 = pick_or(d,"N1"), pct1 = percent1,
      n2 = pick_or(d,"n2"), N2 = pick_or(d,"N2"), pct2 = percent2,
      p_adj_wbh = dplyr::coalesce(
        pick_or(d,"p_adj_rank_wbh"),
        pick_or(d,"padj_wbh"),
        pick_or(d,"p_adj"),
        pick_or(d,"p_raw")
      )
    ) %>%
    dplyr::filter(is.finite(pooled), is.finite(delta)) %>%
    dplyr::mutate(pooled_clip = pmin(pmax(as.numeric(pooled), 1e-6), 1 - 1e-6))
  if (!nrow(w)) stop("No finite rows to plot.")

  # jitter (display)
  if (isTRUE(add_jitter)) {
    if (!is.null(jitter_seed)) { old_seed <- .Random.seed; set.seed(jitter_seed); on.exit({ if (exists("old_seed")) .Random.seed <<- old_seed }, add=TRUE) }
    jx <- stats::runif(nrow(w), -point_jitter_width,  point_jitter_width)
    jy <- stats::runif(nrow(w), -point_jitter_height, point_jitter_height)
    x_jit <- pmin(pmax(w$pooled_clip + jx, 1e-6), 1 - 1e-6)
    y_jit <- w$delta + jy
  } else { x_jit <- w$pooled_clip; y_jit <- w$delta }

  # smooth (opcjonalnie)
  smooth_df <- NULL
  if (isTRUE(add_smooth) && requireNamespace("mgcv", quietly = TRUE)) {
    fit <- mgcv::gam(delta ~ s(pooled_clip, k = smooth_k), data = w)
    xs  <- seq(min(w$pooled_clip), max(w$pooled_clip), length.out = 200)
    smooth_df <- data.frame(pooled_clip = xs,
                            delta = stats::predict(fit, newdata = data.frame(pooled_clip = xs), type = "response"))
  }

  # --- GEOMETRIA STRZAŁEK/ETYKIET ---
  # osie w "data units"
  xr <- c(0, 1)                               # pooled_clip zawsze [0,1]
  yr <- range(w$delta, na.rm = TRUE)
  if (!is.finite(diff(yr)) || diff(yr) == 0) yr <- c(-0.05, 0.05)  # awaryjnie
  xspan <- diff(xr);  yspan <- diff(yr)

  arrow_x_frac <- max(0.002, min(0.995, arrow_x_frac))
  arrow_x   <- xr[1] + xspan * arrow_x_frac
  arrow_len <- max(1e-6, arrow_length_frac * yspan)

  # etykiety: gwarantowanie nad/pod zerem, z minimalnym buforem eps
  eps <- max(yspan * 0.005, 1e-6)
  label_x <- max(xr[1], arrow_x - label_x_gap_frac * xspan)

  y_up   <- min(yr[2] - label_y_gap_frac * yspan,  0 + 0.6 * arrow_len)
  if (y_up <= 0)  y_up <- min(0 + eps, yr[2] - label_y_gap_frac * yspan)

  y_down <- max(yr[1] + label_y_gap_frac * yspan,  0 - 0.6 * arrow_len)
  if (y_down >= 0) y_down <- max(0 - eps, yr[1] + label_y_gap_frac * yspan)

  # hover
  fmt_pct <- function(x) sprintf("%.1f%%", x)
  fmt_p   <- function(p) ifelse(is.na(p), "NA", formatC(p, format="e", digits=2))
  hover_text <- sprintf(
    paste0("<b>%s</b><br>",
           "%s: %s/%s (%s)<br>",
           "%s: %s/%s (%s)<br>",
           "p (Fisher, wBH): %s<br>",
           "pooled: %s<br>",
           "\u0394: %s"),
    w$feature,
    g1_lab, ifelse(is.na(w$n1),"NA",w$n1), ifelse(is.na(w$N1),"NA",w$N1), fmt_pct(w$pct1),
    g2_lab, ifelse(is.na(w$n2),"NA",w$n2), ifelse(is.na(w$N2),"NA",w$N2), fmt_pct(w$pct2),
    fmt_p(w$p_adj_wbh),
    scales::percent(w$pooled_clip, accuracy = 0.1),
    scales::percent(w$delta,       accuracy = 0.1)
  )

  plt <- plotly::plot_ly() |>
    plotly::add_trace(
      type="scatter", mode="markers",
      x=x_jit, y=y_jit, text=hover_text, hoverinfo="text",
      marker=list(size=point_size, opacity=point_alpha)
    ) |>
    plotly::add_trace(
      type="scatter", mode="lines",
      x = xr, y = c(0,0),
      hoverinfo="skip", line=list(dash="dash"),
      showlegend=FALSE
    )

  if (!is.null(smooth_df)) {
    plt <- plt |>
      plotly::add_trace(
        type="scatter", mode="lines",
        x=smooth_df$pooled_clip, y=smooth_df$delta,
        hoverinfo="skip", showlegend=FALSE
      )
  }

  plt <- plt |>
    plotly::layout(
      title = list(text = htmltools::HTML(
        paste0(plot_title %||% sprintf("Per-peptide shift vs pooled prevalence ( %s \u2212 %s )", g2_lab, g1_lab),
               if (!is.null(plot_subtitle)) sprintf("<br><sup>%s</sup>", plot_subtitle) else "")
      )),
      xaxis = list(title=sprintf("Pooled prevalence (%s & %s)", g1_lab, g2_lab),
                   tickformat=".0%%", range=xr),
      yaxis = list(title=sprintf("\u0394 prevalence (%s \u2212 %s)", g2_lab, g1_lab),
                   tickformat=".1%"),
      annotations = list(
        # strzałka w górę
        list(x=arrow_x, y= arrow_len, xref="x", yref="y",
             ax=arrow_x, ay=0, axref="x", ayref="y",
             text="", showarrow=TRUE, arrowhead=2, arrowsize=0.6,
             arrowwidth=2, arrowcolor=arrow_color),
        # strzałka w dół
        list(x=arrow_x, y=-arrow_len, xref="x", yref="y",
             ax=arrow_x, ay=0, axref="x", ayref="y",
             text="", showarrow=TRUE, arrowhead=2, arrowsize=0.6,
             arrowwidth=2, arrowcolor=arrow_color),

        # etykieta nad zerem
        list(x=label_x, y=y_up, xref="x", yref="y",
             text=paste0("More in ", g2_lab), showarrow=FALSE,
             xanchor="right", font=list(color=arrow_color, size=12),
             bgcolor="rgba(255,255,255,0.65)", bordercolor=arrow_color, borderwidth=1),
        # etykieta pod zerem
        list(x=label_x, y=y_down, xref="x", yref="y",
             text=paste0("More in ", g1_lab), showarrow=FALSE,
             xanchor="right", font=list(color=arrow_color, size=12),
             bgcolor="rgba(255,255,255,0.65)", bordercolor=arrow_color, borderwidth=1)
      ),
      margin = list(l=60, r=80, b=60, t=60)
    )

  plt
}

# -------------------------------------------------------------------------
# forest_delta: top/bottom forest plot for DELTA/Stouffer results (RAW T_obs)
# -------------------------------------------------------------------------

#' Forest plot of top/bottom raw Stouffer T by rank (with arrows)
#'
#' @description
#' Builds a forest plot for the most extreme DELTA/Stouffer results within a
#' chosen rank using the **raw** `T_obs` (no division by any SD). Rows can be
#' optionally filtered by a significance column (e.g., `p_adj_rank`, `padj_wbh`)
#' with a user-specified threshold.
#'
#' @param results_tbl A data frame/tibble with at least:
#'   `rank, feature, group1, group2, design, T_obs, p_perm, p_adj_rank, category_rank_bh`.
#' @param rank_of_interest Character scalar, e.g. `"species"`.
#' @param n_each Integer; how many items from negative and positive tails (default 15).
#' @param use_diverging_colors Logical; if `TRUE`, lines/points are shaded blue (left, T<0) to
#'   red (right, T>0) with higher contrast; otherwise monochrome.
#' @param filter_significant character; name of the column to filter on, or `"none"` to disable filtering.
#'   if the named column is numeric, rows are kept where `col <= sig_level`.
#'   if the named column is non-numeric (e.g., `category_rank_bh`), rows are kept where
#'   `col == "significant (BH, per rank)"`. default: `"none"`.
#' @param sig_level numeric; significance threshold used when `filter_significant` is a numeric column (default 0.05).
#' @param add_signed_z_from_p Logical; if `TRUE`, computes a signed Z from
#'   permutation p (`sign(T_obs) * qnorm(1 - p_perm/2)`) and returns it in data.
#' @param left_label Text for the left arrow/side label.
#' @param right_label Text for the right arrow/side label.
#' @param arrow_length_frac Fraction of max |T_obs| used as half-length of arrows.
#' @param label_x_gap_frac Horizontal gap for arrow labels beyond arrow tips
#'   (fraction of max |T_obs|).
#' @param y_pad Vertical padding above the top category for arrows/labels (in y units).
#' @param label_y_offset Additional vertical offset for arrow-end labels (in y units).
#' @param label_vjust Passed to `annotate("text", vjust=...)` (default -0.30).
#' @param arrow_color Arrow/label color.
#' @param arrow_linewidth Arrow line width.
#' @param arrow_head_length_mm Arrow head length in mm.
#' @param statistic_to_plot Which statistic to plot: \code{"T"},
#'   \code{"T_stand"}, or \code{"Z_from_p"}.
#'
#' @return A list with:
#' \itemize{
#'   \item `data` – tibble used for plotting (selected top/bottom items),
#'   \item `plot` – ggplot object.
#' }
#'
#' @details
#' - `T_obs` here is plotted **as is** (raw Stouffer statistic from your shift).
#' - For calibrated inference or cross-figure comparability, prefer permutation
#'   p-values (optionally expose `add_signed_z_from_p=TRUE` for reporting).
#' @export
forest_delta <- function(
    results_tbl,
    rank_of_interest,
    n_each = NULL,
    n_neg_each = 15,
    n_pos_each = 15,
    use_diverging_colors = FALSE,
    filter_significant = "none",
    sig_level = 0.05,
    add_signed_z_from_p = FALSE,
    # arrow/label params
    left_label    = "More in group1",
    right_label   = "More in group2",
    arrow_length_frac = 0.35,
    label_x_gap_frac  = 0.06,
    y_pad         = 0.6,
    label_y_offset = 0.00,
    label_vjust    = -0.30,
    arrow_color   = "red",
    arrow_linewidth = 0.6,
    arrow_head_length_mm = 3.0,
    # global typography & aesthetics
    base_text_pt  = 12,
    font_family   = "Montserrat",
    seg_width     = 1.2,
    point_size    = 3.6,
    show_grid     = FALSE,
    statistic_to_plot = c("T", "T_stand", "Z_from_p")
) {
  if (!is.data.frame(results_tbl)) stop("`results_tbl` must be a data.frame/tibble.")
  statistic_to_plot <- match.arg(statistic_to_plot)

  need_cols <- c("rank","feature","group1","group2","design",
                 "T_obs","p_perm","p_adj_rank","category_rank_bh")
  miss <- setdiff(need_cols, colnames(results_tbl))
  if (length(miss)) stop(paste("`results_tbl` is missing required columns:", paste(miss, collapse = ", ")))
  if (!is.null(n_each)) { n_neg_each <- n_pos_each <- as.integer(n_each) }

  # helper: pt -> ggplot size units
  .pt_to_gg <- function(pt) pt / 2.845

  df_rk <- dplyr::filter(results_tbl, .data$rank == rank_of_interest)

  if (!identical(filter_significant, "none")) {
    if (!filter_significant %in% colnames(df_rk)) {
      stop(paste0("column '", filter_significant, "' not found"))
    }
    col_vals <- df_rk[[filter_significant]]
    if (is.numeric(col_vals)) {
      df_rk <- dplyr::filter(df_rk, .data[[filter_significant]] <= sig_level)
    } else {
      df_rk <- dplyr::filter(
        df_rk,
        .data[[filter_significant]] == "significant (BH, per rank)"
      )
    }
  }

  if (nrow(df_rk) == 0L) {
    return(list(
      data = df_rk,
      plot = .empty_placeholder_plot(
        sprintf("No rows to plot (rank='%s')", rank_of_interest)
      )
    ))
  }

  # choose statistic used for ranking and selection
  if (identical(statistic_to_plot, "T")) {
    df_rk$stat_for_sort <- df_rk$T_obs
    stat_title       <- "Stouffer T (raw)"
    stat_title_short <- "Stouffer T"
  } else if (identical(statistic_to_plot, "T_stand")) {
    if (!"T_obs_stand" %in% names(df_rk)) {
      stop("Column `T_obs_stand` not found in `x`.")
    }
    df_rk$stat_for_sort <- df_rk$T_obs_stand
    stat_title       <- "Stouffer T (permutation-standardized)"
    stat_title_short <- "T (standardized)"
  } else {  # "Z_from_p"
    if (!"Z_from_p" %in% names(df_rk)) {
      stop("Column `Z_from_p` not found; pass ph_prevalence_shift() results or provide Z_from_p yourself.")
    }
    df_rk$stat_for_sort <- df_rk$Z_from_p
    stat_title       <- "Z from permutation p-values"
    stat_title_short <- "Z_from_p"
  }


  # selection based on chosen statistic
  n_pos <- sum(df_rk$stat_for_sort > 0, na.rm = TRUE)
  n_neg <- sum(df_rk$stat_for_sort < 0, na.rm = TRUE)

  top_pos <- df_rk |>
    dplyr::arrange(dplyr::desc(.data$stat_for_sort)) |>
    dplyr::slice_head(n = min(n_pos_each, n_pos))

  top_neg <- df_rk |>
    dplyr::arrange(.data$stat_for_sort) |>
    dplyr::slice_head(n = min(n_neg_each, n_neg))

  df_plot <- dplyr::bind_rows(top_neg, top_pos) |>
    dplyr::mutate(
      species_label = .data$feature,
      species_label = forcats::fct_reorder(.data$species_label, .data$stat_for_sort)
    )

  # this is what goes on the x-axis
  if (identical(statistic_to_plot, "T")) {
    df_plot$stat_val <- df_plot$T_obs
  } else if (identical(statistic_to_plot, "T_stand")) {
    df_plot$stat_val <- df_plot$T_obs_stand
  } else {  # "Z" or "Z_from_p"
    df_plot$stat_val <- df_plot$Z_from_p
  }


  g1  <- if (nrow(df_plot)) df_plot$group1[1] else ""
  g2  <- if (nrow(df_plot)) df_plot$group2[1] else ""
  des <- if (nrow(df_plot)) df_plot$design[1] else ""

  # color scaling based on plotted statistic
  vals    <- df_plot$stat_val
  max_neg <- max(abs(vals[vals < 0]), na.rm = TRUE)
  max_pos <- max(abs(vals[vals > 0]), na.rm = TRUE)
  if (!is.finite(max_neg) || max_neg == 0) max_neg <- max(abs(vals), 1, na.rm = TRUE)
  if (!is.finite(max_pos) || max_pos == 0) max_pos <- max(abs(vals), 1, na.rm = TRUE)

  df_plot$T_col <- dplyr::case_when(
    vals < 0 ~ -abs(vals) / max_neg,
    vals > 0 ~  vals      / max_pos,
    TRUE     ~  0
  )
  gamma <- 0.85
  df_plot$T_col <- sign(df_plot$T_col) * (abs(df_plot$T_col))^gamma

  # ---- base ggplot with global typography ----
  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = .data$stat_val, y = .data$species_label)
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      linewidth  = 0.4,
      alpha      = 0.6
    ) +
    {
      if (isTRUE(use_diverging_colors)) {
        ggplot2::geom_segment(
          ggplot2::aes(
            x = 0, xend = .data$stat_val,
            y = .data$species_label, yend = .data$species_label,
            color = .data$T_col
          ),
          linewidth   = seg_width,
          alpha       = 0.9,
          show.legend = FALSE
        )
      } else {
        ggplot2::geom_segment(
          ggplot2::aes(
            x = 0, xend = .data$stat_val,
            y = .data$species_label, yend = .data$species_label
          ),
          linewidth = seg_width,
          alpha     = 0.85
        )
      }
    } +
    {
      if (isTRUE(use_diverging_colors)) {
        ggplot2::geom_point(
          ggplot2::aes(color = .data$T_col),
          size        = point_size,
          show.legend = FALSE
        )
      } else {
        ggplot2::geom_point(size = point_size)
      }
    } +
    ggplot2::labs(
      title = sprintf("Top/Bottom %s - rank: %s", stat_title_short, rank_of_interest),
      subtitle = sprintf(
        "Contrast: %s vs %s (%s); shown: %d neg + %d pos",
        g1, g2, des, nrow(top_neg), nrow(top_pos)
      ),
      x = stat_title,
      y = NULL
    ) +
    theme_phip() +
    ggplot2::theme(
      text          = ggplot2::element_text(size = base_text_pt, family = font_family),
      plot.title    = ggplot2::element_text(size = base_text_pt * 1.15, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_text_pt * 0.95),
      axis.title    = ggplot2::element_text(size = base_text_pt),
      axis.text     = ggplot2::element_text(size = base_text_pt * 0.90),
      legend.text   = ggplot2::element_text(size = base_text_pt * 0.90),
      panel.grid.major = if (isTRUE(show_grid)) ggplot2::element_line() else ggplot2::element_blank(),
      panel.grid.minor = if (isTRUE(show_grid)) ggplot2::element_line() else ggplot2::element_blank()
    )

  if (isTRUE(use_diverging_colors)) {
    p <- p + ggplot2::scale_color_gradientn(
      colors = c("#1f4e79", "#6f94c2", "#eeeeee", "#e7a39c", "#8e1b10"),
      values = c(0, 0.25, 0.5, 0.75, 1),
      limits = c(-1, 1)
    )
  }

  # ---- red arrows + labels (sizes tied to base_text_pt) ----
  lvl_n  <- length(levels(df_plot$species_label))
  if (lvl_n == 0L) lvl_n <- length(unique(df_plot$species_label))
  y_top  <- lvl_n + y_pad

  maxabs <- max(abs(df_plot$stat_val), na.rm = TRUE)
  if (!is.finite(maxabs) || maxabs == 0) maxabs <- 1
  x_off  <- arrow_length_frac * maxabs
  x_left_text  <- -x_off - label_x_gap_frac * maxabs
  x_right_text <-  x_off + label_x_gap_frac * maxabs

  # update subtitle with left/right labels
  if (!is.null(p$labels$subtitle)) {
    s <- p$labels$subtitle
    s <- sub(
      "^Contrast:\\s*.*?;",
      paste0("Contrast: ", left_label, " vs ", right_label, ";"),
      s
    )
    p$labels$subtitle <- s
  }

  arrow_text_size <- .pt_to_gg(base_text_pt * 1.05)

  p <- p +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(mult = c(0.05, 0.10))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::annotate(
      "segment",
      x = 0, xend = -x_off, y = y_top, yend = y_top,
      colour   = arrow_color,
      linewidth = arrow_linewidth,
      arrow = ggplot2::arrow(
        length = grid::unit(arrow_head_length_mm, "mm"),
        type   = "closed",
        ends   = "last"
      )
    ) +
    ggplot2::annotate(
      "text",
      x = x_left_text, y = y_top + label_y_offset, label = left_label,
      size     = arrow_text_size,
      colour   = arrow_color,
      fontface = "bold",
      family   = font_family,
      hjust    = 1,
      vjust    = label_vjust
    ) +
    ggplot2::annotate(
      "segment",
      x = 0, xend =  x_off, y = y_top, yend = y_top,
      colour   = arrow_color,
      linewidth = arrow_linewidth,
      arrow = ggplot2::arrow(
        length = grid::unit(arrow_head_length_mm, "mm"),
        type   = "closed",
        ends   = "last"
      )
    ) +
    ggplot2::annotate(
      "text",
      x = x_right_text, y = y_top + label_y_offset, label = right_label,
      size     = arrow_text_size,
      colour   = arrow_color,
      fontface = "bold",
      family   = font_family,
      hjust    = 0,
      vjust    = label_vjust
    )

  list(data = df_plot, plot = p)
}

#' @title Interactive forest plot of top/bottom raw Stouffer T by rank
#' @description Plotly version of \code{forest_delta()} with the same selection
#'   logic and styling parameters.
#' @param results_tbl A data frame/tibble with required columns documented in
#'   \code{forest_delta()}.
#' @param rank_of_interest Character scalar, e.g. \code{"species"}.
#' @param n_each Integer; how many items from negative and positive tails.
#' @param n_neg_each,n_pos_each Integers; per-tail counts when \code{n_each} is
#'   \code{NULL}.
#' @param use_diverging_colors Logical; if \code{TRUE}, apply a blue-to-red gradient.
#' @param filter_significant Column name to filter on, or \code{"none"}.
#' @param sig_level Threshold for numeric significance columns. Default 0.05.
#' @param add_signed_z_from_p Logical; include signed Z in output data.
#' @param left_label,right_label Arrow-end labels.
#' @param arrow_length_frac Fraction of max |T_obs| for arrow half-length.
#' @param label_x_gap_frac Horizontal label gap (fraction of max |T_obs|).
#' @param y_pad Vertical padding above top category (y units).
#' @param label_y_offset Additional vertical label offset.
#' @param label_vjust Included for API parity (ggplot-only).
#' @param arrow_color Arrow/label color.
#' @param arrow_linewidth Arrow line width (plotly units).
#' @param arrow_head_length_mm Arrow head size (approximate, plotly units).
#' @param show_grid Logical; show grid lines.
#' @param base_text_pt Base text size.
#' @param font_family Font family name.
#' @param seg_width Segment line width.
#' @param point_size Point size.
#' @param statistic_to_plot Which statistic to plot: \code{"T"},
#'   \code{"T_stand"}, or \code{"Z_from_p"}.
#' @return A list with \code{data} and \code{plot} (plotly object).
#' @export
forest_delta_plotly <- function(
    results_tbl,
    rank_of_interest,
    n_each = NULL,
    n_neg_each = 15,
    n_pos_each = 15,
    use_diverging_colors = FALSE,
    filter_significant = "none",
    sig_level = 0.05,
    add_signed_z_from_p = FALSE,
    left_label    = "More in group1",
    right_label   = "More in group2",
    arrow_length_frac = 0.35,
    label_x_gap_frac  = 0.06,
    y_pad         = 0.6,
    label_y_offset = 0.00,
    label_vjust    = -0.30,
    arrow_color   = "red",
    arrow_linewidth = 0.6,
    arrow_head_length_mm = 3.0,
    show_grid     = FALSE,
    base_text_pt  = 12,
    font_family   = "Montserrat",
    seg_width     = 1.6,
    point_size    = 11,
    statistic_to_plot = c("T", "T_stand", "Z_from_p")
) {
  if (!is.data.frame(results_tbl)) stop("`results_tbl` must be a data.frame/tibble.")
  statistic_to_plot <- match.arg(statistic_to_plot)

  need_cols <- c("rank","feature","group1","group2","design",
                 "T_obs","p_perm","p_adj_rank","category_rank_bh")
  miss <- setdiff(need_cols, colnames(results_tbl))
  if (length(miss)) stop(paste("`results_tbl` is missing required columns:", paste(miss, collapse = ", ")))
  if (!is.null(n_each)) { n_neg_each <- n_pos_each <- as.integer(n_each) }

  df_rk <- dplyr::filter(results_tbl, .data$rank == rank_of_interest)

  if (!identical(filter_significant, "none")) {
    if (!filter_significant %in% colnames(df_rk)) {
      stop(paste0("column '", filter_significant, "' not found"))
    }
    if (is.numeric(df_rk[[filter_significant]])) {
      df_rk <- dplyr::filter(df_rk, .data[[filter_significant]] <= sig_level)
    } else {
      df_rk <- dplyr::filter(
        df_rk,
        .data[[filter_significant]] == "significant (BH, per rank)"
      )
    }
  }

  if (identical(statistic_to_plot, "T")) {
    stat_title       <- "Stouffer T (raw)"
    stat_title_short <- "Stouffer T"
  } else if (identical(statistic_to_plot, "T_stand")) {
    stat_title       <- "Stouffer T (permutation-standardized)"
    stat_title_short <- "T (standardized)"
  } else {  # "Z" or "Z_from_p"
    stat_title       <- "Z from permutation p-values"
    stat_title_short <- "Z_from_p"
  }

  if (nrow(df_rk) == 0L) {
    plt <- plotly::plot_ly(type = "scatter", mode = "text") |>
      plotly::layout(
        title = list(
          text = paste0("Top/Bottom ", stat_title_short, " - rank: ", rank_of_interest)
        ),
        xaxis = list(visible = FALSE, showgrid = FALSE),
        yaxis = list(visible = FALSE, showgrid = FALSE),
        annotations = list(
          list(
            text = sprintf("No rows to plot (rank='%s')", rank_of_interest),
            xref = "paper", yref = "paper", x = 0.5, y = 0.5,
            showarrow = FALSE,
            font = list(size = base_text_pt, family = font_family)
          )
        ),
        font = list(family = font_family, size = base_text_pt),
        margin = list(l = 20, r = 20, t = 60, b = 20)
      )
    return(list(data = df_rk, plot = plt))
  }

  # choose statistic used for ranking and selection
  if (identical(statistic_to_plot, "T")) {
    df_rk$stat_for_sort <- df_rk$T_obs
  } else if (identical(statistic_to_plot, "T_stand")) {
    if (!"T_obs_stand" %in% names(df_rk)) {
      stop("Column `T_obs_stand` not found in `x`.")
    }
    df_rk$stat_for_sort <- df_rk$T_obs_stand
  } else {  # "Z" or "Z_from_p"
    if (!"Z_from_p" %in% names(df_rk)) {
      stop("Column `Z_from_p` not found; pass ph_prevalence_shift() results or provide Z_from_p yourself.")
    }
    df_rk$stat_for_sort <- df_rk$Z_from_p
  }

  # selection based on chosen statistic
  n_pos <- sum(df_rk$stat_for_sort > 0, na.rm = TRUE)
  n_neg <- sum(df_rk$stat_for_sort < 0, na.rm = TRUE)

  top_pos <- df_rk |>
    dplyr::arrange(dplyr::desc(.data$stat_for_sort)) |>
    dplyr::slice_head(n = min(n_pos_each, n_pos))

  top_neg <- df_rk |>
    dplyr::arrange(.data$stat_for_sort) |>
    dplyr::slice_head(n = min(n_neg_each, n_neg))

  df_plot <- dplyr::bind_rows(top_neg, top_pos) |>
    dplyr::mutate(
      species_label = .data$feature,
      species_label = forcats::fct_reorder(.data$species_label, .data$stat_for_sort)
    )

  # statistic on x-axis
  if (identical(statistic_to_plot, "T")) {
    df_plot$stat_val <- df_plot$T_obs
  } else if (identical(statistic_to_plot, "T_stand")) {
    df_plot$stat_val <- df_plot$T_obs_stand
  } else {  # "Z" or "Z_from_p"
    df_plot$stat_val <- df_plot$Z_from_p
  }

  g1  <- if (nrow(df_plot)) df_plot$group1[1] else ""
  g2  <- if (nrow(df_plot)) df_plot$group2[1] else ""
  des <- if (nrow(df_plot)) df_plot$design[1] else ""

  # color scaling based on plotted statistic
  vals    <- df_plot$stat_val
  max_neg <- max(abs(vals[vals < 0]), na.rm = TRUE)
  max_pos <- max(abs(vals[vals > 0]), na.rm = TRUE)
  if (!is.finite(max_neg) || max_neg == 0) max_neg <- max(abs(vals), 1, na.rm = TRUE)
  if (!is.finite(max_pos) || max_pos == 0) max_pos <- max(abs(vals), 1, na.rm = TRUE)

  df_plot$T_col <- dplyr::case_when(
    vals < 0 ~ -abs(vals) / max_neg,
    vals > 0 ~  vals      / max_pos,
    TRUE     ~  0
  )
  gamma <- 0.85
  df_plot$T_col <- sign(df_plot$T_col) * (abs(df_plot$T_col))^gamma

  if (!"n_peptides_used" %in% names(df_plot)) df_plot$n_peptides_used <- NA_integer_
  df_plot$padj_hover <- if ("padj_wbh" %in% names(df_plot)) df_plot$padj_wbh else {
    if ("p_adj_rank" %in% names(df_plot)) df_plot$p_adj_rank else NA_real_
  }
  df_plot$interpretation <- ifelse(
    df_plot$stat_val < 0, paste("More in", df_plot$group1),
    ifelse(df_plot$stat_val > 0, paste("More in", df_plot$group2), "No difference")
  )

  hover_text <- sprintf(
    "Feature: %s<br>n_peptides_used: %s<br>%s: %s<br>Interpretation: %s<br>padj_wbh: %s",
    df_plot$feature,
    ifelse(is.na(df_plot$n_peptides_used), "NA", format(df_plot$n_peptides_used, big.mark = ",")),
    stat_title_short,
    format(round(df_plot$stat_val, 4), nsmall = 4),
    df_plot$interpretation,
    ifelse(is.na(df_plot$padj_hover), "NA", formatC(df_plot$padj_hover, format = "e", digits = 2))
  )

  map_hex <- local({
    cols <- c("#1f4e79", "#6f94c2", "#eeeeee", "#e7a39c", "#8e1b10")
    ramp <- grDevices::colorRampPalette(cols)(256)
    function(v) {
      v <- pmax(-1, pmin(1, v))
      idx <- round(((v + 1) / 2) * 255) + 1
      ramp[idx]
    }
  })
  seg_hex <- if (isTRUE(use_diverging_colors)) map_hex(df_plot$T_col) else rep("rgba(0,0,0,0.90)", nrow(df_plot))
  pt_hex  <- if (isTRUE(use_diverging_colors)) map_hex(df_plot$T_col) else rep("rgba(0,0,0,1)",     nrow(df_plot))

  y_levels <- levels(df_plot$species_label)
  df_plot$species_chr <- as.character(df_plot$species_label)

  plt <- plotly::plot_ly()

  plt <- plotly::layout(
    plt,
    shapes = list(
      list(
        type = "line", x0 = 0, x1 = 0, xref = "x",
        y0 = 0, y1 = 1, yref = "paper",
        line = list(dash = "dash", width = 1, color = "rgba(0,0,0,0.5)")
      )
    )
  )

  if (nrow(df_plot)) {
    for (i in seq_len(nrow(df_plot))) {
      plt <- plotly::add_segments(
        plt,
        x = 0, xend = df_plot$stat_val[i],
        y = df_plot$species_chr[i], yend = df_plot$species_chr[i],
        showlegend = FALSE,
        hoverinfo  = "text",
        text       = hover_text[i],
        line       = list(color = seg_hex[i], width = seg_width),
        opacity    = if (isTRUE(use_diverging_colors)) 0.95 else 0.90
      )
    }
  }

  plt <- plotly::add_markers(
    plt,
    x = df_plot$stat_val,
    y = df_plot$species_chr,
    showlegend = FALSE,
    hoverinfo  = "text",
    text       = hover_text,
    marker     = list(size = point_size, color = pt_hex, opacity = 1)
  )

  maxabs <- max(abs(df_plot$stat_val), na.rm = TRUE)
  if (!is.finite(maxabs) || maxabs == 0) maxabs <- 1
  x_off  <- arrow_length_frac * maxabs
  x_left_text  <- -x_off - label_x_gap_frac * maxabs
  x_right_text <-  x_off + label_x_gap_frac * maxabs

  ann_list <- list(
    list(
      x = -x_off, y = 1.02, ax = 0, ay = 1.02, xref = "x", yref = "paper",
      axref = "x", ayref = "paper", showarrow = TRUE, arrowhead = 3,
      arrowsize = arrow_head_length_mm / 1.875,
      arrowwidth = arrow_linewidth * 4.5, arrowcolor = arrow_color
    ),
    list(
      x =  x_off, y = 1.02, ax = 0, ay = 1.02, xref = "x", yref = "paper",
      axref = "x", ayref = "paper", showarrow = TRUE, arrowhead = 3,
      arrowsize = arrow_head_length_mm / 1.875,
      arrowwidth = arrow_linewidth * 4.5, arrowcolor = arrow_color
    ),
    list(
      x = x_left_text, y = 1.04 + label_y_offset, xref = "x", yref = "paper",
      text = left_label, showarrow = FALSE,
      font = list(size = round(base_text_pt * 1.05), color = arrow_color, family = font_family),
      xanchor = "right", yanchor = "bottom"
    ),
    list(
      x = x_right_text, y = 1.04 + label_y_offset, xref = "x", yref = "paper",
      text = right_label, showarrow = FALSE,
      font = list(size = round(base_text_pt * 1.05), color = arrow_color, family = font_family),
      xanchor = "left", yanchor = "bottom"
    )
  )

  subtitle_txt <- sprintf(
    "Contrast: %s vs %s (%s); shown: %d neg + %d pos",
    g1, g2, des,
    sum(df_plot$stat_val < 0, na.rm = TRUE),
    sum(df_plot$stat_val > 0, na.rm = TRUE)
  )
  subtitle_txt <- sub(
    "^Contrast:\\s*.*?;",
    paste0("Contrast: ", left_label, " vs ", right_label, ";"),
    subtitle_txt
  )

  plt <- plotly::layout(
    plt,
    title = list(
      text = paste0(
        "Top/Bottom ", stat_title_short, " - rank: ", rank_of_interest,
        "<br><sub>", subtitle_txt, "</sub>"
      ),
      font = list(size = round(base_text_pt * 1.15), family = font_family)
    ),
    xaxis = list(
      title    = stat_title,
      zeroline = FALSE,
      showgrid = isTRUE(show_grid),
      titlefont = list(size = base_text_pt, family = font_family),
      tickfont  = list(size = round(base_text_pt * 0.9), family = font_family)
    ),
    yaxis = list(
      title = NULL,
      categoryorder = "array",
      categoryarray = y_levels,
      showgrid = isTRUE(show_grid),
      tickfont = list(size = round(base_text_pt * 0.9), family = font_family)
    ),
    font = list(family = font_family, size = base_text_pt),
    margin = list(
      l = 10 + max(nchar(df_plot$species_chr)) * 5,
      r = 10, t = 90, b = 40
    ),
    annotations = ann_list
  )

  list(data = df_plot, plot = plt)
}

# --- small internal helper: white placeholder with red label ------------------
.empty_placeholder_plot <- function(label = "No significant results") {
  ggplot2::ggplot() +
    ggplot2::geom_text(
      data = data.frame(x = 0.5, y = 0.5),
      ggplot2::aes(x = x, y = y),
      label = label, color = "red", fontface = "bold", size = 6
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = NA))
}

# ==============================================================================
# ECDF of per-peptide prevalence for two groups (static)
# Input: tibble from ph_compute_prevalence() with group1, group2, prop1, prop2
# Optional: feature, n1/N1/n2/N2, percent1/percent2
# ==============================================================================
#' @title ECDF of per-peptide prevalence for two groups (static)
#' @description Plot ECDF curves for \code{prop1} and \code{prop2} from a
#'   prevalence table.
#' @param prev_tbl Data frame with columns \code{group1}, \code{group2},
#'   \code{prop1}, \code{prop2}. Optional \code{feature} columns are ignored.
#' @param group_pair_values Optional length-2 vector \code{c(g1, g2)} selecting a pair.
#' @param group_labels Optional length-2 vector of display labels.
#' @param line_width Line width for ECDF steps.
#' @param line_alpha Line alpha for ECDF steps.
#' @param group1_color,group2_color Line colors for group1 and group2.
#' @param show_median Logical; add median lines. Default \code{TRUE}.
#' @param show_ks Logical; add KS test summary to subtitle. Default \code{TRUE}.
#' @param plot_title,plot_subtitle Optional plot labels.
#' @param x_label,y_label Optional axis labels.
#' @return A `ggplot` object.
#' @export
ecdf_prevalence <- function(
    prev_tbl,
    group_pair_values = NULL,
    group_labels     = NULL,
    # styling
    line_width   = 1.0,
    line_alpha   = 1.0,
    group1_color = "#1f77b4",
    group2_color = "#d62728",
    show_median  = TRUE,
    show_ks      = TRUE,
    plot_title   = NULL,
    plot_subtitle = NULL,
    x_label      = NULL,
    y_label      = NULL
) {
  `%||%` <- function(a,b) if (!is.null(a)) a else b

  d <- prev_tbl
  need <- c("group1","group2","prop1","prop2")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("ecdf_prevalence(): missing columns: ", paste(miss, collapse = ", "))

  # ---- select exactly one pair ----
  if (!is.null(group_pair_values)) {
    if (length(group_pair_values) != 2L) {
      stop("group_pair_values must be length-2 vector: c(g1, g2).")
    }
    d <- d[d$group1 == group_pair_values[1] &
             d$group2 == group_pair_values[2], , drop = FALSE]
    if (!nrow(d)) stop("No rows for group_pair_values = c('",
                       group_pair_values[1], "', '", group_pair_values[2], "').")
    g1_raw <- group_pair_values[1]; g2_raw <- group_pair_values[2]
  } else {
    pairs <- unique(d[, c("group1","group2")])
    if (nrow(pairs) != 1L) {
      stop("Multiple (group1, group2) pairs detected. Pass group_pair_values = c(g1, g2).")
    }
    g1_raw <- pairs$group1[1]; g2_raw <- pairs$group2[1]
  }

  if (is.null(group_labels)) {
    g1_lab <- as.character(g1_raw); g2_lab <- as.character(g2_raw)
  } else {
    if (length(group_labels) != 2L) stop("group_labels must be length-2 vector.")
    g1_lab <- group_labels[1]; g2_lab <- group_labels[2]
  }

  # ---- data vectors ----
  x1 <- as.numeric(d$prop1)
  x2 <- as.numeric(d$prop2)
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (!length(x1) || !length(x2)) stop("No finite prop values to plot.")

  # helper to convert ecdf() to data.frame (step function)
  ecdf_df <- function(x) {
    if (!length(x)) return(data.frame(x = numeric(0), y = numeric(0)))
    xs <- sort(unique(x))
    F  <- stats::ecdf(x)
    data.frame(x = xs, y = F(xs))
  }

  df1 <- ecdf_df(x1); df1$group <- g1_lab
  df2 <- ecdf_df(x2); df2$group <- g2_lab
  ec  <- rbind(df1, df2)

  med1 <- stats::median(x1, na.rm = TRUE)
  med2 <- stats::median(x2, na.rm = TRUE)
  dmed <- med2 - med1

  ks_txt <- NULL
  if (isTRUE(show_ks)) {
    ks <- try(suppressWarnings(stats::ks.test(x1, x2)), silent = TRUE)
    if (!inherits(ks, "try-error")) {
      ks_txt <- sprintf("KS D=%.3f  p=%s", unname(ks$statistic), formatC(ks$p.value, format = "e", digits = 2))
    }
  }

  p <- ggplot2::ggplot(ec, ggplot2::aes(x, y, color = group)) +
    ggplot2::geom_step(linewidth = line_width, alpha = line_alpha, direction = "hv") +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    ggplot2::labs(
      title    = plot_title %||% sprintf("ECDF of per-peptide prevalence (%s vs %s)", g2_lab, g1_lab),
      subtitle = plot_subtitle %||% if (!is.null(ks_txt)) sprintf("%s | Δ median = %s", ks_txt, scales::percent(dmed, accuracy = 0.1)) else NULL,
      x        = x_label %||% "Prevalence",
      y        = y_label %||% "ECDF"
    ) +
    ggplot2::scale_color_manual(values = setNames(c(group1_color, group2_color), c(g1_lab, g2_lab))) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  if (isTRUE(show_median)) {
    p <- p +
      ggplot2::geom_vline(xintercept = med1, color = group1_color, linetype = 3) +
      ggplot2::geom_vline(xintercept = med2, color = group2_color, linetype = 3)
  }

  p
}

# ==============================================================================
# ECDF of per-peptide prevalence for two groups (interactive, plotly)
# Same input/semantics as ecdf_prevalence()
# ==============================================================================
#' @title ECDF of per-peptide prevalence for two groups (interactive)
#' @description Plotly version of \code{ecdf_prevalence()}.
#' @param prev_tbl Data frame with columns \code{group1}, \code{group2},
#'   \code{prop1}, \code{prop2}.
#' @param group_pair_values Optional length-2 vector \code{c(g1, g2)} selecting a pair.
#' @param group_labels Optional length-2 vector of display labels.
#' @param line_width Line width for ECDF steps (plotly units).
#' @param line_alpha Line alpha for ECDF steps.
#' @param group1_color,group2_color Line colors for group1 and group2.
#' @param show_median Logical; add median lines. Default \code{TRUE}.
#' @param show_ks Logical; add KS test summary to subtitle. Default \code{TRUE}.
#' @param plot_title,plot_subtitle Optional plot labels.
#' @return A plotly object.
#' @export
ecdf_prevalence_interactive <- function(
    prev_tbl,
    group_pair_values = NULL,
    group_labels     = NULL,
    # styling
    line_width   = 2.0,      # px
    line_alpha   = 1.0,
    group1_color = "#1f77b4",
    group2_color = "#d62728",
    show_median  = TRUE,
    show_ks      = TRUE,
    plot_title   = NULL,
    plot_subtitle = NULL
) {
  `%||%` <- function(a,b) if (!is.null(a)) a else b

  d <- prev_tbl
  need <- c("group1","group2","prop1","prop2")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("ecdf_prevalence_interactive(): missing columns: ", paste(miss, collapse = ", "))

  # ---- pair ----
  if (!is.null(group_pair_values)) {
    if (length(group_pair_values) != 2L) {
      stop("group_pair_values must be length-2 vector: c(g1, g2).")
    }
    d <- d[d$group1 == group_pair_values[1] &
             d$group2 == group_pair_values[2], , drop = FALSE]
    if (!nrow(d)) stop("No rows for group_pair_values = c('",
                       group_pair_values[1], "', '", group_pair_values[2], "').")
    g1_raw <- group_pair_values[1]; g2_raw <- group_pair_values[2]
  } else {
    pairs <- unique(d[, c("group1","group2")])
    if (nrow(pairs) != 1L) stop("Multiple (group1, group2) pairs detected. Pass group_pair_values = c(g1, g2).")
    g1_raw <- pairs$group1[1]; g2_raw <- pairs$group2[1]
  }

  if (is.null(group_labels)) {
    g1_lab <- as.character(g1_raw); g2_lab <- as.character(g2_raw)
  } else {
    if (length(group_labels) != 2L) stop("group_labels must be length-2 vector.")
    g1_lab <- group_labels[1]; g2_lab <- group_labels[2]
  }

  x1 <- as.numeric(d$prop1); x1 <- x1[is.finite(x1)]
  x2 <- as.numeric(d$prop2); x2 <- x2[is.finite(x2)]
  if (!length(x1) || !length(x2)) stop("No finite prop values to plot.")

  ecdf_df <- function(x) {
    xs <- sort(unique(x))
    F  <- stats::ecdf(x)
    data.frame(x = xs, y = F(xs))
  }

  df1 <- ecdf_df(x1)
  df2 <- ecdf_df(x2)
  med1 <- stats::median(x1, na.rm = TRUE)
  med2 <- stats::median(x2, na.rm = TRUE)
  dmed <- med2 - med1

  ks_txt <- NULL
  if (isTRUE(show_ks)) {
    ks <- try(suppressWarnings(stats::ks.test(x1, x2)), silent = TRUE)
    if (!inherits(ks, "try-error")) {
      ks_txt <- sprintf("KS D=%.3f  p=%s", unname(ks$statistic), formatC(ks$p.value, format = "e", digits = 2))
    }
  }

  # hover formatting
  fmt_pct <- function(x) scales::percent(x, accuracy = 0.1)
  hover1 <- sprintf("<b>%s</b><br>x: %s<br>F(x): %s", g1_lab, fmt_pct(df1$x), fmt_pct(df1$y))
  hover2 <- sprintf("<b>%s</b><br>x: %s<br>F(x): %s", g2_lab, fmt_pct(df2$x), fmt_pct(df2$y))

  plt <- plotly::plot_ly() |>
    plotly::add_trace(
      type = "scatter", mode = "lines",
      x = df1$x, y = df1$y, text = hover1, hoverinfo = "text",
      line = list(width = line_width, color = group1_color, shape = "hv", opacity = line_alpha),
      name = g1_lab
    ) |>
    plotly::add_trace(
      type = "scatter", mode = "lines",
      x = df2$x, y = df2$y, text = hover2, hoverinfo = "text",
      line = list(width = line_width, color = group2_color, shape = "hv", opacity = line_alpha),
      name = g2_lab
    ) |>
    plotly::layout(
      title = list(text = htmltools::HTML(
        paste0(plot_title %||% sprintf("ECDF of per-peptide prevalence ( %s vs %s )", g2_lab, g1_lab),
               if (!is.null(plot_subtitle) || !is.null(ks_txt)) {
                 sub <- plot_subtitle %||% ""
                 km  <- if (!is.null(ks_txt)) sprintf("%s | Δ median = %s", ks_txt, fmt_pct(dmed)) else ""
                 sprintf("<br><sup>%s%s%s</sup>",
                         sub, if (nzchar(sub) && nzchar(km)) " — " else "", km)
               } else "")
      )),
      xaxis = list(title = "Prevalence", tickformat = ".0%", range = c(0,1)),
      yaxis = list(title = "ECDF", tickformat = ".0%", range = c(0,1)),
      legend = list(orientation = "h", x = 0, y = 1.1),
      margin = list(l = 60, r = 40, b = 60, t = 70)
    )

  if (isTRUE(show_median)) {
    plt <- plt |>
      plotly::add_segments(x = med1, xend = med1, y = 0, yend = 1,
                           line = list(dash = "dot", width = 1, color = group1_color),
                           hoverinfo = "skip", showlegend = FALSE) |>
      plotly::add_segments(x = med2, xend = med2, y = 0, yend = 1,
                           line = list(dash = "dot", width = 1, color = group2_color),
                           hoverinfo = "skip", showlegend = FALSE)
  }

  plt
}
