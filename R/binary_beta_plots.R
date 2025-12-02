# -------------------------------------------------------------------
# MAIN: plot_beta_pcoa (prototype, shaded group×time)
# -------------------------------------------------------------------
#' @param x phip_beta_diversity
#' @param view which view (defaults to first)
#' @param subview which subview (defaults to first)
#' @param axes integer pair, e.g. c(1,2) or c(2,3)
#' @param group_col,time_fac,time_cont column names present in scores tibble
#' @param shade "auto" | "time" | "group" | "none"  (kept for compatibility; points use shaded legend automatically when both factors present)
#' @param show_centroids,show_ellipses logical
#' @param centroid_by "auto" | "group" | "time" | "group_time"
#' @param connect_centroids "none" | "group" | "time"
#' @param ellipse_type "t" | "norm" | "euclid"
#' @param ellipse_level numeric confidence for ellipses (default 0.95)
#' @param point_size,point_alpha numeric
#' @param centroid_size numeric (uses shape 21 internally)
#' @param base_colors optional named vector of base colors per group (e.g., c(g1="#...", g2="#..."))
#' @param shade_range numeric length-2 in [0,1], how much to tint toward white for time shades
#' @return a ggplot object
# -------------------------------------------------------------------
# plot_beta_pcoa — shaded group×time points + centroids + flexible ellipses
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# plot_beta_pcoa — shaded group×time points + centroids + flexible ellipses
# (+ optional CAP biplot overlay for categorical time)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# plot_beta_pcoa — shaded group×time points + centroids + flexible ellipses
# + optional CAP biplot overlay
# now supports continuous-time shading (tints within each group)
# -------------------------------------------------------------------
plot_beta_pcoa <- function(x,
                           view = NULL,
                           subview = NULL,
                           axes = c(1, 2),
                           group_col = "group",
                           time_fac = "time_fac",
                           time_cont = "time_cont",
                           # shading ----------------------------------------------------------------
                           shade = c("auto", "time", "group", "none"),
                           shade_norm = c("global", "by_group"), # for continuous time
                           shade_direction = c("light_to_dark", "dark_to_light"),
                           shade_range = c(0.15, 0.85), # min→max tint to white
                           base_colors = NULL, # named by group (optional)
                           # centroids / ellipses ---------------------------------------------------
                           show_centroids = TRUE,
                           centroid_by = c("auto", "group", "time", "group_time"),
                           connect_centroids = c("none", "group", "time"),
                           show_ellipses = TRUE,
                           ellipse_by = c("group"),
                           ellipse_type = c("t", "norm", "euclid"),
                           ellipse_level = 0.95,
                           # geoms -------------------------------------------------------------------
                           point_size = 1.5,
                           point_alpha = 0.55,
                           centroid_size = 3,
                           # CAP biplot --------------------------------------------------------------
                           add_biplot = FALSE,
                           biplot_label_size = 3) {
  # ---- tiny helpers ----------------------------------------------------------
  `%||%` <- function(a, b) if (is.null(a)) b else a
  .has <- function(nm, df) !is.null(nm) && nm %in% names(df)
  .ph_info <- function(msg) {
    if (exists(".ph_log_info", mode = "function")) {
      .ph_log_info(msg, step = "plot_beta_pcoa")
    } else if (exists(".phi", mode = "function")) .phi(msg) else message(sprintf("[INFO] %s", msg))
  }
  .ph_err <- function(msg) {
    if (exists(".ph_abort", mode = "function")) {
      .ph_abort(msg)
    } else if (exists(".phe", mode = "function")) .phe(msg) else stop(msg, call. = FALSE)
  }

  .blend_hex <- function(c1, c2 = "#FFFFFF", w = 0.5) {
    r1 <- grDevices::col2rgb(c1) / 255
    r2 <- grDevices::col2rgb(c2) / 255
    out <- (1 - w) * r1 + w * r2
    grDevices::rgb(out[1], out[2], out[3])
  }
  .make_shades <- function(base_hex, n, range = c(0.15, 0.85)) {
    if (n <= 1) {
      return(base_hex)
    }
    w <- seq(range[1], range[2], length.out = n)
    vapply(w, function(wi) .blend_hex(base_hex, "#FFFFFF", wi), character(1))
  }
  .phip_pal <- function(n) {
    if (exists("phip_palette", inherits = TRUE)) {
      pal <- get("phip_palette", inherits = TRUE)
      if (length(pal) < n) rep(pal, length.out = n) else pal[seq_len(n)]
    } else {
      c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928") |> rep(length.out = n)
    }
  }
  .build_shaded_map <- function(groups, times, base_colors = NULL,
                                default_palette = NULL, shade_range = c(0.15, 0.85),
                                sep = " | ") {
    groups <- as.character(groups)
    times <- as.character(times)
    if (is.null(default_palette)) default_palette <- .phip_pal(length(groups))
    # bases (respect user-provided names)
    if (is.null(base_colors)) {
      base_colors <- stats::setNames(default_palette[seq_along(groups)], groups)
    } else {
      if (is.null(names(base_colors))) names(base_colors) <- groups
      base_colors <- base_colors[groups]
    }
    fill_vals <- outline_vals <- character(0)
    for (g in groups) {
      shades <- .make_shades(base_colors[[g]], length(times), range = shade_range)
      names(shades) <- paste(g, times, sep = sep)
      fill_vals <- c(fill_vals, shades)
      outline_vals <- c(outline_vals, stats::setNames(rep(base_colors[[g]], length(times)), names(shades)))
    }
    list(fill = fill_vals, outline = outline_vals, base = base_colors)
  }

  # ---- parse args / meta -----------------------------------------------------
  shade <- match.arg(shade)
  shade_norm <- match.arg(shade_norm)
  shade_direction <- match.arg(shade_direction)
  centroid_by <- match.arg(centroid_by)
  connect_centroids <- match.arg(connect_centroids)
  ellipse_type <- match.arg(ellipse_type)
  ellipse_by <- match.arg(ellipse_by, c("none", "group", "time", "group_time"), several.ok = TRUE)
  if ("none" %in% ellipse_by) ellipse_by <- character(0)

  method <- attr(x, "method_pcoa") %||% "joint"
  if (is.null(view)) view <- names(x)[1]
  if (is.null(subview)) subview <- names(x[[view]])[1]

  .ph_info(sprintf(
    "method: %s | view: %s | subview: %s | axes: %s:%s",
    method, view, subview, axes[1], axes[2]
  ))

  sv <- x[[view]][[subview]]
  df <- sv$pcoa
  if (is.null(df) || !nrow(df)) .ph_err("PCoA/CAP scores not found.")

  has_group <- .has(group_col, df)
  has_timeF <- .has(time_fac, df) && any(!is.na(df[[time_fac]]))
  has_timeC <- .has(time_cont, df) && any(!is.na(df[[time_cont]]))

  # axis names + labels w/ % explained
  axis_prefix <- if (tolower(method) == "cap") "CAP" else "PCoA"
  ax_names <- paste0(axis_prefix, axes)
  if (!all(ax_names %in% names(df))) {
    .ph_err(sprintf("Requested axes not found: %s", paste(ax_names[!ax_names %in% names(df)], collapse = ", ")))
  }
  ax1 <- ax_names[1]
  ax2 <- ax_names[2]
  ve <- sv$var_explained
  if (!is.null(ve) && nrow(ve) == 1) {
    p1 <- suppressWarnings(as.numeric(ve[[paste0("%", axis_prefix, axes[1])]]))
    p2 <- suppressWarnings(as.numeric(ve[[paste0("%", axis_prefix, axes[2])]]))
  } else {
    p1 <- p2 <- NA_real_
  }
  xlab <- if (is.finite(p1)) sprintf("%s%d (%.1f%%)", axis_prefix, axes[1], p1) else ax1
  ylab <- if (is.finite(p2)) sprintf("%s%d (%.1f%%)", axis_prefix, axes[2], p2) else ax2

  # defaults if missing
  if (!has_group) df[[group_col]] <- "All samples"
  groups <- droplevels(factor(df[[group_col]]))
  group_levels <- levels(groups)

  if (has_timeF) {
    df[[time_fac]] <- droplevels(factor(df[[time_fac]]))
    time_levels <- levels(df[[time_fac]])
  } else {
    time_levels <- character(0)
  }

  .ph_info(sprintf(
    "n=%d points | groups: %d | time_cont: %s | time_fac: %s",
    nrow(df), nlevels(groups), has_timeC, has_timeF
  ))

  # ---- palettes / shading ----------------------------------------------------
  pal_default <- .phip_pal(length(group_levels))
  if (is.null(base_colors)) {
    base_map <- stats::setNames(pal_default[seq_along(group_levels)], group_levels)
  } else {
    if (is.null(names(base_colors))) names(base_colors) <- group_levels
    base_map <- base_colors[group_levels]
  }

  # choose shading mode automatically if requested
  if (shade == "auto") {
    shade <- if (has_timeF) "time" else if (has_group) "group" else "none"
  }

  # For categorical time: build "group × time" palette/legend (fill)
  shade_maps <- NULL
  if (shade == "time" && has_timeF) {
    df$gt <- interaction(df[[group_col]], df[[time_fac]],
      sep = " | ", drop = TRUE, lex.order = TRUE
    )
    shade_maps <- .build_shaded_map(group_levels, time_levels,
      base_colors = base_map,
      default_palette = pal_default,
      shade_range = shade_range
    )
  }

  # For continuous time: compute a per-point hex fill (identity scale)
  if (shade == "time" && !has_timeF && has_timeC) {
    # normalize months across ALL points or per group
    if (shade_norm == "by_group") {
      df$..t_norm <- ave(df[[time_cont]], df[[group_col]], FUN = function(z) {
        z <- as.numeric(z)
        rng <- range(z, na.rm = TRUE)
        if (!is.finite(rng[1]) || rng[1] == rng[2]) {
          return(rep(0.5, length(z)))
        }
        (z - rng[1]) / (rng[2] - rng[1])
      })
    } else {
      z <- as.numeric(df[[time_cont]])
      rng <- range(z, na.rm = TRUE)
      df$..t_norm <- if (!is.finite(rng[1]) || rng[1] == rng[2]) rep(0.5, length(z)) else (z - rng[1]) / (rng[2] - rng[1])
    }
    if (shade_direction == "dark_to_light") df$..t_norm <- 1 - df$..t_norm
    # map to tint weights inside shade_range
    df$..w <- shade_range[1] + df$..t_norm * (shade_range[2] - shade_range[1])
    # hex per point (tint base group color towards white)
    df$..fill_hex <- mapply(function(g, w) .blend_hex(base_map[[as.character(g)]], "#FFFFFF", w),
      df[[group_col]], df$..w,
      USE.NAMES = FALSE
    )
  }

  # ---- centroids (granularity) ----------------------------------------------
  if (centroid_by == "auto") {
    centroid_by <- if (has_group && has_timeF) "group_time" else if (has_group) "group" else if (has_timeF) "time" else "group"
  }

  cent <- NULL
  if (isTRUE(show_centroids)) {
    if (centroid_by == "group") {
      cent <- df |>
        dplyr::group_by(.data[[group_col]]) |>
        dplyr::summarise(
          cx = mean(.data[[ax1]], na.rm = TRUE),
          cy = mean(.data[[ax2]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(group_lab = !!group_col)
    } else if (centroid_by == "time") {
      if (!has_timeF) .ph_err("centroid_by='time' requires a categorical time column.")
      cent <- df |>
        dplyr::group_by(.data[[time_fac]]) |>
        dplyr::summarise(
          cx = mean(.data[[ax1]], na.rm = TRUE),
          cy = mean(.data[[ax2]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(time_lab = !!time_fac)
    } else if (centroid_by == "group_time") {
      if (!has_timeF) .ph_err("centroid_by='group_time' requires a categorical time column.")
      cent <- df |>
        dplyr::group_by(.data[[group_col]], .data[[time_fac]]) |>
        dplyr::summarise(
          cx = mean(.data[[ax1]], na.rm = TRUE),
          cy = mean(.data[[ax2]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(group_lab = !!group_col, time_lab = !!time_fac)
      cent$time_ord <- as.integer(factor(cent$time_lab, levels = time_levels))
      cent$gt <- interaction(cent$group_lab, cent$time_lab, sep = " | ", drop = TRUE)
    }
  }

  # ---- base plot -------------------------------------------------------------
  p <- ggplot2::ggplot()

  if (shade == "time" && has_timeF) {
    # categorical-shaded group×time legend (single legend)
    p <- p +
      ggplot2::geom_point(
        data = df,
        ggplot2::aes(
          x = .data[[ax1]], y = .data[[ax2]],
          fill = .data$gt, color = .data[[group_col]]
        ),
        shape = 21, size = point_size, alpha = point_alpha, stroke = 0.6
      ) +
      ggplot2::scale_fill_manual(values = shade_maps$fill, name = "group × time") +
      ggplot2::scale_color_manual(values = base_map, guide = "none")
  } else if (shade == "time" && !has_timeF && has_timeC) {
    # continuous-time tints within each group (identity fill, no fill legend)
    p <- p +
      ggplot2::geom_point(
        data = df,
        ggplot2::aes(
          x = .data[[ax1]], y = .data[[ax2]],
          fill = .data$..fill_hex, color = .data[[group_col]]
        ),
        shape = 21, size = point_size, alpha = point_alpha, stroke = 0.6
      ) +
      ggplot2::scale_fill_identity(guide = "none") +
      ggplot2::scale_color_manual(values = base_map, name = group_col)
  } else if (has_group) { # shade == "group" or time shading disabled
    p <- p +
      ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[ax1]], y = .data[[ax2]], color = .data[[group_col]]),
        size = point_size, alpha = point_alpha
      ) +
      ggplot2::scale_color_manual(values = base_map, name = group_col)
  } else if (has_timeF) {
    # only time present, categorical
    tcols <- stats::setNames(.phip_pal(length(time_levels))[seq_along(time_levels)], time_levels)
    p <- p +
      ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[ax1]], y = .data[[ax2]], color = .data[[time_fac]]),
        size = point_size, alpha = point_alpha
      ) +
      ggplot2::scale_color_manual(values = tcols, name = time_fac)
  } else {
    p <- p + ggplot2::geom_point(
      data = df, ggplot2::aes(x = .data[[ax1]], y = .data[[ax2]]),
      size = point_size, alpha = point_alpha, color = "grey40"
    )
  }

  # ---- ellipses --------------------------------------------------------------
  if (isTRUE(show_ellipses) && length(ellipse_by)) {
    if ("group" %in% ellipse_by && has_group) {
      for (g in group_levels) {
        p <- p + ggplot2::stat_ellipse(
          data = df[df[[group_col]] == g, , drop = FALSE],
          ggplot2::aes(x = .data[[ax1]], y = .data[[ax2]]),
          type = ellipse_type, level = ellipse_level,
          linewidth = 0.7, alpha = 0.7, color = base_map[[g]],
          show.legend = FALSE
        )
      }
    }
    if ("time" %in% ellipse_by && has_timeF) {
      tcols <- stats::setNames(.phip_pal(length(time_levels))[seq_along(time_levels)], time_levels)
      for (tlev in time_levels) {
        p <- p + ggplot2::stat_ellipse(
          data = df[df[[time_fac]] == tlev, , drop = FALSE],
          ggplot2::aes(x = .data[[ax1]], y = .data[[ax2]]),
          type = ellipse_type, level = ellipse_level,
          linewidth = 0.6, alpha = 0.6, color = tcols[[tlev]],
          show.legend = FALSE
        )
      }
    }
    if ("group_time" %in% ellipse_by && has_group && has_timeF) {
      for (g in group_levels) {
        for (tlev in time_levels) {
          dsub <- df[df[[group_col]] == g & df[[time_fac]] == tlev, , drop = FALSE]
          if (nrow(dsub) < 3) next
          p <- p + ggplot2::stat_ellipse(
            data = dsub,
            ggplot2::aes(x = .data[[ax1]], y = .data[[ax2]]),
            type = ellipse_type, level = ellipse_level,
            linewidth = 0.5, alpha = 0.7, color = base_map[[g]],
            show.legend = FALSE
          )
        }
      }
    }
  }

  # ---- centroids -------------------------------------------------------------
  if (isTRUE(show_centroids) && !is.null(cent) && nrow(cent)) {
    if (has_group && has_timeF && "gt" %in% names(cent)) {
      # same shaded legend as points
      p <- p + ggplot2::geom_point(
        data = cent,
        ggplot2::aes(x = .data$cx, y = .data$cy, fill = .data$gt),
        shape = 21, size = centroid_size, stroke = 0.9,
        color = base_map[as.character(cent$group_lab)]
      )
      if (connect_centroids == "group") {
        for (g in group_levels) {
          df_g <- cent[cent$group_lab == g, , drop = FALSE]
          if (nrow(df_g) > 1) {
            df_g <- df_g[order(df_g$time_ord), , drop = FALSE]
            p <- p + ggplot2::geom_path(
              data = df_g, ggplot2::aes(x = .data$cx, y = .data$cy),
              color = base_map[[g]], linewidth = 0.7, alpha = 0.7, show.legend = FALSE
            )
          }
        }
      } else if (connect_centroids == "time") {
        for (tlev in time_levels) {
          df_t <- cent[cent$time_lab == tlev, , drop = FALSE]
          if (nrow(df_t) > 1) {
            p <- p + ggplot2::geom_path(
              data = df_t, ggplot2::aes(x = .data$cx, y = .data$cy),
              color = "grey50", linewidth = 0.6, alpha = 0.6, show.legend = FALSE
            )
          }
        }
      }
    } else if ("group_lab" %in% names(cent)) {
      p <- p + ggplot2::geom_point(
        data = cent,
        ggplot2::aes(x = .data$cx, y = .data$cy, color = .data$group_lab),
        size = centroid_size
      ) +
        ggplot2::scale_color_manual(values = base_map, name = group_col)
    } else if ("time_lab" %in% names(cent)) {
      tcols <- stats::setNames(.phip_pal(length(time_levels))[seq_along(time_levels)], time_levels)
      p <- p + ggplot2::geom_point(
        data = cent,
        ggplot2::aes(x = .data$cx, y = .data$cy, color = .data$time_lab),
        size = centroid_size
      ) +
        ggplot2::scale_color_manual(values = tcols, name = time_fac)
    }
  }

  # ---- optional CAP biplot overlay ------------------------------------------
  if (isTRUE(add_biplot) && identical(tolower(method), "cap")) {
    rank_nm <- (attr(x, "ranks") %||% names(sv$cap_fits))[1]
    cap_fit <- sv$cap_fits[[rank_nm]]
    if (!is.null(cap_fit)) {
      bp <- as.data.frame(vegan::scores(cap_fit, display = "bp", choices = axes))
      if (nrow(bp)) {
        bp$term <- rownames(bp)
        if (has_timeF) {
          bp <- bp[grepl("^group", bp$term) |
            grepl("^time_", bp$term) |
            grepl(":time_", bp$term), , drop = FALSE]
        } else {
          bp <- bp[grepl("^group", bp$term) |
            grepl("time_cont$", bp$term) |
            grepl(":time_cont$", bp$term), , drop = FALSE]
        }
        if (nrow(bp)) {
          sites <- as.data.frame(vegan::scores(cap_fit, display = "sites", choices = axes))
          max_site <- max(abs(unlist(sites[, 1:2])), na.rm = TRUE)
          max_bp <- max(abs(unlist(bp[, 1:2])), na.rm = TRUE)
          scl <- if (is.finite(max_site) && is.finite(max_bp) && max_bp > 0) 0.95 * max_site / max_bp else 1
          bp$X1s <- bp[[1]] * scl
          bp$X2s <- bp[[2]] * scl

          p <- p +
            ggplot2::geom_segment(
              data = bp,
              ggplot2::aes(x = 0, y = 0, xend = .data$X1s, yend = .data$X2s),
              arrow = grid::arrow(length = grid::unit(0.02, "npc")),
              linewidth = 0.7, color = "grey25", alpha = 0.9, show.legend = FALSE
            ) +
            ggrepel::geom_text_repel(
              data = bp,
              ggplot2::aes(x = .data$X1s, y = .data$X2s, label = .data$term),
              size = biplot_label_size, segment.color = NA, color = "grey20", show.legend = FALSE
            )
        }
      }
    } else {
      .ph_info("CAP fit not found in object; skipping biplot overlay.")
    }
  }

  # ---- finish ----------------------------------------------------------------
  p <- p +
    ggplot2::labs(x = xlab, y = ylab) +
    theme_phip() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  .ph_info("Finished building plot.")
  p
}

# ======================================================================
# CAP (capscale) small-multiples: CAP{1,2} ~ time_cont  (by group)
# ======================================================================

#' Plot CAP axes vs time (capscale; continuous time)
#'
#' @param beta A `phip_beta_diversity` object produced with method_pcoa = "cap".
#' @param axes Character vector of CAP axes to plot (default c("CAP1","CAP2")).
#' @param point_alpha Point transparency (default 0.25).
#' @param point_size  Point size (default 1).
#' @param linewidth   Smooth line width (default 1).
#' @param free_y      Use free y-scales across facets (default TRUE).
#' @return A ggplot object (facetted small multiples).
#' @examples
#' \dontrun{
#' p <- ph_plot_cap_axes_vs_time(beta_cap)
#' print(p)
#' }
ph_plot_cap_axes_vs_time <- function(beta,
                                     axes = c("CAP1", "CAP2"),
                                     point_alpha = 0.25,
                                     point_size = 1,
                                     linewidth = 1,
                                     free_y = TRUE) {
  .ph_with_timing("CAP: axes vs time", {
    # --- sanity checks -------------------------------------------------
    if (!inherits(beta, "phip_beta_diversity")) {
      .ph_abort("`beta` is not a phip_beta_diversity object")
    }
    method <- attr(beta, "method_pcoa") %||% NA_character_
    if (!identical(method, "cap")) {
      .ph_abort("This plot requires method_pcoa = 'cap' (got {method})")
    }

    view_name <- attr(beta, "group_cols") %||% "big_group"
    if (!view_name %in% names(beta)) {
      .ph_abort("View '{view_name}' not found in `beta`")
    }
    tab <- beta[[view_name]][["CAP_global"]]
    if (is.null(tab) || !is.list(tab)) {
      .ph_abort("Missing `CAP_global` in view '{view_name}'")
    }

    pcoa <- tab$pcoa
    needed <- c("sample_id", "group", "time_cont")
    if (!all(needed %in% names(pcoa))) {
      .ph_abort("`pcoa` tibble lacks required columns: {paste(setdiff(needed, names(pcoa)), collapse=', ')}")
    }
    if (!all(axes %in% names(pcoa))) {
      .ph_abort("Requested axes not present in `pcoa`: {paste(setdiff(axes, names(pcoa)), collapse=', ')}")
    }
    if (all(is.na(pcoa$time_cont))) {
      .ph_abort("`time_cont` is all NA; this plot expects continuous time")
    }

    # --- data prep -----------------------------------------------------
    df_long <- pcoa |>
      dplyr::select(dplyr::all_of(c("sample_id", "group", "time_cont", axes))) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(axes),
        names_to = "axis",
        values_to = "score"
      )

    .ph_log_info("Fitting LM smooths per group & axis",
      step = "CAP_axes_vs_time",
      bullets = sprintf(
        "n = %s rows; axes = %s",
        format(nrow(df_long), big.mark = "'"),
        paste(axes, collapse = ", ")
      )
    )

    # --- plot ----------------------------------------------------------
    p <- ggplot2::ggplot(df_long, ggplot2::aes(time_cont, score, color = group)) +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = linewidth) +
      ggplot2::facet_wrap(~axis, scales = if (isTRUE(free_y)) "free_y" else "fixed") +
      scale_color_phip(name = "group") +
      ggplot2::labs(
        title = "CAP (constrained) axes vs time",
        subtitle = sprintf("view = %s | method = capscale | continuous time", view_name),
        x = "time (continuous)", y = "CAP score"
      ) +
      theme_phip()

    .ph_log_ok("Plot ready", step = "CAP_axes_vs_time")
    return(p)
  })
}

# ======================================================================
# Dispersion (betadisper distance to group centroid) ~ time_cont
# ======================================================================

#' Plot dispersion (distance to centroid) vs time (capscale; continuous time)
#'
#' @param beta A `phip_beta_diversity` object produced with method_pcoa = "cap".
#' @param point_alpha Point transparency (default 0.25).
#' @param point_size  Point size (default 1.2).
#' @param linewidth   Smooth line width (default 1).
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' p <- ph_plot_dispersion_vs_time(beta_cap)
#' print(p)
#' }
ph_plot_dispersion_vs_time <- function(beta,
                                       point_alpha = 0.25,
                                       point_size = 1.2,
                                       linewidth = 1) {
  .ph_with_timing("CAP: dispersion vs time", {
    # --- sanity checks -------------------------------------------------
    if (!inherits(beta, "phip_beta_diversity")) {
      .ph_abort("`beta` is not a phip_beta_diversity object")
    }
    method <- attr(beta, "method_pcoa") %||% NA_character_
    if (!identical(method, "cap")) {
      .ph_abort("This plot requires method_pcoa = 'cap' (got {method})")
    }

    view_name <- attr(beta, "group_cols") %||% "big_group"
    if (!view_name %in% names(beta)) {
      .ph_abort("View '{view_name}' not found in `beta`")
    }
    tab <- beta[[view_name]][["CAP_global"]]
    if (is.null(tab) || !is.list(tab)) {
      .ph_abort("Missing `CAP_global` in view '{view_name}'")
    }

    pcoa <- tab$pcoa
    if (!all(c("sample_id", "group", "time_cont") %in% names(pcoa))) {
      .ph_abort("`pcoa` tibble lacks required columns `sample_id`, `group`, `time_cont`")
    }
    if (all(is.na(pcoa$time_cont))) {
      .ph_abort("`time_cont` is all NA; this plot expects continuous time")
    }

    disp <- tab$dispersion
    if (is.null(disp) || !"distance" %in% names(disp)) {
      .ph_abort("`dispersion` table is missing or lacks `distance` column")
    }

    # --- subset to global group dispersion and join time ----------------
    disp2 <- disp |>
      dplyr::filter(.data$scope == "group", .data$contrast == "<global>") |>
      dplyr::left_join(
        pcoa |> dplyr::select(sample_id, group, time_cont),
        by = "sample_id"
      )

    if (!nrow(disp2)) {
      .ph_abort("No dispersion rows for scope='group' & contrast='<global>'")
    }

    .ph_log_info("Fitting LM smooths per group",
      step = "dispersion_vs_time",
      bullets = sprintf("n = %s rows", format(nrow(disp2), big.mark = "'"))
    )

    # --- plot ----------------------------------------------------------
    p <- ggplot2::ggplot(disp2, ggplot2::aes(time_cont, distance, color = group)) +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = linewidth) +
      scale_color_phip(name = "group") +
      ggplot2::labs(
        title = "Dispersion vs time",
        subtitle = sprintf("distance to group centroid | view = %s", view_name),
        x = "time (continuous)", y = "distance"
      ) +
      theme_phip()

    .ph_log_ok("Plot ready", step = "dispersion_vs_time")
    return(p)
  })
}

# ======================================================================
# phiper-style: Dispersion (distance to group centroid) — boxplots + Wilcoxon
# Works with ANY phip_beta_diversity object (capscale or plain PCoA), single function.
# No extra helpers beyond phiper logging + theme_phip()/scale_*_phip().
# ======================================================================

#' Plot dispersion (distance to centroid) by group
#'
#' @param beta            A `phip_beta_diversity` object (works for method_pcoa = "cap" or "pcoa").
#' @param view            Name of the view; by default uses attr(beta, "group_cols") or "big_group".
#' @param scope_filter    Dispersion scope to use (default "group").
#' @param contrast_filter Dispersion contrast to use (default "<global>").
#' @param sig_level       p-value threshold for pairwise Wilcoxon annotations (default 0.05).
#' @param label_format    Passed to ggpubr::stat_compare_means (e.g., "p.format" or "p.signif").
#' @param custom_colors   Optional named vector for group fill colors (names = groups). If NULL, uses scale_fill_phip().
#' @param show_points     Add jittered points on top of boxplots (default TRUE).
#' @param point_alpha     Point alpha (default 0.30).
#' @param point_size      Point size (default 1).
#' @param rotate_x        Rotate x labels by 45° (default TRUE).
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' p <- ph_plot_dispersion_box(beta_cap)
#' print(p)
#' }
ph_plot_dispersion_box <- function(beta,
                                   view = attr(beta, "group_cols") %||% "big_group",
                                   scope_filter = "group",
                                   contrast_filter = "<global>",
                                   sig_level = 0.05,
                                   label_format = "p.format",
                                   custom_colors = NULL,
                                   show_points = TRUE,
                                   point_alpha = 0.30,
                                   point_size = 1,
                                   rotate_x = TRUE) {
  .ph_with_timing("Dispersion boxplots", {
    # --- checks --------------------------------------------------------
    if (!inherits(beta, "phip_beta_diversity")) {
      .ph_abort("`beta` is not a phip_beta_diversity object")
    }
    if (!view %in% names(beta)) {
      .ph_abort("View '{view}' not found in `beta`")
    }
    view_tab <- beta[[view]]
    if (!is.list(view_tab) || length(view_tab) == 0L) {
      .ph_abort("Empty view '{view}' in `beta`")
    }

    # --- locate the sublist holding dispersion (CAP_global / PCoA_global / etc.) ----
    sub_names <- names(view_tab)
    hit <- NULL
    for (nm in sub_names) {
      cand <- view_tab[[nm]]
      if (is.list(cand) && ("dispersion" %in% names(cand))) {
        hit <- nm
        break
      }
    }
    if (is.null(hit)) {
      .ph_abort("Could not find a component with a `dispersion` table inside view '{view}'")
    }
    tab <- view_tab[[hit]]

    disp <- tab$dispersion
    if (is.null(disp) || !"distance" %in% names(disp)) {
      .ph_abort("`dispersion` table is missing or lacks `distance` column")
    }

    # --- filter to desired scope/contrast and prepare data -------------
    disp2 <- disp |>
      dplyr::filter(.data$scope == !!scope_filter, .data$contrast == !!contrast_filter)

    if (!nrow(disp2)) {
      .ph_abort("No dispersion rows for scope='{scope_filter}' & contrast='{contrast_filter}'")
    }

    # Most beta outputs store the grouping column as `group_disp`; fall back to `group` if needed.
    grp_col <- if ("group_disp" %in% names(disp2)) "group_disp" else if ("group" %in% names(disp2)) "group" else NULL
    if (is.null(grp_col)) .ph_abort("Could not find grouping column (`group_disp`/`group`) in dispersion table")
    grp_sym <- rlang::sym(grp_col)

    # Drop NA groups/distances
    disp2 <- disp2 |>
      dplyr::filter(!is.na(!!grp_sym), !is.na(.data$distance))

    if (!nrow(disp2)) .ph_abort("All rows are NA after cleaning group/distance")

    # counts for axis labels
    df_counts <- disp2 |>
      dplyr::count(!!grp_sym, name = "sample_count") |>
      dplyr::ungroup()

    x_labels <- stats::setNames(
      paste0(df_counts[[grp_col]], "\n(n = ", df_counts$sample_count, ")"),
      df_counts[[grp_col]]
    )

    # pairwise comparisons over present groups
    groups_present <- levels(factor(disp2[[grp_col]]))
    pairwise_comparisons <- utils::combn(groups_present, 2, simplify = FALSE)

    # compute Wilcoxon for each pair and keep significant
    sig_comparisons <- purrr::keep(pairwise_comparisons, function(pair) {
      g1 <- pair[1]
      g2 <- pair[2]
      x <- disp2 |>
        dplyr::filter(!!grp_sym == g1) |>
        dplyr::pull(.data$distance)
      y <- disp2 |>
        dplyr::filter(!!grp_sym == g2) |>
        dplyr::pull(.data$distance)
      if (!length(x) || !length(y)) {
        return(FALSE)
      }
      stats::wilcox.test(x, y, exact = FALSE)$p.value < sig_level
    })

    .ph_log_info("Pairwise Wilcoxon tests",
      step = "dispersion_box",
      bullets = c(
        paste("view:", view, "| component:", hit),
        paste("groups:", paste(groups_present, collapse = ", ")),
        paste("n rows:", nrow(disp2)),
        paste("significant pairs:", length(sig_comparisons))
      )
    )

    # --- plot ----------------------------------------------------------
    p <- ggplot2::ggplot(disp2, ggplot2::aes(x = !!grp_sym, y = .data$distance, fill = !!grp_sym)) +
      ggplot2::geom_boxplot(show.legend = FALSE, outlier.shape = NA)

    if (isTRUE(show_points)) {
      p <- p + ggplot2::geom_jitter(
        color = "black", size = point_size, width = 0.2,
        alpha = point_alpha, show.legend = FALSE
      )
    }

    if (!is.null(custom_colors)) {
      p <- p + ggplot2::scale_fill_manual(values = custom_colors)
    } else {
      # falls back to your package palette
      p <- p + scale_fill_phip()
    }

    p <- p +
      ggplot2::scale_x_discrete(labels = x_labels) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
      ggplot2::labs(
        title = "Dispersion by group",
        subtitle = sprintf(
          "Distance to group centroid · view = %s · scope = %s · contrast = %s",
          view, scope_filter, contrast_filter
        ),
        x = "Group",
        y = "Distance to centroid"
      ) +
      theme_phip() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 1, 0, 1, unit = "pt")
      )

    if (isTRUE(rotate_x)) {
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.6, hjust = 0.5))
    }

    # Add significant pair annotations if any
    if (length(sig_comparisons) > 0) {
      p <- p + ggpubr::stat_compare_means(
        method      = "wilcox.test",
        comparisons = sig_comparisons,
        label       = label_format,
        hide.ns     = FALSE,
        size        = 4.5,
        tip.length  = 0.02
      )
    }

    .ph_log_ok("Plot ready", step = "dispersion_box")
    return(p)
  })
}
