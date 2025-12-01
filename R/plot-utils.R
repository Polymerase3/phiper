# keep ggplot2’s replace operator available
"%+replace%" <- ggplot2::"%+replace%"

#' @title PHIP default colour palette
#'
#' @description A concise, publication-ready qualitative palette tailored for
#'   PHIP figures. Colors are stable across sessions to keep styling
#'   reproducible.
#'
#' @format A character vector of hex colors.
#' @export
phip_palette <- c(
  # --- original 12 ------------------------------------------------------------
  "#1f78b4", # dark blue
  "#fb9a99", # light salmon
  "#33a02c", # medium green
  "#fdbf6f", # soft orange
  "#6a3d9a", # purple
  "#b2df8a", # pale green
  "#e31a1c", # red
  "#a6cee3", # pale blue
  "#ff7f00", # vivid orange
  "#cab2d6", # lavender
  "#b15928", # brown
  "#ffff99", # light yellow

  # --- extension (18) ---------------------------------------------------------
  "#66c2a5", # teal green
  "#fc8d62", # salmon orange
  "#8da0cb", # periwinkle
  "#e78ac3", # pink-magenta
  "#a6d854", # lime green
  "#ffd92f", # medium yellow
  "#e5c494", # tan
  "#b3b3b3", # neutral grey

  "#1b9e77", # deep teal
  "#d95f02", # burnt orange
  "#7570b3", # indigo
  "#e7298a", # fuchsia
  "#66a61e", # olive green
  "#e6ab02", # goldenrod
  "#a6761d", # ochre brown
  "#666666", # dark grey

  "#7fc97f", # soft green
  "#386cb0" # strong blue
)

#' @title Discrete colour & fill scales using the PHIP palette
#'
#' @description Thin wrappers around `ggplot2::scale_*_manual()` that apply
#'   `phip_palette`.
#'
#' @inheritParams ggplot2::scale_colour_manual
#' @return A ggplot2 scale
#' @family phip-ggplot
#' @export
scale_colour_phip <- function(...) {
  ggplot2::scale_colour_manual(values = phip_palette, ...)
}

#' @rdname scale_colour_phip
#' @export
scale_color_phip <- scale_colour_phip

#' @rdname scale_colour_phip
#' @export
scale_fill_phip <- function(...) {
  ggplot2::scale_fill_manual(values = phip_palette, ...)
}

# ------------------------------------------------------------------------------
# Fonts via showtext/sysfonts
# ------------------------------------------------------------------------------
#' @title Register and enable the Montserrat font for plotting (showtext)
#'
#' @description Registers the Google font **Montserrat** (via **sysfonts**) and
#'   enables **showtext** so the font renders consistently in all devices (PNG,
#'   PDF, etc.). Safe to call multiple times; silently no-ops if already active.
#'
#' @param family Internal family alias to register. Default `"Montserrat"`.
#' @param enable If `TRUE`, turns on `showtext::showtext_auto(TRUE)`.
#'
#' @details Requires packages **showtext** and **sysfonts**. If the Google
#'   download fails (e.g., offline), the theme will still try to use a locally
#'   installed Montserrat by name; otherwise it falls back to the default device
#'   font.
#'
#' @return (Invisibly) the family name to use in themes.
#' @family phip-ggplot
#' @export
phip_use_montserrat <- function(family = "Montserrat",
                                enable = TRUE) {
  if (!requireNamespace("sysfonts", quietly = TRUE) ||
    !requireNamespace("showtext", quietly = TRUE)) {
    .ph_warn(
      headline = "Missing dependencies for Montserrat registration.",
      step = "phip_use_montserrat()",
      bullets = c(
        "Packages 'sysfonts' and 'showtext' are required.",
        "Install with: install.packages(c('sysfonts','showtext'))"
      )
    )

    return(invisible(family))
  }

  # Try to add from Google; ignore errors if it already exists or offline.
  try(
    {
      sysfonts::font_add_google(name = "Montserrat", family = family)
    },
    silent = TRUE
  )

  if (isTRUE(enable)) {
    showtext::showtext_auto(enable = TRUE)
  }
  invisible(family)
}

# ------------------------------------------------------------------------------
# Theme
# ------------------------------------------------------------------------------
#' @title Theme `theme_phip`
#'
#' @description A clean, publication-ready ggplot2 theme inspired by the
#'   provided `theme_Publication` snippet, tuned for **facetted** plots and
#'   consistent use of the **Montserrat** font (register it with
#'   [phip_use_montserrat()]).
#'
#' @param base_size   Base font size.
#' @param base_family Base font family (default `"Montserrat"`). Call
#'   [phip_use_montserrat()] once per session to register and enable rendering.
#'
#' @return A ggplot2 `theme` object.
#' @family phip-ggplot
#' @examples
#' \dontrun{
#' # Register Montserrat once per session
#' phip_use_montserrat()
#'
#' ggplot2::ggplot(iris, ggplot2::aes(Sepal.Length, Sepal.Width, colour = Species)) +
#'   ggplot2::geom_point() +
#'   scale_colour_phip() +
#'   theme_phip()
#' }
#' @export
theme_phip <- function(base_size = 14,
                       base_family = "Montserrat") {
  # Start from a minimal, unobtrusive foundation
  ggplot2::theme_minimal(
    base_size = base_size,
    base_family = base_family
  ) %+replace%
    ggplot2::theme(
      # titles & text ----------------------------------------------------------
      plot.title = ggplot2::element_text(
        family = base_family, face = "bold",
        size = ggplot2::rel(1.4), hjust = 0.5,
        margin = ggplot2::margin(b = 18)
      ),
      plot.subtitle = ggplot2::element_text(
        family = base_family, margin = ggplot2::margin(b = 8)
      ),
      plot.caption = ggplot2::element_text(
        family = base_family, size = ggplot2::rel(0.9),
        hjust = 1, margin = ggplot2::margin(t = 8)
      ),
      text = ggplot2::element_text(family = base_family),
      axis.title = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1)
      ),
      axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
      axis.title.x = ggplot2::element_text(vjust = -0.2),
      axis.text = ggplot2::element_text(),

      # panel & background -----------------------------------------------------
      panel.background = ggplot2::element_rect(fill = NA, colour = NA),
      plot.background = ggplot2::element_rect(fill = NA, colour = NA),
      panel.border = ggplot2::element_rect(fill = NA, colour = NA),


      # axes & gridlines -------------------------------------------------------
      axis.line = ggplot2::element_line(
        colour = "black",
        linewidth = 0.4
      ),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      # panel.grid.major = ggplot2::element_line(colour = "#f0f0f0",
      #                                          linewidth = 0.4),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # legend -----------------------------------------------------------------
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key = ggplot2::element_rect(colour = NA, fill = NA),
      legend.key.size = grid::unit(0.2, "cm"),
      legend.spacing.x = grid::unit(4, "pt"),
      legend.box.margin = grid::unit(c(0, 0, 0, 0), "cm"),
      legend.title = ggplot2::element_text(face = "italic"),

      # Margins ----------------------------------------------------------------
      plot.margin = grid::unit(c(10, 5, 5, 5), "mm"),

      # Facets ------------------------------------------------------------------
      strip.background = ggplot2::element_rect(
        fill = "#f0f0f0",
        colour = "black",
        linewidth = 0.6
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.1),
        margin = ggplot2::margin(
          t = 4,
          b = 4
        )
      ),
      panel.spacing = grid::unit(6, "pt")
    )
}

# -------------------------------------------------------------------
# internal: per-plot PHIP palette mapping (recycled & named)
# -------------------------------------------------------------------
.phip_palette_map <- function(levels_vec, palette = phip_palette) {
  if (length(levels_vec) == 0L) {
    return(stats::setNames(character(0), character(0)))
  }
  stats::setNames(rep_len(palette, length(levels_vec)), levels_vec)
}

# ---- helpers -----------------------------------------------------------------
# ---------- tiny logging shims (use your phiper ones if available) ----------
.phi <- get0(".ph_log_info", ifnotfound = function(msg, step = "plot_beta_pcoa") {
  cli::cat_line(sprintf(
    "[%s] INFO  %s\n                 -> %s",
    format(Sys.time(), "%H:%M:%S"), msg, step
  ))
})
.phe <- get0(".ph_abort", ifnotfound = function(msg) {
  rlang::abort(msg)
})

# ------------------------------- palette utils ------------------------------
.get_palette <- function(n) {
  pal <- get0("phip_palette", ifnotfound = NULL, inherits = TRUE)
  if (is.null(pal)) {
    return(scales::hue_pal()(n))
  }
  if (is.function(pal)) {
    return(pal(n))
  }
  rep_len(pal, n)
}


# hex <-> rgb (robust, no named indexing; clamps to [0,1])
.hex2rgb <- function(hex) {
  if (is.null(hex) || is.na(hex) || !nzchar(hex)) {
    return(c(NA_real_, NA_real_, NA_real_))
  }
  as.numeric(grDevices::col2rgb(hex)) / 255
}
.rgb2hex <- function(x) {
  if (length(x) < 3 || anyNA(x)) {
    return("#B3B3B3")
  }
  x <- pmin(1, pmax(0, as.numeric(x[1:3])))
  grDevices::rgb(x[1], x[2], x[3])
}
.mix_cols <- function(hex_a, hex_b, t) {
  if (is.na(hex_a) || is.na(hex_b) || is.na(t)) {
    return("#B3B3B3")
  }
  a <- .hex2rgb(hex_a)
  b <- .hex2rgb(hex_b)
  if (anyNA(a) || anyNA(b)) {
    return("#B3B3B3")
  }
  .rgb2hex((1 - t) * a + t * b)
}
.tint <- function(hex, t_light = 0) .mix_cols(hex, "#FFFFFF", t_light)

.ph_cols <- function(n) rep_len(phip_palette, n)

.pick_axes <- function(df, axes = c(1, 2)) {
  stopifnot(length(axes) == 2)
  for (stem in c("CAP", "PCoA", "PC", "MDS")) {
    cols <- paste0(stem, axes)
    if (all(cols %in% names(df))) {
      return(cols)
    }
  }
  .ph_abort(sprintf(
    "Could not find axis columns %s (looked for CAP/PCoA/PC/MDS).",
    paste(axes, collapse = ":")
  ))
}

.axis_labels_with_pct <- function(var_ex, axis_cols) {
  if (is.null(var_ex) || !nrow(var_ex)) {
    setNames(axis_cols, axis_cols)
  } else {
    pct_cols <- paste0("%", axis_cols)
    miss <- setdiff(pct_cols, names(var_ex))
    if (length(miss)) {
      .ph_log_info(sprintf("No %% columns for %s — using raw axis names.", paste(miss, collapse = ", ")),
        step = "plot_beta_pcoa"
      )
      setNames(axis_cols, axis_cols)
    } else {
      pct <- suppressWarnings(as.numeric(var_ex[1, pct_cols, drop = TRUE]))
      nms <- sprintf("%s (%.1f%%)", axis_cols, pct)
      setNames(nms, axis_cols)
    }
  }
}

.centroids_df <- function(df, axes_cols, by) {
  if (is.null(by) || !(by %in% names(df))) {
    return(NULL)
  }
  stats::aggregate(df[axes_cols], by = list(df[[by]]), FUN = mean, na.rm = TRUE) |>
    `names<-`(c(by, axes_cols))
}

# Shade points by time (continuous or fac) on top of base group color.
.shaded_colors <- function(df,
                           group_col = "group",
                           time_cont = "time_cont",
                           time_fac = "time_fac",
                           shade = c("auto", "time", "group", "none")) {
  shade <- match.arg(shade)

  # base color per group (robust to missing)
  grp_chr <- as.character(df[[group_col]])
  grp_chr[is.na(grp_chr)] <- "<missing>"
  groups <- droplevels(factor(grp_chr, levels = unique(grp_chr)))
  base_map <- setNames(.ph_cols(nlevels(groups)), levels(groups))
  cols <- unname(base_map[as.character(groups)])

  has_tc <- time_cont %in% names(df) && any(!is.na(df[[time_cont]]))
  has_tf <- time_fac %in% names(df) && any(!is.na(df[[time_fac]]))

  if (shade %in% c("auto", "time")) {
    if (has_tc) {
      # normalize within group for nicer gradients; default mid if degenerate
      norm_t <- ave(df[[time_cont]], groups, FUN = function(z) {
        z <- as.numeric(z)
        rng <- range(z, na.rm = TRUE)
        if (!is.finite(diff(rng)) || diff(rng) == 0) {
          return(rep(0.5, length(z)))
        }
        (z - rng[1]) / diff(rng)
      })
      norm_t[is.na(norm_t)] <- 0.5
      cols <- mapply(function(hex, t) .tint(hex, 0.35 * t), cols, norm_t, USE.NAMES = FALSE)
    } else if (has_tf) {
      lev <- levels(droplevels(factor(df[[time_fac]])))
      rnk <- as.numeric(factor(df[[time_fac]], levels = lev))
      tt <- if (max(rnk, na.rm = TRUE) > 1) (rnk - 1) / (max(rnk, na.rm = TRUE) - 1) else 0.5
      tt[is.na(tt)] <- 0.5
      cols <- mapply(function(hex, t) .tint(hex, 0.35 * t), cols, tt, USE.NAMES = FALSE)
    }
  }

  # final guard
  cols[is.na(cols)] <- "#B3B3B3"
  cols
}

.first_subview_name <- function(x_view) {
  nms <- names(x_view)
  if (length(nms) == 0) .ph_abort("No subviews found under the selected view.")
  nms[1]
}

# shaded fills for points; keeps group legend by mapping color=group and fill=identity
.make_point_fills <- function(df, base_map, group_col = "group",
                              time_fac = "time_fac", time_cont = "time_cont",
                              mode = c("auto", "time", "group", "none"), max_tint = 0.35) {
  mode <- match.arg(mode)
  grp <- as.character(df[[group_col]])
  base <- unname(base_map[grp])

  has_tc <- time_cont %in% names(df) && any(!is.na(df[[time_cont]]))
  has_tf <- time_fac %in% names(df) && any(!is.na(df[[time_fac]]))

  if (mode == "auto") mode <- if (has_tc || has_tf) "time" else "group"

  if (mode == "time" && has_tc) {
    # normalise time within group for nicer within-group gradients
    t_norm <- ave(df[[time_cont]], grp, FUN = function(z) {
      z <- as.numeric(z)
      r <- range(z, na.rm = TRUE)
      if (!is.finite(diff(r)) || diff(r) == 0) {
        return(rep(0.5, length(z)))
      }
      (z - r[1]) / diff(r)
    })
    t_norm[is.na(t_norm)] <- 0.5
    return(mapply(function(hex, t) .tint(hex, max_tint * t), base, t_norm, USE.NAMES = FALSE))
  }

  if (mode == "time" && has_tf) {
    lev <- levels(droplevels(factor(df[[time_fac]])))
    rnk <- as.numeric(factor(df[[time_fac]], levels = lev))
    t_norm <- if (length(lev) > 1) (rnk - 1) / (length(lev) - 1) else rep(0.5, length(rnk))
    t_norm[is.na(t_norm)] <- 0.5
    return(mapply(function(hex, t) .tint(hex, max_tint * t), base, t_norm, USE.NAMES = FALSE))
  }

  # fallback: just base colors (no shading)
  base
}

# ---- color helpers -----------------------------------------------------------
.blend_hex <- function(c1, c2 = "#FFFFFF", w = 0.5) {
  # simple RGB blend: (1 - w) * c1 + w * c2
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

# Build fill/outline mappings for interaction(group, time)
.build_shaded_map <- function(groups, times, base_colors = NULL,
                              default_palette = phip_palette,
                              shade_range = c(0.15, 0.85),
                              sep = " | ") {
  groups <- as.character(groups)
  times <- as.character(times)

  # default bases from phip_palette if not provided
  if (is.null(base_colors)) {
    base_colors <- stats::setNames(default_palette[seq_along(groups)], groups)
  } else {
    # recycle / align names
    if (is.null(names(base_colors))) names(base_colors) <- groups
    base_colors <- base_colors[groups]
  }

  fill_vals <- character(0)
  outline_vals <- character(0)

  for (g in groups) {
    shades <- .make_shades(base_colors[[g]], length(times), range = shade_range)
    names(shades) <- paste(g, times, sep = sep)
    fill_vals <- c(fill_vals, shades)
    outline_vals <- c(outline_vals, stats::setNames(rep(base_colors[[g]], length(times)), names(shades)))
  }
  list(fill = fill_vals, outline = outline_vals)
}
