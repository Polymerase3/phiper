#' @importFrom ggplot2 %+replace%
NULL

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

#' @title Discrete colour scale using the PHIP palette
#'
#' @description A thin wrapper around `ggplot2::scale_colour_manual()` that
#'   applies `phip_palette`. Set as the session default on `library(phiper)`.
#'
#' @inheritParams ggplot2::scale_colour_manual
#' @return A ggplot2 scale.
#' @family phip-ggplot
#' @examples
#' ggplot2::ggplot(iris,
#'   ggplot2::aes(Sepal.Length, Sepal.Width, colour = Species)) +
#'   ggplot2::geom_point() +
#'   scale_colour_phip()
#' @export
scale_colour_phip <- function(...) {
  ggplot2::scale_colour_manual(values = phip_palette, ...)
}

#' @rdname scale_colour_phip
#' @export
scale_color_phip <- scale_colour_phip

#' @title Discrete fill scale using the PHIP palette
#'
#' @description A thin wrapper around `ggplot2::scale_fill_manual()` that
#'   applies `phip_palette`. Set as the session default on `library(phiper)`.
#'
#' @inheritParams ggplot2::scale_fill_manual
#' @return A ggplot2 scale.
#' @family phip-ggplot
#' @examples
#' ggplot2::ggplot(iris,
#'   ggplot2::aes(Species, Sepal.Length, fill = Species)) +
#'   ggplot2::geom_boxplot() +
#'   scale_fill_phip()
#' @export
scale_fill_phip <- function(...) {
  ggplot2::scale_fill_manual(values = phip_palette, ...)
}

# ------------------------------------------------------------------------------
# Theme
# ------------------------------------------------------------------------------
#' @title Theme `theme_phip`
#'
#' @description A clean, publication-ready ggplot2 theme tuned for **facetted**
#'   plots, with \code{phip_palette} as the default discrete colour and
#'   fill palette.
#'
#' @param base_size   Base font size.
#' @param base_family Base font family (default `"sans"`).
#'
#' @return A ggplot2 `theme` object.
#' @family phip-ggplot
#' @examples
#' \dontrun{
#' ggplot2::ggplot(iris, ggplot2::aes(Sepal.Length, Sepal.Width, colour = Species)) +
#'   ggplot2::geom_point() +
#'   scale_colour_phip() +
#'   theme_phip()
#' }
#' @export
theme_phip <- function(base_size = 14,
                       base_family = "sans") {
  # Start from a minimal, unobtrusive foundation
  ggplot2::theme_minimal(
    base_size = base_size,
    base_family = base_family
  ) %+replace%
    ggplot2::theme(
      # default discrete palettes ----------------------------------------------
      palette.colour.discrete = phip_palette,
      palette.fill.discrete   = phip_palette,
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

# -------------------------------------------------------------------
# internal: axis / label formatters
# -------------------------------------------------------------------
# proportion -> "12.3%"; `accuracy` is the rounding step (1, 0.1, ...)
.ph_percent <- function(x, accuracy = 1) {
  digits <- max(0L, -floor(log10(accuracy)))
  # "+ 0" turns -0 into 0 so tiny negatives print as "0%", not "-0%"
  out <- paste0(formatC(round(x * 100 / accuracy) * accuracy + 0,
                        format = "f", digits = digits), "%")
  out[is.na(x)] <- NA_character_
  out
}

.ph_percent_format <- function(accuracy = 1) {
  function(x) .ph_percent(x, accuracy = accuracy)
}

# 1000 -> expression(10^3), for log10 axes
.ph_log10_labels <- function(x) {
  parse(text = paste0("10^", log10(x)))
}

# Longest chain of overlapping segments [start, end] (inclusive positions),
# where each segment overlaps the region covered by the chain so far. Ties go
# to the first chain along the protein. Returns c(n, start, end) of the chain.
.ph_overlap_run <- function(start, end) {
  o <- order(start)
  start <- start[o]
  end <- end[o]
  best <- c(1L, 1L)
  first <- 1L
  reach <- end[1]
  for (i in seq_along(start)[-1]) {
    if (start[i] <= reach) {
      reach <- max(reach, end[i])
    } else {
      first <- i
      reach <- end[i]
    }
    if (i - first > best[2] - best[1]) best <- c(first, i)
  }
  run <- best[1]:best[2]
  c(n = length(run), start = start[best[1]], end = max(end[run]))
}
