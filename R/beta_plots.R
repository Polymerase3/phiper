#' @title Plot Principal Coordinates Analysis (PCoA)
#' @description
#' Creates a 2D PCoA scatterplot from a `"beta_pcoa"` object (the output of
#' [compute_pcoa()]), with optional grouping by a categorical variable,
#' optional time factor, centroids, centroid trajectories and ellipses.
#'
#' @param pcoa_res A `"beta_pcoa"` object as returned by [compute_pcoa()].
#'   Must contain at least a `sample_coords` component with columns
#'   `sample_id` and `PCoA1`, `PCoA2`, ..., and a `var_explained` component
#'   (one-row tibble) with columns such as `"%PCoA1"`, `"%PCoA2"`, etc.
#' @param axes Integer vector of length 2 giving the PCoA axes to plot
#'   (e.g., `c(1, 2)` for PCoA1 vs PCoA2). Defaults to `c(1, 2)`.
#' @param group_col Optional name of the grouping column in
#'   `pcoa_res$sample_coords` (e.g., `"group"`, `"type_person"`). Used for
#'   colouring points and for group-level centroids / ellipses. If `NULL`,
#'   the function attempts to auto-detect `"group"` if present.
#' @param time_col Optional name of a **categorical** time column in
#'   `pcoa_res$sample_coords` (e.g., `"time"` or `"timepoint"`). Used for
#'   shaping points and for time/interaction centroids / ellipses. If `NULL`,
#'   the function attempts to auto-detect `"time"` if present. Continuous
#'   time is **not** supported here.
#' @param show_centroids Logical; if `TRUE` (default), centroids of the
#'   chosen grouping are shown.
#' @param centroid_by Granularity of centroids. One of:
#'   \itemize{
#'     \item `"auto"` (default): if both `group_col` and `time_col` are
#'       present, centroids are computed for each group-time combination
#'       (`"group_time"`); otherwise by `"group"` or `"time"` depending on
#'       which factor is available.
#'     \item `"group"`: centroids per group (requires `group_col`).
#'     \item `"time"`: centroids per time level (requires `time_col`).
#'     \item `"group_time"`: centroids per group*time combination
#'       (requires both `group_col` and `time_col`).
#'   }
#' @param connect_centroids Character; if centroids are shown and
#'   `centroid_by = "group_time"`, this controls how they are joined by
#'   paths:
#'   \itemize{
#'     \item `"none"` (default): no connecting paths;
#'     \item `"group"`: connect centroids along time within each group;
#'     \item `"time"`: connect centroids across groups within each time
#'       level.
#'   }
#'   For other `centroid_by` values, connections are ignored.
#' @param show_ellipses Logical; if `TRUE` (default), confidence ellipses
#'   (via [ggplot2::stat_ellipse()]) are drawn for selected factors.
#' @param ellipse_by Character vector specifying which factor(s) to use for
#'   ellipses. Possible values are `"group"`, `"time"`, `"group_time"`.
#'   Default is `"group"`. Values that require a missing factor are ignored.
#' @param ellipse_type Type of ellipse passed to
#'   [ggplot2::stat_ellipse()], one of `"t"`, `"norm"`, `"euclid"`.
#'   Defaults to `"t"`.
#' @param ellipse_level Numeric confidence level for ellipses (default
#'   `0.95`).
#' @param point_size Numeric size of points. Default `1.5`.
#' @param point_alpha Numeric alpha (opacity) of points in `[0, 1]`.
#'   Default `0.6`.
#' @param centroid_size Numeric size of centroid points. Default `3`.
#'
#' @details
#' The function expects that any grouping / time variables you wish to use
#' (e.g., `"type_person"`, `"age_group"`, `"timepoint"`) have already been
#' merged into `pcoa_res$sample_coords` (typically by joining with a metadata
#' table via `sample_id`).
#'
#' Axis labels are automatically annotated with the percentage of variance
#' explained, using the `var_explained` component of `pcoa_res` if available.
#'
#' Styled with [theme_phip()].

#'
#' @return A [ggplot2::ggplot] object representing the PCoA scatter plot.
#'
#' @examples
#' \dontrun{
#'   # pcoa_res <- compute_pcoa(dist_bc)
#'   # Suppose you have metadata with columns sample_id, type_person, timepoint:
#'   meta_sc <- dplyr::left_join(
#'     pcoa_res$sample_coords,
#'     meta,
#'     by = "sample_id"
#'   )
#'   pcoa_res$sample_coords <- meta_sc
#'
#'   # Basic PCoA coloured by type_person with group-level ellipses:
#'   plot_pcoa(
#'     pcoa_res,
#'     axes        = c(1, 2),
#'     group_col   = "type_person",
#'     time_col    = "timepoint",
#'     ellipse_by  = "group",
#'     show_centroids = TRUE,
#'     centroid_by    = "group_time",
#'     connect_centroids = "group"
#'   )
#' }
#' @export
plot_pcoa <- function(pcoa_res,
                      axes = c(1, 2),
                      group_col = NULL,
                      time_col = NULL,
                      show_centroids = TRUE,
                      centroid_by = c("auto", "group", "time", "group_time"),
                      connect_centroids = c("none", "group", "time"),
                      show_ellipses = TRUE,
                      ellipse_by = c("group"),
                      ellipse_type = c("t", "norm", "euclid"),
                      ellipse_level = 0.95,
                      point_size = 1.5,
                      point_alpha = 0.6,
                      centroid_size = 3) {

  # --- small helpers ---------------------------------------------------------
  .abort <- function(msg) {
    if (exists(".ph_abort", mode = "function")) {
      .ph_abort(msg, step = "plot_pcoa")
    } else {
      stop(msg, call. = FALSE)
    }
  }
  .warn <- function(msg) {
    if (exists(".ph_warn", mode = "function")) {
      .ph_warn(msg, step = "plot_pcoa")
    } else {
      warning(msg, call. = FALSE)
    }
  }
  .info <- function(msg) {
    if (exists(".ph_log_info", mode = "function")) {
      .ph_log_info(msg, step = "plot_pcoa")
    }
  }

  # --- checks on input -------------------------------------------------------
  if (is.null(pcoa_res) || !is.list(pcoa_res)) {
    .abort("`pcoa_res` must be a list (typically an object of class 'beta_pcoa').")
  }
  if (!"sample_coords" %in% names(pcoa_res)) {
    .abort("`pcoa_res` must contain a `sample_coords` component.")
  }

  centroid_by      <- match.arg(centroid_by)
  connect_centroids <- match.arg(connect_centroids)
  ellipse_type     <- match.arg(ellipse_type)

  if (length(axes) != 2L || any(!is.finite(axes))) {
    .abort("`axes` must be an integer vector of length 2, e.g. c(1, 2).")
  }
  axes <- as.integer(axes)

  df <- pcoa_res$sample_coords
  df <- tibble::as_tibble(df)

  # --- axis names + labels ---------------------------------------------------
  ax_names <- paste0("PCoA", axes)
  missing_ax <- setdiff(ax_names, names(df))
  if (length(missing_ax) > 0L) {
    .abort(
      paste0(
        "Requested axes not found in `sample_coords`: ",
        paste(missing_ax, collapse = ", ")
      )
    )
  }

  # variance explained (if available)
  ve <- pcoa_res$var_explained
  pct1 <- pct2 <- NA_real_
  if (!is.null(ve) && nrow(ve) == 1L) {
    nm1 <- paste0("%PCoA", axes[1])
    nm2 <- paste0("%PCoA", axes[2])
    if (nm1 %in% names(ve)) pct1 <- suppressWarnings(as.numeric(ve[[nm1]]))
    if (nm2 %in% names(ve)) pct2 <- suppressWarnings(as.numeric(ve[[nm2]]))
  }

  xlab <- if (is.finite(pct1)) {
    sprintf("PCoA%d (%.1f%%)", axes[1], pct1)
  } else {
    ax_names[1]
  }
  ylab <- if (is.finite(pct2)) {
    sprintf("PCoA%d (%.1f%%)", axes[2], pct2)
  } else {
    ax_names[2]
  }

  # --- auto-detect group/time columns if not given --------------------------
  if (is.null(group_col) && "group" %in% names(df)) {
    group_col <- "group"
  }
  if (is.null(time_col) && "time" %in% names(df)) {
    time_col <- "time"
  }

  has_group <- !is.null(group_col) && group_col %in% names(df)
  has_time  <- !is.null(time_col)  && time_col  %in% names(df)

  if (has_group) {
    df[[group_col]] <- droplevels(as.factor(df[[group_col]]))
  }
  if (has_time) {
    df[[time_col]] <- droplevels(as.factor(df[[time_col]]))
  }

  # Reconcile centroid_by "auto"
  if (centroid_by == "auto") {
    centroid_by <- if (has_group && has_time) {
      "group_time"
    } else if (has_group) {
      "group"
    } else if (has_time) {
      "time"
    } else {
      "group"
    }
  }

  # Sanity: if user chose centroid_by requiring data that aren't there,
  # we log and disable centroids.
  if (centroid_by == "group" && !has_group) {
    .warn("`centroid_by = 'group'` requested but `group_col` not available; centroids will be disabled.")
    show_centroids <- FALSE
  }
  if (centroid_by == "time" && !has_time) {
    .warn("`centroid_by = 'time'` requested but `time_col` not available; centroids will be disabled.")
    show_centroids <- FALSE
  }
  if (centroid_by == "group_time" && (!has_group || !has_time)) {
    .warn("`centroid_by = 'group_time'` requested but `group_col` or `time_col` not available; centroids will be disabled.")
    show_centroids <- FALSE
  }

  # Ellipse targets
  ellipse_by <- match.arg(
    ellipse_by,
    choices = c("group", "time", "group_time"),
    several.ok = TRUE
  )

  # --- centroids table (if requested) ----------------------------------------
  cent <- NULL
  if (isTRUE(show_centroids)) {
    if (centroid_by == "group" && has_group) {
      cent <- df |>
        dplyr::group_by(.data[[group_col]]) |>
        dplyr::summarise(
          cx = mean(.data[[ax_names[1]]], na.rm = TRUE),
          cy = mean(.data[[ax_names[2]]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(group_lab = !!rlang::sym(group_col))
    } else if (centroid_by == "time" && has_time) {
      cent <- df |>
        dplyr::group_by(.data[[time_col]]) |>
        dplyr::summarise(
          cx = mean(.data[[ax_names[1]]], na.rm = TRUE),
          cy = mean(.data[[ax_names[2]]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(time_lab = !!rlang::sym(time_col))
    } else if (centroid_by == "group_time" && has_group && has_time) {
      cent <- df |>
        dplyr::group_by(.data[[group_col]], .data[[time_col]]) |>
        dplyr::summarise(
          cx = mean(.data[[ax_names[1]]], na.rm = TRUE),
          cy = mean(.data[[ax_names[2]]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(
          group_lab = !!rlang::sym(group_col),
          time_lab  = !!rlang::sym(time_col)
        )
      cent$time_ord <- as.integer(cent$time_lab)
    }
  }

  .info(sprintf(
    "Plotting PCoA: n=%d | group_col=%s | time_col=%s | centroid_by=%s",
    nrow(df),
    if (has_group) group_col else "<none>",
    if (has_time)  time_col  else "<none>",
    centroid_by
  ))

  # --- base plot: points ------------------------------------------------------
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[ax_names[1]]],
                                        y = .data[[ax_names[2]]]))

  if (has_group && has_time) {
    # colour by group, shape by time
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[group_col]], shape = .data[[time_col]]),
        size  = point_size,
        alpha = point_alpha
      )
  } else if (has_group) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[group_col]]),
        size  = point_size,
        alpha = point_alpha
      )
  } else if (has_time) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[time_col]]),
        size  = point_size,
        alpha = point_alpha
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        size  = point_size,
        alpha = point_alpha,
        color = "grey40"
      )
  }

  # --- ellipses --------------------------------------------------------------
  if (isTRUE(show_ellipses) && nrow(df) >= 3L) {
    # group-level ellipses
    if ("group" %in% ellipse_by && has_group) {
      for (g in stats::na.omit(unique(df[[group_col]]))) {
        dsub <- df[df[[group_col]] == g, , drop = FALSE]
        if (nrow(dsub) < 3L) next
        p <- p +
          ggplot2::stat_ellipse(
            data  = dsub,
            ggplot2::aes(x = .data[[ax_names[1]]],
                         y = .data[[ax_names[2]]]),
            type  = ellipse_type,
            level = ellipse_level,
            linewidth = 0.7,
            alpha     = 0.7,
            show.legend = FALSE
          )
      }
    }

    # time-level ellipses
    if ("time" %in% ellipse_by && has_time) {
      for (tlev in stats::na.omit(unique(df[[time_col]]))) {
        dsub <- df[df[[time_col]] == tlev, , drop = FALSE]
        if (nrow(dsub) < 3L) next
        p <- p +
          ggplot2::stat_ellipse(
            data  = dsub,
            ggplot2::aes(x = .data[[ax_names[1]]],
                         y = .data[[ax_names[2]]]),
            type  = ellipse_type,
            level = ellipse_level,
            linewidth = 0.6,
            alpha     = 0.6,
            show.legend = FALSE
          )
      }
    }

    # group*time ellipses
    if ("group_time" %in% ellipse_by && has_group && has_time) {
      for (g in stats::na.omit(unique(df[[group_col]]))) {
        for (tlev in stats::na.omit(unique(df[[time_col]]))) {
          dsub <- df[df[[group_col]] == g & df[[time_col]] == tlev, , drop = FALSE]
          if (nrow(dsub) < 3L) next
          p <- p +
            ggplot2::stat_ellipse(
              data  = dsub,
              ggplot2::aes(x = .data[[ax_names[1]]],
                           y = .data[[ax_names[2]]]),
              type  = ellipse_type,
              level = ellipse_level,
              linewidth = 0.5,
              alpha     = 0.7,
              show.legend = FALSE
            )
        }
      }
    }
  }

  # --- centroids + optional connections --------------------------------------
  if (isTRUE(show_centroids) && !is.null(cent) && nrow(cent) > 0L) {
    if ("group_lab" %in% names(cent) && "time_lab" %in% names(cent)) {
      # group*time centroids
      p <- p +
        ggplot2::geom_point(
          data = cent,
          ggplot2::aes(x = .data$cx, y = .data$cy, color = .data$group_lab),
          size  = centroid_size,
          shape = 4,
          stroke = 1
        )

      if (connect_centroids == "group") {
        for (g in stats::na.omit(unique(cent$group_lab))) {
          dsub <- cent[cent$group_lab == g, , drop = FALSE]
          if (nrow(dsub) < 2L) next
          dsub <- dsub[order(dsub$time_ord), , drop = FALSE]
          p <- p +
            ggplot2::geom_path(
              data = dsub,
              ggplot2::aes(x = .data$cx, y = .data$cy),
              linewidth = 0.7,
              alpha     = 0.7,
              show.legend = FALSE
            )
        }
      } else if (connect_centroids == "time") {
        for (tlev in stats::na.omit(unique(cent$time_lab))) {
          dsub <- cent[cent$time_lab == tlev, , drop = FALSE]
          if (nrow(dsub) < 2L) next
          p <- p +
            ggplot2::geom_path(
              data = dsub,
              ggplot2::aes(x = .data$cx, y = .data$cy),
              linewidth = 0.7,
              alpha     = 0.7,
              show.legend = FALSE
            )
        }
      }

    } else if ("group_lab" %in% names(cent)) {
      # group-only centroids
      p <- p +
        ggplot2::geom_point(
          data = cent,
          ggplot2::aes(x = .data$cx, y = .data$cy, color = .data$group_lab),
          size  = centroid_size,
          shape = 4,
          stroke = 1
        )

    } else if ("time_lab" %in% names(cent)) {
      # time-only centroids
      p <- p +
        ggplot2::geom_point(
          data = cent,
          ggplot2::aes(x = .data$cx, y = .data$cy, color = .data$time_lab),
          size  = centroid_size,
          shape = 4,
          stroke = 1
        )
    }
  }

  # --- theme + labels --------------------------------------------------------
  base_theme <- theme_phip()

  p <- p +
    ggplot2::labs(x = xlab, y = ylab) +
    base_theme +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  p
}

#' @title Plot CAP/db-RDA Results (Constrained Ordination)
#' @description
#' Creates a 2D scatter plot from a `"beta_capscale"` object (output of
#' [compute_capscale()]), showing sample scores on CAP axes. Points can be
#' coloured by group and shaped by (categorical) time. Optional ellipses and
#' centroids can be overlaid for group/time/group*time summaries.
#'
#' @param cap_res A `"beta_capscale"` object as returned by [compute_capscale()].
#'   Must contain at least a `sample_coords` tibble with columns `sample_id`
#'   and `CAP1`, `CAP2`, ..., and a numeric `eigenvalues` vector for constrained
#'   axes.
#' @param axes Integer vector of length 2 giving the CAP axes to plot
#'   (e.g., `c(1, 2)` for CAP1 vs CAP2). Default is `c(1, 2)`.
#' @param group_col Optional name of a grouping column in
#'   `cap_res$sample_coords` used to colour points (e.g., `"group"`,
#'   `"type_person"`). If `NULL`, the function auto-detects `"group"`
#'   if present.
#' @param time_col Optional name of a **categorical** time column in
#'   `cap_res$sample_coords` used to shape points (e.g., `"time"` or
#'   `"timepoint"`). If `NULL`, the function auto-detects `"time"` if present.
#'   Continuous time is not supported in this plotting function.
#' @param show_centroids Logical; if `TRUE` (default), plots centroids for
#'   groups / time levels / group*time combinations depending on `centroid_by`.
#' @param centroid_by How to define centroids. One of
#'   `"auto"`, `"group"`, `"time"`, `"group_time"`. `"auto"` uses:
#'   group*time if both factors are available, otherwise group if present,
#'   otherwise time.
#' @param connect_centroids How to connect centroids with lines. One of
#'   `"none"` (default), `"group"`, `"time"`. The logic is:
#'   - if centroids are by `group_time` and `connect_centroids = "group"`,
#'     centroids of each group are connected along time levels;
#'   - if centroids are by `group_time` and `connect_centroids = "time"`,
#'     centroids of each time level are connected across groups;
#'   - for other `centroid_by` values, connections are ignored.
#' @param show_ellipses Logical; if `TRUE` (default), draws ellipses for
#'   selected factors.
#' @param ellipse_by Character vector describing which ellipses to draw.
#'   Any combination of `"group"`, `"time"`, `"group_time"`, or `"none"`
#'   (to disable). Default is `"group"`.
#' @param ellipse_type Ellipse type passed to [ggplot2::stat_ellipse()],
#'   one of `"t"`, `"norm"`, `"euclid"`. Default `"t"`.
#' @param ellipse_level Numeric confidence level for ellipses (default `0.95`).
#' @param point_size Numeric size of sample points. Default `1.5`.
#' @param point_alpha Numeric opacity of sample points in `[0, 1]`.
#'   Default `0.6`.
#' @param centroid_size Numeric size for centroid points. Default `3`.
#'
#' @details
#' Axis labels are annotated with the percentage of constrained variance
#' explained by each CAP axis, computed as
#' \eqn{100 * \lambda_k / \sum_j \lambda_j}, where \eqn{\lambda_k} are the
#' constrained eigenvalues from `cap_res$eigenvalues`. Only positive parts
#' of eigenvalues are used when computing percentages.
#'
#' Styled with [theme_phip()].

#'
#' Typically you will want to join sample-level metadata into
#' `cap_res$sample_coords` before plotting, e.g.:
#'
#' ```r
#' cap_res$sample_coords <- dplyr::left_join(
#'   cap_res$sample_coords,
#'   meta,             # data frame with sample_id, group, time, etc.
#'   by = "sample_id"
#' )
#' ```
#'
#' @return A [ggplot2::ggplot] object representing the CAP ordination.
#'
#' @examples
#' \dontrun{
#'   # cap_res <- compute_capscale(dist_bc, ps = ps, formula = ~ type_person + age)
#'   cap_res$sample_coords <- dplyr::left_join(
#'     cap_res$sample_coords,
#'     meta,  # contains type_person, timepoint, etc.
#'     by = "sample_id"
#'   )
#'
#'   plot_cap(
#'     cap_res,
#'     axes            = c(1, 2),
#'     group_col       = "type_person",
#'     time_col        = "timepoint",
#'     show_centroids  = TRUE,
#'     centroid_by     = "group_time",
#'     connect_centroids = "group",
#'     show_ellipses   = TRUE,
#'     ellipse_by      = c("group", "group_time")
#'   )
#' }
#' @export
plot_cap <- function(cap_res,
                     axes = c(1, 2),
                     group_col = NULL,
                     time_col = NULL,
                     show_centroids = TRUE,
                     centroid_by = c("auto", "group", "time", "group_time"),
                     connect_centroids = c("none", "group", "time"),
                     show_ellipses = TRUE,
                     ellipse_by = c("group"),
                     ellipse_type = c("t", "norm", "euclid"),
                     ellipse_level = 0.95,
                     point_size = 1.5,
                     point_alpha = 0.6,
                     centroid_size = 3) {
  # --- small helpers ---------------------------------------------------------
  .abort <- function(msg) {
    if (exists(".ph_abort", mode = "function")) {
      .ph_abort(msg, step = "plot_cap")
    } else {
      stop(msg, call. = FALSE)
    }
  }
  .warn <- function(msg) {
    if (exists(".ph_warn", mode = "function")) {
      .ph_warn(msg, step = "plot_cap")
    } else {
      warning(msg, call. = FALSE)
    }
  }
  .info <- function(msg) {
    if (exists(".ph_log_info", mode = "function")) {
      .ph_log_info(msg, step = "plot_cap")
    }
  }
  .phip_pal <- function(n) {
    if (exists("phip_palette", inherits = TRUE)) {
      pal <- get("phip_palette", inherits = TRUE)
      if (length(pal) < n) rep(pal, length.out = n) else pal[seq_len(n)]
    } else {
      c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00",
        "#6a3d9a", "#b15928") |> rep(length.out = n)
    }
  }
  `%||%` <- function(a, b) if (is.null(a)) b else a

  centroid_by       <- match.arg(centroid_by)
  connect_centroids <- match.arg(connect_centroids)
  ellipse_type      <- match.arg(ellipse_type)
  ellipse_by <- match.arg(
    ellipse_by,
    c("none", "group", "time", "group_time"),
    several.ok = TRUE
  )
  if ("none" %in% ellipse_by) ellipse_by <- character(0)

  # --- basic checks ----------------------------------------------------------
  if (!is.list(cap_res) || !("sample_coords" %in% names(cap_res))) {
    .abort("`cap_res` must be a list containing `sample_coords`.")
  }
  df <- tibble::as_tibble(cap_res$sample_coords)

  if (length(axes) != 2L || any(!is.finite(axes))) {
    .abort("`axes` must be an integer vector of length 2, e.g. c(1, 2).")
  }
  axes <- as.integer(axes)
  ax_names <- paste0("CAP", axes)

  missing_ax <- setdiff(ax_names, names(df))
  if (length(missing_ax) > 0L) {
    .abort(
      paste0(
        "Requested axes not found in `sample_coords`: ",
        paste(missing_ax, collapse = ", ")
      )
    )
  }

  # --- variance explained for labels ----------------------------------------
  eig <- cap_res$eigenvalues %||% numeric(0L)
  pct1 <- pct2 <- NA_real_
  if (length(eig) >= max(axes)) {
    pos_eig <- pmax(eig, 0)
    total   <- sum(pos_eig, na.rm = TRUE)
    if (total > 0) {
      pct1 <- 100 * pos_eig[axes[1]] / total
      pct2 <- 100 * pos_eig[axes[2]] / total
    }
  }

  xlab <- if (is.finite(pct1)) sprintf("CAP%d (%.1f%%)", axes[1], pct1) else ax_names[1]
  ylab <- if (is.finite(pct2)) sprintf("CAP%d (%.1f%%)", axes[2], pct2) else ax_names[2]

  # --- auto-detect group / time ---------------------------------------------
  if (is.null(group_col) && "group" %in% names(df)) {
    group_col <- "group"
  }
  if (is.null(time_col) && "time" %in% names(df)) {
    time_col <- "time"
  }

  has_group <- !is.null(group_col) && group_col %in% names(df)
  has_time  <- !is.null(time_col)  && time_col  %in% names(df)

  if (has_group) df[[group_col]] <- droplevels(as.factor(df[[group_col]]))
  if (has_time)  df[[time_col]]  <- droplevels(as.factor(df[[time_col]]))

  group_levels <- if (has_group) levels(df[[group_col]]) else character(0)
  time_levels  <- if (has_time)  levels(df[[time_col]])  else character(0)

  .info(sprintf(
    "CAP plot: n=%d samples | groups=%d | times=%d",
    nrow(df), length(group_levels), length(time_levels)
  ))

  # palettes
  group_cols <- if (has_group) {
    cols <- .phip_pal(length(group_levels))
    stats::setNames(cols, group_levels)
  } else NULL

  time_cols <- if (has_time) {
    cols <- .phip_pal(length(time_levels))
    stats::setNames(cols, time_levels)
  } else NULL

  # --- base ggplot: samples --------------------------------------------------
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[ax_names[1]]],
                 y = .data[[ax_names[2]]])
  )

  if (has_group && has_time) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[group_col]],
                     shape = .data[[time_col]]),
        size = point_size,
        alpha = point_alpha
      ) +
      ggplot2::scale_color_manual(values = group_cols, name = group_col)
  } else if (has_group) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[group_col]]),
        size = point_size,
        alpha = point_alpha
      ) +
      ggplot2::scale_color_manual(values = group_cols, name = group_col)
  } else if (has_time) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[time_col]]),
        size = point_size,
        alpha = point_alpha
      ) +
      ggplot2::scale_color_manual(values = time_cols, name = time_col)
  } else {
    p <- p +
      ggplot2::geom_point(
        size = point_size,
        alpha = point_alpha,
        color = "grey40"
      )
  }

  # --- ellipses --------------------------------------------------------------
  if (isTRUE(show_ellipses) && length(ellipse_by) > 0L) {
    # group ellipses
    if ("group" %in% ellipse_by && has_group) {
      for (g in group_levels) {
        dsub <- df[df[[group_col]] == g, , drop = FALSE]
        if (nrow(dsub) < 3) next
        p <- p +
          ggplot2::stat_ellipse(
            data = dsub,
            ggplot2::aes(
              x = .data[[ax_names[1]]],
              y = .data[[ax_names[2]]]
            ),
            type = ellipse_type,
            level = ellipse_level,
            linewidth = 0.7,
            alpha = 0.6,
            color = group_cols[[g]],
            show.legend = FALSE
          )
      }
    }
    # time ellipses
    if ("time" %in% ellipse_by && has_time) {
      for (tlev in time_levels) {
        dsub <- df[df[[time_col]] == tlev, , drop = FALSE]
        if (nrow(dsub) < 3) next
        p <- p +
          ggplot2::stat_ellipse(
            data = dsub,
            ggplot2::aes(
              x = .data[[ax_names[1]]],
              y = .data[[ax_names[2]]]
            ),
            type = ellipse_type,
            level = ellipse_level,
            linewidth = 0.6,
            alpha = 0.6,
            color = time_cols[[tlev]],
            show.legend = FALSE
          )
      }
    }
    # group*time ellipses
    if ("group_time" %in% ellipse_by && has_group && has_time) {
      for (g in group_levels) {
        for (tlev in time_levels) {
          dsub <- df[df[[group_col]] == g & df[[time_col]] == tlev, , drop = FALSE]
          if (nrow(dsub) < 3) next
          p <- p +
            ggplot2::stat_ellipse(
              data = dsub,
              ggplot2::aes(
                x = .data[[ax_names[1]]],
                y = .data[[ax_names[2]]]
              ),
              type = ellipse_type,
              level = ellipse_level,
              linewidth = 0.5,
              alpha = 0.7,
              color = group_cols[[g]],
              show.legend = FALSE
            )
        }
      }
    }
  }

  # --- centroids -------------------------------------------------------------
  cent <- NULL

  if (isTRUE(show_centroids)) {
    if (centroid_by == "auto") {
      centroid_by <- if (has_group && has_time) {
        "group_time"
      } else if (has_group) {
        "group"
      } else if (has_time) {
        "time"
      } else {
        "none"
      }
    }

    if (centroid_by == "group" && has_group) {
      cent <- df |>
        dplyr::group_by(.data[[group_col]]) |>
        dplyr::summarise(
          cx = mean(.data[[ax_names[1]]], na.rm = TRUE),
          cy = mean(.data[[ax_names[2]]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(group_lab = !!rlang::sym(group_col))
      cent$group_lab <- factor(cent$group_lab, levels = group_levels)
    } else if (centroid_by == "time" && has_time) {
      cent <- df |>
        dplyr::group_by(.data[[time_col]]) |>
        dplyr::summarise(
          cx = mean(.data[[ax_names[1]]], na.rm = TRUE),
          cy = mean(.data[[ax_names[2]]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(time_lab = !!rlang::sym(time_col))
      cent$time_lab <- factor(cent$time_lab, levels = time_levels)
    } else if (centroid_by == "group_time" && has_group && has_time) {
      cent <- df |>
        dplyr::group_by(.data[[group_col]], .data[[time_col]]) |>
        dplyr::summarise(
          cx = mean(.data[[ax_names[1]]], na.rm = TRUE),
          cy = mean(.data[[ax_names[2]]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::rename(
          group_lab = !!rlang::sym(group_col),
          time_lab  = !!rlang::sym(time_col)
        )
      cent$group_lab <- factor(cent$group_lab, levels = group_levels)
      cent$time_lab  <- factor(cent$time_lab,  levels = time_levels)
      cent$time_ord  <- as.integer(cent$time_lab)
    }

    if (!is.null(cent) && nrow(cent) > 0L) {
      # plot centroids
      if ("group_lab" %in% names(cent) && "time_lab" %in% names(cent)) {
        p <- p +
          ggplot2::geom_point(
            data = cent,
            ggplot2::aes(
              x = .data$cx,
              y = .data$cy,
              color = .data$group_lab
            ),
            size = centroid_size
          )
        # connect centroids if requested
        if (connect_centroids == "group") {
          for (g in group_levels) {
            cg <- cent[cent$group_lab == g, , drop = FALSE]
            if (nrow(cg) > 1L) {
              cg <- cg[order(cg$time_ord), , drop = FALSE]
              p <- p +
                ggplot2::geom_path(
                  data = cg,
                  ggplot2::aes(x = .data$cx, y = .data$cy),
                  color = group_cols[[g]],
                  linewidth = 0.7,
                  alpha = 0.7,
                  show.legend = FALSE
                )
            }
          }
        } else if (connect_centroids == "time") {
          for (tlev in time_levels) {
            ct <- cent[cent$time_lab == tlev, , drop = FALSE]
            if (nrow(ct) > 1L) {
              p <- p +
                ggplot2::geom_path(
                  data = ct,
                  ggplot2::aes(x = .data$cx, y = .data$cy),
                  color = "grey50",
                  linewidth = 0.6,
                  alpha = 0.6,
                  show.legend = FALSE
                )
            }
          }
        }
      } else if ("group_lab" %in% names(cent)) {
        p <- p +
          ggplot2::geom_point(
            data = cent,
            ggplot2::aes(
              x = .data$cx,
              y = .data$cy,
              color = .data$group_lab
            ),
            size = centroid_size
          )
      } else if ("time_lab" %in% names(cent)) {
        p <- p +
          ggplot2::geom_point(
            data = cent,
            ggplot2::aes(
              x = .data$cx,
              y = .data$cy,
              color = .data$time_lab
            ),
            size = centroid_size
          )
      }
    }
  }

  # --- finish ----------------------------------------------------------------
  base_theme <- theme_phip()

  p <- p +
    ggplot2::labs(x = xlab, y = ylab) +
    base_theme +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  p
}
#' @title Scree Plot for PCoA Eigenvalues
#' @description
#' Creates a scree plot from a `"beta_pcoa"` object (output of
#' [compute_pcoa()]), showing the percentage of variance explained
#' for the first `n_axes` axes.
#'
#' @param pcoa_res A `"beta_pcoa"` object as returned by [compute_pcoa()].
#'   Must contain an `eigenvalues` component (numeric vector of eigenvalues).
#' @param n_axes Integer giving the number of axes to display. If `NULL`
#'   (default), the function uses `min(10, length(eigenvalues))`.
#' @param type Type of scree plot to draw: `"bar"` (default) for a bar plot
#'   of variance explained, or `"line"` for a line/point plot.
#'
#' @details
#' Percentages are computed from the positive part of the eigenvalues,
#' i.e. `pmax(eigenvalues, 0)`, to handle possible small negative
#' eigenvalues from non-perfectly Euclidean distance matrices.
#'
#' Styled with [theme_phip()].

#'
#' @return A [ggplot2::ggplot] object representing the scree plot.
#'
#' @examples
#' \dontrun{
#'   # pcoa_res <- compute_pcoa(dist_bc)
#'   plot_scree(pcoa_res)
#'
#'   # More axes as line plot:
#'   plot_scree(pcoa_res, n_axes = 15, type = "line")
#' }
#' @export
plot_scree <- function(pcoa_res,
                       n_axes = NULL,
                       type = c("bar", "line")) {
  # --- helpers ---------------------------------------------------------------
  .abort <- function(msg) {
    if (exists(".ph_abort", mode = "function")) {
      .ph_abort(msg, step = "plot_scree")
    } else {
      stop(msg, call. = FALSE)
    }
  }

  type <- match.arg(type)

  # --- checks ----------------------------------------------------------------
  if (!is.list(pcoa_res) || !("eigenvalues" %in% names(pcoa_res))) {
    .abort("`pcoa_res` must be a list containing an `eigenvalues` component.")
  }

  eig <- as.numeric(pcoa_res$eigenvalues)
  if (length(eig) == 0L) {
    .abort("`pcoa_res$eigenvalues` is empty; cannot produce scree plot.")
  }

  pos_eig <- pmax(eig, 0)
  tot <- sum(pos_eig, na.rm = TRUE)
  k_max <- length(eig)

  if (is.null(n_axes)) {
    n_axes <- min(10L, k_max)
  } else {
    n_axes <- max(1L, min(as.integer(n_axes), k_max))
  }

  eig_use <- eig[seq_len(n_axes)]
  pos_use <- pmax(eig_use, 0)

  pct <- if (tot > 0) 100 * pos_use / tot else rep(NA_real_, n_axes)

  dat <- tibble::tibble(
    axis     = seq_len(n_axes),
    eigenval = eig_use,
    pct      = pct
  )

  base_theme <- theme_phip()

  # --- ggplot ----------------------------------------------------------------
  if (type == "bar") {
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$axis, y = .data$pct)) +
      ggplot2::geom_col(fill = "grey60")
  } else {
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$axis, y = .data$pct)) +
      ggplot2::geom_line() +
      ggplot2::geom_point()
  }

  p +
    ggplot2::labs(
      x = "PCoA axis",
      y = "Variance explained [%]"
    ) +
    base_theme
}


#' @title Plot beta dispersion (distance to centroid)
#' @description
#' Plot distances to group/time centroids from a \code{beta_dispersion} object
#' (output of \code{compute_dispersion()}) as violins, hollow boxplots and/or
#' jittered points, by level on the x-axis.
#'
#' @param x A \code{beta_dispersion} object returned by
#'   \code{compute_dispersion()}.
#' @param scope Character scalar indicating which dispersion scope to plot.
#'   Typically one of \code{"group"}, \code{"time"} or \code{"group:time"}.
#'   Must match values in \code{x$distances$scope}. Default \code{"group"}.
#' @param contrast Character scalar indicating which contrast within
#'   \code{scope} to plot, e.g. \code{"<global>"} or \code{"A vs B"}.
#'   Must match values in \code{x$distances$contrast}. Default \code{"<global>"}.
#' @param show_violin Logical; if \code{TRUE} (default) draw a violin layer.
#' @param show_box Logical; if \code{TRUE} (default) draw hollow boxplots
#'   (\code{fill = NA}) per level.
#' @param show_points Logical; if \code{TRUE} (default) draw jittered points.
#' @param point_size Numeric point size for jitter layer. Default 1.5.
#' @param point_alpha Numeric alpha for jitter layer. Default 0.6.
#' @param width_violin Width of violin layer. Default 0.9.
#' @param width_box Width of boxplot layer. Default 0.5.
#' @param jitter_width,jitter_height Jitter width/height for points.
#'   Defaults: \code{0.1} and \code{0}.
#' @param add_pvalue Logical; if \code{TRUE} (default) and a matching row is
#'   found in \code{x$tests}, the p-value is shown in the subtitle.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#'   disp_res <- compute_dispersion(dist_bc, ps, group_col = "group_char")
#'
#'   # Global group dispersion:
#'   plot_dispersion(disp_res, scope = "group", contrast = "<global>")
#'
#'   # Pairwise contrast, with only boxplots + points:
#'   plot_dispersion(
#'     disp_res,
#'     scope        = "group",
#'     contrast     = "dementia vs control",
#'     show_violin  = FALSE,
#'     show_box     = TRUE,
#'     show_points  = TRUE
#'   )
#' }
#' @export
plot_dispersion <- function(x,
                            scope          = "group",
                            contrast       = "<global>",
                            show_violin    = TRUE,
                            show_box       = TRUE,
                            show_points    = TRUE,
                            point_size     = 1.5,
                            point_alpha    = 0.6,
                            width_violin   = 0.9,
                            width_box      = 0.5,
                            jitter_width   = 0.1,
                            jitter_height  = 0,
                            add_pvalue     = TRUE) {
  if (!inherits(x, "beta_dispersion")) {
    .ph_abort("`x` must be a <beta_dispersion> object (from `compute_dispersion()`).",
              step = "plot_dispersion")
  }

  dist_df <- x$distances
  tests   <- x$tests %||% tibble::tibble()

  if (!all(c("sample_id", "distance", "level", "scope", "contrast") %in% names(dist_df))) {
    .ph_abort("`x$distances` does not have the expected columns.",
              step = "plot_dispersion")
  }

  # Filter to requested scope/contrast ----------------------------------------
  dist_sub <- dist_df |>
    dplyr::filter(.data$scope == !!scope, .data$contrast == !!contrast)

  if (!nrow(dist_sub)) {
    .ph_abort(
      paste0(
        "No rows in `x$distances` for scope = '", scope,
        "' and contrast = '", contrast, "'."
      ),
      step = "plot_dispersion"
    )
  }

  dist_sub$level <- factor(dist_sub$level)

  .ph_log_info(
    paste0(
      "Plotting dispersion for scope = '", scope,
      "', contrast = '", contrast, "' (n = ", nrow(dist_sub), ")."
    ),
    step = "plot_dispersion"
  )

  # Optional p-value for subtitle ---------------------------------------------
  subtitle <- NULL
  if (isTRUE(add_pvalue) && nrow(tests)) {
    test_row <- tests |>
      dplyr::filter(.data$scope == !!scope, .data$contrast == !!contrast) |>
      dplyr::slice(1)
    if (nrow(test_row) == 1L && is.finite(test_row$p_value[1])) {
      subtitle <- sprintf(
        "%s (%s): p = %.3g (n_perm = %d)",
        scope, contrast,
        test_row$p_value[1],
        test_row$n_perm[1] %||% NA_integer_
      )
    }
  }

  # Base ggplot ---------------------------------------------------------------
  p <- ggplot2::ggplot(dist_sub, ggplot2::aes(x = .data$level, y = .data$distance))

  if (isTRUE(show_violin)) {
    p <- p +
      ggplot2::geom_violin(
        ggplot2::aes(fill = .data$level),
        width = width_violin,
        alpha = 0.25,
        color = NA,
        trim  = FALSE
      )
  }

  if (isTRUE(show_box)) {
    p <- p +
      ggplot2::geom_boxplot(
        ggplot2::aes(color = .data$level),
        width = width_box,
        fill  = NA,
        outlier.shape = NA
      )
  }

  if (isTRUE(show_points)) {
    p <- p +
      ggplot2::geom_jitter(
        ggplot2::aes(color = .data$level),
        width  = jitter_width,
        height = jitter_height,
        size   = point_size,
        alpha  = point_alpha
      )
  }

  p <- p +
    ggplot2::labs(
      x        = "Level",
      y        = "Distance to centroid",
      title    = paste("Beta-dispersion:", scope, "-", contrast),
      subtitle = subtitle,
      color    = "Level",
      fill     = "Level"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
    )

  p <- p + theme_phip()

  p
}

#' Plot t-SNE embeddings
#'
#' @description
#' Plot t-SNE embeddings computed by [compute_tsne()]. Supports both
#' 2D (ggplot2) and 3D (plotly) views.
#'
#' @param tsne_res A `phip_tsne` object returned by [compute_tsne()].
#'   The function also works with a plain tibble/data frame that contains
#'   at least columns `tSNE1` and `tSNE2` (and `tSNE3` for 3D).
#' @param view Character, either `"2d"` or `"3d"`. `"2d"` returns a
#'   `ggplot` object; `"3d"` returns a `plotly` HTML widget.
#' @param colour Optional name of a column in `tsne_res` to map to point
#'   colour. If `NULL`, the function uses the first metadata column stored
#'   in `attr(tsne_res, "meta_cols")`, if available.
#' @param size Numeric point size for scatter plots. Defaults to `1.5`.
#' @param alpha Numeric transparency for points (0-1). Defaults to `0.8`.
#' @param palette Optional vector of colour values passed to
#'   `scale_color_manual()` (2D) or `colors` (3D plotly).
#' @param ... Currently ignored; reserved for future extensions.
#'
#' @return
#' For `view = "2d"`, a `ggplot` object.
#' For `view = "3d"`, a `plotly` object (`htmlwidget`).
#'
#' @examples
#' \dontrun{
#' tsne_res <- compute_tsne(ps, compute_distance(ps))
#'
#' # 2D plot
#' p2d <- plot_tsne(tsne_res, view = "2d", colour = "type_person")
#'
#' # 3D interactive plot
#' p3d <- plot_tsne(tsne_res, view = "3d", colour = "type_person")
#' }
#'
#' @export
plot_tsne <- function(tsne_res,
                      view = c("2d", "3d"),
                      colour = NULL,
                      size = 1.5,
                      alpha = 0.8,
                      palette = NULL,
                      ...) {
  view <- match.arg(view)

  # Basic sanity checks -------------------------------------------------------
  if (!("tSNE1" %in% names(tsne_res)) ||
      !("tSNE2" %in% names(tsne_res))) {
    .ph_abort(
      "Input `tsne_res` must contain columns `tSNE1` and `tSNE2`.",
      step = "plot_tsne"
    )
  }

  # Default colour column: first attached metadata column (if any)
  if (is.null(colour)) {
    meta_cols_attr <- attr(tsne_res, "meta_cols")
    if (!is.null(meta_cols_attr) && length(meta_cols_attr) > 0L) {
      colour <- meta_cols_attr[[1L]]
    }
  }

  # 2D plot -------------------------------------------------------------------
  if (identical(view, "2d")) {
    .ph_log_info("Creating 2D t-SNE plot (ggplot2).",
                 step = "plot_tsne")

    df <- tsne_res

    # Build ggplot with optional colour mapping
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$tSNE1, y = .data$tSNE2))

    if (!is.null(colour) && colour %in% names(df)) {
      col_sym <- rlang::sym(colour)
      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(color = !!col_sym),
          size = size,
          alpha = alpha
        )
    } else {
      p <- p + ggplot2::geom_point(size = size, alpha = alpha)
    }

    # Palette (if given and we actually mapped colour)
    if (!is.null(palette) && !is.null(colour) && colour %in% names(df)) {
      p <- p + ggplot2::scale_color_manual(values = palette)
    }

    p <- p +
      ggplot2::labs(
        x = "t-SNE 1",
        y = "t-SNE 2",
        colour = if (!is.null(colour)) colour else NULL
      )

    p <- p + theme_phip()

    return(p)
  }

  # 3D plot -------------------------------------------------------------------
  if (!("tSNE3" %in% names(tsne_res))) {
    .ph_abort(
      "3D view requested but column `tSNE3` is missing. ",
      "Re-run `compute_tsne()` with `dims = 3L`.",
      step = "plot_tsne"
    )
  }

  if (!rlang::is_installed("plotly")) {
    .ph_abort(
      "Package 'plotly' is required for 3D t-SNE plots. Please install it.",
      step = "plot_tsne"
    )
  }

  .ph_log_info("Creating 3D t-SNE plot (plotly).",
               step = "plot_tsne")

  df <- tsne_res

  # Dynamic formula for colour mapping (if requested and available)
  color_formula <- NULL
  if (!is.null(colour) && colour %in% names(df)) {
    color_formula <- stats::as.formula(paste0("~`", colour, "`"))
  }

  # Build 3D scatter plot
  p3d <- plotly::plot_ly(
    data  = df,
    x     = ~tSNE1,
    y     = ~tSNE2,
    z     = ~tSNE3,
    type  = "scatter3d",
    mode  = "markers",
    color = color_formula,
    colors = palette,
    marker = list(size = size * 2, opacity = alpha)
  ) |>
    plotly::layout(
      scene = list(
        xaxis = list(title = "t-SNE 1"),
        yaxis = list(title = "t-SNE 2"),
        zaxis = list(title = "t-SNE 3"),
        aspectmode = "data"
      ),
      legend = list(title = list(text = if (!is.null(colour)) colour else ""))
    )

  p3d
}
