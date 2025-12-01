# ==============================================================================
# plot_enrichment_counts() — single exported function for <phip_data>
# ==============================================================================

#' Plot enrichment counts per group (and optional interaction)
#'
#' @description
#' Visualizes per-sample peptide enrichment counts across groups in a
#' `<phip_data>` object, with optional interaction of multiple grouping
#' variables. Presence is defined by `exist > 0`. The plot(s) are produced by
#' an internal helper `.plot_enrichment_counts_one()`.
#'
#' @details
#' - If `group_cols = NULL`, a single plot is returned for all samples.
#' - If `group_cols` is a character vector, a list of plots is returned (one per
#'   grouping column) unless `interaction_only = TRUE`.
#' - If `group_interaction = TRUE` and at least two `group_cols` are supplied,
#'   an additional interaction plot is created whose label joins groups using
#'   `interaction_sep`.
#' - If `interaction_only = TRUE`, only the interaction plot is returned (this
#'   requires `group_interaction = TRUE` and at least two `group_cols`).
#'
#' @param phip_data A `<phip_data>` object with `data_long` containing at least
#'   `sample_id`, `peptide_id`, and `exist`.
#' @param group_cols Character vector of grouping columns in `data_long`, or
#'   `NULL` to plot all samples together.
#' @param prevalence_threshold Numeric in `[0,1]`; minimum prevalence used by
#'   `.plot_enrichment_counts_one()` to filter/annotate bins (default `0.05`).
#' @param custom_colors Optional named vector for group colors passed through to
#'   `.plot_enrichment_counts_one()` (default `NULL`).
#' @param binwidth Numeric bin width for histograms (default `1`).
#' @param group_interaction Logical; also compute a plot for the interaction of
#'   all `group_cols` (default `FALSE`).
#' @param interaction_only Logical; if `TRUE`, return only the interaction plot
#'   (requires `group_interaction = TRUE` and at least two `group_cols`).
#' @param interaction_sep Character separator for interaction labels (default `" * "`).
#' @param annotation_size Numeric; size of the in-plot threshold annotations
#'   (passed to `geom_text(size = ...)`). Typical range 3–6. Default `4`.
#' @param ... Reserved for future extensions; ignored.
#'
#' @return
#' - A single plot object (when `group_cols = NULL`, or when only one plot is
#'   produced), or
#' - A named list of plot objects (when multiple plots are produced).
#'
#' @examples
#' # per-group plots
#' pd <- phip_load_example_data()
#' p <- plot_enrichment_counts(pd, group_cols = c("group","timepoint"))
#'
#' # add interaction plot
#' p2 <- plot_enrichment_counts(pd,
#'   group_cols = c("group","timepoint"),
#'   group_interaction = TRUE
#' )
#'
#' # interaction only
#' p3 <- plot_enrichment_counts(pd,
#'   group_cols = c("group","timepoint"),
#'   group_interaction = TRUE,
#'   interaction_only = TRUE
#' )
#' @export
plot_enrichment_counts <- function(phip_data,
                                   group_cols = NULL,
                                   prevalence_threshold = 0.05,
                                   custom_colors = NULL,
                                   binwidth = 1,
                                   group_interaction = FALSE,
                                   interaction_only = FALSE,
                                   interaction_sep = " * ",
                                   annotation_size = 4,
                                   ...) {
  # -- basic checks -------------------------------------------------------------
  x <- phip_data
  stopifnot(inherits(x, "phip_data"))
  .data <- rlang::.data

  # -- enforce interaction_only preconditions ----------------------------------
  if (isTRUE(interaction_only)) {
    if (!isTRUE(group_interaction) || is.null(group_cols) || length(group_cols) < 2L) {
      .ph_abort(
        headline = "`interaction_only = TRUE` requires an interaction.",
        step = "argument validation",
        bullets = c(
          "set group_interaction = TRUE",
          "provide at least two grouping columns in group_cols"
        )
      )
    }
  }

  .ph_with_timing(
    headline = "Plotting enrichment counts (<phip_data>)",
    step = if (is.null(group_cols)) {
      "group_cols: <none>"
    } else {
      sprintf("group_cols: %s", paste(add_quotes(group_cols, 1L), collapse = ", "))
    },
    expr = {
      tbl <- x$data_long

      # -- prune full-cross rows to speed up plotting --------------------------
      if (isTRUE(x$meta$full_cross) && ("exist" %in% colnames(tbl))) {
        red_txt <- tryCatch({
          ep <- as.numeric(x$meta$exist_prop)
          if (is.finite(ep) && ep > 0) sprintf("~%.1fx", 1 / ep) else "<unknown>"
        }, error = function(e) "<unknown>")
        .ph_log_info(
          "Full-cross detected; pruning non-existing rows before plotting",
          bullets = c("rule: keep exist == 1", sprintf("estimated reduction: %s", red_txt))
        )
        tbl <- dplyr::filter(tbl, .data$exist == 1L)
      }

      # -- validate requested grouping columns ---------------------------------
      if (!is.null(group_cols) && length(group_cols)) {
        missing_gcs <- setdiff(group_cols, colnames(tbl))
        if (length(missing_gcs)) {
          .ph_abort(
            headline = "Missing grouping columns in data_long.",
            step = "input validation",
            bullets = sprintf("missing: %s", paste(add_quotes(missing_gcs, 1L), collapse = ", "))
          )
        }
      }

      # -- branch: no groups -> single plot ------------------------------------
      if (is.null(group_cols) || !length(group_cols)) {
        return(
          .plot_enrichment_counts_one(
            tbl,
            group_col = NULL,
            prevalence_threshold = prevalence_threshold,
            custom_colors = custom_colors,
            binwidth = binwidth
          )
        )
      }

      # -- groups present -------------------------------------------------------
      # interaction-only shortcut: compute and return just that plot
      if (isTRUE(interaction_only)) {
        combo_nm  <- paste(group_cols, collapse = interaction_sep)
        inter_col <- "..phip_interaction.."
        tbl_inter <- dplyr::mutate(
          tbl, !!rlang::sym(inter_col) := paste(!!!rlang::syms(group_cols), sep = interaction_sep)
        )
        return(
          .plot_enrichment_counts_one(
            tbl_inter,
            group_col = inter_col,
            prevalence_threshold = prevalence_threshold,
            custom_colors = custom_colors,
            binwidth = binwidth,
            title_label = combo_nm,
            annotation_size = annotation_size
          )
        )
      }

      # normal mode: per-group plots (+ optional interaction)
      plots <- lapply(group_cols, function(gc) {
        .plot_enrichment_counts_one(
          tbl,
          group_col = gc,
          prevalence_threshold = prevalence_threshold,
          custom_colors = custom_colors,
          binwidth = binwidth,
          title_label = gc,
          annotation_size = annotation_size
        )
      })
      names(plots) <- group_cols

      if (isTRUE(group_interaction) && length(group_cols) >= 2L) {
        combo_nm  <- paste(group_cols, collapse = interaction_sep)
        inter_col <- "..phip_interaction.."
        tbl_inter <- dplyr::mutate(
          tbl, !!rlang::sym(inter_col) := paste(!!!rlang::syms(group_cols), sep = interaction_sep)
        )
        plots[[combo_nm]] <- .plot_enrichment_counts_one(
          tbl_inter,
          group_col = inter_col,
          prevalence_threshold = prevalence_threshold,
          custom_colors = custom_colors,
          binwidth = binwidth,
          title_label = combo_nm,
          annotation_size = annotation_size
        )
      } else if (isTRUE(group_interaction) && length(group_cols) < 2L) {
        .ph_warn(
          headline = "group_interaction requested but fewer than 2 group_cols supplied.",
          step     = "interaction plot",
          bullets  = "interaction plot skipped."
        )
      }

      if (length(plots) == 1L) plots[[1L]] else plots
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ==============================================================================
# plot_alpha_diversity() — single exported function for precomputed alpha data
# (output of compute_alpha_diversity)
# ==============================================================================

#' Plot alpha diversity (richness/Shannon/Simpson) from precomputed results
#'
#' @description
#' Plot alpha diversity metrics from the precomputed output of
#' [compute_alpha_diversity()]. Supports filtering groups/ranks and optional
#' faceting by rank. If the input contains an interaction table, you can request
#' only that via `interaction_only = TRUE`.
#'
#' @details
#' - `x` can be:
#'   - the named list returned by [compute_alpha_diversity()] (class
#'     `"phip_alpha_diversity"`), or
#'   - a single data frame taken from that list.
#' - when `interaction_only = TRUE`, the function tries to select the element
#'   named by `paste(attr(x, "group_cols"), collapse = interaction_sep)`. if not
#'   found, it falls back to the first element whose name contains the separator.
#'
#' @param x a `"phip_alpha_diversity"` list (output of
#'   [compute_alpha_diversity()]) or a single alpha-diversity data frame.
#' @param metric one of `"richness"`, `"shannon_diversity"`, `"simpson_diversity"`.
#' @param group_col name of the grouping column in the alpha table
#'   (default `"group"` when `group_cols = NULL` in the computation step).
#' @param rank_col name of the rank column (default `"rank"`).
#' @param filter_groups optional character vector; keep only these group levels.
#' @param filter_ranks optional character vector; keep only these ranks.
#' @param custom_colors optional named vector of colors for groups.
#' @param facet_by_rank logical; facet by `rank_col` if multiple ranks present.
#' @param ncol integer; number of columns in facet wrap (default `2`).
#' @param facet_scales `"fixed"`, `"free_x"`, `"free_y"`, or `"free"`.
#' @param interaction_only logical; if `TRUE`, plot only the interaction table
#'   (when available in `x`). useful when `compute_alpha_diversity()` was run
#'   with `group_interaction = TRUE`.
#' @param interaction_sep character; separator used to join interaction labels.
#'   default is taken from `attr(x, "interaction_sep")` if present, otherwise `" * "`.
#'
#' @return a `ggplot` object.
#'
#' @examples
#' \dontrun{
#' # precomputed alpha (list) -> boxplot per group
#' p <- plot_alpha_diversity(alpha_list, metric = "richness", group_col = "Cohort")
#'
#' # only the interaction table (if available)
#' p_int <- plot_alpha_diversity(
#'   alpha_list,
#'   metric = "shannon_diversity",
#'   group_col = "Cohort * timepoint",
#'   interaction_only = TRUE
#' )
#' }
#' @export
plot_alpha_diversity <- function(
    x,
    metric = c("richness", "shannon_diversity", "simpson_diversity"),
    group_col = "group",
    rank_col = "rank",
    filter_groups = NULL,
    filter_ranks  = NULL,
    custom_colors = NULL,
    facet_by_rank = TRUE,
    ncol = 2,
    facet_scales = "fixed",
    interaction_only = FALSE,
    interaction_sep = NULL,
    # --- new UI-like controls (parity with interactive) ---
    jitter_width = 0.25,
    point_size   = 1.8,
    point_alpha  = 0.70,
    text_size    = 12,
    font_family  = "Montserrat",
    show_grids   = TRUE,
    x_order      = NULL,
    x_labels     = NULL,
    y_range      = NULL,
    x_tickangle  = 0
) {
  .data  <- rlang::.data
  metric <- tolower(match.arg(metric))

  .ph_with_timing(
    headline = "plotting alpha diversity (precomputed)",
    step     = sprintf("metric: %s", metric),
    expr = {
      # ---- normalize to a single data frame -----------------------------------
      alpha_df <- NULL
      if (is.list(x) && any(vapply(x, inherits, logical(1), "data.frame"))) {
        if (isTRUE(interaction_only)) {
          interaction_sep <- interaction_sep %||% attr(x, "interaction_sep") %||% " * "
          gc <- attr(x, "group_cols", exact = TRUE)
          expected_name <- if (!is.null(gc) && length(gc) >= 2L) paste(gc, collapse = interaction_sep) else NULL

          pick_name <- NULL
          if (!is.null(expected_name) && !is.null(names(x)) && expected_name %in% names(x)) {
            pick_name <- expected_name
          } else if (!is.null(names(x)) && any(grepl(interaction_sep, names(x), fixed = TRUE))) {
            pick_name <- names(x)[which(grepl(interaction_sep, names(x), fixed = TRUE))[1L]]
          }

          if (!is.null(pick_name)) {
            alpha_df <- tibble::as_tibble(x[[pick_name]])
            .ph_log_info("selected interaction table", bullets = pick_name)
          } else {
            .ph_abort("interaction_only requested but interaction table not found.", step = "input selection")
          }
        } else {
          inner    <- x[vapply(x, inherits, logical(1), "data.frame")]
          alpha_df <- suppressWarnings(dplyr::bind_rows(inner))
        }
      } else {
        alpha_df <- tibble::as_tibble(x)
      }

      # ---- checks --------------------------------------------------------------
      need_metric <- metric
      if (!need_metric %in% names(alpha_df)) {
        .ph_abort("metric column not found in alpha data.", step = "input validation",
                  bullets = sprintf("missing: %s", add_quotes(need_metric, 1L)))
      }
      if (!is.null(group_col) && !group_col %in% names(alpha_df)) {
        .ph_abort("grouping column not found in alpha data.", step = "input validation",
                  bullets = sprintf("group_col: %s", add_quotes(group_col, 1L)))
      }
      if (!is.null(rank_col) && !rank_col %in% names(alpha_df)) {
        .ph_warn("rank column not found; disabling faceting by rank.", step = "input validation",
                 bullets = sprintf("rank_col: %s", add_quotes(rank_col, 1L)))
        rank_col <- NULL
      }

      # ---- filtering -----------------------------------------------------------
      if (!is.null(group_col) && !is.null(filter_groups)) {
        alpha_df <- dplyr::filter(alpha_df, .data[[group_col]] %in% !!filter_groups)
      }
      if (!is.null(rank_col) && !is.null(filter_ranks)) {
        alpha_df <- dplyr::filter(alpha_df, .data[[rank_col]] %in% !!filter_ranks)
      }

      # ---- labels --------------------------------------------------------------
      ylab <- switch(
        need_metric,
        richness          = "Richness",
        shannon_diversity = "Shannon diversity",
        simpson_diversity = "Simpson diversity (1 - \u03A3 p^2)"
      )

      # ---- small theme helper (non-invasive to your theme_phip) ---------------
      .theme_overrides <- ggplot2::theme(
        text = ggplot2::element_text(family = font_family, size = text_size),
        axis.text.x = ggplot2::element_text(angle = x_tickangle, hjust = if (x_tickangle == 0) 0.5 else 1),
        panel.grid.major = if (isTRUE(show_grids)) ggplot2::element_line(colour = "grey85", linewidth = 0.3) else ggplot2::element_blank(),
        panel.grid.minor = if (isTRUE(show_grids)) ggplot2::element_line(colour = "grey92", linewidth = 0.2) else ggplot2::element_blank()
      )

      # ======================= no grouping: single box ==========================
      if (is.null(group_col)) {
        p <- ggplot2::ggplot(alpha_df, ggplot2::aes(y = .data[[need_metric]], x = "")) +
          ggplot2::geom_boxplot(outlier.shape = NA, width = 0.7, linewidth = 0.5, fill = "grey70", colour = "black") +
          ggplot2::geom_jitter(width = jitter_width, height = 0, size = point_size, alpha = point_alpha) +
          ggplot2::labs(x = NULL, y = ylab) +
          theme_phip() + .theme_overrides

        if (!is.null(y_range) && length(y_range) == 2) {
          p <- p + ggplot2::coord_cartesian(ylim = y_range)
        }

        if (isTRUE(facet_by_rank) && !is.null(rank_col) && length(unique(alpha_df[[rank_col]])) > 1) {
          p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", rank_col)), ncol = ncol, scales = facet_scales)
        }
        return(p)
      }

      # ======================= grouped: boxes per group =========================
      gsym <- rlang::sym(group_col)
      msym <- rlang::sym(need_metric)

      # establish stable order of x levels
      gvals <- as.character(alpha_df[[group_col]])
      group_levels <- {
        if (!is.null(x_order)) {
          unique(intersect(x_order, gvals))
        } else if (is.factor(alpha_df[[group_col]])) {
          intersect(levels(alpha_df[[group_col]]), unique(gvals))
        } else {
          sort(unique(gvals))
        }
      }
      alpha_df[[group_col]] <- factor(alpha_df[[group_col]], levels = group_levels)

      # counts and labels
      df_counts <- alpha_df |>
        dplyr::group_by(!!gsym) |>
        dplyr::summarise(sample_count = dplyr::n(), .groups = "drop")

      default_lab <- setNames(
        paste0(as.character(df_counts[[group_col]]), "\n(n = ", df_counts$sample_count, ")"),
        as.character(df_counts[[group_col]])
      )
      if (!is.null(x_labels) && length(x_labels)) {
        # keep order and only replace where provided
        xlab_map <- default_lab
        hits <- intersect(names(x_labels), names(default_lab))
        xlab_map[hits] <- x_labels[hits]
      } else {
        xlab_map <- default_lab
      }

      # build plot
      p <- ggplot2::ggplot(alpha_df, ggplot2::aes(x = !!gsym, y = !!msym, fill = !!gsym, colour = !!gsym)) +
        ggplot2::geom_boxplot(outlier.shape = NA, width = 0.7, linewidth = 0.5, show.legend = FALSE, color = "black") +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = jitter_width, height = 0), shape = 21, color = "black",
                            size = point_size, alpha = point_alpha, show.legend = FALSE) +
        ggplot2::scale_x_discrete(limits = group_levels, labels = xlab_map) +
        ggplot2::labs(x = "Group", y = ylab, fill = group_col, colour = group_col) +
        theme_phip() + .theme_overrides

      # y range (without dropping data)
      if (!is.null(y_range) && length(y_range) == 2) {
        p <- p + ggplot2::coord_cartesian(ylim = y_range)
      }

      # palette
      if (is.null(custom_colors)) {
        p <- p + scale_color_phip(limits = group_levels) + scale_fill_phip(limits = group_levels)
      } else {
        # allow named vector; recycle if unnamed
        if (!is.null(names(custom_colors)) && all(group_levels %in% names(custom_colors))) {
          vals <- unname(custom_colors[group_levels])
        } else {
          vals <- rep(custom_colors, length.out = length(group_levels))
        }
        p <- p + ggplot2::scale_color_manual(values = vals, limits = group_levels) +
          ggplot2::scale_fill_manual(values = vals, limits = group_levels)
      }

      # faceting
      if (isTRUE(facet_by_rank) && !is.null(rank_col) && length(unique(alpha_df[[rank_col]])) > 1) {
        p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", rank_col)),
                                     ncol = ncol, scales = facet_scales)
      }

      p
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

#' Plot alpha diversity (precomputed) — interactive (plotly)
#'
#' @description
#' native plotly version of `plot_alpha_diversity()`. mirrors the ggplot look:
#' per-group boxplots with jittered points and optional faceting by rank.
#'
#' @inheritParams plot_alpha_diversity
#' @param custom_colors named character vector of hex colors for groups (like ggplot scale_fill_manual()).
#' @param x_order optional character vector with desired group order (levels).
#' @param x_labels optional named character vector mapping group -> label (used on x axis).
#' @param y_range optional numeric length-2; y axis range (e.g., c(0, 2300)).
#' @param x_tickangle numeric; tick label rotation in degrees (default 0; e.g., 25).
#' @param quartile_method one of c("exclusive","inclusive","linear"); passed to plotly box (default "exclusive" ~ ggplot).
#' @return A plotly htmlwidget.
#' @examples
#' \dontrun{
#' plot_alpha_diversity_interactive(df, metric = "richness", group_col = "group")
#' }
#' @export
plot_alpha_diversity_interactive <- function(
    x,
    metric = c("richness", "shannon_diversity", "simpson_diversity"),
    group_col = "group",
    rank_col = "rank",
    filter_groups = NULL,
    filter_ranks  = NULL,
    custom_colors = NULL,
    facet_by_rank = TRUE,
    ncol = 2,
    facet_scales = "fixed",
    interaction_only = FALSE,
    interaction_sep = NULL,
    # --- 1:1 parity with static ---
    x_order = NULL,            # user can define the order here
    x_labels = NULL,
    y_range = NULL,
    x_tickangle = 0,
    quartile_method = c("exclusive","inclusive","linear"),
    jitter_width = 0.25,
    point_size   = 6,
    point_alpha  = 0.85,
    text_size    = 12,
    font_family  = "Montserrat",
    show_grids   = TRUE
) {
  .data  <- rlang::.data
  metric <- tolower(match.arg(metric))
  quartile_method <- match.arg(quartile_method)

  .ph_with_timing(
    headline = "plotting alpha diversity (precomputed, interactive)",
    step     = sprintf("metric: %s", metric),
    expr = {
      # ------------ normalize ------------
      alpha_df <- NULL
      if (is.list(x) && any(vapply(x, inherits, logical(1), "data.frame"))) {
        if (isTRUE(interaction_only)) {
          interaction_sep <- interaction_sep %||% attr(x, "interaction_sep") %||% " * "
          gc <- attr(x, "group_cols", exact = TRUE)
          expected_name <- if (!is.null(gc) && length(gc) >= 2L) paste(gc, collapse = interaction_sep) else NULL
          pick_name <- NULL
          if (!is.null(expected_name) && !is.null(names(x)) && expected_name %in% names(x)) {
            pick_name <- expected_name
          } else if (!is.null(names(x)) && any(grepl(interaction_sep, names(x), fixed = TRUE))) {
            pick_name <- names(x)[which(grepl(interaction_sep, names(x), fixed = TRUE))[1L]]
          }
          if (!is.null(pick_name)) {
            alpha_df <- tibble::as_tibble(x[[pick_name]])
            .ph_log_info("selected interaction table", bullets = pick_name)
          } else {
            .ph_abort("interaction_only requested but interaction table not found.", step = "input selection")
          }
        } else {
          inner    <- x[vapply(x, inherits, logical(1), "data.frame")]
          alpha_df <- suppressWarnings(dplyr::bind_rows(inner))
        }
      } else {
        alpha_df <- tibble::as_tibble(x)
      }

      # ------------ checks ------------
      need_metric <- metric
      if (!need_metric %in% names(alpha_df)) {
        .ph_abort("metric column not found in alpha data.", step = "input validation",
                  bullets = sprintf("missing: %s", add_quotes(need_metric, 1L)))
      }
      if (!is.null(group_col) && !group_col %in% names(alpha_df)) {
        .ph_abort("grouping column not found in alpha data.", step = "input validation",
                  bullets = sprintf("group_col: %s", add_quotes(group_col, 1L)))
      }
      if (!is.null(rank_col) && !rank_col %in% names(alpha_df)) {
        .ph_warn("rank column not found; disabling faceting by rank.", step = "input validation",
                 bullets = sprintf("rank_col: %s", add_quotes(rank_col, 1L)))
        rank_col <- NULL
      }

      # ------------ filters ------------
      if (!is.null(group_col) && !is.null(filter_groups)) {
        alpha_df <- dplyr::filter(alpha_df, .data[[group_col]] %in% !!filter_groups)
      }
      if (!is.null(rank_col) && !is.null(filter_ranks)) {
        alpha_df <- dplyr::filter(alpha_df, .data[[rank_col]] %in% !!filter_ranks)
      }

      # ------------ helpers ------------
      ylab <- switch(
        need_metric,
        richness          = "Richness",
        shannon_diversity = "Shannon diversity",
        simpson_diversity = "Simpson diversity (1 - \u03A3 p^2)"
      )
      .as_rgba <- function(col, alpha = 0.6) {
        if (is.null(col)) return(NULL)
        rgb <- grDevices::col2rgb(col) / 255
        sprintf("rgba(%d,%d,%d,%.3f)", round(rgb[1]*255), round(rgb[2]*255), round(rgb[3]*255), alpha)
      }
      .color_map <- function(levels_chr, custom_colors) {
        if (is.null(levels_chr) || !length(levels_chr)) return(NULL)
        if (is.null(custom_colors)) return(NULL)
        if (!is.null(names(custom_colors)) && all(levels_chr %in% names(custom_colors))) {
          unname(custom_colors[levels_chr])
        } else {
          rep(custom_colors, length.out = length(levels_chr))
        }
      }

      .build_panel <- function(df_panel, panel_title = NULL) {
        # safe hover
        subj <- if ("subject_id" %in% names(df_panel)) df_panel$subject_id else rep(NA_character_, nrow(df_panel))
        samp <- if ("sample_id"  %in% names(df_panel)) df_panel$sample_id  else rep(NA, nrow(df_panel))
        shan <- if ("shannon_diversity" %in% names(df_panel)) df_panel$shannon_diversity else rep(NA_real_, nrow(df_panel))

        fmt <- function(x) ifelse(is.na(x), "NA", as.character(x))
        hovertext_vec <- paste0(
          "<b>subject_id</b>: ", fmt(subj),
          "<br><b>sample_id</b>: ", fmt(samp),
          "<br><b>", need_metric, "</b>: ", fmt(df_panel[[need_metric]]),
          "<br><b>shannon_diversity</b>: ", fmt(shan)
        )

        if (is.null(group_col)) {
          base_x <- rep(1, nrow(df_panel))
          set.seed(1L)
          xjit   <- base_x + stats::runif(nrow(df_panel), -jitter_width, jitter_width)

          p <- plotly::plot_ly()
          p <- plotly::add_boxplot(
            p,
            x = base_x, y = df_panel[[need_metric]],
            name = "",
            boxpoints = "outliers",
            marker = list(opacity = 0),
            hoverinfo = "skip",
            line = list(width = 1),
            showlegend = FALSE,
            quartilemethod = quartile_method
          )
          p <- plotly::add_markers(
            p,
            x = xjit, y = df_panel[[need_metric]],
            marker = list(size = point_size, opacity = point_alpha),
            text = hovertext_vec, hovertemplate = "%{text}<extra></extra>",
            showlegend = FALSE
          )
          return(plotly::layout(
            p,
            xaxis = list(
              title = NULL, tickmode = "array", tickvals = 1, ticktext = "",
              range = c(0.5, 1.5), showgrid = show_grids, tickfont = list(size = text_size, family = font_family)
            ),
            yaxis = c(list(title = ylab, showgrid = show_grids, tickfont = list(size = text_size, family = font_family)),
                      if (!is.null(y_range)) list(range = y_range) else NULL),
            title = panel_title %||% NULL,
            hoverlabel = list(font = list(family = font_family, size = text_size)),
            font = list(family = font_family, size = text_size),
            margin = list(t = 40, r = 10, l = 60, b = 40)
          ))
        }

        gvals_raw <- as.character(df_panel[[group_col]])

        # --- ORDER: honor x_order if provided (keeps user order, drops missing) ---
        group_levels <- {
          if (!is.null(x_order)) {
            intersect(x_order, unique(gvals_raw))
          } else if (is.factor(df_panel[[group_col]])) {
            intersect(levels(df_panel[[group_col]]), unique(gvals_raw))
          } else {
            unique(gvals_raw)
          }
        }

        pos    <- seq_along(group_levels)
        pos_of <- setNames(pos, group_levels)
        base_x <- unname(pos_of[gvals_raw])

        if (!is.null(x_labels) && !is.null(names(x_labels))) {
          ticktext <- unname(vapply(group_levels, function(lev) x_labels[[lev]] %||% lev, ""))
        } else {
          counts <- df_panel |>
            dplyr::group_by(.data[[group_col]]) |>
            dplyr::summarise(n = dplyr::n(), .groups = "drop")
          n_map <- setNames(counts$n, as.character(counts[[group_col]]))
          ticktext <- unname(vapply(group_levels, function(lev) sprintf("%s<br>(n = %s)", lev, n_map[[lev]] %||% "NA"), ""))
        }

        # --- COLORS: boxes and per-point colors mapped from custom_colors ---
        cols_hex  <- .color_map(group_levels, custom_colors)
        cols_fill <- if (!is.null(cols_hex)) vapply(cols_hex, .as_rgba, "", alpha = 0.60) else NULL
        col_per_point <- NULL
        if (!is.null(cols_hex)) {
          map_idx <- match(gvals_raw, group_levels)
          col_per_point <- cols_hex[map_idx]
        }

        set.seed(1L)
        xjit <- base_x + stats::runif(length(base_x), -jitter_width, jitter_width)

        p <- plotly::plot_ly()
        for (i in seq_along(group_levels)) {
          lev <- group_levels[i]
          idx <- which(gvals_raw == lev)
          p <- plotly::add_boxplot(
            p,
            x = rep(pos[i], length(idx)),
            y = df_panel[[need_metric]][idx],
            name = lev,
            boxpoints = "outliers",
            marker    = list(opacity = 0),
            hoverinfo = "skip",
            line      = list(color = cols_hex[i] %||% NULL, width = 1),
            fillcolor = cols_fill[i] %||% NULL,
            showlegend = FALSE,
            quartilemethod = quartile_method
          )
        }

        # points NOW get per-point colors matching their group's box
        p <- plotly::add_markers(
          p,
          x = xjit,
          y = df_panel[[need_metric]],
          marker = list(size = point_size, opacity = point_alpha, color = col_per_point %||% NULL),
          text = hovertext_vec, hovertemplate = "%{text}<extra></extra>",
          showlegend = FALSE
        )

        plotly::layout(
          p,
          xaxis = list(
            title = "Group",
            tickmode = "array", tickvals = pos, ticktext = ticktext,
            tickangle = x_tickangle %||% 0,
            range = c(min(pos) - 0.5, max(pos) + 0.5),
            showgrid = show_grids,
            tickfont = list(size = text_size, family = font_family)
          ),
          yaxis = c(
            list(title = ylab, showgrid = show_grids, tickfont = list(size = text_size, family = font_family)),
            if (!is.null(y_range)) list(range = y_range) else NULL
          ),
          title = panel_title %||% NULL,
          hoverlabel = list(font = list(family = font_family, size = text_size)),
          font = list(family = font_family, size = text_size),
          margin = list(t = 40, r = 10, l = 60, b = 70)
        )
      }

      # ------------ faceting or single ------------
      if (isTRUE(facet_by_rank) && !is.null(rank_col) && length(unique(alpha_df[[rank_col]])) > 1) {
        ranks <- unique(alpha_df[[rank_col]])
        panels <- lapply(
          ranks,
          function(rk) .build_panel(alpha_df[alpha_df[[rank_col]] == rk, , drop = FALSE],
                                    panel_title = paste0(rank_col, ": ", rk))
        )
        p_all <- plotly::subplot(
          panels,
          nrows  = ceiling(length(panels)/max(1L, ncol)),
          shareX = TRUE,
          shareY = !identical(facet_scales, "free"),
          margin = 0.02,
          titleX = TRUE, titleY = TRUE
        )
        return(plotly::layout(p_all, showlegend = FALSE,
                              font = list(family = font_family, size = text_size),
                              hoverlabel = list(font = list(family = font_family, size = text_size))))
      }
      .build_panel(alpha_df, panel_title = NULL)
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ------------------------------------------------------------------------------
# internal helper: build one enrichment-count plot (ggplot)
# - tbl: long table with sample_id, peptide_id, exist, and optional group_col
# - group_col: grouping variable to facet by; NULL => aggregate over all samples
# - prevalence_threshold: fraction of samples used for the dashed threshold line
# - custom_colors: optional named vector for cohort colors
# - binwidth: histogram bin width
# - title_label: optional title label (defaults to group_col or "all samples")
# ------------------------------------------------------------------------------

.plot_enrichment_counts_one <- function(
    tbl,
    group_col = NULL,
    prevalence_threshold = 0.05,
    custom_colors = NULL,
    binwidth = 1,
    title_label = if (is.null(group_col)) { "all samples" } else { group_col },
    annotation_size = annotation_size
) {
  .data  <- rlang::.data
  grp_sym <- if (!is.null(group_col)) rlang::sym(group_col) else NULL

  .ph_with_timing(
    headline = "building enrichment count plot",
    step = if (is.null(group_col)) "no grouping (aggregate)" else sprintf("grouping variable: '%s'", group_col),
    expr = {

      # -- validate required columns -------------------------------------------
      need <- c("sample_id", "peptide_id", "exist")
      if (!is.null(group_col)) need <- c(need, group_col)
      miss <- setdiff(need, colnames(tbl))
      if (length(miss)) {
        .ph_abort(
          headline = "missing required columns for plotting.",
          step     = "input validation",
          bullets  = sprintf("missing: %s", paste(add_quotes(miss, 1L), collapse = ", "))
        )
      }

      # -- build a unified cohort column (character) ---------------------------
      tbl <- if (is.null(group_col)) {
        dplyr::mutate(tbl, Cohort = "All samples")
      } else {
        dplyr::mutate(tbl, Cohort = !!grp_sym)
      }

      # -- keep only present peptide calls ------------------------------------
      tbl_present <- dplyr::filter(tbl, .data$exist > 0)

      # -- per-cohort sample sizes (for thresholds and labels) -----------------
      cohort_sizes <- tbl |>
        dplyr::distinct(.data$sample_id, .data$Cohort) |>
        dplyr::count(.data$Cohort, name = "n_samples") |>
        dplyr::collect()

      if (nrow(cohort_sizes) == 0L || any(cohort_sizes$n_samples <= 0)) {
        .ph_warn(
          headline = "no samples for at least one cohort (facet may be empty).",
          step     = "cohort sizes"
        )
        cohort_sizes <- dplyr::filter(cohort_sizes, .data$n_samples > 0)
      }

      # -- count present calls per peptide within cohort -----------------------
      pep_counts <- tbl_present |>
        dplyr::group_by(.data$Cohort, .data$peptide_id) |>
        dplyr::summarise(
          n_present = dplyr::n_distinct(.data$sample_id),
          .groups   = "drop"
        ) |>
        dplyr::filter(.data$n_present > 0) |>
        dplyr::collect()

      # -- compute cohort-specific prevalence thresholds -----------------------
      thresholds <- cohort_sizes |>
        dplyr::mutate(thresh = ceiling(.data$n_samples * prevalence_threshold))

      n_thresh_tbl <- pep_counts |>
        dplyr::inner_join(dplyr::select(thresholds, .data$Cohort, .data$thresh), by = "Cohort") |>
        dplyr::group_by(.data$Cohort) |>
        dplyr::summarise(
          n_peptides_thresh = sum(.data$n_present >= .data$thresh),
          .groups = "drop"
        )

      overall_tbl <- pep_counts |>
        dplyr::group_by(.data$Cohort) |>
        dplyr::summarise(
          n_overall = dplyr::n_distinct(.data$peptide_id),
          .groups   = "drop"
        )

      thresholds_df <- thresholds |>
        dplyr::left_join(n_thresh_tbl, by = "Cohort") |>
        dplyr::left_join(overall_tbl,  by = "Cohort") |>
        dplyr::mutate(
          n_peptides_thresh = dplyr::coalesce(.data$n_peptides_thresh, 0L),
          n_overall         = dplyr::coalesce(.data$n_overall, 0L),
          y_line            = pmax(.data$n_peptides_thresh, 1L),  # keep line on-log visible
          x_mid             = (.data$thresh + .data$n_samples) / 2  # centered label pos
        )

      # -- order cohorts deterministically for facets and colors ---------------
      real_order <- thresholds_df$Cohort

      count_df <- pep_counts |>
        dplyr::mutate(Cohort = factor(.data$Cohort, levels = real_order))

      # -- palette mapping (use phip palette unless custom provided) -----------
      lvl <- levels(count_df$Cohort)
      pal_map <- if (is.null(custom_colors)) {
        .phip_palette_map(lvl)
      } else {
        if (is.null(names(custom_colors))) {
          stats::setNames(rep_len(custom_colors, length(lvl)), lvl)
        } else {
          stats::setNames(unname(custom_colors[lvl]), lvl)
        }
      }

      count_df <- dplyr::mutate(
        count_df,
        fill_col = unname(pal_map[as.character(.data$Cohort)])
      )

      # -- plot ----------------------------------------------------------------
      p <- ggplot2::ggplot(
        count_df,
        ggplot2::aes(x = .data$n_present, fill = .data$fill_col)
      ) +
        ggplot2::geom_histogram(
          binwidth = binwidth,
          position = "identity",
          alpha = 1,
          colour = NA
        ) +
        ggplot2::scale_y_log10(
          breaks = 10^(0:6),
          labels = scales::trans_format("log10", scales::math_format(10^.x)),
          expand = ggplot2::expansion(mult = c(0, .15))
        ) +
        ggplot2::annotation_logticks(sides = "l", scaled = TRUE) +
        ggplot2::scale_fill_identity(guide = "none") +
        ggplot2::labs(
          x = "# of observations",
          y = expression("# of significantly bound peptides (" * log[10] * ")")
        ) +
        ggplot2::geom_segment(
          data = thresholds_df,
          ggplot2::aes(
            x = .data$thresh, xend = .data$n_samples,
            y = .data$y_line, yend = .data$y_line
          ),
          inherit.aes = FALSE,
          linetype = "dashed",
          color = "black",
          linewidth = 0.4,
          arrow = ggplot2::arrow(length = grid::unit(0.1, "cm"), ends = "both")
        ) +
        ggplot2::geom_text(
          data = thresholds_df,
          ggplot2::aes(
            x = .data$x_mid, y = .data$y_line,
            label = paste0(
              .data$n_peptides_thresh, " peptides in ≥",
              round(prevalence_threshold * 100), "%"
            )
          ),
          inherit.aes = FALSE,
          vjust = -0.6,
          size  = annotation_size
        ) +
        ggplot2::geom_text(
          data = thresholds_df,
          ggplot2::aes(
            x = .data$x_mid, y = .data$y_line,
            label = paste0(.data$n_overall, " peptides overall")
          ),
          inherit.aes = FALSE,
          vjust = 1.6,
          size  = annotation_size
        ) +
        theme_phip() +
        ggplot2::ggtitle(sprintf("Enrichment counts by %s", title_label))

      # facet by cohort when grouping is used
      if (!is.null(group_col)) {
        p <- p + ggplot2::facet_wrap(~Cohort, ncol = 2, scales = "free_x")
      }

      .ph_log_ok("plot built")
      p
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}
