# ==============================================================================
#' Static scatterplot of percent1 vs percent2 from `ph_prevalence_compare()`
#'
#' @description
#' A ggplot2 clone of `scatter_interactive()` with the same filtering behavior
#' and default coloring rules:
#' - If `color_by = NULL`, points are colored by significance category.
#'   If `category_rank_wbh` exists, it is used directly.
#'   Otherwise, categories are derived with priority:
#'     p_adj_rank_wbh -> p_adj_rank -> p_raw (threshold = `alpha`).
#'   If only one category remains, it auto-falls back to p-value bins for
#'   informative coloring.
#' - If `color_by` is provided, it joins peptide metadata (like interactive).
#'
#' @param df A `ph_prev_result` object or a data frame with prevalence results.
#' @param pair optional group pair (character length-2).
#' @param rank optional single rank (character) to keep.
#' @param universe optional `group_col` value or regex (if `universe_regex = TRUE`).
#' @param features optional character vector or regex patterns (if `features_regex = TRUE`).
#' @param features_regex logical; treat `features` as regex patterns.
#' @param universe_regex logical; treat `universe` as regex pattern(s).
#' @param xlab,ylab axis labels; if missing and `pair` is provided, they default
#'   to `pair[1]` and `pair[2]`.
#' @param alpha numeric in (0,1]; used only for nominal labels.
#' @param prefer_flags logical; reserved for future use (kept for back-compat).
#' @param color_by optional peptide-level meta column name to color by.
#' @param color_title optional legend title for `color_by`.
#' @param point_size Numeric; marker size for points.
#' @param point_alpha Numeric in (0,1); marker opacity.
#' @param jitter_width_pp Numeric; jitter width in percentage points.
#' @param jitter_height_pp Numeric; jitter height in percentage points.
#' @param font_family Character; font family for plot text.
#' @param font_size Numeric; font size for plot text.
#'
#' @return A ggplot object.
#' @export
scatter_static <- function(df,
                           pair = NULL,
                           rank = NULL,
                           universe = NULL,
                           features = NULL,
                           features_regex = FALSE,
                           universe_regex = FALSE,
                           xlab = NULL,
                           ylab = NULL,
                           alpha = 0.05,
                           prefer_flags = TRUE,
                           color_by = NULL,
                           color_title = NULL,
                           # ---- NEW ----
                           point_size      = 2.0,
                           point_alpha     = 0.85,
                           jitter_width_pp = 0.0,   # jitter in percentage points (x)
                           jitter_height_pp= 0.0,   # jitter in percentage points (y)
                           font_family     = NULL,
                           font_size       = 12) {
  clamp01 <- function(x) pmin(1, pmax(0, as.numeric(x)))
  point_alpha <- clamp01(point_alpha)

  df_original <- df

  # same filtering behavior as before -----------------------------------------
  if (inherits(df, "ph_prev_result")) {
    if (!is.null(pair)) {
      df <- prev_filter_pairs(
        df,
        gA = pair[1], gB = pair[2],
        ranks = rank,
        features = features, features_regex = features_regex,
        group_universe = universe, universe_regex = universe_regex,
        drop_na = TRUE
      )
    } else {
      df <- as.data.frame(df)
      if (!is.null(rank) && "rank" %in% names(df)) {
        df <- df[df$rank %in% rank, , drop = FALSE]
      }
      if (!is.null(universe) && "group_col" %in% names(df)) {
        gvec <- as.character(df$group_col)
        gvec[is.na(gvec)] <- ""
        if (isTRUE(universe_regex)) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(universe)) {
            keep <- keep | grepl(p, gvec, ignore.case = TRUE, perl = TRUE)
          }
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(gvec) %in% tolower(as.character(universe)), , drop = FALSE]
        }
      }
      if (!is.null(features) && "feature" %in% names(df)) {
        fvec <- as.character(df$feature)
        fvec[is.na(fvec)] <- ""
        if (isTRUE(features_regex)) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(features)) {
            keep <- keep | grepl(p, fvec, ignore.case = TRUE, perl = TRUE)
          }
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(fvec) %in% tolower(as.character(features)), , drop = FALSE]
        }
      }
    }
  } else {
    df <- as.data.frame(df)
  }

  if (is.null(df) || !nrow(df)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "no data for this contrast") +
        ggplot2::theme_void()
    )
  }

  # axis labels
  if (!is.null(pair)) {
    if (is.null(xlab)) xlab <- pair[1]
    if (is.null(ylab)) ylab <- pair[2]
  }
  if (is.null(xlab)) xlab <- "group a"
  if (is.null(ylab)) ylab <- "group b"

  # significance category if color_by is NULL ---------------------------------
  color_var <- NULL
  color_val_label <- NULL
  if (is.null(color_by)) {
    has_cat <- "category_rank_wbh" %in% names(df) &&
      any(nzchar(trimws(as.character(df$category_rank_wbh %||% ""))))
    if (has_cat) {
      df <- df %>%
        dplyr::mutate(
          category_raw = tolower(trimws(as.character(.data$category_rank_wbh %||% NA_character_))),
          category = dplyr::case_when(
            grepl("wbh", category_raw) & grepl("significant", category_raw) ~ "significant (wBH, per rank)",
            grepl("nominal", category_raw) ~ "nominal only",
            grepl("not significant", category_raw) ~ "not significant",
            TRUE ~ "not significant"
          )
        )
    } else {
      df <- df %>%
        dplyr::mutate(
          category = dplyr::case_when(
            !is.na(.data$p_adj_rank_wbh) & .data$p_adj_rank_wbh <= alpha ~ "significant (wBH, per rank)",
            !is.na(.data$p_raw) & .data$p_raw <= alpha ~ "nominal only",
            TRUE ~ "not significant"
          )
        )
    }
    cat_levels <- c("significant (wBH, per rank)",
                    "significant (BH, per rank)",
                    "nominal only",
                    "not significant")
    df$category <- factor(df$category, levels = cat_levels)
    if (length(stats::na.omit(unique(df$category))) <= 1L) {
      base_p <- if ("p_adj_rank_wbh" %in% names(df)) df$p_adj_rank_wbh else df$p_raw
      df$p_bin <- cut(base_p,
        breaks = c(0, 1e-3, 1e-2, 5e-2, 1, Inf),
        labels = c("<+1e-3", "(1e-3,1e-2]", "(1e-2,0.05]", ">0.05", "NA"),
        include.lowest = TRUE, right = TRUE
      )
      color_var <- "p_bin"
      color_val_label <- "p-value"
    } else {
      color_var <- "category"
      color_val_label <- NULL
    }
  }

  # join peptide-level meta if color_by ---------------------------------------
  if (!is.null(color_by)) {
    if ("peptide_id" %in% names(df)) {
      df <- dplyr::mutate(df, pep_key = as.character(.data$peptide_id))
    } else if (all(c("feature","rank") %in% names(df))) {
      df <- dplyr::mutate(df, pep_key = dplyr::if_else(.data$rank=="peptide_id",
                                                       as.character(.data$feature), NA_character_))
    } else if ("feature" %in% names(df)) {
      df <- dplyr::mutate(df, pep_key = as.character(.data$feature))
    } else {
      df <- dplyr::mutate(df, pep_key = NA_character_)
    }

    has_any_keys <- any(!is.na(df$pep_key))
    if (has_any_keys) {
      lib_handle <- tryCatch(
        {
          if (inherits(df_original, "ph_prev_result")) {
            attr(df_original, "prev_meta")$peptide_library
          } else {
            NULL
          }
        },
        error = function(...) NULL
      )

      pm <- if (!is.null(lib_handle)) {
        lib_handle %>%
          dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
          dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
          dplyr::collect()
      } else {
        get_peptide_library() %>%
          dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
          dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
          dplyr::collect()
      }

      df <- df %>% dplyr::left_join(pm, by = dplyr::join_by(pep_key == peptide_id))
      if (color_by %in% names(df) && is.list(df[[color_by]])) {
        df[[color_by]] <- vapply(df[[color_by]], function(z){
          if (length(z)==0) return(NA_character_)
          if (is.atomic(z)) return(as.character(z)[1])
          as.character(z[[1]])
        }, character(1))
      }
    } else {
      df[[color_by]] <- NA
    }
    if (is.logical(df[[color_by]])) {
      df[[color_by]] <- dplyr::case_when(
        df[[color_by]] %in% TRUE ~ "yes",
        df[[color_by]] %in% FALSE ~ "no",
        is.na(df[[color_by]]) ~ "na"
      )
    }
    df[[color_by]] <- as.factor(as.character(df[[color_by]]))
    color_var <- color_by
    if (is.null(color_title)) color_title <- color_by
    color_val_label <- color_title
  }

  # precompute jittered aesthetics (percent scale, 0..100) --------------------
  jw <- as.numeric(jitter_width_pp  %||% 0)
  jh <- as.numeric(jitter_height_pp %||% 0)
  pdat <- df
  if (jw > 0 || jh > 0) {
    set.seed(1L)  # deterministic jitter in the app
    n <- nrow(pdat)
    if (!"percent1" %in% names(pdat) && "prop1" %in% names(pdat)) pdat$percent1 <- pdat$prop1 * 100
    if (!"percent2" %in% names(pdat) && "prop2" %in% names(pdat)) pdat$percent2 <- pdat$prop2 * 100
    pdat$percent1 <- pmin(100, pmax(0, pdat$percent1 + stats::rnorm(n, 0, jw)))
    pdat$percent2 <- pmin(100, pmax(0, pdat$percent2 + stats::rnorm(n, 0, jh)))
  }

  # build plot -----------------------------------------------------------------
  p <- ggplot2::ggplot(pdat, ggplot2::aes(x = percent1, y = percent2))

  if (identical(color_var, "category")) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = .data$category),
                          alpha = point_alpha, size = point_size) +
      ggplot2::scale_color_manual(
        values = c(
          "significant (wBH, per rank)" = "#e31a1c",
          "nominal only"                = "#1b9e77",
          "not significant"             = "#386cb0"
        ),
        name = NULL, drop = FALSE
      )
  } else if (!is.null(color_var)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = .data[[color_var]]),
                          alpha = point_alpha, size = point_size) +
      ggplot2::labs(color = color_val_label)
  } else {
    p <- p + ggplot2::geom_point(alpha = point_alpha, size = point_size, color = "grey35")
  }

  rng <- range(c(pdat$percent1, pdat$percent2), na.rm = TRUE)
  p <- p +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         linewidth = 0.4, color = "grey60") +
    ggplot2::expand_limits(x = rng, y = rng) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(x = xlab, y = ylab)

  p <- p +
    theme_phip() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      text = ggplot2::element_text(family = font_family %||% ggplot2::theme_get()$text$family,
                                   size   = font_size)
    )

  p
}

# ==============================================================================
#' Interactive prevalence scatter for `ph_prev_result`
#'
#' @description
#' creates an interactive scatter (plotly) comparing prevalence in **group a**
#' vs **group b** for per-feature results produced by `ph_prevalence_compare()`.
#' accepts either a `ph_prev_result` (recommended) or a plain data.frame with the
#' same columns (`percent1`, `percent2`, `feature`, `group1`, `group2`, etc.).
#'
#' when a `ph_prev_result` is passed, you may specify a concrete pair of group
#' levels via `pair = c("A","B")`, and optionally restrict to a `rank`, a
#' `universe` (`group_col`) and/or `features`. internally, subsetting uses
#' `prev_filter_pairs()`, preserving metadata.
#'
#' color mapping:
#' - by default (when `color_by = NULL`), points are colored by `category_rank_wbh`
#'   with three levels: "significant (wBH, per rank)", "nominal only",
#'   "not significant".
#' - when `color_by` is a peptide-level meta column name (e.g. `"is_flagellum"`),
#'   the function tries to join peptide metadata using `peptide_id` (or `feature`
#'   if `rank == "peptide_id"`), and color by that variable.
#'
#' @param df a `ph_prev_result` or a `data.frame` with columns:
#'   `percent1`, `percent2`, `feature`, `group1`, `group2`,
#'   `p_raw`, `p_adj_rank`, `p_adj_rank_wbh`, `category_rank_wbh`,
#'   optionally `n_peptides`, `rank`, `peptide_id`, `group_col`.
#' @param pair optional length-2 character, e.g. `c("kid_serum::T2","kid_serum::T8")`.
#'   only used when `df` is a `ph_prev_result`; filters rows to that exact pair
#'   (order-agnostic).
#' @param rank optional single rank (character) to keep (only used with `ph_prev_result`).
#' @param universe optional `group_col` value or regex (if `universe_regex = TRUE`);
#'   only used with `ph_prev_result`.
#' @param features optional character vector or regex patterns (if `features_regex = TRUE`);
#'   only used with `ph_prev_result`.
#' @param features_regex logical; treat `features` as regex patterns (or).
#' @param universe_regex logical; treat `universe` as regex pattern(s) (or).
#' @param xlab,ylab axis labels; if missing and `pair` is provided, they default
#'   to `pair[1]` and `pair[2]`.
#' @param alpha numeric in (0,1]; used only for nominal labels; not the plotly alpha.
#' @param prefer_flags logical; reserved for future use (kept for back-compat).
#' @param color_by optional peptide-level meta column name to color by.
#' @param color_title optional legend title for `color_by`.
#' @param peplib Optional peptide metadata table used to color points by
#'   `color_by` when not present in `df`.
#' @param category_colors Optional named vector mapping categories to colors.
#' @param show_background Logical; plot background peptides if provided.
#' @param background_df Optional data frame of background points.
#' @param background_name Character; legend name for background points.
#' @param background_color Character; color for background points.
#' @param background_size Numeric; size for background points.
#' @param background_alpha Numeric in (0,1); opacity for background points.
#' @param background_max_n Integer; maximum number of background points to plot.
#' @param background_seed Integer; RNG seed for background downsampling.
#' @param point_line_width Numeric; outline width for points.
#' @param point_line_color Character; outline color for points.
#' @param point_size Numeric; marker size for points.
#' @param point_alpha Numeric in (0,1); marker opacity.
#' @param jitter_width_pp Numeric; jitter width in percentage points.
#' @param jitter_height_pp Numeric; jitter height in percentage points.
#' @param font_family Character; font family for plot text.
#' @param font_size Numeric; font size for plot text.
#'
#' @return a `plotly` object.
#'
#' @examples
#' \dontrun{
#' # typical usage with ph_prev_result:
#' p <- scatter_interactive(scatters,
#'   pair     = c("kid_serum::T2","kid_serum::T8"),
#'   rank     = "peptide_id",
#'   color_by = "is_flagellum",
#'   color_title = "Flagellum"
#' )
#' p
#' }
#'
#' @export
scatter_interactive <- function(df,
                                pair = NULL,
                                rank = NULL,
                                universe = NULL,
                                features = NULL,
                                features_regex = FALSE,
                                universe_regex = FALSE,
                                xlab = NULL,
                                ylab = NULL,
                                alpha = 0.05,
                                prefer_flags = TRUE,
                                color_by = NULL,
                                color_title = NULL,

                                # ---- NEW ----
                                peplib = NULL,

                                # Colors for the default category mapping (when color_by is NULL)
                                category_colors = c(
                                  "significant (wBH, per rank)" = "#009E73",
                                  "nominal only"                = "#E69F00",
                                  "not significant"             = "#0072B2",
                                  "other/NA"                    = "#999999"
                                ),

                                # Background overlay (separate trace, toggleable in legend)
                                show_background   = FALSE,
                                background_df     = NULL,
                                background_name   = "background peptides",
                                background_color  = "#808080",
                                background_size   = 6,
                                background_alpha  = 0.30,
                                background_max_n  = 3000,
                                background_seed   = 1L,

                                # Marker styling (main points)
                                point_size       = 7,
                                point_alpha      = 0.85,
                                point_line_width = 0.7,
                                point_line_color = "rgba(0,0,0,0.6)",

                                jitter_width_pp  = 0.0,
                                jitter_height_pp = 0.0,
                                font_family      = NULL,
                                font_size        = 12) {

  clamp01 <- function(x) pmin(1, pmax(0, as.numeric(x)))
  point_alpha <- clamp01(point_alpha)

  escape_html <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;",  x, fixed = TRUE)
    x <- gsub(">", "&gt;",  x, fixed = TRUE)
    x <- gsub('"', "&quot;", x, fixed = TRUE)
    x
  }
  fmt_na <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- "NA"
    x
  }

  df_original <- df

  # optional filtering like before --------------------------------------------
  if (inherits(df, "ph_prev_result")) {
    if (!is.null(pair)) {
      df <- prev_filter_pairs(
        df,
        gA = pair[1], gB = pair[2],
        ranks = rank,
        features = features, features_regex = features_regex,
        group_universe = universe, universe_regex = universe_regex,
        drop_na = TRUE
      )
    } else {
      df <- as.data.frame(df)
      if (!is.null(rank) && "rank" %in% names(df)) df <- df[df$rank %in% rank, , drop = FALSE]
      if (!is.null(universe) && "group_col" %in% names(df)) {
        gvec <- as.character(df$group_col); gvec[is.na(gvec)] <- ""
        if (universe_regex) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(universe)) keep <- keep | grepl(p, gvec, ignore.case = TRUE, perl = TRUE)
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(gvec) %in% tolower(as.character(universe)), , drop = FALSE]
        }
      }
      if (!is.null(features) && "feature" %in% names(df)) {
        fvec <- as.character(df$feature); fvec[is.na(fvec)] <- ""
        if (features_regex) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(features)) keep <- keep | grepl(p, fvec, ignore.case = TRUE, perl = TRUE)
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(fvec) %in% tolower(as.character(features)), , drop = FALSE]
        }
      }
    }
  } else {
    df <- as.data.frame(df)
  }

  if (is.null(df) || !nrow(df)) {
    p <- plotly::plot_ly()
    return(plotly::layout(
      p,
      annotations = list(
        x = 0.5, y = 0.5, text = "no data for this contrast",
        showarrow = FALSE, xref = "paper", yref = "paper"
      ),
      xaxis = list(title = if (is.null(xlab)) "group a" else xlab, zeroline = FALSE),
      yaxis = list(title = if (is.null(ylab)) "group b" else ylab, zeroline = FALSE),
      legend = list(orientation = "h", y = -0.15)
    ))
  }

  # axis labels
  if (!is.null(pair)) {
    if (is.null(xlab)) xlab <- pair[1]
    if (is.null(ylab)) ylab <- pair[2]
  }
  if (is.null(xlab)) xlab <- "group a"
  if (is.null(ylab)) ylab <- "group b"

  # categories if color_by is NULL --------------------------------------------
  cat_levels <- c("significant (wBH, per rank)", "nominal only", "not significant")
  if (is.null(color_by)) {
    df <- df %>%
      dplyr::mutate(
        category_raw = tolower(trimws(as.character(.data$category_rank_wbh %||% NA_character_))),
        category = dplyr::case_when(
          grepl("wbh", category_raw) & grepl("significant", category_raw) ~ cat_levels[1],
          grepl("nominal", category_raw)                                  ~ cat_levels[2],
          grepl("not significant", category_raw)                          ~ cat_levels[3],
          TRUE ~ "other/NA"
        ),
        category = factor(category, levels = c(cat_levels, "other/NA"))
      )
  }

  # join peptide-level meta for color_by (keep old behaviour) -----------------
  color_var <- NULL
  color_val_label <- NULL
  if (!is.null(color_by)) {
    if ("peptide_id" %in% names(df)) {
      df <- dplyr::mutate(df, pep_key = as.character(.data$peptide_id))
    } else if (all(c("feature","rank") %in% names(df))) {
      df <- dplyr::mutate(df, pep_key = dplyr::if_else(.data$rank == "peptide_id",
                                                       as.character(.data$feature), NA_character_))
    } else if ("feature" %in% names(df)) {
      df <- dplyr::mutate(df, pep_key = as.character(.data$feature))
    } else {
      df <- dplyr::mutate(df, pep_key = NA_character_)
    }

    has_any_keys <- any(!is.na(df$pep_key))
    if (has_any_keys) {
      pm <- NULL
      if (!is.null(peplib) && ("peptide_id" %in% names(peplib)) && (color_by %in% names(peplib))) {
        pm <- peplib %>%
          dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
          dplyr::distinct(peptide_id, .keep_all = TRUE)
      } else {
        lib_handle <- NULL
        if (inherits(df_original, "ph_prev_result")) {
          meta <- attr(df_original, "prev_meta")
          if (!is.null(meta) && "peptide_library" %in% names(meta)) lib_handle <- meta$peptide_library
        }
        pm <- if (!is.null(lib_handle)) {
          lib_handle %>%
            dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
            dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
            dplyr::collect()
        } else {
          get_peptide_library() %>%
            dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
            dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
            dplyr::collect()
        }
      }

      df <- df %>% dplyr::left_join(pm, by = dplyr::join_by(pep_key == peptide_id))

      if (color_by %in% names(df) && is.list(df[[color_by]])) {
        df[[color_by]] <- vapply(df[[color_by]], function(z) {
          if (length(z) == 0) return(NA_character_)
          if (is.atomic(z)) return(as.character(z)[1])
          as.character(z[[1]])
        }, character(1))
      }
    } else {
      df[[color_by]] <- NA
    }

    if (is.logical(df[[color_by]])) {
      df[[color_by]] <- dplyr::case_when(
        df[[color_by]] %in% TRUE  ~ "yes",
        df[[color_by]] %in% FALSE ~ "no",
        is.na(df[[color_by]])     ~ "na"
      )
    }

    df[[color_by]] <- as.factor(as.character(df[[color_by]]))
    color_var <- color_by
    if (is.null(color_title)) color_title <- color_by
    color_val_label <- color_title
  }

  # Create jittered copies in percent scale -----------------------------------
  jw <- as.numeric(jitter_width_pp  %||% 0)
  jh <- as.numeric(jitter_height_pp %||% 0)

  pdat <- df
  if (!"percent1" %in% names(pdat) && "prop1" %in% names(pdat)) pdat$percent1 <- pdat$prop1 * 100
  if (!"percent2" %in% names(pdat) && "prop2" %in% names(pdat)) pdat$percent2 <- pdat$prop2 * 100

  if (jw > 0 || jh > 0) {
    set.seed(1L)
    n <- nrow(pdat)
    pdat$percent1 <- pmin(100, pmax(0, pdat$percent1 + stats::rnorm(n, 0, jw)))
    pdat$percent2 <- pmin(100, pmax(0, pdat$percent2 + stats::rnorm(n, 0, jh)))
  }

  # ---- peptide_id + peplib join (ONLY for peptide points) -------------------
  pdat$pep_id <- NA_character_
  if (all(c("rank", "feature") %in% names(pdat))) {
    pdat$pep_id <- ifelse(pdat$rank == "peptide_id", as.character(pdat$feature), NA_character_)
  } else if ("peptide_id" %in% names(pdat)) {
    pdat$pep_id <- as.character(pdat$peptide_id)
  }

  if (!is.null(peplib) && any(!is.na(pdat$pep_id))) {
    pm2 <- peplib %>%
      dplyr::select(dplyr::any_of(c("peptide_id", "species", "common", "Fullname", "Description"))) %>%
      dplyr::distinct(peptide_id, .keep_all = TRUE)
    pdat <- dplyr::left_join(pdat, pm2, by = c("pep_id" = "peptide_id"))
  } else {
    if (!"species" %in% names(pdat))      pdat$species <- NA_character_
    if (!"common" %in% names(pdat))       pdat$common <- NA_character_
    if (!"Fullname" %in% names(pdat))     pdat$Fullname <- NA_character_
    if (!"Description" %in% names(pdat))  pdat$Description <- NA_character_
  }

  # enriched hover text (KEEP OLD PREVALENCE REPORTING) -----------------------
  pdat <- pdat %>%
    dplyr::mutate(
      p1r  = round(.data$percent1, 2),
      p2r  = round(.data$percent2, 2),
      rr   = dplyr::coalesce(round(
        if (!"ratio" %in% names(.) & all(c("prop1","prop2") %in% names(.)))
          (pmax(prop1,1e-12)/pmax(prop2,1e-12)) else ratio, 2), NA_real_),
      praw = if ("p_raw" %in% names(.)) round(.data$p_raw, 3) else NA_real_,
      padj_wbh = if ("p_adj_rank_wbh" %in% names(.)) round(.data$p_adj_rank_wbh, 3) else NA_real_,
      padj_bh  = if ("p_adj_rank" %in% names(.)) round(.data$p_adj_rank, 3) else NA_real_,
      padj_bh_str = ifelse(is.na(padj_bh), "NA", sprintf("%.3f", padj_bh)),
      npep = dplyr::coalesce(
        if ("n_peptides" %in% names(.)) .data[["n_peptides"]] else NULL,
        if ("n_peptide"  %in% names(.)) .data[["n_peptide" ]] else NULL
      ),
      color_info = if (!is.null(color_var)) as.character(.data[[color_var]]) else NA_character_,
      meta_block = dplyr::if_else(
        !is.na(.data$pep_id),
        paste0(
          "<br>species: ",     escape_html(fmt_na(.data$species)),
          "<br>common: ",      escape_html(fmt_na(.data$common)),
          "<br>Fullname: ",    escape_html(fmt_na(.data$Fullname)),
          "<br>Description: ", escape_html(fmt_na(.data$Description))
        ),
        ""
      ),
      text = sprintf(
        "<b>%s</b>%s<br>%s: %d/%d (%.2f%%)<br>%s: %d/%d (%.2f%%)%s<br>peptides: %s<br>p: %.3f<br>p_adj (bh): %s<br>p_adj (wbh): %.3f",
        .data$feature,
        .data$meta_block,
        xlab, .data$n1, .data$N1, .data$p1r,
        ylab, .data$n2, .data$N2, .data$p2r,
        if (!is.null(color_var)) sprintf("<br>%s: %s", color_val_label, .data$color_info) else "",
        format(.data$npep, big.mark = "\u2009", scientific = FALSE),
        .data$praw, .data$padj_bh_str, .data$padj_wbh
      )
    )

  pdat$text <- as.character(pdat$text)

  # ---- background overlay prep ----------------------------------------------
  bg <- NULL
  if (isTRUE(show_background) && !is.null(background_df) && nrow(background_df) > 0) {
    bg <- as.data.frame(background_df)

    if (!"percent1" %in% names(bg) && "prop1" %in% names(bg)) bg$percent1 <- bg$prop1 * 100
    if (!"percent2" %in% names(bg) && "prop2" %in% names(bg)) bg$percent2 <- bg$prop2 * 100

    bg <- bg[, intersect(c("percent1", "percent2"), names(bg)), drop = FALSE]
    bg$percent1 <- as.numeric(bg$percent1)
    bg$percent2 <- as.numeric(bg$percent2)
    bg <- bg[stats::complete.cases(bg), , drop = FALSE]

    if (is.finite(background_max_n) && nrow(bg) > background_max_n) {
      set.seed(as.integer(background_seed))
      bg <- bg[sample.int(nrow(bg), size = background_max_n), , drop = FALSE]
    }
  }

  # ---- build plot: background first (bottom), then diagonal, then main -------
  p <- plotly::plot_ly()

  if (!is.null(bg) && nrow(bg) > 0) {
    p <- plotly::add_trace(
      p,
      data = bg,
      x = ~percent1, y = ~percent2,
      type = "scatter", mode = "markers",
      name = background_name,
      showlegend = TRUE,
      hoverinfo = "skip",
      inherit = FALSE,
      legendrank = 999,
      marker = list(
        size    = background_size,
        color   = background_color,
        opacity = clamp01(background_alpha)
      )
    )
  }

  rng <- range(c(pdat$percent1, pdat$percent2), na.rm = TRUE)
  p <- plotly::add_trace(
    p,
    x = rng, y = rng,
    type = "scatter", mode = "lines",
    inherit = FALSE,
    showlegend = FALSE,
    hoverinfo = "skip"
  )

  if (is.null(color_var)) {
    # Explicit colors per category = always works (no plotly palette surprises)
    pdat$category <- droplevels(pdat$category)
    levs <- levels(pdat$category)

    default_pal <- c(
      "significant (wBH, per rank)" = "#009E73",
      "nominal only"                = "#E69F00",
      "not significant"             = "#0072B2",
      "other/NA"                    = "#999999"
    )
    pal_map <- default_pal
    if (!is.null(category_colors) && length(category_colors) > 0) {
      pal_map[names(category_colors)] <- category_colors
    }
    # ensure every present level has a color
    missing_cols <- setdiff(levs, names(pal_map))
    if (length(missing_cols) > 0) pal_map[missing_cols] <- "#999999"

    for (lv in levs) {
      dd <- pdat[pdat$category == lv, , drop = FALSE]
      if (!nrow(dd)) next
      p <- plotly::add_trace(
        p,
        data = dd,
        x = ~percent1, y = ~percent2,
        type = "scatter", mode = "markers",
        name = lv,
        showlegend = TRUE,
        inherit = FALSE,
        text = ~text,
        hovertemplate = "%{text}<extra></extra>",
        marker = list(
          size    = point_size,
          opacity = point_alpha,
          color   = unname(pal_map[lv]),
          line    = list(width = point_line_width, color = point_line_color)
        )
      )
    }

  } else {
    # keep old color_by behaviour
    p <- plotly::add_trace(
      p,
      data = pdat,
      x = ~percent1, y = ~percent2,
      color = stats::as.formula(paste0("~`", color_var, "`")),
      type = "scatter", mode = "markers",
      name = if (!is.null(color_title)) color_title else color_var,
      inherit = FALSE,
      text = ~text,
      hovertemplate = "%{text}<extra></extra>",
      marker = list(
        size    = point_size,
        opacity = point_alpha,
        line    = list(width = point_line_width, color = point_line_color)
      )
    )
  }

  plotly::layout(
    p,
    font  = list(family = font_family, size = font_size),
    xaxis = list(title = xlab, zeroline = FALSE, range = c(0, 100)),
    yaxis = list(title = ylab, zeroline = FALSE, range = c(0, 100)),
    legend = list(orientation = "h", y = -0.15)
  )
}

# ------------------------------------------------------------------------------
# Internal helper to prepare volcano data (shared by static/interactive)
# ------------------------------------------------------------------------------
.volcano_prepare <- function(df,
                             pair = NULL,
                             rank = NULL,
                             universe = NULL,
                             features = NULL,
                             features_regex = FALSE,
                             universe_regex = FALSE,
                             color_by = NULL,
                             color_title = NULL,
                             fc_cut = 1,
                             p_cut = 0.05,
                             p_mode = c("raw","bh","wbh"),
                             significant_colors = c(
                               "not significant"                 = "#386cb0",
                               "significant prior correction"    = "#1b9e77",
                               "significant post fdr correction" = "#e31a1c"
                             )) {
  # --- choose p mode ----------------------------------------------------------
  p_mode <- match.arg(p_mode)

  # Keep original (for peptide library handle)
  df_original <- df

  # --- Filtering identical to scatter_* --------------------------------------
  if (inherits(df, "ph_prev_result")) {
    if (!is.null(pair)) {
      df <- prev_filter_pairs(
        df,
        gA = pair[1], gB = pair[2],
        ranks = rank,
        features = features, features_regex = features_regex,
        group_universe = universe, universe_regex = universe_regex,
        drop_na = TRUE
      )
    } else {
      df <- as.data.frame(df)
      if (!is.null(rank) && "rank" %in% names(df)) {
        df <- df[df$rank %in% rank, , drop = FALSE]
      }
      if (!is.null(universe) && "group_col" %in% names(df)) {
        gvec <- as.character(df$group_col); gvec[is.na(gvec)] <- ""
        if (isTRUE(universe_regex)) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(universe)) keep <- keep | grepl(p, gvec, ignore.case = TRUE, perl = TRUE)
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(gvec) %in% tolower(as.character(universe)), , drop = FALSE]
        }
      }
      if (!is.null(features) && "feature" %in% names(df)) {
        fvec <- as.character(df$feature); fvec[is.na(fvec)] <- ""
        if (isTRUE(features_regex)) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(features)) keep <- keep | grepl(p, fvec, ignore.case = TRUE, perl = TRUE)
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(fvec) %in% tolower(as.character(features)), , drop = FALSE]
        }
      }
    }
  } else {
    df <- as.data.frame(df)
    if (!is.null(rank) && "rank" %in% names(df)) {
      df <- df[df$rank %in% rank, , drop = FALSE]
    }
  }

  # Empty-guard
  if (is.null(df) || !nrow(df)) {
    return(list(empty = TRUE))
  }

  # Required columns
  req <- c("rank","feature","group1","group2")
  miss <- setdiff(req, names(df))
  if (length(miss)) {
    stop("volcano_*(): missing required columns: ", paste(miss, collapse = ", "))
  }

  # Coercions & pair label
  df$rank    <- as.character(df$rank)
  df$feature <- as.character(df$feature)
  df$group1  <- as.character(df$group1)
  df$group2  <- as.character(df$group2)
  df$pair_label <- paste0(df$group1, " vs ", df$group2)

  # Ensure ratio (fallback from prop1/prop2)
  if (!"ratio" %in% names(df)) {
    if (!all(c("prop1","prop2") %in% names(df))) {
      stop("volcano_*(): need either 'ratio' or both 'prop1' and 'prop2'.")
    }
    p1 <- ifelse(is.na(df$prop1) | df$prop1 <= 0, df$prop1 + 1e-12, df$prop1)
    p2 <- ifelse(is.na(df$prop2) | df$prop2 <= 0, df$prop2 + 1e-12, df$prop2)
    df$ratio <- p1 / p2
  }

  # Pick p column
  p_col <- switch(p_mode,
                  raw = "p_raw",
                  bh  = "p_adj_rank",
                  wbh = "p_adj_rank_wbh")
  if (!p_col %in% names(df)) {
    stop("volcano_*(): requested p_mode='", p_mode, "' but column '", p_col, "' not found.")
  }

  # Axes
  df$log2ratio <- log2(df$ratio)
  p_use <- pmax(df[[p_col]], .Machine$double.xmin)
  df$nlog10p <- -log10(p_use)

  if (!any(is.finite(df$log2ratio)) || !any(is.finite(df$nlog10p))) {
    stop("volcano_*(): Non-finite values for log2(ratio) or -log10(p).")
  }

  # Default categories if color_by is NULL (mirror scatter logic)
  df$category_use <- NA_character_
  if ("category_rank_wbh" %in% names(df)) {
    raw <- tolower(trimws(as.character(df$category_rank_wbh)))
    df$category_use <- ifelse(grepl("significant", raw) & grepl("wbh", raw),
                              "significant post fdr correction",
                              ifelse(grepl("nominal", raw),
                                     "significant prior correction",
                                     "not significant"))
  } else {
    rawp  <- if ("p_raw" %in% names(df)) df$p_raw else p_use
    sigfc <- abs(df$log2ratio) >= fc_cut
    post  <- (p_use <= p_cut) & sigfc
    prior <- (rawp  <= p_cut) & sigfc & !post
    df$category_use <- ifelse(post,  "significant post fdr correction",
                              ifelse(prior,"significant prior correction","not significant"))
  }
  df$category_use <- factor(df$category_use,
                            levels = c("significant post fdr correction",
                                       "significant prior correction",
                                       "not significant"))

  # Join peptide metadata if color_by requested (use saved library)
  color_var <- NULL
  color_label <- NULL
  if (!is.null(color_by)) {
    if ("peptide_id" %in% names(df)) {
      df$pep_key <- as.character(df$peptide_id)
    } else if (all(c("feature", "rank") %in% names(df))) {
      df$pep_key <- ifelse(df$rank == "peptide_id", as.character(df$feature), NA_character_)
    } else if ("feature" %in% names(df)) {
      df$pep_key <- as.character(df$feature)
    } else {
      df$pep_key <- NA_character_
    }

    has_keys <- any(!is.na(df$pep_key))
    if (has_keys) {
      lib_handle <- tryCatch({
        if (inherits(df_original, "ph_prev_result"))
          attr(df_original, "prev_meta")$peptide_library
        else NULL
      }, error = function(...) NULL)

      if (!is.null(lib_handle)) {
        pm <- lib_handle %>%
          dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
          dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
          dplyr::collect()
      } else {
        pm <- get_peptide_library() %>%
          dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
          dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
          dplyr::collect()
      }

      df <- dplyr::left_join(df, pm, by = dplyr::join_by(pep_key == peptide_id))

      if (color_by %in% names(df) && is.list(df[[color_by]])) {
        df[[color_by]] <- vapply(
          df[[color_by]],
          function(z) {
            if (length(z) == 0) return(NA_character_)
            if (is.atomic(z)) return(as.character(z)[1])
            as.character(z[[1]])
          },
          character(1)
        )
      }
      if (is.logical(df[[color_by]])) {
        df[[color_by]] <- dplyr::case_when(
          df[[color_by]] %in% TRUE  ~ "yes",
          df[[color_by]] %in% FALSE ~ "no",
          is.na(df[[color_by]])     ~ "na"
        )
      }
    } else {
      df[[color_by]] <- NA
    }

    df[[color_by]] <- as.factor(as.character(df[[color_by]]))
    color_var   <- color_by
    color_label <- color_title %||% color_by
  }

  list(
    empty = FALSE,
    data = df,
    color_var = color_var,
    color_label = color_label,
    significant_colors = significant_colors,
    fc_cut = fc_cut,
    p_cut = p_cut,
    p_mode = p_mode
  )
}

# ------------------------------------------------------------------------------
# Static volcano (ggplot2)
# ------------------------------------------------------------------------------
#' Static volcano plot (log2 ratio vs -log10 p)
#'
#' @param df A `ph_prev_result` object or a data frame with prevalence results.
#' @param pair optional group pair (character length-2).
#' @param rank optional single rank (character) to keep.
#' @param universe optional `group_col` value or regex (if `universe_regex = TRUE`).
#' @param features optional character vector or regex patterns (if `features_regex = TRUE`).
#' @param features_regex logical; treat `features` as regex patterns.
#' @param universe_regex logical; treat `universe` as regex pattern(s).
#' @param color_by optional peptide-level meta column name to color by.
#' @param color_title optional legend title for `color_by`.
#' @param fc_cut Numeric; absolute log2 fold-change cutoff.
#' @param p_cut Numeric; p-value cutoff.
#' @param p_mode One of `c("raw","bh","wbh")` controlling which p-values to use.
#' @param significant_colors Named vector of colors for significance categories.
#'
#' @return A `ggplot` object.
#' @export
volcano_static <- function(df,
                           pair = NULL,
                           rank = NULL,
                           universe = NULL,
                           features = NULL,
                           features_regex = FALSE,
                           universe_regex = FALSE,
                           color_by = NULL,
                           color_title = NULL,
                           fc_cut = 1,
                           p_cut = 0.05,
                           p_mode = c("raw","bh","wbh"),
                           significant_colors = c(
                             "not significant"                 = "#386cb0",
                             "significant prior correction"    = "#1b9e77",
                             "significant post fdr correction" = "#e31a1c"
                           )) {
  prep <- .volcano_prepare(
    df, pair, rank, universe, features, features_regex, universe_regex,
    color_by, color_title, fc_cut, p_cut, p_mode, significant_colors
  )
  if (isTRUE(prep$empty)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "no data for this contrast") +
        ggplot2::theme_void()
    )
  }
  dd <- prep$data

  aes_base <- ggplot2::aes(x = log2ratio, y = nlog10p)
  p <- if (!is.null(prep$color_var)) {
    ggplot2::ggplot(dd, aes_base + ggplot2::aes(color = .data[[prep$color_var]])) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::labs(color = prep$color_label)
  } else {
    ggplot2::ggplot(dd, aes_base + ggplot2::aes(color = category_use)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::scale_color_manual(values = prep$significant_colors, name = NULL)
  }

  p +
    ggplot2::geom_hline(yintercept = -log10(prep$p_cut), linetype = "dashed", color = "gray50") +
    ggplot2::geom_vline(xintercept = c(-prep$fc_cut, 0, prep$fc_cut), linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      x = "log\u2082 ratio (group1 / group2)",
      y = "-log\u2081\u2080(p)"
    ) +
    theme_phip(base_size = 12) +
    ggplot2::theme(legend.position = if (!is.null(prep$color_var)) "bottom" else "none")
}

# ------------------------------------------------------------------------------
# Interactive volcano (plotly)
# ------------------------------------------------------------------------------
#' Interactive volcano plot (log2 ratio vs -log10 p)
#'
#' @inheritParams volcano_static
#' @return A `plotly` htmlwidget.
#' @export
volcano_interactive <- function(df,
                                pair = NULL,
                                rank = NULL,
                                universe = NULL,
                                features = NULL,
                                features_regex = FALSE,
                                universe_regex = FALSE,
                                color_by = NULL,
                                color_title = NULL,
                                fc_cut = 1,
                                p_cut = 0.05,
                                p_mode = c("raw","bh","wbh"),
                                significant_colors = c(
                                  "not significant"                 = "#386cb0",
                                  "significant prior correction"    = "#1b9e77",
                                  "significant post fdr correction" = "#e31a1c"
                                )) {
  prep <- .volcano_prepare(
    df, pair, rank, universe, features, features_regex, universe_regex,
    color_by, color_title, fc_cut, p_cut, p_mode, significant_colors
  )
  if (isTRUE(prep$empty)) {
    p <- plotly::plot_ly()
    return(plotly::layout(
      p,
      annotations = list(x = 0.5, y = 0.5, text = "no data for this contrast",
                         showarrow = FALSE, xref = "paper", yref = "paper"),
      xaxis = list(title = "log\u2082 ratio (group1 / group2)", zeroline = FALSE),
      yaxis = list(title = "-log\u2081\u2080(p)", zeroline = FALSE),
      legend = list(orientation = "h", y = -0.15)
    ))
  }

  dd <- prep$data

  # Hover text (concise)
  hov <- sprintf(
    "<b>%s</b><br>rank: %s<br>pair: %s vs %s<br>log2(ratio): %.3f<br>-log10(p): %.3f%s",
    dd$feature, dd$rank, dd$group1, dd$group2, dd$log2ratio, dd$nlog10p,
    if (!is.null(prep$color_var))
      sprintf("<br>%s: %s", prep$color_label, as.character(dd[[prep$color_var]]))
    else ""
  )

  # Colors
  if (!is.null(prep$color_var)) {
    pal <- grDevices::hcl.colors(
      n = max(3L, length(levels(dd[[prep$color_var]]))),
      palette = "Dark 3"
    )
    pal <- stats::setNames(pal[seq_along(levels(dd[[prep$color_var]]))],
                    levels(dd[[prep$color_var]]))
    col_vec <- pal[as.character(dd[[prep$color_var]])]
  } else {
    pal <- prep$significant_colors
    cat_vec <- as.character(dd$category_use)
    col_vec <- pal[cat_vec]; col_vec[is.na(col_vec)] <- pal[["not significant"]]
  }

  p <- plotly::plot_ly(
    dd,
    x = ~log2ratio, y = ~nlog10p,
    type = "scatter", mode = "markers",
    text = hov,
    hovertemplate = "%{text}<extra></extra>",
    marker = list(size = 7, opacity = 0.85, color = col_vec)
  )

  # Threshold lines
  p <- plotly::add_segments(p, x = -Inf, xend = Inf,
                            y = -log10(prep$p_cut), yend = -log10(prep$p_cut),
                            inherit = FALSE, showlegend = FALSE)
  p <- plotly::add_segments(p, x =  prep$fc_cut, xend =  prep$fc_cut,
                            y = -Inf, yend = Inf, inherit = FALSE, showlegend = FALSE)
  p <- plotly::add_segments(p, x = -prep$fc_cut, xend = -prep$fc_cut,
                            y = -Inf, yend = Inf, inherit = FALSE, showlegend = FALSE)
  p <- plotly::add_segments(p, x = 0, xend = 0,
                            y = -Inf, yend = Inf, inherit = FALSE, showlegend = FALSE)

  leg <- list(orientation = "h", y = -0.15)
  if (!is.null(prep$color_var) && !is.null(prep$color_label)) {
    leg$title <- list(text = prep$color_label)
  }

  plotly::layout(
    p,
    xaxis = list(title = "log\u2082 ratio (group1 / group2)", zeroline = FALSE),
    yaxis = list(title = "-log\u2081\u2080(p)", zeroline = FALSE),
    legend = leg,
    showlegend = !is.null(prep$color_var)
  )
}
