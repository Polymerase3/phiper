# ==============================================================================
# Internal helper: join peptide metadata and create a grouping column from a
# named color_by vector, e.g. c("is_flagellum" = TRUE) or
# c("species" = "Staphylococcus aureus"). Multiple entries are supported;
# last-match wins for points satisfying more than one condition.
# Returns list(df = df_with_.color_group, group_labels = character vector).
.build_color_group <- function(df, color_by, peplib = NULL) {
  col_names <- unique(names(color_by))

  # derive peptide key
  if ("peptide_id" %in% names(df)) {
    df$pep_key <- as.character(df$peptide_id)
  } else if (all(c("feature", "rank") %in% names(df))) {
    df$pep_key <- ifelse(df$rank == "peptide_id", as.character(df$feature), NA_character_)
  } else if ("feature" %in% names(df)) {
    df$pep_key <- as.character(df$feature)
  } else {
    df$pep_key <- NA_character_
  }

  if (any(!is.na(df$pep_key))) {
    if (!is.null(peplib) &&
        all(c("peptide_id", col_names) %in% names(peplib))) {
      pm <- peplib %>%
        dplyr::select("peptide_id", tidyselect::all_of(col_names)) %>%
        dplyr::distinct(peptide_id, .keep_all = TRUE)
    } else {
      pm <- get_peptide_library() %>%
        dplyr::select("peptide_id", tidyselect::all_of(col_names)) %>%
        dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
        dplyr::collect()
    }
    df <- dplyr::left_join(df, pm, by = dplyr::join_by(pep_key == peptide_id))

    # fix list-columns that may arise from Arrow/duckdb
    for (cn in col_names) {
      if (cn %in% names(df) && is.list(df[[cn]])) {
        df[[cn]] <- vapply(df[[cn]], function(z) {
          if (length(z) == 0L) return(NA_character_)
          if (is.atomic(z))    return(as.character(z)[1L])
          as.character(z[[1L]])
        }, character(1L))
      }
    }
  }

  # build grouping column (last matching condition wins)
  group_labels    <- character(length(color_by))
  df$.color_group <- "other"
  for (i in seq_along(color_by)) {
    cn  <- names(color_by)[i]
    val <- color_by[[i]]
    lbl <- paste0(cn, ": ", as.character(val))
    group_labels[i] <- lbl
    if (cn %in% names(df)) {
      matches <- !is.na(df[[cn]]) &
        (as.character(df[[cn]]) == as.character(val))
      df$.color_group[matches] <- lbl
    }
  }
  present <- group_labels[group_labels %in% df$.color_group]
  df$.color_group <- factor(df$.color_group,
                            levels = c(present, "other"))
  list(df = df, group_labels = group_labels)
}

# ==============================================================================
#' Static scatterplot of percent1 vs percent2 from `ph_prevalence_compare()`
#'
#' @description
#' A ggplot2 scatterplot comparing prevalence in **group a** vs **group b**.
#' Default coloring uses BH-corrected p-values computed per-plot from `p_raw`:
#' "significant (BH)", "nominal only", "not significant".
#' If only one category is present the plot falls back to p-value bins.
#'
#' When `color_by` is supplied as a named vector, peptide metadata is joined and
#' points matching the specified values are highlighted.  Multiple groups may be
#' given simultaneously:
#' \preformatted{
#'   color_by = c("is_flagellum" = TRUE)
#'   color_by = c("species" = "Staphylococcus aureus")
#'   color_by = c("is_flagellum" = TRUE, "species" = "Staphylococcus aureus")
#' }
#'
#' @param df A data frame with prevalence results.
#' @param pair optional group pair (character length-2).
#' @param rank optional single rank (character) to keep.
#' @param xlab,ylab axis labels; defaults to `pair[1]`/`pair[2]` when `pair` is given.
#' @param alpha numeric in (0,1]; significance threshold for category labels.
#' @param color_by optional named vector identifying peptide-library values to
#'   highlight, e.g. `c("is_flagellum" = TRUE)` or
#'   `c("species" = "Staphylococcus aureus")`.
#' @param color_title optional legend title when `color_by` is used.
#' @param ... graphical parameters: `point_size` (default 2), `point_alpha`
#'   (default 0.85), `jitter_width_pp` (default 0), `jitter_height_pp`
#'   (default 0), `font_family`, `font_size` (default 12).
#'
#' @return A ggplot object.
#'
#' @examples
#' set.seed(1)
#' prev <- data.frame(
#'   rank       = "peptide_id",
#'   feature    = paste0("pep", 1:30),
#'   group1     = "A",
#'   group2     = "B",
#'   prop1      = runif(30),
#'   prop2      = runif(30),
#'   percent1   = runif(30, 0, 100),
#'   percent2   = runif(30, 0, 100),
#'   ratio      = runif(30, 0.1, 10),
#'   p_raw      = runif(30),
#'   n_peptides = 1L
#' )
#'
#' # basic plot
#' scatter_static(prev)
#'
#' # filter to a specific pair and set axis labels
#' scatter_static(prev,
#'   pair  = c("A", "B"),
#'   xlab  = "Group A (%)",
#'   ylab  = "Group B (%)",
#'   alpha = 0.05
#' )
#' @export
scatter_static <- function(df,
                           pair        = NULL,
                           rank        = NULL,
                           xlab        = NULL,
                           ylab        = NULL,
                           alpha       = 0.05,
                           color_by    = NULL,
                           color_title = NULL,
                           ...) {
  dots             <- list(...)
  point_size       <- dots[["point_size"]]       %||% 2.0
  point_alpha      <- dots[["point_alpha"]]      %||% 0.85
  jitter_width_pp  <- dots[["jitter_width_pp"]]  %||% 0.0
  jitter_height_pp <- dots[["jitter_height_pp"]] %||% 0.0
  font_family      <- dots[["font_family"]]      %||% NULL
  font_size        <- dots[["font_size"]]        %||% 12

  clamp01     <- function(x) pmin(1, pmax(0, as.numeric(x)))
  point_alpha <- clamp01(point_alpha)

  # filtering ----------------------------------------------------------------
  df <- as.data.frame(df)
  if (!is.null(pair)) {
    df <- .ph_filter_pairs(
      df,
      gA = pair[1], gB = pair[2],
      ranks   = rank,
      drop_na = TRUE
    )
  } else if (!is.null(rank) && "rank" %in% names(df)) {
    df <- df[df$rank %in% rank, , drop = FALSE]
  }

  if (is.null(df) || !nrow(df)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "no data for this contrast") +
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

  # default significance coloring (BH per-plot) ------------------------------
  color_var       <- NULL
  color_val_label <- NULL
  if (is.null(color_by)) {
    df$p_bh    <- stats::p.adjust(df$p_raw, method = "BH")
    cat_levels <- c("significant (BH)", "nominal only", "not significant")
    df <- df %>%
      dplyr::mutate(
        category = dplyr::case_when(
          !is.na(.data$p_bh)  & .data$p_bh  <= alpha ~ "significant (BH)",
          !is.na(.data$p_raw) & .data$p_raw <= alpha  ~ "nominal only",
          TRUE                                         ~ "not significant"
        )
      )
    df$category <- factor(df$category, levels = cat_levels)
    if (length(stats::na.omit(unique(df$category))) <= 1L) {
      df$p_bin <- cut(df$p_bh,
        breaks         = c(0, 1e-3, 1e-2, 5e-2, 1, Inf),
        labels         = c("<1e-3", "(1e-3,1e-2]", "(1e-2,0.05]", ">0.05", "NA"),
        include.lowest = TRUE, right = TRUE
      )
      color_var       <- "p_bin"
      color_val_label <- "p-value (BH)"
    } else {
      color_var <- "category"
    }
  }

  # peptide-level coloring via color_by --------------------------------------
  if (!is.null(color_by)) {
    df        <- .build_color_group(df, color_by)$df
    color_var <- ".color_group"
    if (is.null(color_title)) color_title <- paste(names(color_by), collapse = " / ")
    color_val_label <- color_title
  }

  # jitter -------------------------------------------------------------------
  jw   <- as.numeric(jitter_width_pp  %||% 0)
  jh   <- as.numeric(jitter_height_pp %||% 0)
  pdat <- df
  if (!"percent1" %in% names(pdat) && "prop1" %in% names(pdat))
    pdat$percent1 <- pdat$prop1 * 100
  if (!"percent2" %in% names(pdat) && "prop2" %in% names(pdat))
    pdat$percent2 <- pdat$prop2 * 100
  if (jw > 0 || jh > 0) {
    set.seed(1L)
    n             <- nrow(pdat)
    pdat$percent1 <- pmin(100, pmax(0, pdat$percent1 + stats::rnorm(n, 0, jw)))
    pdat$percent2 <- pmin(100, pmax(0, pdat$percent2 + stats::rnorm(n, 0, jh)))
  }

  # build plot ---------------------------------------------------------------
  p <- ggplot2::ggplot(pdat, ggplot2::aes(x = percent1, y = percent2))

  if (identical(color_var, "category")) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = .data$category),
                          alpha = point_alpha, size = point_size) +
      ggplot2::scale_color_manual(
        values = c(
          "significant (BH)" = "#e31a1c",
          "nominal only"     = "#1b9e77",
          "not significant"  = "#386cb0"
        ),
        name = NULL, drop = FALSE
      )
  } else if (!is.null(color_var)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = .data[[color_var]]),
                          alpha = point_alpha, size = point_size) +
      ggplot2::labs(color = color_val_label)
  } else {
    p <- p +
      ggplot2::geom_point(alpha = point_alpha, size = point_size, color = "grey35")
  }

  rng <- range(c(pdat$percent1, pdat$percent2), na.rm = TRUE)
  p <- p +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         linewidth = 0.4, color = "grey60") +
    ggplot2::expand_limits(x = rng, y = rng) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(x = xlab, y = ylab) +
    theme_phip() +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      text = ggplot2::element_text(
        family = font_family %||% ggplot2::theme_get()$text$family,
        size   = font_size
      )
    )

  p
}

# ==============================================================================
#' Interactive prevalence scatter for prevalence results
#'
#' @description
#' Creates an interactive scatter (plotly) comparing prevalence in **group a**
#' vs **group b** for per-feature results produced by `ph_prevalence_compare()`.
#' Accepts a plain data.frame with columns `percent1`, `percent2`, `feature`,
#' `group1`, `group2`, `p_raw`, etc.
#'
#' When `pair` is given, subsetting uses `.ph_filter_pairs()`.
#'
#' Color mapping:
#' - Default (`color_by = NULL`): BH-corrected p-values computed per-plot from
#'   `p_raw` ("significant (BH)", "nominal only", "not significant").
#' - `color_by` as a named vector highlights points matching the specified
#'   peptide-library values:
#' \preformatted{
#'   color_by = c("is_flagellum" = TRUE)
#'   color_by = c("species" = "Staphylococcus aureus")
#'   color_by = c("is_flagellum" = TRUE, "species" = "Staphylococcus aureus")
#' }
#'
#' @param df A `data.frame` with columns `percent1`, `percent2`, `feature`,
#'   `group1`, `group2`, `p_raw`, optionally `n_peptides`, `rank`, `peptide_id`.
#' @param pair optional length-2 character, e.g. `c("kid_serum::T2","kid_serum::T8")`.
#' @param rank optional single rank (character) to keep.
#' @param xlab,ylab axis labels; defaults to `pair[1]`/`pair[2]` when `pair` is given.
#' @param alpha numeric in (0,1]; significance threshold; not the plotly alpha.
#' @param color_by optional named vector identifying peptide-library values to
#'   highlight, e.g. `c("is_flagellum" = TRUE)` or
#'   `c("species" = "Staphylococcus aureus")`.
#' @param color_title optional legend title when `color_by` is used.
#' @param peplib Optional peptide metadata table used to resolve `color_by`
#'   when not available via the global library.
#' @param background_df Optional data frame of background points.
#' @param ... graphical parameters: `category_colors`, `show_background`,
#'   `background_name`, `background_color`, `background_size`, `background_alpha`,
#'   `background_max_n`, `background_seed`, `point_line_width`, `point_line_color`,
#'   `point_size`, `point_alpha`, `jitter_width_pp`, `jitter_height_pp`,
#'   `font_family`, `font_size`.
#'
#' @return a `plotly` object.
#'
#' @examples
#' \dontrun{
#' p <- scatter_interactive(scatters,
#'   pair        = c("kid_serum::T2", "kid_serum::T8"),
#'   rank        = "peptide_id",
#'   color_by    = c("is_flagellum" = TRUE),
#'   color_title = "Flagellum"
#' )
#' p
#' }
#'
#' @export
scatter_interactive <- function(df,
                                pair          = NULL,
                                rank          = NULL,
                                xlab          = NULL,
                                ylab          = NULL,
                                alpha         = 0.05,
                                color_by      = NULL,
                                color_title   = NULL,
                                peplib        = NULL,
                                background_df = NULL,
                                ...) {
  dots             <- list(...)
  category_colors  <- dots[["category_colors"]] %||% c(
    "significant (BH)" = "#009E73",
    "nominal only"     = "#E69F00",
    "not significant"  = "#0072B2",
    "other/NA"         = "#999999"
  )
  show_background   <- dots[["show_background"]]   %||% FALSE
  background_name   <- dots[["background_name"]]   %||% "background peptides"
  background_color  <- dots[["background_color"]]  %||% "#808080"
  background_size   <- dots[["background_size"]]   %||% 6
  background_alpha  <- dots[["background_alpha"]]  %||% 0.30
  background_max_n  <- dots[["background_max_n"]]  %||% 3000
  background_seed   <- dots[["background_seed"]]   %||% 1L
  point_size        <- dots[["point_size"]]        %||% 7
  point_alpha       <- dots[["point_alpha"]]       %||% 0.85
  point_line_width  <- dots[["point_line_width"]]  %||% 0.7
  point_line_color  <- dots[["point_line_color"]]  %||% "rgba(0,0,0,0.6)"
  jitter_width_pp   <- dots[["jitter_width_pp"]]   %||% 0.0
  jitter_height_pp  <- dots[["jitter_height_pp"]]  %||% 0.0
  font_family       <- dots[["font_family"]]       %||% NULL
  font_size         <- dots[["font_size"]]         %||% 12

  clamp01     <- function(x) pmin(1, pmax(0, as.numeric(x)))
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

  # filtering ----------------------------------------------------------------
  df <- as.data.frame(df)
  if (!is.null(pair)) {
    df <- .ph_filter_pairs(
      df,
      gA = pair[1], gB = pair[2],
      ranks   = rank,
      drop_na = TRUE
    )
  } else if (!is.null(rank) && "rank" %in% names(df)) {
    df <- df[df$rank %in% rank, , drop = FALSE]
  }

  if (is.null(df) || !nrow(df)) {
    p <- plotly::plot_ly()
    return(plotly::layout(
      p,
      annotations = list(
        x = 0.5, y = 0.5, text = "no data for this contrast",
        showarrow = FALSE, xref = "paper", yref = "paper"
      ),
      xaxis  = list(title = if (is.null(xlab)) "group a" else xlab, zeroline = FALSE),
      yaxis  = list(title = if (is.null(ylab)) "group b" else ylab, zeroline = FALSE),
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

  # default significance coloring (BH per-plot) ------------------------------
  cat_levels <- c("significant (BH)", "nominal only", "not significant")
  color_var       <- NULL
  color_val_label <- NULL
  if (is.null(color_by)) {
    df$p_bh <- stats::p.adjust(df$p_raw, method = "BH")
    df <- df %>%
      dplyr::mutate(
        category = dplyr::case_when(
          !is.na(.data$p_bh)  & .data$p_bh  <= alpha ~ cat_levels[1],
          !is.na(.data$p_raw) & .data$p_raw <= alpha  ~ cat_levels[2],
          TRUE                                         ~ cat_levels[3]
        ),
        category = factor(category, levels = c(cat_levels, "other/NA"))
      )
  }

  # peptide-level coloring via color_by --------------------------------------
  if (!is.null(color_by)) {
    df        <- .build_color_group(df, color_by, peplib = peplib)$df
    color_var <- ".color_group"
    if (is.null(color_title)) color_title <- paste(names(color_by), collapse = " / ")
    color_val_label <- color_title
  }

  # jitter + percent scale ---------------------------------------------------
  jw   <- as.numeric(jitter_width_pp  %||% 0)
  jh   <- as.numeric(jitter_height_pp %||% 0)
  pdat <- df
  if (!"percent1" %in% names(pdat) && "prop1" %in% names(pdat))
    pdat$percent1 <- pdat$prop1 * 100
  if (!"percent2" %in% names(pdat) && "prop2" %in% names(pdat))
    pdat$percent2 <- pdat$prop2 * 100
  if (jw > 0 || jh > 0) {
    set.seed(1L)
    n             <- nrow(pdat)
    pdat$percent1 <- pmin(100, pmax(0, pdat$percent1 + stats::rnorm(n, 0, jw)))
    pdat$percent2 <- pmin(100, pmax(0, pdat$percent2 + stats::rnorm(n, 0, jh)))
  }

  # join peplib for hover metadata (species, common, Fullname, Description) ---
  pdat$pep_id <- NA_character_
  if (all(c("rank", "feature") %in% names(pdat))) {
    pdat$pep_id <- ifelse(pdat$rank == "peptide_id",
                          as.character(pdat$feature), NA_character_)
  } else if ("peptide_id" %in% names(pdat)) {
    pdat$pep_id <- as.character(pdat$peptide_id)
  }

  # ensure hover metadata columns always exist (filled in by peplib join if available)
  if (!"species"      %in% names(pdat)) pdat$species      <- NA_character_
  if (!"common"       %in% names(pdat)) pdat$common       <- NA_character_
  if (!"Fullname"     %in% names(pdat)) pdat$Fullname     <- NA_character_
  if (!"Description"  %in% names(pdat)) pdat$Description  <- NA_character_

  if (any(!is.na(pdat$pep_id))) {
    meta_src <- if (!is.null(peplib)) peplib else
      tryCatch(get_peptide_library(), error = function(...) NULL)
    if (!is.null(meta_src)) {
      pm2 <- meta_src %>%
        dplyr::select(dplyr::any_of(
          c("peptide_id", "species", "common", "Fullname", "Description"))) %>%
        dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
        dplyr::collect()
      meta_cols <- intersect(
        c("species", "common", "Fullname", "Description"),
        setdiff(names(pm2), "peptide_id")
      )
      if (length(meta_cols) > 0L) {
        pdat <- pdat %>%
          dplyr::select(-tidyselect::all_of(meta_cols)) %>%
          dplyr::left_join(pm2, by = c("pep_id" = "peptide_id"))
        for (col in setdiff(c("species", "common", "Fullname", "Description"), names(pdat)))
          pdat[[col]] <- NA_character_
      }
    }
  }

  # hover text ---------------------------------------------------------------
  pdat <- pdat %>%
    dplyr::mutate(
      p1r        = round(.data$percent1, 2),
      p2r        = round(.data$percent2, 2),
      rr         = dplyr::coalesce(round(
        if (!"ratio" %in% names(.) & all(c("prop1", "prop2") %in% names(.)))
          (pmax(prop1, 1e-12) / pmax(prop2, 1e-12)) else ratio, 2), NA_real_),
      praw       = if ("p_raw" %in% names(.)) round(.data$p_raw, 3) else NA_real_,
      pbh        = if ("p_bh"  %in% names(.)) round(.data$p_bh,  3) else NA_real_,
      npep       = dplyr::coalesce(
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
        "<b>%s</b>%s<br>%s: %d/%d (%.2f%%)<br>%s: %d/%d (%.2f%%)%s<br>peptides: %s<br>p: %.3f<br>p_adj (BH): %.3f",
        .data$feature,
        .data$meta_block,
        xlab, .data$n1, .data$N1, .data$p1r,
        ylab, .data$n2, .data$N2, .data$p2r,
        if (!is.null(color_var))
          sprintf("<br>%s: %s", color_val_label, .data$color_info)
        else "",
        format(.data$npep, big.mark = "\u2009", scientific = FALSE),
        .data$praw, .data$pbh
      )
    )

  pdat$text <- as.character(pdat$text)

  # background overlay prep --------------------------------------------------
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

  # build plot: background → diagonal → main points -------------------------
  p <- plotly::plot_ly()

  if (!is.null(bg) && nrow(bg) > 0) {
    p <- plotly::add_trace(
      p,
      data        = bg,
      x           = ~percent1, y = ~percent2,
      type        = "scatter", mode = "markers",
      name        = background_name,
      showlegend  = TRUE,
      hoverinfo   = "skip",
      inherit     = FALSE,
      legendrank  = 999,
      marker      = list(
        size    = background_size,
        color   = background_color,
        opacity = clamp01(background_alpha)
      )
    )
  }

  rng <- range(c(pdat$percent1, pdat$percent2), na.rm = TRUE)
  p   <- plotly::add_trace(
    p,
    x = rng, y = rng,
    type       = "scatter", mode = "lines",
    inherit    = FALSE,
    showlegend = FALSE,
    hoverinfo  = "skip"
  )

  if (is.null(color_var)) {
    pdat$category <- droplevels(pdat$category)
    levs          <- levels(pdat$category)

    default_pal <- c(
      "significant (BH)" = "#009E73",
      "nominal only"     = "#E69F00",
      "not significant"  = "#0072B2",
      "other/NA"         = "#999999"
    )
    pal_map <- default_pal
    if (!is.null(category_colors) && length(category_colors) > 0)
      pal_map[names(category_colors)] <- category_colors
    missing_cols <- setdiff(levs, names(pal_map))
    if (length(missing_cols) > 0) pal_map[missing_cols] <- "#999999"

    for (lv in levs) {
      dd <- pdat[pdat$category == lv, , drop = FALSE]
      if (!nrow(dd)) next
      p <- plotly::add_trace(
        p,
        data             = dd,
        x                = ~percent1, y = ~percent2,
        type             = "scatter", mode = "markers",
        name             = lv,
        showlegend       = TRUE,
        inherit          = FALSE,
        text             = ~text,
        hovertemplate    = "%{text}<extra></extra>",
        marker           = list(
          size    = point_size,
          opacity = point_alpha,
          color   = unname(pal_map[lv]),
          line    = list(width = point_line_width, color = point_line_color)
        )
      )
    }
  } else {
    # Per-level loop (same pattern as category case) to ensure ~text resolves
    # correctly inside each trace — using color = ~formula causes plotly to
    # split traces internally and the %{text} reference becomes unreliable.
    levs_cb  <- levels(pdat[[color_var]])
    n_groups <- sum(levs_cb != "other")
    cb_pal   <- if (n_groups > 0L)
      stats::setNames(
        grDevices::hcl.colors(max(3L, n_groups), palette = "Dark 3")[seq_len(n_groups)],
        levs_cb[levs_cb != "other"]
      ) else character(0)
    if ("other" %in% levs_cb) cb_pal["other"] <- "#aaaaaa"

    for (lv in levs_cb) {
      dd <- pdat[as.character(pdat[[color_var]]) == lv, , drop = FALSE]
      if (!nrow(dd)) next
      p <- plotly::add_trace(
        p,
        data          = dd,
        x             = ~percent1, y = ~percent2,
        type          = "scatter", mode = "markers",
        name          = lv,
        showlegend    = TRUE,
        inherit       = FALSE,
        text          = ~text,
        hovertemplate = "%{text}<extra></extra>",
        marker        = list(
          size    = point_size,
          opacity = point_alpha,
          color   = unname(cb_pal[lv]),
          line    = list(width = point_line_width, color = point_line_color)
        )
      )
    }
  }

  plotly::layout(
    p,
    font   = list(family = font_family, size = font_size),
    xaxis  = list(title = xlab, zeroline = FALSE, range = c(0, 100)),
    yaxis  = list(title = ylab, zeroline = FALSE, range = c(0, 100)),
    legend = list(orientation = "h", y = -0.15)
  )
}

# ==============================================================================
# Internal helper to prepare volcano data (shared by static/interactive)
# ==============================================================================
.volcano_prepare <- function(df,
                             pair = NULL,
                             rank = NULL,
                             color_by = NULL,
                             color_title = NULL,
                             fc_cut = 1,
                             p_cut = 0.05,
                             p_mode = c("raw", "bh"),
                             significant_colors = c(
                               "not significant"                 = "#386cb0",
                               "significant prior correction"    = "#1b9e77",
                               "significant post fdr correction" = "#e31a1c"
                             ),
                             peplib = NULL) {
  p_mode <- match.arg(p_mode)

  # filtering ----------------------------------------------------------------
  df <- as.data.frame(df)
  if (!is.null(pair)) {
    df <- .ph_filter_pairs(
      df,
      gA = pair[1], gB = pair[2],
      ranks   = rank,
      drop_na = TRUE
    )
  } else if (!is.null(rank) && "rank" %in% names(df)) {
    df <- df[df$rank %in% rank, , drop = FALSE]
  }

  if (is.null(df) || !nrow(df)) return(list(empty = TRUE))

  req  <- c("rank", "feature", "group1", "group2")
  miss <- setdiff(req, names(df))
  if (length(miss))
    stop("volcano_*(): missing required columns: ", paste(miss, collapse = ", "))

  df$rank      <- as.character(df$rank)
  df$feature   <- as.character(df$feature)
  df$group1    <- as.character(df$group1)
  df$group2    <- as.character(df$group2)
  df$pair_label <- paste0(df$group1, " vs ", df$group2)

  # ensure ratio
  if (!"ratio" %in% names(df)) {
    if (!all(c("prop1", "prop2") %in% names(df)))
      stop("volcano_*(): need either 'ratio' or both 'prop1' and 'prop2'.")
    p1       <- ifelse(is.na(df$prop1) | df$prop1 <= 0, df$prop1 + 1e-12, df$prop1)
    p2       <- ifelse(is.na(df$prop2) | df$prop2 <= 0, df$prop2 + 1e-12, df$prop2)
    df$ratio <- p1 / p2
  }

  # p column: raw or BH computed per-plot
  if (p_mode == "bh") {
    if (!"p_raw" %in% names(df))
      stop("volcano_*(): 'p_raw' column required when p_mode='bh'.")
    df$p_use <- stats::p.adjust(df$p_raw, method = "BH")
  } else {
    if (!"p_raw" %in% names(df))
      stop("volcano_*(): 'p_raw' column not found.")
    df$p_use <- df$p_raw
  }

  df$log2ratio <- log2(df$ratio)
  p_use_vec    <- pmax(df$p_use, .Machine$double.xmin)
  df$nlog10p   <- -log10(p_use_vec)

  if (!any(is.finite(df$log2ratio)) || !any(is.finite(df$nlog10p)))
    stop("volcano_*(): Non-finite values for log2(ratio) or -log10(p).")

  # categories
  rawp   <- if ("p_raw" %in% names(df)) df$p_raw else df$p_use
  sigfc  <- abs(df$log2ratio) >= fc_cut
  post   <- (df$p_use <= p_cut) & sigfc
  prior  <- (rawp     <= p_cut) & sigfc & !post
  df$category_use <- factor(
    ifelse(post,  "significant post fdr correction",
    ifelse(prior, "significant prior correction", "not significant")),
    levels = c("significant post fdr correction",
               "significant prior correction",
               "not significant")
  )

  # peptide-level coloring via color_by
  color_var   <- NULL
  color_label <- NULL
  if (!is.null(color_by)) {
    df          <- .build_color_group(df, color_by, peplib = peplib)$df
    color_var   <- ".color_group"
    color_label <- color_title %||% paste(names(color_by), collapse = " / ")
  }

  list(
    empty              = FALSE,
    data               = df,
    color_var          = color_var,
    color_label        = color_label,
    significant_colors = significant_colors,
    fc_cut             = fc_cut,
    p_cut              = p_cut,
    p_mode             = p_mode
  )
}

# ==============================================================================
#' Static volcano plot (log2 ratio vs -log10 p)
#'
#' @param df A data frame with prevalence results.
#' @param pair optional group pair (character length-2).
#' @param rank optional single rank (character) to keep.
#' @param color_by optional named vector identifying peptide-library values to
#'   highlight, e.g. `c("is_flagellum" = TRUE)`.
#' @param color_title optional legend title when `color_by` is used.
#' @param fc_cut Numeric; absolute log2 fold-change cutoff.
#' @param p_cut Numeric; p-value cutoff.
#' @param p_mode One of `c("raw","bh")`; `"bh"` applies BH correction per-plot
#'   from `p_raw`.
#' @param significant_colors Named vector of colors for significance categories.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' set.seed(2)
#' prev <- data.frame(
#'   rank       = "peptide_id",
#'   feature    = paste0("pep", 1:40),
#'   group1     = "A",
#'   group2     = "B",
#'   prop1      = runif(40),
#'   prop2      = runif(40),
#'   percent1   = runif(40, 0, 100),
#'   percent2   = runif(40, 0, 100),
#'   ratio      = runif(40, 0.1, 10),
#'   p_raw      = c(runif(10, 0, 0.01), runif(30, 0.1, 1)),
#'   n_peptides = 1L
#' )
#'
#' # basic volcano
#' volcano_static(prev)
#'
#' # BH correction, custom cutoffs
#' volcano_static(prev, fc_cut = 1.5, p_cut = 0.01, p_mode = "bh")
#'
#' # filter to one pair
#' volcano_static(prev, pair = c("A", "B"), rank = "peptide_id")
#' @export
volcano_static <- function(df,
                           pair               = NULL,
                           rank               = NULL,
                           color_by           = NULL,
                           color_title        = NULL,
                           fc_cut             = 1,
                           p_cut              = 0.05,
                           p_mode             = c("raw", "bh"),
                           significant_colors = c(
                             "not significant"                 = "#386cb0",
                             "significant prior correction"    = "#1b9e77",
                             "significant post fdr correction" = "#e31a1c"
                           )) {
  prep <- .volcano_prepare(
    df, pair, rank,
    color_by, color_title, fc_cut, p_cut, p_mode, significant_colors
  )
  if (isTRUE(prep$empty)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "no data for this contrast") +
        ggplot2::theme_void()
    )
  }
  dd <- prep$data

  p <- if (!is.null(prep$color_var)) {
    ggplot2::ggplot(dd,
      ggplot2::aes(x = log2ratio, y = nlog10p,
                   color = .data[[prep$color_var]])) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::labs(color = prep$color_label)
  } else {
    ggplot2::ggplot(dd,
      ggplot2::aes(x = log2ratio, y = nlog10p, color = category_use)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::scale_color_manual(values = prep$significant_colors, name = NULL)
  }

  p +
    ggplot2::geom_hline(yintercept = -log10(prep$p_cut),
                        linetype = "dashed", color = "gray50") +
    ggplot2::geom_vline(xintercept = c(-prep$fc_cut, 0, prep$fc_cut),
                        linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      x = "log\u2082 ratio (group1 / group2)",
      y = "-log\u2081\u2080(p)"
    ) +
    theme_phip(base_size = 12) +
    ggplot2::theme(
      legend.position = if (!is.null(prep$color_var)) "bottom" else "none"
    )
}

# ==============================================================================
#' Interactive volcano plot (log2 ratio vs -log10 p)
#'
#' @inheritParams volcano_static
#' @return A `plotly` htmlwidget.
#'
#' @examples
#' \dontrun{
#' set.seed(3)
#' prev <- data.frame(
#'   rank       = "peptide_id",
#'   feature    = paste0("pep", 1:40),
#'   group1     = "A",
#'   group2     = "B",
#'   prop1      = runif(40),
#'   prop2      = runif(40),
#'   percent1   = runif(40, 0, 100),
#'   percent2   = runif(40, 0, 100),
#'   ratio      = runif(40, 0.1, 10),
#'   p_raw      = c(runif(10, 0, 0.01), runif(30, 0.1, 1)),
#'   n_peptides = 1L
#' )
#'
#' # interactive volcano — hover to inspect each peptide
#' volcano_interactive(prev)
#'
#' # BH correction
#' volcano_interactive(prev, p_mode = "bh", fc_cut = 1.5, p_cut = 0.01)
#' }
#' @export
volcano_interactive <- function(df,
                                pair               = NULL,
                                rank               = NULL,
                                color_by           = NULL,
                                color_title        = NULL,
                                fc_cut             = 1,
                                p_cut              = 0.05,
                                p_mode             = c("raw", "bh"),
                                significant_colors = c(
                                  "not significant"                 = "#386cb0",
                                  "significant prior correction"    = "#1b9e77",
                                  "significant post fdr correction" = "#e31a1c"
                                )) {
  prep <- .volcano_prepare(
    df, pair, rank,
    color_by, color_title, fc_cut, p_cut, p_mode, significant_colors
  )
  if (isTRUE(prep$empty)) {
    p <- plotly::plot_ly()
    return(plotly::layout(
      p,
      annotations = list(x = 0.5, y = 0.5, text = "no data for this contrast",
                         showarrow = FALSE, xref = "paper", yref = "paper"),
      xaxis  = list(title = "log\u2082 ratio (group1 / group2)", zeroline = FALSE),
      yaxis  = list(title = "-log\u2081\u2080(p)", zeroline = FALSE),
      legend = list(orientation = "h", y = -0.15)
    ))
  }

  dd  <- prep$data
  hov <- sprintf(
    "<b>%s</b><br>rank: %s<br>pair: %s vs %s<br>log2(ratio): %.3f<br>-log10(p): %.3f%s",
    dd$feature, dd$rank, dd$group1, dd$group2, dd$log2ratio, dd$nlog10p,
    if (!is.null(prep$color_var))
      sprintf("<br>%s: %s", prep$color_label, as.character(dd[[prep$color_var]]))
    else ""
  )

  if (!is.null(prep$color_var)) {
    pal     <- grDevices::hcl.colors(
      n       = max(3L, length(levels(dd[[prep$color_var]]))),
      palette = "Dark 3"
    )
    pal     <- stats::setNames(pal[seq_along(levels(dd[[prep$color_var]]))],
                               levels(dd[[prep$color_var]]))
    col_vec <- pal[as.character(dd[[prep$color_var]])]
  } else {
    pal     <- prep$significant_colors
    cat_vec <- as.character(dd$category_use)
    col_vec <- pal[cat_vec]
    col_vec[is.na(col_vec)] <- pal[["not significant"]]
  }

  p <- plotly::plot_ly(
    dd,
    x = ~log2ratio, y = ~nlog10p,
    type          = "scatter", mode = "markers",
    text          = hov,
    hovertemplate = "%{text}<extra></extra>",
    marker        = list(size = 7, opacity = 0.85, color = col_vec)
  )

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
  if (!is.null(prep$color_var) && !is.null(prep$color_label))
    leg$title <- list(text = prep$color_label)

  plotly::layout(
    p,
    xaxis      = list(title = "log\u2082 ratio (group1 / group2)", zeroline = FALSE),
    yaxis      = list(title = "-log\u2081\u2080(p)", zeroline = FALSE),
    legend     = leg,
    showlegend = !is.null(prep$color_var)
  )
}
