#' @title Fold‑change Q–Q plot for `phip_data`
#'
#' @description A convenience wrapper around **ggplot2**'s `stat_qq()` / `stat_qq_line()`
#' that produces a ready‑to‑modify Q–Q plot of the *fold_change* column in a
#' `phip_data` object.
#'
#' @param pd        A `phip_data` object.
#' @param by        Optional character vector of column names used for grouping.
#'                  The first element is mapped to `colour` / `fill`; users can
#'                  facet on any columns afterwards via the usual `+ facet_*()`.
#' @param prop      Proportion of rows (per group) to sample **before** plotting
#'                  (numeric, `0 < prop <= 1`).  Defaults to `1` (all rows).
#' @param distribution,dparams Passed straight to `stat_qq()` and
#'                  `stat_qq_line()` for alternative theoretical quantiles.
#' @param xlim      Two‑element numeric vector passed to `coord_cartesian(xlim
#'                  = xlim)`.  Leave `NULL` to auto‑scale.
#'
#' @return A `ggplot` object ready for further layers / faceting.
#' @examples
#' library(ggplot2)
#' library(ggpubr)
#' pd <- phip_load_example_data()
#' p = fold_qq_plot(pd, by = c("group", "timepoint"), prop = 0.1) +
#'   facet_wrap(~group, nrow = 2) +
#'   ggpubr::theme_pubr()
#' @export
fold_qq_plot <- function(pd,
                         by = NULL,
                         prop = 1,
                         distribution = stats::qnorm,
                         dparams = list(),
                         xlim = NULL) {
  # ––– sanity checks ---------------------------------------------------------
  .chk_cond(
    !inherits(pd, "phip_data"),
    "`pd` has to be an object of class `phip_data`!"
  )

  df <- pd$data_long
  vars <- dplyr::tbl_vars(df)

  .chk_cond(
    !"fold_change" %in% vars,
    "Column 'fold_change' not found in pd$data_long."
  )

  .chk_cond(
    !is.null(by) && !all(by %in% vars),
    sprintf(
      "Grouping column(s) not found: %s",
      paste(setdiff(by, vars), collapse = ", ")
    )
  )

  .chk_cond(
    !is.numeric(prop) || prop <= 0 || prop > 1,
    "`prop` must be a number in (0, 1]."
  )

  # ––– grouping --------------------------------------------------------------
  if (!is.null(by) && length(by)) {
    df <- dplyr::group_by(df, dplyr::across(dplyr::all_of(by)))
  }

  # ––– optional down‑sampling ------------------------------------------------
  if (prop < 1) {
    is_remote <- inherits(df, "tbl_sql") || inherits(df, "tbl_arrow")

    if (is_remote) {
      # Collect minimal columns locally, then sample -------------------------
      need_cols <- unique(c(by, "fold_change"))
      df <- df %>%
        dplyr::select(dplyr::all_of(need_cols)) %>%
        dplyr::collect()

      if (!is.null(by) && length(by)) {
        df <- dplyr::group_by(df, dplyr::across(dplyr::all_of(by)))
      }

      df <- dplyr::slice_sample(df, prop = prop)
    } else {
      df <- dplyr::slice_sample(df, prop = prop)
    }
  }

  df <- dplyr::ungroup(df)

  # ––– build aes mapping -----------------------------------------------------
  aes_list <- list(sample = rlang::expr(fold_change))
  if (!is.null(by) && length(by)) {
    aes_list$colour <- rlang::sym(by[[1]])
    aes_list$fill <- rlang::sym(by[[1]])
  }

  #--- calculate meaningful breaks -------------------------------------------

  # ––– plot ------------------------------------------------------------------
  p <- ggplot2::ggplot(df, ggplot2::aes(!!!aes_list)) +
    ggplot2::stat_qq(distribution = distribution, dparams = dparams) +
    ggplot2::stat_qq_line(
      distribution = distribution, dparams = dparams,
      size = 1.3
    ) +
    ggplot2::xlab("Theoretical Quantiles") +
    ggplot2::ylab("Sample fold_change") +
    ggplot2::scale_y_continuous(
      breaks = function(lims) pretty(lims, n = 10), # 10-ish “round” breaks
      labels = function(x) format(x, scientific = FALSE) # plain numbers
    ) +
    theme_phip()

  # apply PHIP palette automatically when grouping was provided -------------
  if (!is.null(by) && length(by)) {
    p <- p + scale_colour_phip()
  }

  if (!is.null(xlim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }

  return(p)
}


#' @title Box‑Cox‑transform the *fold_change* column of a `phip_data`
#'
#' @description Applies a Box‑Cox power transformation to `fold_change`, either
#'   globally or separately within groups.  The transformed values overwrite the
#'   original column **in place** (use `confirm = FALSE` to suppress the
#'   interactive prompt).
#'
#' @section How it works:
#' * **Backend‑agnostic** – works on in‑memory, DuckDB, Arrow, … sources.
#' * For remote tables the function pulls just the grouping cols +
#'   `fold_change` into R to estimate an optimal \eqn{\lambda} via
#'   `MASS::boxcox()`, then joins the lambdas back and mutates remotely.
#' * If \eqn{\lambda = 0} the transformation falls back to `log(x)`.
#'
#' @param pd      A `phip_data` object.
#' @param by      Optional character vector of grouping columns.  If `NULL`
#'                (default) the whole data set shares one Box‑Cox \eqn{\lambda}.
#' @param confirm Logical.  When `TRUE` (default) *and* interactive, prints a
#'                warning that `fold_change` will be overwritten and asks the
#'                user to continue.
#'
#' @return The same `phip_data`, with `fold_change` replaced by its Box‑Cox
#'         transform and the attribute `boxcox_lambda` containing the per‑group
#'         lambdas (tibble).
#'
#' @examples
#' library(tidyr)
#' pd <- phip_load_example_data()
#' pd2 <- transform_fold_boxcox(pd, by = "group", confirm=FALSE)
#' @export
#' @importFrom MASS boxcox
transform_fold_boxcox <- function(pd, by = NULL, confirm = interactive()) {
  .chk_cond(
    !inherits(pd, "phip_data"),
    "`pd` has to be an object of class `phip_data`!"
  )

  df <- pd$data_long
  vars <- dplyr::tbl_vars(df)

  .chk_cond(
    !"fold_change" %in% vars,
    "Column 'fold_change' not found in pd$data_long."
  )

  .chk_cond(
    !is.null(by) && !all(by %in% vars),
    sprintf(
      "Grouping column(s) not found: %s",
      paste(setdiff(by, vars), collapse = ", ")
    )
  )

  # ––– user confirmation -----------------------------------------------------
  if (confirm && interactive()) {
    warn <- if (requireNamespace("cli", quietly = TRUE)) cli::cli_alert_warning else warning
    warn("`fold_change` will be **overwritten in place** by its Box‑Cox transform.")
    ok <- .cli_yesno("Proceed? (y/n) ")
    if (!isTRUE(ok)) {
      chk::chk_abort("Operation aborted by user.", call. = FALSE)
    }
  }

  grp_cols <- by %||% character(0)

  # ––– estimate lambdas locally ---------------------------------------------
  rlang::check_installed("MASS")

  lambda_tbl <- df %>%
    dplyr::select(dplyr::all_of(c(grp_cols, "fold_change"))) %>%
    dplyr::collect() %>%
    {
      if (length(grp_cols)) dplyr::group_by(., dplyr::across(dplyr::all_of(grp_cols))) else .
    } %>%
    dplyr::summarise(
      lambda = {
        bc <- MASS::boxcox(fold_change ~ 1, plot = FALSE)
        bc$x[which.max(bc$y)]
      },
      .groups = "drop"
    )

  # ––– apply transformation --------------------------------------------------
  df2 <- df %>%
    dplyr::left_join(lambda_tbl, by = grp_cols, copy = TRUE) %>%
    dplyr::mutate(
      fold_change = dplyr::if_else(
        lambda == 0,
        log(fold_change),
        (fold_change^lambda - 1) / lambda
      )
    ) %>%
    dplyr::select(-lambda)

  out <- .modify_pd(pd, df2)
  attr(out, "boxcox_lambda") <- lambda_tbl
  out
}

#' Histogram of *fold_change* for `phip_data`
#'
#' Generates a histogram of the `fold_change` column, with optional per‑group
#' colouring/filling, sampling, and remote‑table friendliness.
#'
#' @inheritParams fold_qq_plot
#' @param bins Number of histogram bins (default = 100).
#' @param xlim Optional numeric vector of length‑2 for x‑axis limits.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' library(ggplot2)
#' pd <- phip_load_example_data()
#' p = fold_hist_plot(pd, by = "group", prop = 0.1, bins = 100) +
#'   facet_wrap(~group, nrow = 2)
#' }
fold_hist_plot <- function(pd,
                           by = NULL,
                           prop = 1,
                           bins = 50,
                           xlim = NULL) {
  stopifnot(inherits(pd, "phip_data"))

  df <- pd$data_long
  vars <- dplyr::tbl_vars(df)

  .chk_cond(
    !"fold_change" %in% vars,
    "Column 'fold_change' not found in pd$data_long."
  )

  .chk_cond(
    !is.null(by) && !all(by %in% vars),
    sprintf(
      "Grouping column(s) not found: %s",
      paste(setdiff(by, vars), collapse = ", ")
    )
  )

  .chk_cond(
    !is.numeric(prop) || prop <= 0 || prop > 1,
    "`prop` must be a number in (0, 1]."
  )

  # group remote sampling similar to fold_qq_plot ---------------------------------
  if (!is.null(by) && length(by)) {
    df <- dplyr::group_by(df, dplyr::across(dplyr::all_of(by)))
  }

  if (prop < 1) {
    is_remote <- inherits(df, "tbl_sql") || inherits(df, "tbl_arrow")
    if (is_remote) {
      need_cols <- unique(c(by, "fold_change"))
      df <- df %>%
        dplyr::select(dplyr::all_of(need_cols)) %>%
        dplyr::collect()
      if (!is.null(by) && length(by)) {
        df <- dplyr::group_by(df, dplyr::across(dplyr::all_of(by)))
      }
      df <- dplyr::slice_sample(df, prop = prop)
    } else {
      df <- dplyr::slice_sample(df, prop = prop)
    }
  }

  df <- dplyr::ungroup(df)

  aes_list <- list(x = rlang::expr(fold_change))
  if (!is.null(by) && length(by)) {
    aes_list$fill <- rlang::sym(by[[1]])
    aes_list$colour <- rlang::sym(by[[1]])
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(!!!aes_list)) +
    ggplot2::geom_histogram(
      bins = bins,
      alpha = 0.4,
      colour = "black",
      position = "identity"
    ) +
    # ggplot2::xlab("Theoretical Quantiles") +
    # ggplot2::ylab("Sample fold_change") +
    ggplot2::scale_x_continuous(
      breaks = function(lims) pretty(lims, n = 7), # 10-ish “round” breaks
      labels = function(x) format(x, scientific = FALSE) # plain numbers
    ) +
    theme_phip()

  if (!is.null(by) && length(by)) {
    p <- p + scale_fill_phip()
  }

  if (!is.null(xlim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }

  p
}
