# ==============================================================================
# compute_delta2(): input adapters for the delta plotters
# ==============================================================================
# deltaplot() / ecdf_plot() read a prevalence table (group1, group2, prop1,
# prop2) and forestplot() a results table (rank, feature, group1, group2,
# design, T_obs, p_perm). These turn a delta2_result into those shapes, so
# the plotters themselves stay as they are.

# The group_col of a delta2_result to plot: the only one, or the one named
.ph_delta2_pick_group_col <- function(x, group_col = NULL, available) {
  if (!length(available)) {
    .ph_abort("The delta2_result has no rows to plot.")
  }
  if (is.null(group_col)) {
    if (length(available) > 1L) {
      .ph_abort(
        "The delta2_result has several group_cols; pick one with group_col.",
        bullets = paste("-", available)
      )
    }
    return(available)
  }
  chk::chk_string(group_col)
  if (!group_col %in% available) {
    .ph_abort(sprintf("group_col %s is not in the delta2_result.", group_col),
      bullets = paste("-", available)
    )
  }
  group_col
}

#' @title Internal helper: .ph_delta2_prev_tbl
#' @description Per-peptide prevalence of a `delta2_result` in the shape of
#'   [compute_pop()]'s table: one row per peptide and pair of levels, with
#'   `prop1` / `prop2` the prevalences at `group1` / `group2`.
#' @param x A `delta2_result` with `$peptides`.
#' @param term `"group"` (levels of the group column) or `"time"`
#'   (timepoints).
#' @param group_col The group column; may be `NULL` when there is one.
#' @param adjusted Logical; use `prevalence_adj` (inverse-propensity
#'   weighted group levels).
#' @return A tibble with `feature`, `peptide_id`, `group1`, `group2`,
#'   `prop1`, `prop2`, `n1`, `n2`.
#' @keywords internal
.ph_delta2_prev_tbl <- function(x, term = c("group", "time"),
                                group_col = NULL, adjusted = FALSE) {
  term <- match.arg(term)
  chk::chk_flag(adjusted)
  pep <- x$peptides
  if (is.null(pep)) {
    .ph_abort("The delta2_result has no $peptides table; rerun
      compute_delta2() with return_peptides = TRUE.")
  }
  gc <- .ph_delta2_pick_group_col(x, group_col, unique(pep$group_col))
  pep <- pep[pep$group_col == gc & pep$factor == term, , drop = FALSE]
  if (!nrow(pep)) {
    .ph_abort(sprintf("No %s levels to plot for group_col %s.", term, gc))
  }
  meta <- x$meta$group_cols[[gc]]$design
  levels <- as.character(if (term == "group") {
    meta$group_levels
  } else {
    meta$time_levels
  })
  levels <- intersect(levels, unique(pep$level))
  if (length(levels) < 2L) {
    .ph_abort(sprintf("Need two %s levels to plot.", term))
  }
  val <- if (adjusted) pep$prevalence_adj else pep$prevalence
  pairs <- utils::combn(levels, 2L)
  out <- lapply(seq_len(ncol(pairs)), function(k) {
    a <- pep$level == pairs[1, k]
    b <- pep$level == pairs[2, k]
    ib <- match(pep$peptide_id[a], pep$peptide_id[b])
    tibble::tibble(
      feature = pep$peptide_id[a],
      peptide_id = pep$peptide_id[a],
      group1 = pairs[1, k],
      group2 = pairs[2, k],
      prop1 = val[a],
      prop2 = val[b][ib],
      n1 = pep$n[a],
      n2 = pep$n[b][ib]
    )
  })
  dplyr::bind_rows(out)
}

#' @title Internal helper: .ph_delta2_forest_tbl
#' @description One term of a `delta2_result` in the shape of
#'   [compute_delta()]'s result table for [forestplot()].
#' @param x A `delta2_result`.
#' @param term `"group"`, `"time"`, `"interaction"` or `"pairwise"`. The
#'   three terms must have one degree of freedom (two levels or scores):
#'   an omnibus `T` is a sum of squares and has no sign to plot.
#' @param group_col The group column; may be `NULL` when there is one.
#' @param group_pair_values For `term = "pairwise"`, the pair `c(group1,
#'   group2)`; may be `NULL` when there is one pair.
#' @return The results table (`rank`, `feature`, `group1`, `group2`,
#'   `design`, `T_obs`, `p_perm`, `T_obs_stand`, `Z_from_p`, plus
#'   `n_peptides_used` and, for pairs, `p_maxT`), with attribute `labels`:
#'   the default left and right arrow labels.
#' @keywords internal
.ph_delta2_forest_tbl <- function(x,
                                  term = c(
                                    "group", "time", "interaction",
                                    "pairwise"
                                  ),
                                  group_col = NULL,
                                  group_pair_values = NULL) {
  term <- match.arg(term)
  if (term == "pairwise") {
    pw <- x$pairwise
    gc <- .ph_delta2_pick_group_col(x, group_col, unique(pw$group_col))
    pw <- pw[pw$group_col == gc, , drop = FALSE]
    pairs <- unique(pw[, c("group1", "group2")])
    if (is.null(group_pair_values)) {
      if (nrow(pairs) > 1L) {
        .ph_abort(
          "Several pairs; pick one with group_pair_values = c(g1, g2).",
          bullets = paste("-", pairs$group1, "vs", pairs$group2)
        )
      }
      group_pair_values <- c(pairs$group1, pairs$group2)
    }
    chk::chk_character(group_pair_values)
    chk::chk_length(group_pair_values, 2L)
    pw <- pw[pw$group1 == group_pair_values[1] &
      pw$group2 == group_pair_values[2], , drop = FALSE]
    if (!nrow(pw)) {
      .ph_abort(sprintf(
        "No pairwise rows for %s vs %s (group1 vs group2).",
        group_pair_values[1], group_pair_values[2]
      ))
    }
    out <- tibble::tibble(
      rank = pw$rank, feature = pw$feature,
      group1 = pw$group1, group2 = pw$group2,
      design = "pairwise",
      T_obs = pw$T_obs, p_perm = pw$p_perm, T_obs_stand = NA_real_,
      Z_from_p = sign(pw$T_obs) * stats::qnorm(pw$p_perm / 2,
        lower.tail = FALSE
      ),
      p_maxT = pw$p_maxT
    )
    attr(out, "labels") <- paste("More in", group_pair_values)
    return(out)
  }

  tt <- x$terms
  col <- function(nm) paste0(nm, "_", term)
  if (!col("T") %in% names(tt)) {
    .ph_abort(sprintf("The %s term was not tested.", term))
  }
  tested <- unique(tt$group_col[!is.na(tt[[col("T")]])])
  gc <- .ph_delta2_pick_group_col(x, group_col, tested)
  tt <- tt[tt$group_col == gc & !is.na(tt[[col("T")]]), , drop = FALSE]
  if (any(tt[[col("df")]] != 1L)) {
    .ph_abort(sprintf(
      "The %s term of %s is an omnibus test (df > 1) with no sign to plot;
      plot term = \"pairwise\", or pass scores for a 1-df trend.",
      term, gc
    ))
  }

  des <- x$meta$group_cols[[gc]]$design
  ends <- function(levels, scores) {
    levels <- as.character(levels)
    if (length(scores)) {
      levels[c(which.min(scores), which.max(scores))]
    } else {
      levels[1:2]
    }
  }
  g <- ends(des$group_levels, des$group_scores)
  t <- ends(des$time_levels, des$time_scores)
  lv <- switch(term, group = g, time = t, interaction = g)
  labels <- switch(term,
    group = paste("More in", g),
    time = paste("Higher at", t),
    interaction = paste("Rises more in", g)
  )
  out <- tibble::tibble(
    rank = tt$rank, feature = tt$feature,
    group1 = lv[1], group2 = lv[2],
    design = paste(term, "term"),
    T_obs = tt[[col("T")]], p_perm = tt[[col("p")]],
    T_obs_stand = tt[[col("T_stand")]], Z_from_p = tt[[col("Z_from_p")]],
    n_peptides_used = tt$n_peptides_used
  )
  attr(out, "labels") <- labels
  out
}
