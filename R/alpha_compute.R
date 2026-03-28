#' @title Compute alpha diversity per sample / group across ranks
#'
#' @description
#' Computes **richness**, **Shannon**, **Simpson**, **Pielou's evenness**, and
#' **Berger-Parker dominance** per sample and per grouping variable at one or
#' more **ranks** (columns describing peptides).
#'
#' @details
#' ## Ranks
#' Ranks are the peptide identities or characteristics you aggregate by. They
#' must be **exact column names**:
#' - For `<phip_data>`: columns in the PHIPER peptide library (e.g. `peptide_id`,
#'   lineage/taxa fields).
#' - For `data.frame`: columns present on your long count table.
#'
#' ## Modes
#' The `mode` parameter controls what `n` represents in every diversity formula:
#' - `"binary"` (default): presence = `exist > 0`. `n` = count of distinct
#'   enriched peptides per `(sample, rank_val)`.
#' - `"threshold"`: presence = `abundance_col > threshold`. `n` = count of
#'   peptides passing the threshold. `threshold` is required; `abundance_col`
#'   defaults to `"fold_change"` when `NULL`.
#' - `"abundance"`: no binarisation. `abundance_col` is required. At
#'   `rank = "peptide_id"`, `n` = the raw abundance value per peptide. At
#'   higher ranks, peptides are aggregated within each `(sample, rank_val)`
#'   using `abundance_agg` (`"mean"`, `"sum"`, or `"max"`).
#'
#' ## Grouping, interactions, and interaction-only mode
#' - `group_cols` can be a character vector; the return value is a **named list**
#'   of data frames, one per `group_col`.
#' - If `group_cols = NULL`, a single non-facetted table is returned under the
#'   name `"all_samples"`.
#' - If `group_interaction = TRUE` (and you supplied >= 2 `group_cols`), an
#'   additional element is computed for the interaction of all group columns,
#'   with labels joined by `interaction_sep`.
#' - If `interaction_only = TRUE`, you get **only** that interaction element
#'   (requires `group_interaction = TRUE` and at least 2 `group_cols`).
#'
#' @param x A `<phip_data>` object or a long `data.frame`.
#' @param group_cols Character vector of grouping columns, or `NULL` for a single
#'   aggregate (non-facetted) table. All columns must be present on the input.
#' @param ranks Character vector of **exact column names** to aggregate by.
#'   Typical values: `"peptide_id"` or taxonomy/lineage columns.
#' @param mode One of `"binary"` (default), `"threshold"`, or `"abundance"`.
#'   Controls how per-category counts are computed; see **Modes** in Details.
#' @param abundance_col Character scalar naming the abundance column. Required
#'   for `mode = "abundance"`; optional for `mode = "threshold"` (defaults to
#'   `"fold_change"`); ignored for `mode = "binary"`.
#' @param threshold Numeric scalar. Required for `mode = "threshold"`. Ignored
#'   for other modes.
#' @param abundance_agg One of `"mean"` (default), `"sum"`, or `"max"`. Used
#'   only in `mode = "abundance"` when `rank != "peptide_id"` to aggregate
#'   peptide-level abundances within each rank category.
#' @param metrics Character vector of diversity metrics to compute. Any subset
#'   of: `"richness"`, `"shannon"`, `"simpson"`, `"pielou_evenness"`,
#'   `"berger_parker"`. Defaults to all five. Output columns follow the
#'   canonical names: `richness`, `shannon_diversity`, `simpson_diversity`,
#'   `pielou_evenness`, `berger_parker_dominance`.
#' @param shannon_base One of `"ln"`, `"log2"`, `"log10"`; reporting base for
#'   the Shannon index (via base change from natural log).
#' @param carry_cols Optional character vector of extra columns to carry forward
#'   into the output if present (e.g. sample metadata in your table).
#' @param group_interaction Logical; also compute the interaction of all
#'   `group_cols` (default `FALSE`).
#' @param interaction_only Logical; if `TRUE`, return only the interaction table
#'   (requires `group_interaction = TRUE` and at least two `group_cols`).
#' @param interaction_sep Separator used for the interaction label (default
#'   `" * "`). The resulting string (e.g. `"GroupA * T1"`) becomes both a list
#'   name and a value in the output column, so choose a separator that does not
#'   appear in your group levels.
#' @param shannon_log Deprecated. Use `shannon_base` instead.
#'
#' @return A **named list** of data frames with S3 class
#'   `"phip_alpha_diversity"`. Each element (per `group_col`, plus optional
#'   interaction or `"all_samples"`) contains: `rank`, `sample_id`, the grouping
#'   column (or `group` when `group_cols = NULL`), any `carry_cols`, and one
#'   column per requested metric (see `metrics`). `pielou_evenness` is `NA`
#'   when richness <= 1; `berger_parker_dominance` is `NA` when richness == 0.
#'
#' @examples
#' pd <- load_example_data()
#' # phip_data input: peptide-level diversity by group
#' out <- compute_alpha(
#'   pd, group_cols = "group", ranks = "peptide_id"
#' )
#'
#' # include interaction of multiple grouping variables
#' out2 <- compute_alpha(
#'   pd,
#'   group_cols = c("group", "timepoint"),
#'   ranks = c("peptide_id", "family", "genus"),
#'   group_interaction = TRUE
#' )
#'
#' # interaction only (returns a single element named "group * timepoint")
#' out3 <- compute_alpha(
#'   pd,
#'   group_cols = c("group", "timepoint"),
#'   ranks = "peptide_id",
#'   group_interaction = TRUE,
#'   interaction_only = TRUE
#' )
#'
#' \dontrun{
#' # data.frame input: ranks must be columns in the data
#' out_df <- compute_alpha(
#'   df_long, group_cols = NULL, ranks = "peptide_id"
#' )
#' }
#'
#' # threshold mode: presence = fold_change > 1.5
#' out_thr <- compute_alpha(
#'   pd, group_cols = "group", ranks = "peptide_id",
#'   mode = "threshold", threshold = 1.5
#' )
#'
#' @export
compute_alpha <- function(x,
                                    group_cols = NULL,
                                    ranks = "peptide_id",
                                    mode = c("binary", "threshold", "abundance"),
                                    abundance_col = NULL,
                                    threshold = NULL,
                                    abundance_agg = c("mean", "sum", "max"),
                                    metrics = .alpha_metric_names,
                                    shannon_base = c("ln", "log2", "log10"),
                                    carry_cols = NULL,
                                    group_interaction = FALSE,
                                    interaction_only = FALSE,
                                    interaction_sep = " * ",
                                    shannon_log = NULL) {
  .data <- rlang::.data

  # -- deprecation shim for shannon_log -----------------------------------------
  if (!is.null(shannon_log)) {
    .ph_warn(
      headline = "`shannon_log` is deprecated; use `shannon_base` instead.",
      step = "argument check"
    )
    shannon_base <- shannon_log
  }
  shannon_base  <- match.arg(shannon_base)
  mode          <- match.arg(mode)
  abundance_agg <- match.arg(abundance_agg)
  metrics       <- match.arg(metrics, .alpha_metric_names, several.ok = TRUE)

  # -- mode-specific validation --------------------------------------------------
  .ph_check_cond(mode == "binary" && !is.null(abundance_col),
    "`abundance_col` is ignored in mode = 'binary'.", error = FALSE, step = "argument check")
  .ph_check_cond(mode == "binary" && !is.null(threshold),
    "`threshold` is ignored in mode = 'binary'.", error = FALSE, step = "argument check")

  .ph_check_cond(
    mode == "threshold" && (is.null(threshold) || !is.numeric(threshold) || length(threshold) != 1L),
    "`threshold` must be a numeric scalar when mode = 'threshold'.", step = "argument validation"
  )
  .ph_check_cond(
    mode == "threshold" && is.numeric(threshold) && length(threshold) == 1L && !is.finite(threshold),
    "`threshold` must be finite (not NA or Inf) when mode = 'threshold'.", step = "argument validation"
  )
  .ph_check_cond(
    mode == "threshold" && !is.null(abundance_col) &&
      (!is.character(abundance_col) || length(abundance_col) != 1L),
    "`abundance_col` must be a character scalar when mode = 'threshold'.", step = "argument validation"
  )

  .ph_check_cond(
    mode == "abundance" && (is.null(abundance_col) || !is.character(abundance_col) || length(abundance_col) != 1L),
    "`abundance_col` must be a character scalar when mode = 'abundance'.", step = "argument validation"
  )
  .ph_check_cond(mode == "abundance" && !is.null(threshold),
    "`threshold` is ignored in mode = 'abundance'.", error = FALSE, step = "argument check")

  # -- choose data table and mapping provider depending on input class ----------
  if (inherits(x, "phip_data")) {
    tbl <- x$data_long
    # roster kept on unpruned source so samples with all-zero exist are retained
    tbl_roster <- tbl

    # log + prune full-cross rows (exist == 0) early to reduce volume
    if (isTRUE(x$meta$full_cross) && ("exist" %in% colnames(tbl))) {
      red_txt <- tryCatch({
        ep <- as.numeric(x$meta$exist_prop)
        if (is.finite(ep) && ep > 0) sprintf("~%.1fx", 1 / ep) else "<unknown>"
      }, error = function(e) "<unknown>")
      .ph_log_info(
        "Full-cross detected; pruning non-existing rows before alpha calc",
        bullets = c("rule: keep exist == 1", sprintf("estimated reduction: %s", red_txt))
      )
      tbl <- dplyr::filter(tbl, .data$exist == 1L)
    }

    # validate grouping columns if provided
    if (!is.null(group_cols) && length(group_cols)) {
      miss_gc <- setdiff(group_cols, colnames(tbl))
      .ph_check_cond(
        length(miss_gc) > 0L,
        "Grouping columns not found in data_long.",
        step = "input validation",
        bullets = sprintf("missing: %s", paste(.ph_add_quotes(miss_gc, 1L), collapse = ", "))
      )
    }

    .ph_check_cond(
      !is.null(abundance_col) && !(abundance_col %in% colnames(tbl)),
      sprintf("`abundance_col` '%s' not found in data_long.", abundance_col),
      step = "input validation"
    )

    # attach peplib on main connection and prepare a mapper peptide_id -> rank_val
    peplib_main <- .ph_peplib_on_main(x)
    peplib_cols <- colnames(peplib_main)

    map_provider <- function(rank_name) {
      # return two columns: peptide_id, rank_val (or NULL to skip with warning)
      if (!(rank_name %in% peplib_cols)) {
        .ph_warn(
          headline = "Rank not found in peptide_library (skipping).",
          step     = "rank mapping",
          bullets  = sprintf("rank: %s", .ph_add_quotes(rank_name, 1L))
        )
        return(NULL)
      }
      # select before distinct: deduplicates on (peptide_id, rank_val) only,
      # not on extra peplib columns that would otherwise survive distinct()
      peplib_main |>
        dplyr::select(peptide_id, rank_val = .data[[rank_name]]) |>
        dplyr::distinct()
    }

  } else if (is.data.frame(x)) {
    tbl <- x
    tbl_roster <- tbl

    # validate grouping columns if provided
    if (!is.null(group_cols) && length(group_cols)) {
      miss_gc <- setdiff(group_cols, colnames(tbl))
      .ph_check_cond(
        length(miss_gc) > 0L,
        "Grouping columns not found in data.frame.",
        step = "input validation",
        bullets = sprintf("missing: %s", paste(.ph_add_quotes(miss_gc, 1L), collapse = ", "))
      )
    }

    .ph_check_cond(
      !is.null(abundance_col) && !(abundance_col %in% colnames(tbl)),
      sprintf("`abundance_col` '%s' not found in data.frame.", abundance_col),
      step = "input validation"
    )

    # mapping for data.frame uses columns already in tbl
    map_provider <- function(rank_name) {
      if (!(rank_name %in% colnames(tbl))) {
        .ph_warn(
          headline = "Rank not found in data.frame (skipping).",
          step     = "rank mapping",
          bullets  = sprintf("rank: %s", .ph_add_quotes(rank_name, 1L))
        )
        return(NULL)
      }
      # select two cols before distinct: guards against extra columns
      # in tbl that would otherwise prevent deduplication on (peptide_id, rank_val)
      tibble::tibble(
        peptide_id = tbl$peptide_id,
        rank_val   = tbl[[rank_name]]
      ) |>
        dplyr::distinct()
    }

  } else {
    .ph_abort("`x` must be either a <phip_data> or a data.frame.")
  }

  # -- enforce interaction_only preconditions -----------------------------------
  .ph_check_cond(
    isTRUE(interaction_only) &&
      (!isTRUE(group_interaction) || is.null(group_cols) || length(group_cols) < 2L),
    "`interaction_only = TRUE` requires an interaction.",
    step = "argument validation",
    bullets = c(
      "set group_interaction = TRUE",
      "provide at least two grouping columns in group_cols"
    )
  )

  # -- compute ------------------------------------------------------------------
  .ph_with_timing(
    headline = sprintf("Computing alpha diversity (%s)",
                       if (inherits(x, "phip_data")) "<phip_data>" else "data.frame"),
    step = paste0(
      "group_cols: ", if (is.null(group_cols)) "<none>" else paste(.ph_add_quotes(group_cols, 1L), collapse = ", "),
      "; ranks: ", paste(.ph_add_quotes(ranks, 1L), collapse = ", ")
    ),
    expr = {
      out_list <- list()

      if (is.null(group_cols) || !length(group_cols)) {
        # no grouping -> single table "all_samples"
        out_list[["all_samples"]] <-
          .compute_alpha_for_group(
            tbl,
            group_col = NULL, ranks = ranks,
            mode = mode, abundance_col = abundance_col,
            threshold = threshold, abundance_agg = abundance_agg,
            metrics = metrics, shannon_base = shannon_base,
            carry_cols = carry_cols, map_provider = map_provider,
            tbl_roster = tbl_roster
          )

      } else {
        # grouping case
        if (!isTRUE(interaction_only)) {
          for (gc in group_cols) {
            out_list[[gc]] <-
              .compute_alpha_for_group(
                tbl,
                group_col = gc, ranks = ranks,
                mode = mode, abundance_col = abundance_col,
                threshold = threshold, abundance_agg = abundance_agg,
                metrics = metrics, shannon_base = shannon_base,
                carry_cols = carry_cols, map_provider = map_provider,
                tbl_roster = tbl_roster
              )
          }
        }

        # optional interaction (either in addition, or as the only output)
        if (isTRUE(group_interaction) && length(group_cols) >= 2L) {
          inter_col <- "phip_interaction"
          combo_nm  <- paste(group_cols, collapse = interaction_sep)
          tbl_inter <- dplyr::mutate(
            tbl,
            !!rlang::sym(inter_col) := paste(!!!rlang::syms(group_cols), sep = interaction_sep)
          )
          tbl_roster_inter <- dplyr::mutate(
            tbl_roster,
            !!rlang::sym(inter_col) := paste(!!!rlang::syms(group_cols), sep = interaction_sep)
          )

          out_list[[combo_nm]] <-
            .compute_alpha_for_group(
              tbl_inter,
              group_col = inter_col, ranks = ranks,
              mode = mode, abundance_col = abundance_col,
              threshold = threshold, abundance_agg = abundance_agg,
              metrics = metrics, shannon_base = shannon_base,
              carry_cols = carry_cols, map_provider = map_provider,
              tbl_roster = tbl_roster_inter
            )
        }
      }

      # normalize to tibble + add class/attributes
      out_list <- lapply(out_list, tibble::as_tibble)
      class(out_list) <- c("phip_alpha_diversity", class(out_list))
      attr(out_list, "group_cols")       <- group_cols
      attr(out_list, "ranks")            <- unique(ranks)
      attr(out_list, "mode")             <- mode
      attr(out_list, "abundance_col")    <- abundance_col
      attr(out_list, "threshold")        <- threshold
      attr(out_list, "abundance_agg")    <- abundance_agg
      attr(out_list, "metrics")          <- metrics
      attr(out_list, "shannon_base")     <- shannon_base
      attr(out_list, "interaction")      <- isTRUE(group_interaction)
      attr(out_list, "interaction_only") <- isTRUE(interaction_only)
      attr(out_list, "interaction_sep")  <- interaction_sep
      attr(out_list, "n_samples")        <- tbl_roster |>
        dplyr::distinct(sample_id) |>
        dplyr::collect() |>
        nrow()

      out_list
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ------------------------------------------------------------------------------
# helpers (internal)
# ------------------------------------------------------------------------------

# natural->base-b change factor for shannon (H_b = H_ln / ln(b))
.phip_ln_base <- function(shannon_base = c("ln", "log2", "log10")) {
  shannon_base <- match.arg(shannon_base)
  switch(shannon_base,
    ln    = 1.0,
    log2  = log(2),
    log10 = log(10)
  )
}

# engine used by the exported function
# - tbl: counts table (lazy or local) containing sample_id, peptide_id, exist,
#        optional fold_change/abundance_col, and an optional single "group_col"
# - map_provider(rank_name): returns a two-column table (peptide_id, rank_val),
#        lazy on same con (phip_data) or local tibble (data.frame); NULL -> skip
.compute_alpha_for_group <- function(tbl,
                                     group_col = NULL,
                                     ranks,
                                     mode = c("binary", "threshold", "abundance"),
                                     abundance_col = NULL,
                                     threshold = NULL,
                                     abundance_agg = c("mean", "sum", "max"),
                                     metrics = .alpha_metric_names,
                                     shannon_base = c("ln", "log2", "log10"),
                                     carry_cols = NULL,
                                     map_provider,
                                     tbl_roster = tbl) {
  .data <- rlang::.data
  mode         <- match.arg(mode)
  abundance_agg <- match.arg(abundance_agg)
  shannon_base <- match.arg(shannon_base)
  ln_base      <- .phip_ln_base(shannon_base)

  # validate required columns
  need <- c("sample_id", "peptide_id", "exist")
  if (mode == "threshold") {
    need <- union(need, abundance_col %||% "fold_change")
  } else if (mode == "abundance") {
    need <- union(need, abundance_col)
  }
  if (!is.null(group_col)) need <- union(need, group_col)
  miss <- setdiff(need, colnames(tbl))
  .ph_check_cond(
    length(miss) > 0L,
    "Missing required columns.",
    step = "input validation",
    bullets = sprintf("missing: %s", paste(.ph_add_quotes(miss, 1L), collapse = ", "))
  )

  # unified cohort column (character) — applied to both tbl and tbl_roster
  add_cohort <- function(t) {
    if (is.null(group_col)) dplyr::mutate(t, cohort = "All samples")
    else                    dplyr::mutate(t, cohort = .data[[group_col]])
  }
  tbl        <- add_cohort(tbl)
  tbl_roster <- add_cohort(tbl_roster)

  # presence rule / filtering
  pres_tbl <- switch(mode,
    binary = dplyr::filter(tbl, .data$exist > 0),
    threshold = {
      col <- abundance_col %||% "fold_change"
      dplyr::filter(tbl, .data[[col]] > !!threshold)
    },
    abundance = dplyr::filter(tbl, .data[[abundance_col]] > 0)
  )

  # distinct present tuples (carry .abundance for abundance mode)
  pres_min <- if (mode == "abundance") {
    pres_tbl |>
      dplyr::transmute(
        sample_id  = .data$sample_id,
        peptide_id = .data$peptide_id,
        cohort     = .data$cohort,
        .abundance = .data[[abundance_col]]
      ) |>
      dplyr::distinct()
  } else {
    pres_tbl |>
      dplyr::transmute(
        sample_id  = .data$sample_id,
        peptide_id = .data$peptide_id,
        cohort     = .data$cohort
      ) |>
      dplyr::distinct()
  }

  # validate ranks
  ranks <- unique(ranks)
  .ph_check_cond(
    !length(ranks) || !is.character(ranks),
    "`ranks` must be a non-empty character vector of exact column names."
  )

  # light normalization of column names — defined once, reused per rank
  .norm_names <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))

  # all-samples roster: built from tbl_roster (unpruned) so samples with all
  # exist==0 are retained; collected once before the rank loop
  keep_cols <- intersect(c("sample_id", "cohort", carry_cols), colnames(tbl_roster))
  all_samples_local <- tbl_roster |>
    dplyr::distinct(dplyr::across(dplyr::all_of(keep_cols))) |>
    dplyr::collect()

  # aggregation function for abundance mode at higher ranks
  agg_fn <- switch(abundance_agg, sum = sum, mean = mean, max = max)

  # compute per-rank
  compute_one_rank <- function(rank_name) {
    # map peptide_id -> rank_val; for "peptide_id" it's identity
    if (identical(rank_name, "peptide_id")) {
      ranked <- if (mode == "abundance") {
        pres_min |>
          dplyr::transmute(sample_id, cohort, rank_val = peptide_id, .abundance)
      } else {
        pres_min |>
          dplyr::transmute(sample_id, cohort, rank_val = peptide_id)
      }
    } else {
      map_tbl <- map_provider(rank_name)
      if (is.null(map_tbl)) return(NULL)
      ranked <- pres_min |>
        dplyr::inner_join(map_tbl, by = "peptide_id") |>
        dplyr::filter(!is.na(.data$rank_val))
    }

    # counts per (sample, cohort, rank_val)
    per_cat <- if (mode != "abundance") {
      ranked |>
        dplyr::group_by(sample_id, cohort, rank_val) |>
        dplyr::summarise(n = dplyr::n(), .groups = "drop")
    } else {
      ranked |>
        dplyr::group_by(sample_id, cohort, rank_val) |>
        dplyr::summarise(n = agg_fn(.data$.abundance, na.rm = TRUE), .groups = "drop")
    }

    pc <- dplyr::collect(per_cat)

    # diversity by sample — driven by the metric registry
    if (nrow(pc)) {
      by_sample <- pc |>
        dplyr::group_by(sample_id, cohort) |>
        dplyr::group_modify(~ {
          n_vec <- .x$n
          vals  <- vapply(metrics, function(m) {
            .alpha_metric_fns[[m]](n_vec, ln_base = ln_base)
          }, numeric(1))
          tibble::as_tibble_row(vals)
        }) |>
        dplyr::ungroup()
    } else {
      empty_cols <- stats::setNames(rep(list(numeric(0)), length(metrics)), metrics)
      by_sample  <- rlang::inject(tibble::tibble(
        sample_id = character(0),
        cohort    = character(0),
        !!!empty_cols
      ))
    }

    out <- all_samples_local |>
      dplyr::left_join(by_sample, by = c("sample_id", "cohort")) |>
      dplyr::mutate(rank = rank_name, .before = 1)

    names(out) <- .norm_names(names(out))

    # fill zeros for metrics where S=0 means 0, keep NA for undefined ones
    for (m in intersect(metrics, .zero_fill_metrics)) {
      out[[m]] <- tidyr::replace_na(out[[m]], if (m == "richness") 0L else 0)
    }

    # apply canonical output column names (e.g. shannon -> shannon_diversity)
    to_rename <- .alpha_metric_output_names[metrics]
    to_rename <- to_rename[to_rename != names(to_rename)]
    if (length(to_rename)) {
      out <- dplyr::rename(out, !!!stats::setNames(names(to_rename), to_rename))
    }

    out
  }

  rank_results <- lapply(ranks, compute_one_rank)
  .ph_check_cond(
    all(vapply(rank_results, is.null, logical(1))),
    "No valid ranks found.",
    step    = "rank validation",
    bullets = sprintf("requested: %s", paste(.ph_add_quotes(ranks, 1L), collapse = ", "))
  )
  res <- do.call(dplyr::bind_rows, rank_results)

  # rename cohort back to either "group" (no grouping) or the original column
  if (is.null(group_col)) {
    res <- dplyr::rename(res, group = cohort)
  } else {
    res <- dplyr::rename(res, !!rlang::sym(group_col) := cohort)
  }

  res
}

# ==============================================================================
# Alpha diversity metric registry
#
# Each entry in .alpha_metric_fns is a function(n, ln_base, ...) -> numeric(1).
#   n        : numeric vector of per-rank-val counts for one sample
#   ln_base  : log base conversion factor (from .phip_ln_base); functions that
#              don't use it accept it silently via ...
# ==============================================================================

.alpha_metric_fns <- list(

  richness = function(n, ...) length(n),

  shannon = function(n, ln_base = 1, ...) {
    p <- n / sum(n)
    -sum(p * log(p)) / ln_base
  },

  simpson = function(n, ...) {
    p <- n / sum(n)
    1 - sum(p^2)
  },

  pielou_evenness = function(n, ...) {
    S <- length(n)
    if (S <= 1L) return(NA_real_)
    p <- n / sum(n)
    H <- -sum(p * log(p))
    H / log(S)   # base cancels in H / log(S); ln_base not applied
  },

  berger_parker = function(n, ...) {
    max(n) / sum(n)
  }
)

# canonical names available for the `metrics` parameter
.alpha_metric_names <- names(.alpha_metric_fns)

# registry name -> output column name
.alpha_metric_output_names <- c(
  richness        = "richness",
  shannon         = "shannon_diversity",
  simpson         = "simpson_diversity",
  pielou_evenness = "pielou_evenness",
  berger_parker   = "berger_parker_dominance"
)

# metrics filled with 0 when a sample has no enriched peptides (S = 0)
.zero_fill_metrics <- c("richness", "shannon", "simpson")

# metrics kept as NA when a sample has no enriched peptides (undefined, not 0)
.na_keep_metrics <- c("pielou_evenness", "berger_parker")
