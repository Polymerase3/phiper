# ==============================================================================
# Alpha diversity – single exported function (no S3)
# ==============================================================================

#' @title Compute alpha diversity per sample / group across ranks
#'
#' @description
#' Computes **richness**, **Shannon**, and **Simpson** diversity per sample and
#' per grouping variable at one or more **ranks** (columns describing peptides).
#'
#' @details
#' ## Ranks
#' Ranks are the peptide identities or characteristics you aggregate by. They
#' must be **exact column names**:
#' - For `<phip_data>`: columns in the PHIPER peptide library (e.g. `peptide_id`,
#'   lineage/taxa fields).
#' - For `data.frame`: columns present on your long count table.
#'
#' ## Presence rule
#' - Default: `exist > 0`.
#' - If `fc_threshold` is numeric, presence is `fold_change > fc_threshold`.
#'
#' ## Grouping, interactions, and interaction-only mode
#' - `group_cols` can be a character vector; the return value is a **named list**
#'   of data frames, one per `group_col`.
#' - If `group_cols = NULL`, a single non-facetted table is returned under the
#'   name `"all_samples"`.
#' - If `group_interaction = TRUE` (and you supplied ≥ 2 `group_cols`), an
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
#' @param fc_threshold Numeric or `NULL`. If `NULL` (default), presence is
#'   `exist > 0`. If numeric, presence is `fold_change > fc_threshold`.
#' @param shannon_log One of `"ln"`, `"log2"`, `"log10"`; reporting base for the
#'   Shannon index (via base change from natural log).
#' @param carry_cols Optional character vector of extra columns to carry forward
#'   into the output if present (e.g. sample metadata in your table).
#' @param group_interaction Logical; also compute the interaction of all
#'   `group_cols` (default `FALSE`).
#' @param interaction_only Logical; if `TRUE`, return only the interaction table
#'   (requires `group_interaction = TRUE` and at least two `group_cols`).
#' @param interaction_sep Separator used for the interaction label (default `" * "`).
#'
#' @return A **named list** of data frames with S3 class
#'   `"phip_alpha_diversity"`. Each element (per `group_col`, plus optional
#'   interaction or `"all_samples"`) contains: `rank`, `sample_id`, the grouping
#'   column (or `group` when `group_cols = NULL`), any `carry_cols`, and the
#'   metrics: `richness`, `shannon_diversity`, `simpson_diversity`.
#'
#' @examples
#' pd <- phip_load_example_data()
#' # phip_data input — peptide-level diversity by group
#' out <- compute_alpha_diversity(
#'   pd, group_cols = "group", ranks = "peptide_id"
#' )
#'
#' # include interaction of multiple grouping variables
#' out2 <- compute_alpha_diversity(
#'   pd,
#'   group_cols = c("group", "timepoint"),
#'   ranks = c("peptide_id", "family", "genus"),
#'   group_interaction = TRUE
#' )
#'
#' # interaction only (returns a single element named "group * timepoint")
#' out3 <- compute_alpha_diversity(
#'   pd,
#'   group_cols = c("group", "timepoint"),
#'   ranks = "peptide_id",
#'   group_interaction = TRUE,
#'   interaction_only = TRUE
#' )
#'
#' \dontrun{
#' # data.frame input — ranks must be columns in the data
#' out_df <- compute_alpha_diversity(
#'   df_long, group_cols = NULL, ranks = "peptide_id"
#' )
#' }
#'
#' # presence via fold-change
#' out_fc <- compute_alpha_diversity(
#'   pd, group_cols = "group", ranks = "peptide_id", fc_threshold = 1.5
#' )
#'
#' @export
compute_alpha_diversity <- function(x,
                                    group_cols = NULL,
                                    ranks = "peptide_id",
                                    fc_threshold = NULL,
                                    shannon_log = c("ln", "log2", "log10"),
                                    carry_cols = NULL,
                                    group_interaction = FALSE,
                                    interaction_only = FALSE,
                                    interaction_sep = " * ") {
  .data <- rlang::.data
  shannon_log <- match.arg(shannon_log)

  # -- choose data table and mapping provider depending on input class ----------
  if (inherits(x, "phip_data")) {
    tbl <- x$data_long

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
      if (length(miss_gc)) {
        .ph_abort(
          headline = "Grouping columns not found in data_long.",
          step = "input validation",
          bullets = sprintf("missing: %s", paste(add_quotes(miss_gc, 1L), collapse = ", "))
        )
      }
    }

    # attach peplib on main and prepare a mapper from peptide_id -> rank_val
    peplib_main <- .ensure_peplib_on_main(x)
    peplib_cols <- colnames(peplib_main)
    .ph_log_info(
      "Peptide library attached on main connection",
      bullets = c(
        sprintf(
          "available columns: %s%s",
          paste(utils::head(peplib_cols, 8), collapse = ", "),
          if (length(peplib_cols) > 8) sprintf(" …(+%d)", length(peplib_cols) - 8) else ""
        )
      )
    )

    map_provider <- function(rank_name) {
      # return two columns: peptide_id, rank_val (or NULL to skip with warning)
      if (!(rank_name %in% peplib_cols)) {
        .ph_warn(
          headline = "Rank not found in peptide_library (skipping).",
          step     = "rank mapping",
          bullets  = sprintf("rank: %s", add_quotes(rank_name, 1L))
        )
        return(NULL)
      }
      peplib_main |>
        dplyr::select(peptide_id, rank_val = .data[[rank_name]]) |>
        dplyr::distinct()
    }

  } else if (is.data.frame(x)) {
    tbl <- x

    # validate grouping columns if provided
    if (!is.null(group_cols) && length(group_cols)) {
      miss_gc <- setdiff(group_cols, colnames(tbl))
      if (length(miss_gc)) {
        .ph_abort(
          headline = "Grouping columns not found in data.frame.",
          step = "input validation",
          bullets = sprintf("missing: %s", paste(add_quotes(miss_gc, 1L), collapse = ", "))
        )
      }
    }

    # mapping for data.frame uses columns already in tbl
    map_provider <- function(rank_name) {
      if (!(rank_name %in% colnames(tbl))) {
        .ph_warn(
          headline = "Rank not found in data.frame (skipping).",
          step     = "rank mapping",
          bullets  = sprintf("rank: %s", add_quotes(rank_name, 1L))
        )
        return(NULL)
      }
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

  # -- compute ------------------------------------------------------------------
  .ph_with_timing(
    headline = sprintf("Computing alpha diversity (%s)",
                       if (inherits(x, "phip_data")) "<phip_data>" else "data.frame"),
    step = paste0(
      "group_cols: ", if (is.null(group_cols)) "<none>" else paste(add_quotes(group_cols, 1L), collapse = ", "),
      "; ranks: ", paste(add_quotes(ranks, 1L), collapse = ", ")
    ),
    expr = {
      out_list <- list()

      if (is.null(group_cols) || !length(group_cols)) {
        # no grouping -> single table "all_samples"
        out_list[["all_samples"]] <-
          .compute_alpha_for_group(
            tbl,
            group_col = NULL, ranks = ranks,
            fc_threshold = fc_threshold, shannon_log = shannon_log,
            carry_cols = carry_cols, map_provider = map_provider
          )

      } else {
        # grouping case
        if (!isTRUE(interaction_only)) {
          # per-group tables unless interaction-only is requested
          for (gc in group_cols) {
            out_list[[gc]] <-
              .compute_alpha_for_group(
                tbl,
                group_col = gc, ranks = ranks,
                fc_threshold = fc_threshold, shannon_log = shannon_log,
                carry_cols = carry_cols, map_provider = map_provider
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

          # if interaction_only = TRUE -> we only populate this element
          out_list[[combo_nm]] <-
            .compute_alpha_for_group(
              tbl_inter,
              group_col = inter_col, ranks = ranks,
              fc_threshold = fc_threshold, shannon_log = shannon_log,
              carry_cols = carry_cols, map_provider = map_provider
            )
        }
      }

      # normalize to tibble + add class/attributes
      out_list <- lapply(out_list, tibble::as_tibble)
      class(out_list) <- c("phip_alpha_diversity", class(out_list))
      attr(out_list, "group_cols")     <- group_cols
      attr(out_list, "ranks")          <- unique(ranks)
      attr(out_list, "fc_threshold")   <- fc_threshold
      attr(out_list, "shannon_log")    <- shannon_log
      attr(out_list, "interaction")    <- isTRUE(group_interaction)
      attr(out_list, "interaction_only") <- isTRUE(interaction_only)
      attr(out_list, "interaction_sep")  <- interaction_sep

      out_list
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ------------------------------------------------------------------------------
# helpers (internal)
# ------------------------------------------------------------------------------

# natural->base-b change factor for shannon (H_b = H_ln / ln(b))
.phip_ln_base <- function(shannon_log = c("ln", "log2", "log10")) {
  shannon_log <- match.arg(shannon_log)
  switch(shannon_log,
    ln    = 1.0,
    log2  = log(2),
    log10 = log(10)
  )
}

# engine used by the exported function
# - tbl: counts table (lazy or local) containing sample_id, peptide_id, exist,
#        optional fold_change, and an optional single "group_col"
# - map_provider(rank_name): returns a two-column table (peptide_id, rank_val),
#        lazy on same con (phip_data) or local tibble (data.frame); NULL -> skip
.compute_alpha_for_group <- function(tbl,
                                     group_col = NULL,
                                     ranks,
                                     fc_threshold = NULL,
                                     shannon_log = c("ln", "log2", "log10"),
                                     carry_cols = NULL,
                                     map_provider) {
  .data <- rlang::.data
  shannon_log <- match.arg(shannon_log)
  ln_base <- .phip_ln_base(shannon_log)

  # validate required columns
  need <- c("sample_id", "peptide_id", "exist")
  if (!is.null(fc_threshold)) need <- union(need, "fold_change")
  if (!is.null(group_col))    need <- union(need, group_col)
  miss <- setdiff(need, colnames(tbl))
  if (length(miss)) {
    .ph_abort(
      headline = "Missing required columns.",
      step = "input validation",
      bullets = sprintf("missing: %s", paste(add_quotes(miss, 1L), collapse = ", "))
    )
  }

  # unified cohort column (character)
  tbl <- if (is.null(group_col)) {
    dplyr::mutate(tbl, cohort = "All samples")
  } else {
    dplyr::mutate(tbl, cohort = .data[[group_col]])
  }

  # presence rule
  pres_tbl <- if (is.null(fc_threshold)) {
    dplyr::filter(tbl, .data$exist > 0)
  } else {
    dplyr::filter(tbl, .data$fold_change > !!fc_threshold)
  }

  # distinct present tuples
  pres_min <- pres_tbl |>
    dplyr::transmute(sample_id = .data$sample_id,
                     peptide_id = .data$peptide_id,
                     cohort = .data$cohort) |>
    dplyr::distinct()

  # validate ranks
  ranks <- unique(ranks)
  if (!length(ranks) || !is.character(ranks)) {
    .ph_abort("`ranks` must be a non-empty character vector of exact column names.")
  }

  # compute per-rank
  compute_one_rank <- function(rank_name) {
    # map peptide_id -> rank_val; for "peptide_id" it's identity
    if (identical(rank_name, "peptide_id")) {
      ranked <- pres_min |>
        dplyr::transmute(sample_id, cohort, rank_val = peptide_id)
    } else {
      map_tbl <- map_provider(rank_name)
      if (is.null(map_tbl)) return(NULL)
      ranked <- pres_min |>
        dplyr::inner_join(map_tbl, by = "peptide_id") |>
        dplyr::filter(!is.na(.data$rank_val))
    }

    # counts per (sample, cohort, rank_val)
    per_cat <- ranked |>
      dplyr::group_by(sample_id, cohort, rank_val) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop")

    pc <- dplyr::collect(per_cat)

    # diversity by sample
    if (nrow(pc)) {
      by_sample <- pc |>
        dplyr::group_by(sample_id, cohort) |>
        dplyr::summarise(
          richness = dplyr::n_distinct(rank_val),
          H_ln = { p <- n / sum(n); -sum(p * log(p)) },
          simpson = { p <- n / sum(n); 1 - sum(p * p) },
          .groups = "drop"
        ) |>
        dplyr::mutate(shannon = H_ln / ln_base) |>
        dplyr::select(-H_ln)
    } else {
      by_sample <- tibble::tibble(
        sample_id = character(0), cohort = character(0),
        richness = integer(0), shannon = numeric(0), simpson = numeric(0)
      )
    }

    # carry through requested columns (if present) for completeness
    keep_cols <- c("sample_id", "cohort", carry_cols)
    keep_cols <- intersect(keep_cols, colnames(tbl))
    all_samples <- tbl |>
      dplyr::distinct(dplyr::across(dplyr::all_of(keep_cols))) |>
      dplyr::collect()

    out <- all_samples |>
      dplyr::left_join(by_sample, by = c("sample_id", "cohort")) |>
      dplyr::mutate(
        richness = tidyr::replace_na(richness, 0L),
        shannon  = tidyr::replace_na(shannon, 0),
        simpson  = tidyr::replace_na(simpson, 0)
      ) |>
      dplyr::mutate(rank = rank_name, .before = 1)

    # light normalization of names
    norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
    names(out) <- norm(names(out))

    out |>
      dplyr::rename(
        shannon_diversity = shannon,
        simpson_diversity = simpson
      )
  }

  res <- do.call(dplyr::bind_rows, lapply(ranks, compute_one_rank))

  # rename cohort back to either "group" (no grouping) or the original column
  if (is.null(group_col)) {
    res <- dplyr::rename(res, group = cohort)
  } else {
    res <- dplyr::rename(res, !!rlang::sym(group_col) := cohort)
  }

  res
}
