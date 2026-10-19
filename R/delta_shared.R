# ==============================================================================
# Helpers shared by compute_delta() and compute_delta2()
# ==============================================================================
# Steps that do not depend on the statistic: the rank/feature map of the
# peptides, the rank_feature_keep filter, and the (optionally parallel)
# runner over the master plan.

#' @title Internal helper: .ph_delta_rank_map
#' @description Long map of the retained peptides to `(rank, feature)`,
#'   resolving the peptide library for ranks other than `"peptide_id"` from
#'   `peptide_library`, the `phip_data` object, or a known library.
#' @return Tibble with columns `peptide_id`, `rank`, `feature`.
#' @keywords internal
.ph_delta_rank_map <- function(x, rank_cols, peptides_order,
                               peptide_library = NULL) {
  ranks_need <- setdiff(rank_cols, "peptide_id")

  if (length(ranks_need) == 0L) {
    # Only peptide_id requested --> trivial long map
    return(
      tibble::tibble(
        peptide_id = peptides_order, rank = "peptide_id",
        feature = peptides_order
      )
    )
  }

  lib_src <- NULL
  if (!is.null(peptide_library)) {
    lib_src <- peptide_library
  } else if (inherits(x, "phip_data") && !is.null(x$peptide_library)) {
    lib_src <- x$peptide_library
  } else {
    lib_src <- .ph_library_for_peptides(peptides_order)
    if (is.null(lib_src)) {
      .ph_abort(
        "No known peptide library matches the peptide IDs.",
        bullets = c(
          paste("- rank needing a library:", ranks_need),
          "- supply `peptide_library` with these rank columns"
        )
      )
    }
  }

  # Select needed columns and collect to R (handles DuckDB/lazy)
  lib_needed <- c("peptide_id", ranks_need)
  miss_tax <- setdiff(ranks_need, colnames(lib_src))
  if (length(miss_tax)) {
    .ph_abort("Peptide library missing required columns:",
      bullets = paste("-", miss_tax)
    )
  }

  lib_small <- lib_src |>
    dplyr::select(dplyr::any_of(lib_needed)) |>
    dplyr::distinct()

  if (inherits(lib_small, "tbl_sql") || inherits(lib_small, "tbl_lazy")) {
    lib_small <- lib_small |> dplyr::collect()
  }

  # Keep only peptides we actually retained (peptides_order), then pivot to
  # long
  lib_small <- lib_small |>
    dplyr::filter(.data$peptide_id %in% peptides_order)

  rank_map_long <- lib_small |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(ranks_need),
      names_to  = "rank",
      values_to = "feature",
      values_transform = list(feature = as.character)
    ) |>
    dplyr::filter(!is.na(.data$feature)) |>
    dplyr::select(peptide_id, rank, feature) |>
    dplyr::distinct()

  # Also add the direct peptide_id rank if requested
  if ("peptide_id" %in% rank_cols) {
    rank_map_long <- dplyr::bind_rows(
      rank_map_long,
      tibble::tibble(
        peptide_id = peptides_order, rank = "peptide_id",
        feature = peptides_order
      )
    )
  }

  rank_map_long
}

#' @title Internal helper: .ph_delta_keep_features
#' @description Restrict the rank map to the features listed per rank in
#'   `rank_feature_keep` (ranks not listed, or listed as `NULL`, keep all).
#' @keywords internal
.ph_delta_keep_features <- function(rank_map_long, rank_feature_keep = NULL) {
  if (!is.null(rank_feature_keep) && length(rank_feature_keep)) {
    for (rk in names(rank_feature_keep)) {
      vals <- rank_feature_keep[[rk]]
      if (is.null(vals)) {
        next
      } else {
        keep_vals <- as.character(vals)
        rank_map_long <- rank_map_long |>
          dplyr::filter(!(rank == rk) | feature %in% keep_vals)
      }
    }
  }
  rank_map_long
}

#' @title Internal helper: .ph_delta_run_tasks
#' @description Run `run_one(i)` over the master-plan rows, sequentially or
#'   with `future.apply` when a parallel `future` plan is set (heaviest rows
#'   first, dynamic scheduling), with inner threading pinned to one thread.
#' @param work_weight Numeric work estimate per row.
#' @return `list(results, n_workers)`; `results` may contain `NULL`s.
#' @keywords internal
.ph_delta_run_tasks <- function(work_weight, run_one) {
  n_tasks <- length(work_weight)

  # decide workers based on current future backend
  n_workers <- 1L
  if (rlang::is_installed("future")) {
    n_workers <- tryCatch(future::nbrOfWorkers(), error = function(...) 1L)
  }
  n_workers <- max(1L, min(as.integer(n_workers), n_tasks))

  # ---- dynamic scheduling for future.apply ----
  op_old <- options(
    future.scheduling   = Inf,
    future.chunk.size   = 1
  )
  on.exit(options(op_old), add = TRUE)

  # ---- avoid inner oversubscription (TBB inside RcppParallel) ----
  if (rlang::is_installed("RcppParallel")) {
    RcppParallel::setThreadOptions(numThreads = 1) # fairness > peak single-task
  } #            speed
  env_threads <- c(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1"
  )
  results <- withr::with_envvar(env_threads, {
    if (n_workers == 1L) {
      # --- sequential
      out <- list()
      for (i in seq_len(n_tasks)) {
        row <- run_one(i)
        if (!is.null(row)) out[[length(out) + 1L]] <- row
      }
      out
    } else {
      # --- parallel via future.apply (respect existing future plan)
      idx <- order(work_weight, decreasing = TRUE)
      future.apply::future_lapply(
        X = idx, FUN = run_one,
        future.seed = TRUE,
        future.scheduling = Inf,
        future.chunk.size = 1,
        future.packages = c("dplyr", "tibble", "tidyr")
      )
    }
  })

  list(results = results, n_workers = n_workers)
}
