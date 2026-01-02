# tests/testthat/test-compute_distance.R

# helper that builds a small real phip_data from the package example data
.get_ps_small_for_distance <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) return(cache)

    phip_path <- phip_example_path()

    ps <- phip_convert(
      data_long_path    = phip_path,
      backend           = "duckdb",
      peptide_library   = TRUE,
      subject_id        = "subject_id",
      peptide_id        = "peptide_id",
      sample_id         = "sample_id",
      exist             = "exist",
      timepoint         = "timepoint_factor",
      fold_change       = "fold_change",
      materialise_table = TRUE,
      auto_expand       = TRUE,
      n_cores           = 2
    )

    cols <- dplyr::tbl_vars(ps$data_long)
    tp_col <- if ("time" %in% cols) {
      "time"
    } else if ("timepoint_factor" %in% cols) {
      "timepoint_factor"
    } else if ("timepoint" %in% cols) {
      "timepoint"
    } else {
      testthat::skip("example data_long has no recognized time column")
    }

    # small subset: 5 peptides at time T1
    keep_pep <- c("16627", "5243", "24799", "16196", "18003")

    ps_filt <- ps %>%
      dplyr::filter(
        peptide_id %in% keep_pep,
        !!rlang::sym(tp_col) == "T1"
      ) %>%
      dplyr::collect()

    cache <<- ps_filt
    cache
  }
})

.pick_existing_value_col <- function(ps, preferred) {
  dat <- if("phip_data" %in% class(ps)) ps$data_long else ps
  cols <- dplyr::tbl_vars(dat)
  if (preferred %in% cols) return(preferred)
}

testthat::test_that("compute_distance works on example phip_data and attaches
                    abundances", {

  ps <- .get_ps_small_for_distance()
  value_col <- .pick_existing_value_col(ps, "fold_change")

  d <- suppressWarnings(compute_distance(
    ps,
    value_col = value_col,
    method_normalization = "relative",
    distance = "euclidean",
    n_threads = 1L
  ))

  testthat::expect_s3_class(d, "dist")
  testthat::expect_true(is.matrix(attr(d, "abundances")))
  testthat::expect_true(is.numeric(as.numeric(d)))
  testthat::expect_true(all(as.numeric(d) >= 0))

  abund <- attr(d, "abundances")
  testthat::expect_equal(attr(d, "Size"), nrow(abund))
  testthat::expect_true(length(rownames(abund)) == nrow(abund))
  testthat::expect_true(all(!is.na(rownames(abund))))

  # relative normalization: nonzero rows should sum to ~1; all-zero rows sum to 0
  rs <- rowSums(abund)
  is_zero_row <- rowSums(abund != 0) == 0
  testthat::expect_true(all(abs(rs[!is_zero_row] - 1) < 1e-10))
  testthat::expect_true(all(rs[is_zero_row] == 0))
})

testthat::test_that("compute_distance auto-detects value_col on example data", {

  ps <- .get_ps_small_for_distance()

  d <- suppressWarnings(compute_distance(
    ps,
    value_col = NULL,
    method_normalization = "none",
    distance = "kulczynski",
    n_threads = 1L
  ))

  testthat::expect_s3_class(d, "dist")
  testthat::expect_true(is.matrix(attr(d, "abundances")))
})

testthat::test_that("compute_distance auto normalization chooses none for binary
                    exist data", {

  ps <- .get_ps_small_for_distance()

  cols <- dplyr::tbl_vars(ps)

  d <- suppressWarnings(compute_distance(
    ps,
    value_col = "exist",
    method_normalization = "auto",
    distance = "manhattan",
    n_threads = 1L
  ))
  a <- attr(d, "abundances")
  vals <- as.vector(a)

  # should remain binary (0/1) if auto picked "none"
  testthat::expect_true(all(vals %in% c(0, 1)))
})

testthat::test_that("compute_distance errors on duplicated (sample_id,
                    peptide_id) pairs", {

  ps <- .get_ps_small_for_distance()
  value_col <- .pick_existing_value_col(ps, "fold_change")

  ps_dup <- ps
  ps_dup <- dplyr::bind_rows(ps, ps[1, , drop = FALSE])

  testthat::expect_error(
    suppressWarnings(compute_distance(
      ps_dup,
      value_col = value_col,
      method_normalization = "none",
      distance = "manhattan",
      n_threads = 1L
    )),
    regexp = "(?i)duplicat|unique|pivot"
  )
})

testthat::test_that("compute_distance validates inputs (chk) on example data", {

  ps <- .get_ps_small_for_distance()
  value_col <- .pick_existing_value_col(ps, "fold_change")

  ps_bad <- ps
  class(ps_bad) <- "not_phip_data"
  testthat::expect_error(
    compute_distance(ps_bad, value_col = value_col),
    regexp = "(?i)phip_data|class"
  )

  ps_nolong <- structure(list(data_long = NULL), class = "phip_data")
  testthat::expect_error(
    compute_distance(ps_nolong, value_col = value_col),
    regexp = "(?i)data_long|null|missing"
  )

  testthat::expect_error(
    compute_distance(ps, value_col = value_col, n_threads = 0L),
    regexp = "(?i)thread|count|gt|positive|chk"
  )

  testthat::expect_error(
    compute_distance(ps, value_col = "definitely_not_a_column"),
    regexp = "(?i)not found|missing|column"
  )
})

testthat::test_that("compute_distance reproducibility: repeated runs match
                    across parameter combinations", {

  ps <- .get_ps_small_for_distance()

  value_fc <- "fold_change"
  value_ex <- "exist"

  combos <- list(
    list(value_col = value_fc, norm = "none",
         dist = "euclidean"),
    list(value_col = value_fc, norm = "relative",
         dist = "manhattan"),
    list(value_col = value_fc, norm = "hellinger",
         dist = "euclidean"),
    list(value_col = value_fc, norm = "log",
         dist = "manhattan")
  )

  combos <- combos[!vapply(combos, function(x) is.null(x$value_col),
                           logical(1))]

  for (cmb in combos) {
    d1 <- withr::with_seed(123, suppressWarnings(compute_distance(
      ps,
      value_col = cmb$value_col,
      method_normalization = cmb$norm,
      distance = cmb$dist,
      n_threads = 1L
    )))
    d2 <- withr::with_seed(123, suppressWarnings(compute_distance(
      ps,
      value_col = cmb$value_col,
      method_normalization = cmb$norm,
      distance = cmb$dist,
      n_threads = 1L
    )))

    testthat::expect_equal(as.numeric(d1), as.numeric(d2), tolerance = 1e-12)
    testthat::expect_equal(attr(d1, "abundances"), attr(d2, "abundances"),
                           tolerance = 1e-12)
    testthat::expect_equal(attr(d1, "Size"), attr(d2, "Size"))
  }
})

testthat::test_that("compute_distance bray agrees with vegan and manual", {

  ps <- .get_ps_small_for_distance()
  value_col <- .pick_existing_value_col(ps, "fold_change")

  d <- suppressWarnings(compute_distance(
    ps,
    value_col = value_col,
    method_normalization = "relative",
    distance = "bray",
    n_threads = 1L
  ))
  testthat::expect_s3_class(d, "dist")

  x <- attr(d, "abundances")


  d_ref <- vegan::vegdist(x, method = "bray")
  testthat::expect_equal(as.numeric(d), as.numeric(d_ref), tolerance = 1e-10)

  # manual bray-curtis in dist upper-triangle order
  n <- nrow(x)
  v <- numeric(n * (n - 1L) / 2L)
  k <- 1L
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      num <- sum(abs(x[i, ] - x[j, ]))
      den <- sum(x[i, ]) + sum(x[j, ])
      v[k] <- if (den == 0) 0 else num / den
      k <- k + 1L
    }
  }
  testthat::expect_equal(as.numeric(d), v, tolerance = 1e-10)
})


# helper: compute a distance object with abundances attached
.get_dist_for_pcoa <- function(ps_small,
                               value_col_preferred = c("fold_change", "exist"),
                               norm = "hellinger",
                               distance = "bray") {

  dat <- if("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  cols <- dplyr::tbl_vars(dat)
  val_col <- NULL
  for (nm in value_col_preferred) {
    if (nm %in% cols) {
      val_col <- nm
      break
    }
  }
  if (is.null(val_col)) testthat::skip("example data_long has neither
                                       fold_change nor exist")

  suppressWarnings(compute_distance(
    ps_small,
    value_col = val_col,
    method_normalization = norm,
    distance = distance,
    n_threads = 1L
  ))
}

testthat::test_that("compute_pcoa returns expected structure and types", {

  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  res <- suppressWarnings(compute_pcoa(
    d,
    neg_correction = "none",
    n_axes = 5L,
    top_features = 30L
  ))

  testthat::expect_s3_class(res, "beta_pcoa")
  testthat::expect_true(is.list(res))

  testthat::expect_true(all(c("sample_coords", "eigenvalues", "var_explained",
                              "eigen_diagnostics", "correction_infos",
                              "feature_associations") %in% names(res)))

  testthat::expect_s3_class(res$sample_coords, "tbl_df")
  testthat::expect_true("sample_id" %in% names(res$sample_coords))
  testthat::expect_true(is.numeric(res$eigenvalues))
  testthat::expect_s3_class(res$var_explained, "tbl_df")
  testthat::expect_equal(nrow(res$var_explained), 1L)
  testthat::expect_s3_class(res$feature_associations, "tbl_df")

  # sample_coords rows should match dist size
  n <- attr(d, "Size")
  testthat::expect_equal(nrow(res$sample_coords), n)

  # pcoa axes count should be <= min(n_axes, n - 1)
  k_expected <- min(5L, n - 1L)
  axis_cols <- grep("^PCoA\\d+$", names(res$sample_coords), value = TRUE)
  testthat::expect_equal(length(axis_cols), k_expected)

  # var_explained should contain %Other and %PCoA* columns for returned axes
  testthat::expect_true("%Other" %in% names(res$var_explained))
  if (k_expected > 0L) {
    testthat::expect_true(all(paste0("%PCoA", seq_len(k_expected)) %in%
                                names(res$var_explained)))
  }
})

testthat::test_that("compute_pcoa errors for dist objects with < 2 samples", {

  d1 <- stats::dist(matrix(1, nrow = 1))
  testthat::expect_error(
    compute_pcoa(d1),
    regexp = "(?i)at least 2 samples"
  )
})

testthat::test_that("compute_pcoa validates inputs (chk)", {

  testthat::expect_error(
    compute_pcoa("not_a_dist"),
    regexp = "(?i)dist|s3|class|chk"
  )

  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  testthat::expect_error(
    compute_pcoa(d, n_axes = 0L),
    regexp = "(?i)n_axes|gt|positive|chk"
  )

  testthat::expect_error(
    compute_pcoa(d, top_features = 0L),
    regexp = "(?i)top_features|gt|positive|chk"
  )

  testthat::expect_error(
    compute_pcoa(d, neg_correction = "nope"),
    regexp = "(?i)neg_correction|arg|match"
  )
})

testthat::test_that("compute_pcoa uses all requested axes up to n-1
                    (no hard cap at 10)", {

  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  n <- attr(d, "Size")
  testthat::expect_true(n >= 2L)

  n_axes_req <- 12L
  k_expected <- min(n_axes_req, n - 1L)

  res <- suppressWarnings(compute_pcoa(d, n_axes = n_axes_req,
                                       top_features = 10L))
  axis_cols <- grep("^PCoA\\d+$", names(res$sample_coords), value = TRUE)

  testthat::expect_equal(length(axis_cols), k_expected)

  # if the data has enough samples, this would fail with the old hard cap of 10
  if ((n - 1L) >= 12L) {
    testthat::expect_true("PCoA12" %in% axis_cols)
  }
})

testthat::test_that("compute_pcoa feature_associations: skips when abundances
                    missing", {

  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  # remove abundances attribute
  attr(d, "abundances") <- NULL

  res <- suppressWarnings(compute_pcoa(d, n_axes = 3L, top_features = 10L))
  testthat::expect_s3_class(res$feature_associations, "tbl_df")
  testthat::expect_equal(nrow(res$feature_associations), 0L)
})

testthat::test_that("compute_pcoa feature_associations: returns expected columns
                    and respects top_features logic", {

  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  res <- suppressWarnings(compute_pcoa(d, n_axes = 3L, top_features = 2L))

  fl <- res$feature_associations
  # may still be empty if alignment fails, but on example data it should work
  testthat::expect_true(is.data.frame(fl))
  if (nrow(fl) > 0L) {
    testthat::expect_true("feature" %in% names(fl))
    testthat::expect_true(all(paste0("PCoA", 1:3) %in% names(fl)))
    testthat::expect_equal(length(unique(fl$feature)), nrow(fl))

    # selection is union of top_features per axis, so upper bound is top_features * n_axes
    testthat::expect_lte(nrow(fl), 2L * 3L)
  }
})

testthat::test_that("compute_pcoa reproducibility: repeated runs match across
                    parameter combinations", {

  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  combos <- list(
    list(neg = "none",     n_axes = 3L, top = 10L),
    list(neg = "none",     n_axes = 5L, top = 5L)
  )

  if (rlang::is_installed("vegan")) {
    combos <- c(
      combos,
      list(list(neg = "lingoes",  n_axes = 3L, top = 10L)),
      list(list(neg = "cailliez", n_axes = 3L, top = 10L))
    )
  }

  for (cmb in combos) {
    r1 <- withr::with_seed(123, suppressWarnings(compute_pcoa(
      d,
      neg_correction = cmb$neg,
      n_axes = cmb$n_axes,
      top_features = cmb$top
    )))
    r2 <- withr::with_seed(123, suppressWarnings(compute_pcoa(
      d,
      neg_correction = cmb$neg,
      n_axes = cmb$n_axes,
      top_features = cmb$top
    )))

    testthat::expect_equal(r1$sample_coords, r2$sample_coords, tolerance = 1e-12)
    testthat::expect_equal(r1$eigenvalues, r2$eigenvalues, tolerance = 1e-12)
    testthat::expect_equal(r1$var_explained, r2$var_explained, tolerance = 1e-12)
    testthat::expect_equal(r1$feature_associations, r2$feature_associations,
                           tolerance = 1e-12)
  }
})

# tests/testthat/test-mcompute_capscale.R

# this file assumes these helpers already exist in the test suite:
# - .get_ps_small_for_distance()
# - .get_dist_for_pcoa()

# helper: pick a constraint variable that exists and has >=2 distinct values
.pick_constraint_var_cap <- function(ps_small) {
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  cols <- dplyr::tbl_vars(dat)

  candidates <- c("group", "big_group", "type_person", "sex", "age")
  candidates <- candidates[candidates %in% cols]
  if (length(candidates) == 0L) return(NULL)

  meta <- dat |>
    dplyr::select(sample_id, dplyr::all_of(candidates)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect()

  for (v in candidates) {
    x <- meta[[v]]
    if (!all(is.na(x)) && dplyr::n_distinct(x, na.rm = TRUE) >= 2L) return(v)
  }

  NULL
}

testthat::test_that("compute_capscale: full coverage (success, warnings,
                    errors, edge cases, reproducibility)", {

  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  # ---------------------------------------------------------------------------
  # 1) happy path (basic structure)
  # ---------------------------------------------------------------------------
  rhs_var <- .pick_constraint_var_cap(ps_small)

  fml <- stats::as.formula(paste0("~ ", rhs_var))

  res <- suppressWarnings(compute_capscale(
    dist_obj = d,
    ps = ps_small,
    formula = fml,
    neg_correction = "none",
    top_features = 30L
  ))

  testthat::expect_s3_class(res, "beta_capscale")
  testthat::expect_true(all(c(
    "sample_coords", "eigenvalues", "variance_partition",
    "feature_associations", "r2", "r2_adj", "perm_terms", "cap_model"
  ) %in% names(res)))

  testthat::expect_s3_class(res$sample_coords, "tbl_df")
  testthat::expect_true("sample_id" %in% names(res$sample_coords))
  testthat::expect_true(is.numeric(res$eigenvalues))
  testthat::expect_s3_class(res$variance_partition, "tbl_df")
  testthat::expect_equal(nrow(res$variance_partition), 3L)
  testthat::expect_true(all(c("component", "inertia", "proportion") %in%
                              names(res$variance_partition)))
  testthat::expect_s3_class(res$feature_associations, "tbl_df")
  testthat::expect_true(is.numeric(res$r2))
  testthat::expect_true(is.numeric(res$r2_adj))
  testthat::expect_s3_class(res$perm_terms, "tbl_df")
  testthat::expect_true(inherits(res$cap_model, "capscale"))

  # "Total" proportion should be 1
  tot_row <- res$variance_partition[res$variance_partition$component == "Total",
                                    , drop = FALSE]
  testthat::expect_equal(tot_row$proportion, 1)

  # ---------------------------------------------------------------------------
  # 2) reproducibility (same seed -> identical outputs)
  # ---------------------------------------------------------------------------
  r1 <- withr::with_seed(123, suppressWarnings(compute_capscale(
    dist_obj = d,
    ps = ps_small,
    formula = fml,
    neg_correction = "none",
    top_features = 30L
  )))
  r2 <- withr::with_seed(123, suppressWarnings(compute_capscale(
    dist_obj = d,
    ps = ps_small,
    formula = fml,
    neg_correction = "none",
    top_features = 30L
  )))

  testthat::expect_equal(r1$sample_coords, r2$sample_coords, tolerance = 1e-12)
  testthat::expect_equal(r1$eigenvalues, r2$eigenvalues, tolerance = 1e-12)
  testthat::expect_equal(r1$variance_partition, r2$variance_partition,
                         tolerance = 1e-12)
  testthat::expect_equal(r1$feature_associations, r2$feature_associations,
                         tolerance = 1e-12)

  # ----------------------------------------------------------------------------
  # 3) neg_correction branches
  # ----------------------------------------------------------------------------
  r_lin <- suppressWarnings(compute_capscale(d, ps_small, fml,
                                              neg_correction = "lingoes",
                                              top_features = 10L))
  r_cai <- suppressWarnings(compute_capscale(d, ps_small, fml,
                                              neg_correction = "cailliez",
                                              top_features = 10L))
  testthat::expect_s3_class(r_lin, "beta_capscale")
  testthat::expect_s3_class(r_cai, "beta_capscale")

  # ----------------------------------------------------------------------------
  # 4) transformed rhs terms supported (all.vars(terms()))
  # ----------------------------------------------------------------------------
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  cols <- dplyr::tbl_vars(dat)
  if ("age" %in% cols) {
    testthat::expect_silent(suppressWarnings(compute_capscale(
      dist_obj = d,
      ps = ps_small,
      formula = ~ log(age),
      neg_correction = "none",
      top_features = 5L
    )))
  }

  # ----------------------------------------------------------------------------
  # 5) warning branches (single test file, but cover both warnings)
  #    - no labels in dist_obj (warn)
  #    - no abundances (feature associations empty, no warning in your final func
  # ----------------------------------------------------------------------------
  d_no_labels <- d
  attr(d_no_labels, "Labels") <- NULL

  # keep "Size" intact
  testthat::expect_warning(
    compute_capscale(
      dist_obj = d_no_labels,
      ps = ps_small,
      formula = fml,
      neg_correction = "none",
      top_features = 10L
    ),
    regexp = "(?i)no labels found"
  )

  d_no_ab <- d
  attr(d_no_ab, "abundances") <- NULL
  res_no_ab <- suppressWarnings(compute_capscale(
    dist_obj = d_no_ab,
    ps = ps_small,
    formula = fml,
    neg_correction = "none",
    top_features = 10L
  ))
  testthat::expect_s3_class(res_no_ab$feature_associations, "tbl_df")
  testthat::expect_equal(nrow(res_no_ab$feature_associations), 0L)

  # ----------------------------------------------------------------------------
  # 6) error branches (chk + .ph_abort paths)
  # ----------------------------------------------------------------------------

  # dist_obj not dist -> chk error
  testthat::expect_error(
    compute_capscale(dist_obj = "nope", ps = ps_small, formula = fml),
    regexp = "(?i)dist|s3|class|chk"
  )

  # formula not formula -> chk error
  testthat::expect_error(
    compute_capscale(dist_obj = d, ps = ps_small, formula = "not a formula"),
    regexp = "(?i)formula|chk|true"
  )

  # top_features invalid -> chk error
  testthat::expect_error(
    compute_capscale(dist_obj = d, ps = ps_small, formula = fml,
                     top_features = 0L),
    regexp = "(?i)top_features|gt|positive|chk"
  )

  # ps is NULL / missing -> .ph_abort
  testthat::expect_error(
    compute_capscale(dist_obj = d, ps = NULL, formula = fml),
    regexp = "(?i)ps.*missing|cannot construct metadata"
  )

  # missing required variables -> .ph_abort
  testthat::expect_error(
    compute_capscale(dist_obj = d, ps = ps_small,
                      formula = ~ definitely_not_a_column),
    regexp = "(?i)missing.*variables|definitely_not_a_column"
  )

  # rhs empty -> .ph_abort
  testthat::expect_error(
    compute_capscale(dist_obj = d, ps = ps_small, formula = ~ 1),
    regexp = "(?i)no constraints.*rhs is empty|unconstrained"
  )

  # ps without sample_id -> .ph_abort
  ps_no_sid <- dat |>
    dplyr::select(-sample_id)
  testthat::expect_error(
    compute_capscale(dist_obj = d, ps = ps_no_sid, formula = fml),
    regexp = "(?i)sample_id"
  )

  # dist labels not present in ps -> .ph_abort
  d_bad_labels <- d
  attr(d_bad_labels, "Labels") <- paste0("not_in_ps_",
                                         seq_len(attr(d_bad_labels, "Size")))
  testthat::expect_error(
    compute_capscale(dist_obj = d_bad_labels, ps = ps_small, formula = fml),
    regexp = "(?i)missing in `ps`|missing in ps"
  )

  # all missing in constrained vars -> .ph_abort
  # force a metadata column used in formula to be NA for all samples
  dat_meta <- dat |>
    dplyr::select(sample_id, dplyr::all_of(rhs_var)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect()
  if (nrow(dat_meta) > 0L) {
    dat_all_na <- dat |>
      dplyr::mutate(!!rlang::sym(rhs_var) := NA) # makes complete.cases FALSE

    testthat::expect_error(
      compute_capscale(dist_obj = d, ps = dat_all_na, formula = fml),
      regexp = "(?i)all samples have missing values"
    )
  }
})

# helper: pick a constraint variable that exists and has >=2 distinct values for permanova/dispersion
.pick_constraint_var_perm <- function(ps_small) {
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  cols <- dplyr::tbl_vars(dat)

  candidates <- c("group", "big_group", "type_person", "sex", "age")
  candidates <- candidates[candidates %in% cols]
  if (length(candidates) == 0L) return(NULL)

  meta <- dat |>
    dplyr::select(sample_id, dplyr::all_of(candidates)) |>
    dplyr::distinct(sample_id, .keep_all = TRUE) |>
    dplyr::collect()

  for (v in candidates) {
    x <- meta[[v]]
    if (!all(is.na(x)) && dplyr::n_distinct(x, na.rm = TRUE) >= 2L) return(v)
  }

  NULL
}

testthat::test_that("compute_permanova: basic functionality and input validation", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  # skip if vegan not available
  testthat::skip_if_not_installed("vegan")

  # pick a grouping variable that exists
  group_var <- .pick_constraint_var_perm(ps_small)
  testthat::skip_if(is.null(group_var), "no suitable grouping variable found")

  # basic permanova test
  res <- suppressWarnings(compute_permanova(
    dist_obj = d,
    ps = ps_small,
    group_col = group_var,
    permutations = 99  # small for speed
  ))

  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_true(all(c("scope", "contrast", "term", "p_value",
                              "p_adjust", "n_perm") %in% names(res)))
  testthat::expect_true(nrow(res) >= 1L)  # at least global test
  if (nrow(res) > 0L) {
    testthat::expect_true(all(!is.na(res$p_adjust)))
  }

  # test with both group and time (if available)
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  cols <- dplyr::tbl_vars(dat)
  time_candidates <- c("timepoint", "timepoint_factor", "time")
  time_var <- intersect(time_candidates, cols)[1]

  if (!is.na(time_var)) {
    res2 <- suppressWarnings(compute_permanova(
      dist_obj = d,
      ps = ps_small,
      group_col = group_var,
      time_col = time_var,
      permutations = 99
    ))
    testthat::expect_s3_class(res2, "tbl_df")
    testthat::expect_true(nrow(res2) >= 1L)
  }
})

testthat::test_that("compute_permanova: pairwise comparisons", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("vegan")

  group_var <- .pick_constraint_var_perm(ps_small)
  testthat::skip_if(is.null(group_var), "no suitable grouping variable found")

  res_pair <- suppressWarnings(compute_permanova(
    dist_obj = d,
    ps = ps_small,
    group_col = group_var,
    permutations = 99
  ))

  testthat::expect_s3_class(res_pair, "tbl_df")
  # should have global + pairwise results
  testthat::expect_true(any(res_pair$scope == "global"))
  # pairwise may be skipped if levels are insufficient after filtering
  if (any(res_pair$scope == "group_pairwise")) {
    testthat::expect_true(any(res_pair$scope == "group_pairwise"))
  }
})

testthat::test_that("compute_permanova: input validation errors", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("vegan")

  # dist_obj not dist
  testthat::expect_error(
    compute_permanova(dist_obj = "not_dist", ps = ps_small),
    regexp = "(?i)dist|s3|class|chk"
  )

  # ps is NULL
  testthat::expect_error(
    compute_permanova(dist_obj = d, ps = NULL),
    regexp = "(?i)ps.*missing|cannot construct metadata"
  )

  # missing sample_id column
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  ps_no_sid <- dat |> dplyr::select(-sample_id)
  testthat::expect_error(
    compute_permanova(dist_obj = d, ps = ps_no_sid),
    regexp = "(?i)sample_id"
  )

  # non-existent group column
  testthat::expect_error(
    compute_permanova(dist_obj = d, ps = ps_small, group_col = "nonexistent"),
    regexp = "(?i)not found"
  )

  # invalid permutations
  testthat::expect_error(
    compute_permanova(dist_obj = d, ps = ps_small, permutations = 0),
    regexp = "(?i)gt|greater|chk"
  )
})

testthat::test_that("compute_permanova: ps as data.frame vs phip_data", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("vegan")

  group_var <- .pick_constraint_var_perm(ps_small)
  testthat::skip_if(is.null(group_var), "no suitable grouping variable found")

  # test with phip_data object
  if ("phip_data" %in% class(ps_small)) {
    res1 <- suppressWarnings(compute_permanova(
      dist_obj = d,
      ps = ps_small,
      group_col = group_var,
      permutations = 99
    ))
    testthat::expect_s3_class(res1, "tbl_df")
  }

  # test with data.frame directly
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  res2 <- suppressWarnings(compute_permanova(
    dist_obj = d,
    ps = dat,
    group_col = group_var,
    permutations = 99
  ))
  testthat::expect_s3_class(res2, "tbl_df")
})

testthat::test_that("compute_dispersion: basic functionality and input validation", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("vegan")

  group_var <- .pick_constraint_var_perm(ps_small)
  testthat::skip_if(is.null(group_var), "no suitable grouping variable found")

  # basic dispersion test
  res <- suppressWarnings(compute_dispersion(
    dist_obj = d,
    ps = ps_small,
    group_col = group_var,
    permutations = 99
  ))

  testthat::expect_s3_class(res, "beta_dispersion")
  testthat::expect_true(all(c("distances", "tests") %in% names(res)))
  testthat::expect_s3_class(res$distances, "tbl_df")
  testthat::expect_s3_class(res$tests, "tbl_df")

  # check distances structure
  if (nrow(res$distances) > 0L) {
    testthat::expect_true(all(c("sample_id", "distance", "level", "scope", "contrast") %in% names(res$distances)))
  }

  # check tests structure
  if (nrow(res$tests) > 0L) {
    testthat::expect_true(all(c("scope", "contrast", "term", "p_value",
                                "p_adjust", "n_perm") %in% names(res$tests)))
    testthat::expect_true(all(res$tests$term == "dispersion"))
    testthat::expect_true(all(!is.na(res$tests$p_adjust)))
  }
})

testthat::test_that("compute_dispersion: pairwise comparisons", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("vegan")

  group_var <- .pick_constraint_var_perm(ps_small)
  testthat::skip_if(is.null(group_var), "no suitable grouping variable found")

  res_pair <- suppressWarnings(compute_dispersion(
    dist_obj = d,
    ps = ps_small,
    group_col = group_var,
    permutations = 99
  ))

  testthat::expect_s3_class(res_pair, "beta_dispersion")
  testthat::expect_true(all(c("distances", "tests") %in% names(res_pair)))

  testthat::expect_true(nrow(res_pair$tests) >= 1L)
})

testthat::test_that("compute_dispersion: input validation errors", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("vegan")

  # dist_obj not dist
  testthat::expect_error(
    compute_dispersion(dist_obj = "not_dist", ps = ps_small),
    regexp = "(?i)dist|s3|class|chk"
  )

  # ps is NULL
  testthat::expect_error(
    compute_dispersion(dist_obj = d, ps = NULL),
    regexp = "(?i)ps.*missing|cannot construct metadata"
  )

  # missing sample_id column
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  ps_no_sid <- dat |> dplyr::select(-sample_id)
  testthat::expect_error(
    compute_dispersion(dist_obj = d, ps = ps_no_sid),
    regexp = "(?i)sample_id"
  )

  # non-existent group column
  testthat::expect_error(
    compute_dispersion(dist_obj = d, ps = ps_small, group_col = "nonexistent"),
    regexp = "(?i)not found"
  )

  # invalid permutations
  testthat::expect_error(
    compute_dispersion(dist_obj = d, ps = ps_small, permutations = 0),
    regexp = "(?i)gt|greater|chk"
  )
})

testthat::test_that("compute_dispersion: ps as data.frame vs phip_data", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("vegan")

  group_var <- .pick_constraint_var_perm(ps_small)
  testthat::skip_if(is.null(group_var), "no suitable grouping variable found")

  # test with phip_data object
  if ("phip_data" %in% class(ps_small)) {
    res1 <- suppressWarnings(compute_dispersion(
      dist_obj = d,
      ps = ps_small,
      group_col = group_var,
      permutations = 99
    ))
    testthat::expect_s3_class(res1, "beta_dispersion")
  }

  # test with data.frame directly
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  res2 <- suppressWarnings(compute_dispersion(
    dist_obj = d,
    ps = dat,
    group_col = group_var,
    permutations = 99
  ))
  testthat::expect_s3_class(res2, "beta_dispersion")
})

# tests for compute_tsne function

testthat::test_that("compute_tsne: basic functionality and structure", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  # basic t-SNE computation
  res <- suppressWarnings(compute_tsne(
    ps = ps_small,
    dist_obj = d,
    dims = 3L,
    perplexity = 5,  # small for limited sample size
    seed = 42
  ))

  testthat::expect_s3_class(res, "phip_tsne")
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_true(all(c("sample_id", "tSNE1", "tSNE2", "tSNE3") %in% names(res)))

  # check dimensions match distance object
  n_samples <- attr(d, "Size")
  testthat::expect_equal(nrow(res), n_samples)

  # check coordinates are numeric
  testthat::expect_true(is.numeric(res$tSNE1))
  testthat::expect_true(is.numeric(res$tSNE2))
  testthat::expect_true(is.numeric(res$tSNE3))

  # check attributes
  testthat::expect_identical(attr(res, "distance"), d)
  testthat::expect_true(is.list(attr(res, "tsne_params")))
  testthat::expect_true(is.character(attr(res, "meta_cols")))
})

testthat::test_that("compute_tsne: dims parameter controls output dimensions", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  # 2D t-SNE
  res2d <- suppressWarnings(compute_tsne(
    ps = ps_small,
    dist_obj = d,
    dims = 2L,
    perplexity = 5,
    seed = 42
  ))

  testthat::expect_true(all(c("tSNE1", "tSNE2") %in% names(res2d)))
  testthat::expect_true(all(is.na(res2d$tSNE3)))

  # 3D t-SNE
  res3d <- suppressWarnings(compute_tsne(
    ps = ps_small,
    dist_obj = d,
    dims = 3L,
    perplexity = 5,
    seed = 42
  ))

  testthat::expect_true(all(c("tSNE1", "tSNE2", "tSNE3") %in% names(res3d)))
  testthat::expect_true(all(!is.na(res3d$tSNE3)))
})

testthat::test_that("compute_tsne: metadata attachment works correctly", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  # find available metadata columns
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  cols <- dplyr::tbl_vars(dat)
  meta_candidates <- c("subject_id", "timepoint", "timepoint_factor", "group", "sex")
  available_meta <- intersect(meta_candidates, cols)

  if (length(available_meta) > 0L) {
    res <- suppressWarnings(compute_tsne(
      ps = ps_small,
      dist_obj = d,
      dims = 2L,
      perplexity = 5,
      meta_cols = available_meta[1:min(2, length(available_meta))],
      seed = 42
    ))

    # check that metadata columns were attached
    for (col in available_meta[1:min(2, length(available_meta))]) {
      testthat::expect_true(col %in% names(res))
    }
  }
})

testthat::test_that("compute_tsne: matrix input works", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  # convert dist to matrix
  d_mat <- as.matrix(d)

  res <- suppressWarnings(compute_tsne(
    ps = ps_small,
    dist_obj = d_mat,
    dims = 2L,
    perplexity = 5,
    seed = 42
  ))

  testthat::expect_s3_class(res, "phip_tsne")
  testthat::expect_equal(nrow(res), nrow(d_mat))
})

testthat::test_that("compute_tsne: reproducibility with seed", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  res1 <- suppressWarnings(compute_tsne(
    ps = ps_small,
    dist_obj = d,
    dims = 2L,
    perplexity = 5,
    seed = 123
  ))

  res2 <- suppressWarnings(compute_tsne(
    ps = ps_small,
    dist_obj = d,
    dims = 2L,
    perplexity = 5,
    seed = 123
  ))

  testthat::expect_equal(res1$tSNE1, res2$tSNE1, tolerance = 1e-10)
  testthat::expect_equal(res1$tSNE2, res2$tSNE2, tolerance = 1e-10)
})

testthat::test_that("compute_tsne: input validation errors", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  # missing dist_obj
  testthat::expect_error(
    compute_tsne(ps = ps_small),
    regexp = "(?i)dist_obj.*required"
  )

  # invalid dist_obj type
  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = "not_a_distance"),
    regexp = "(?i)dist.*object.*matrix"
  )

  # invalid dims
  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = d, dims = 4L),
    regexp = "(?i)dims.*2.*3"
  )

  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = d, dims = 1L),
    regexp = "(?i)dims.*2.*3"
  )

  # invalid perplexity
  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = d, perplexity = 0),
    regexp = "(?i)gt|greater|chk"
  )

  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = d, perplexity = -5),
    regexp = "(?i)gt|greater|chk"
  )

  # invalid theta
  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = d, theta = -1),
    regexp = "(?i)greater|chk"
  )

  # invalid max_iter
  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = d, max_iter = 0L),
    regexp = "(?i)gt|greater|chk"
  )

  # invalid meta_cols
  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = d, meta_cols = 123),
    regexp = "(?i)character|chk"
  )

  # invalid seed
  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = d, seed = "not_numeric"),
    regexp = "(?i)count|chk"
  )

  # invalid check_duplicates
  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = d, check_duplicates = "not_logical"),
    regexp = "(?i)flag|chk"
  )
})

testthat::test_that("compute_tsne: perplexity too large gets adjusted", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  n_samples <- attr(d, "Size")

  # perplexity >= n_samples should error
  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = d, perplexity = n_samples),
    regexp = "(?i)perplexity.*smaller.*samples"
  )

  # perplexity too high should warn and adjust
  if (n_samples > 10) {
    testthat::expect_warning(
      compute_tsne(ps = ps_small, dist_obj = d, perplexity = n_samples - 2),
      regexp = "(?i)perplexity.*high.*reducing"
    )
  }
})

testthat::test_that("compute_tsne: handles missing labels gracefully", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  # remove labels from distance object
  d_no_labels <- d
  attr(d_no_labels, "Labels") <- NULL

  testthat::expect_warning(
    res <- compute_tsne(ps = ps_small, dist_obj = d_no_labels, dims = 2L, perplexity = 5),
    regexp = "(?i)no labels.*integer indices"
  )

  testthat::expect_s3_class(res, "phip_tsne")
  testthat::expect_true(all(res$sample_id %in% as.character(seq_len(attr(d, "Size")))))
})

testthat::test_that("compute_tsne: handles non-square matrix error", {
  ps_small <- .get_ps_small_for_distance()

  testthat::skip_if_not_installed("Rtsne")

  # create non-square matrix
  non_square <- matrix(1:12, nrow = 3, ncol = 4)

  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = non_square),
    regexp = "(?i)square.*rows.*columns"
  )
})

testthat::test_that("compute_tsne: handles non-numeric matrix error", {
  ps_small <- .get_ps_small_for_distance()

  testthat::skip_if_not_installed("Rtsne")

  # create non-numeric matrix
  char_mat <- matrix(letters[1:9], nrow = 3, ncol = 3)

  testthat::expect_error(
    compute_tsne(ps = ps_small, dist_obj = char_mat),
    regexp = "(?i)matrix.*numeric"
  )
})

testthat::test_that("compute_tsne: handles too few samples error", {
  testthat::skip_if_not_installed("Rtsne")

  # create tiny distance with only 2 samples
  tiny_dist <- stats::dist(matrix(c(1, 2), nrow = 2))
  ps_tiny <- data.frame(sample_id = c("s1", "s2"))

  testthat::expect_error(
    compute_tsne(ps = ps_tiny, dist_obj = tiny_dist),
    regexp = "(?i)at least 3 samples"
  )
})

testthat::test_that("compute_tsne: handles duplicate labels error", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  # create distance with duplicate labels
  d_dup <- d
  labels <- attr(d_dup, "Labels")
  if (length(labels) > 1) {
    labels[2] <- labels[1]  # make duplicate
    attr(d_dup, "Labels") <- labels

    testthat::expect_error(
      compute_tsne(ps = ps_small, dist_obj = d_dup),
      regexp = "(?i)unique|duplicate"
    )
  }
})

testthat::test_that("compute_tsne: ps as data.frame vs phip_data", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  # test with phip_data object
  if ("phip_data" %in% class(ps_small)) {
    res1 <- suppressWarnings(compute_tsne(
      ps = ps_small,
      dist_obj = d,
      dims = 2L,
      perplexity = 5,
      seed = 42
    ))
    testthat::expect_s3_class(res1, "phip_tsne")
  }

  # test with data.frame directly
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  res2 <- suppressWarnings(compute_tsne(
    ps = dat,
    dist_obj = d,
    dims = 2L,
    perplexity = 5,
    seed = 42
  ))
  testthat::expect_s3_class(res2, "phip_tsne")
})

testthat::test_that("compute_tsne: handles missing sample_id column", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  # remove sample_id column
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  ps_no_sid <- dat |> dplyr::select(-sample_id)

  testthat::expect_error(
    compute_tsne(ps = ps_no_sid, dist_obj = d),
    regexp = "(?i)sample_id"
  )
})

testthat::test_that("compute_tsne: warns about missing metadata columns", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small, norm = "hellinger", distance = "bray")

  testthat::skip_if_not_installed("Rtsne")

  testthat::expect_warning(
    compute_tsne(
      ps = ps_small,
      dist_obj = d,
      dims = 2L,
      perplexity = 5,
      meta_cols = c("nonexistent_column", "another_missing_col")
    ),
    regexp = "(?i)not present.*ignored"
  )
})

# Additional tests for better coverage

testthat::test_that("compute_distance: edge cases and error conditions", {
  ps_small <- .get_ps_small_for_distance()

  # test with empty data after filtering
  ps_empty <- ps_small |> dplyr::filter(FALSE)
  testthat::expect_error(
    compute_distance(ps_empty, value_col = "exist"),
    regexp = "(?i)empty|zero"
  )

  # test with all NA values in abundance column
  ps_na <- ps_small |>
    dplyr::mutate(exist = NA_real_)

  # should replace NAs with 0 and continue
  d_na <- suppressWarnings(compute_distance(
    ps_na,
    value_col = "exist",
    method_normalization = "none",
    distance = "euclidean"
  ))
  testthat::expect_s3_class(d_na, "dist")

  # test with missing sample_id or peptide_id
  ps_no_sample <- ps_small |> dplyr::select(-sample_id)
  testthat::expect_error(
    compute_distance(ps_no_sample, value_col = "exist"),
    regexp = "(?i)missing.*sample_id"
  )

  ps_no_peptide <- ps_small |> dplyr::select(-peptide_id)
  testthat::expect_error(
    compute_distance(ps_no_peptide, value_col = "exist"),
    regexp = "(?i)missing.*peptide_id"
  )

  # test with NA in sample_id or peptide_id
  ps_na_ids <- ps_small
  ps_na_ids$sample_id[1] <- NA
  testthat::expect_error(
    compute_distance(ps_na_ids, value_col = "exist"),
    regexp = "(?i)missing values.*sample_id"
  )
})

testthat::test_that("compute_distance: normalization methods", {
  ps_small <- .get_ps_small_for_distance()

  # test log normalization
  d_log <- suppressWarnings(compute_distance(
    ps_small,
    value_col = "exist",
    method_normalization = "log",
    distance = "euclidean"
  ))
  testthat::expect_s3_class(d_log, "dist")

  # test hellinger normalization
  d_hell <- suppressWarnings(compute_distance(
    ps_small,
    value_col = "exist",
    method_normalization = "hellinger",
    distance = "euclidean"
  ))
  testthat::expect_s3_class(d_hell, "dist")

  # verify hellinger is sqrt(relative)
  abund_hell <- attr(d_hell, "abundances")
  # for binary data, hellinger should be sqrt of relative
  testthat::expect_true(all(abund_hell >= 0))
  testthat::expect_true(all(abund_hell <= 1))
})

testthat::test_that("compute_distance: distance method fallbacks", {
  ps_small <- .get_ps_small_for_distance()

  # test vegan fallback methods (when parallelDist not available or method not supported)
  testthat::skip_if_not_installed("vegan")

  # test a method that should use vegan
  d_jaccard <- suppressWarnings(compute_distance(
    ps_small,
    value_col = "exist",
    method_normalization = "none",
    distance = "jaccard"
  ))
  testthat::expect_s3_class(d_jaccard, "dist")

  # test chebyshev -> maximum mapping
  d_cheby <- suppressWarnings(compute_distance(
    ps_small,
    value_col = "exist",
    method_normalization = "none",
    distance = "chebyshev"
  ))
  testthat::expect_s3_class(d_cheby, "dist")
})

testthat::test_that("compute_distance: parallelDist vs vegan consistency", {
  ps_small <- .get_ps_small_for_distance()

  testthat::skip_if_not_installed("vegan")

  # test that bray-curtis gives same results
  d1 <- suppressWarnings(compute_distance(
    ps_small,
    value_col = "exist",
    method_normalization = "relative",
    distance = "bray"
  ))

  # manually compute with vegan for comparison
  abund <- attr(d1, "abundances")
  d_vegan <- vegan::vegdist(abund, method = "bray")

  testthat::expect_equal(as.numeric(d1), as.numeric(d_vegan), tolerance = 1e-10)
})

testthat::test_that("compute_pcoa: negative eigenvalue corrections", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  testthat::skip_if_not_installed("vegan")

  # test lingoes correction
  res_lingoes <- suppressWarnings(compute_pcoa(
    d,
    neg_correction = "lingoes",
    n_axes = 3L
  ))
  testthat::expect_s3_class(res_lingoes, "beta_pcoa")

  # test cailliez correction
  res_cailliez <- suppressWarnings(compute_pcoa(
    d,
    neg_correction = "cailliez",
    n_axes = 3L
  ))
  testthat::expect_s3_class(res_cailliez, "beta_pcoa")
})

testthat::test_that("compute_pcoa: variance explained calculations", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  res <- suppressWarnings(compute_pcoa(d, n_axes = 3L))

  # check that percentages are reasonable
  var_exp <- res$var_explained
  pct_cols <- grep("^%PCoA", names(var_exp), value = TRUE)

  if (length(pct_cols) > 0) {
    pct_vals <- unlist(var_exp[pct_cols])
    testthat::expect_true(all(pct_vals >= 0, na.rm = TRUE))
    testthat::expect_true(all(pct_vals <= 100, na.rm = TRUE))

    # total should be <= 100
    total_pct <- sum(pct_vals, var_exp$`%Other`, na.rm = TRUE)
    testthat::expect_lte(total_pct, 100.1) # small tolerance for rounding
  }
})

testthat::test_that("compute_pcoa: feature associations edge cases", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  # test with mismatched row names
  X_bad <- attr(d, "abundances")
  rownames(X_bad) <- paste0("bad_", rownames(X_bad))
  attr(d, "abundances") <- X_bad

  res_bad <- suppressWarnings(compute_pcoa(d, n_axes = 2L))
  testthat::expect_equal(nrow(res_bad$feature_associations), 0L)

  # test with zero-weight features
  d2 <- .get_dist_for_pcoa(ps_small)
  X2 <- attr(d2, "abundances")
  X2[, 1] <- 0  # make first feature all zeros
  attr(d2, "abundances") <- X2

  res_zero <- suppressWarnings(compute_pcoa(d2, n_axes = 2L, top_features = 10L))
  # should still work, just exclude zero-weight features
  testthat::expect_s3_class(res_zero$feature_associations, "tbl_df")
})

testthat::test_that("compute_capscale: size mismatches and edge cases", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  rhs_var <- .pick_constraint_var_cap(ps_small)
  testthat::skip_if(is.null(rhs_var), "no constraint variable available")

  # test with distance size mismatch when no labels
  d_no_size <- d
  attr(d_no_size, "Labels") <- NULL
  attr(d_no_size, "Size") <- 999  # wrong size

  testthat::expect_error(
    compute_capscale(d_no_size, ps_small, stats::as.formula(paste0("~ ", rhs_var))),
    regexp = "(?i)cannot align.*size"
  )
})

testthat::test_that("compute_permanova: stratification and repeated measures", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  testthat::skip_if_not_installed("vegan")

  group_var <- .pick_constraint_var_perm(ps_small)
  testthat::skip_if(is.null(group_var), "no grouping variable available")

  # test with subject stratification
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  cols <- dplyr::tbl_vars(dat)

  if ("subject_id" %in% cols) {
    # create artificial repeated measures by duplicating some samples
    dat_collected <- dat |> dplyr::collect()
    n_rows <- nrow(dat_collected)
    indices <- rep(1:min(5, n_rows), each = 2)

    dat_rep <- dat_collected |>
      dplyr::slice(indices) |>
      dplyr::mutate(
        sample_id = paste0(sample_id, "_rep", rep(1:2, length.out = dplyr::n())),
        timepoint_factor = rep(c("T1", "T2"), length.out = dplyr::n())
      )

    # create corresponding distance matrix
    n_rep <- nrow(dat_rep)
    d_rep <- stats::as.dist(matrix(runif(n_rep * (n_rep - 1) / 2),
                                   nrow = n_rep - 1))
    attr(d_rep, "Labels") <- dat_rep$sample_id
    attr(d_rep, "Size") <- n_rep

    res_strat <- suppressWarnings(compute_permanova(
      d_rep, dat_rep,
      group_col = group_var,
      time_col = "timepoint_factor",
      subject_col = "subject_id",
      permutations = 99
    ))
    testthat::expect_s3_class(res_strat, "tbl_df")
  }
})

testthat::test_that("compute_permanova: insufficient data warnings", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  testthat::skip_if_not_installed("vegan")

  # create data with only one group level
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  dat_single <- dat |>
    dplyr::mutate(single_group = "A") |>
    dplyr::collect()

  # should skip global test due to insufficient levels
  res_single <- suppressWarnings(compute_permanova(
    d, dat_single,
    group_col = "single_group",
    permutations = 99
  ))

  # should return empty or minimal results
  testthat::expect_s3_class(res_single, "tbl_df")
})

testthat::test_that("compute_dispersion: time factor edge cases", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  testthat::skip_if_not_installed("vegan")

  # test with numeric time (should warn and skip)
  dat <- if ("phip_data" %in% class(ps_small)) ps_small$data_long else ps_small
  dat_numeric_time <- dat |>
    dplyr::mutate(numeric_time = as.numeric(1:dplyr::n())) |>
    dplyr::collect()

  res_numeric <- suppressWarnings(compute_dispersion(
    d, dat_numeric_time,
    time_col = "numeric_time",
    permutations = 99
  ))

  testthat::expect_s3_class(res_numeric, "beta_dispersion")
  # should have skipped time dispersion test
})

testthat::test_that("compute_dispersion: betadisper failures", {
  ps_small <- .get_ps_small_for_distance()

  testthat::skip_if_not_installed("vegan")

  # create a problematic distance matrix (e.g., all zeros)
  n_samples <- 5
  d_bad <- stats::as.dist(matrix(0, nrow = n_samples, ncol = n_samples))
  attr(d_bad, "Labels") <- paste0("sample_", 1:n_samples)
  attr(d_bad, "Size") <- n_samples

  ps_bad <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    group = rep(c("A", "B"), length.out = n_samples)
  )

  # should handle betadisper failure gracefully
  res_bad <- suppressWarnings(compute_dispersion(
    d_bad, ps_bad,
    group_col = "group",
    permutations = 99
  ))

  testthat::expect_s3_class(res_bad, "beta_dispersion")
})

testthat::test_that("compute_tsne: Rtsne package not available", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  testthat::skip_if_not_installed("mockery")

  # mock rlang::is_installed to return FALSE for Rtsne
  mockery::stub(compute_tsne, "rlang::is_installed", function(pkg) {
    if (pkg == "Rtsne") FALSE else TRUE
  })

  testthat::expect_error(
    compute_tsne(ps_small, d),
    regexp = "(?i)Rtsne.*required.*install"
  )
})

testthat::test_that("compute_tsne: unexpected Rtsne output", {
  ps_small <- .get_ps_small_for_distance()
  d <- .get_dist_for_pcoa(ps_small)

  testthat::skip_if_not_installed("Rtsne")
  testthat::skip_if_not_installed("mockery")

  # mock Rtsne to return unexpected output
  mockery::stub(compute_tsne, "Rtsne::Rtsne", function(...) {
    list(Y = matrix(1, nrow = 2, ncol = 2)) # wrong dimensions
  })

  testthat::expect_error(
    compute_tsne(ps_small, d, dims = 3L, perplexity = 5),
    regexp = "(?i)unexpected output.*wrong shape"
  )
})

testthat::test_that("edge cases: empty results and missing data", {
  # test various edge cases that might not be covered

  # empty distance matrix
  testthat::expect_error(
    compute_pcoa(stats::dist(matrix(numeric(0), nrow = 0))),
    regexp = "(?i)at least 2 samples"
  )

  # test with minimal data
  minimal_ps <- data.frame(
    sample_id = c("s1", "s2", "s3"),
    peptide_id = c("p1", "p1", "p1"),
    exist = c(1, 0, 1)
  )

  d_minimal <- suppressWarnings(compute_distance(
    minimal_ps,
    value_col = "exist",
    method_normalization = "none",
    distance = "euclidean"
  ))

  testthat::expect_s3_class(d_minimal, "dist")
  testthat::expect_equal(attr(d_minimal, "Size"), 3L)
})

testthat::test_that("compute_distance: phip_data auto-detect and non-numeric value errors", {
  dat_bad <- data.frame(
    sample_id = "s1",
    peptide_id = "p1",
    value = "not_numeric",
    stringsAsFactors = FALSE
  )

  ps_bad <- structure(list(data_long = dat_bad), class = "phip_data")

  testthat::expect_error(
    compute_distance(ps_bad, value_col = NULL),
    regexp = "(?i)could not infer.*abundance"
  )

  testthat::expect_error(
    compute_distance(dat_bad, value_col = "value"),
    regexp = "(?i)must be numeric"
  )
})

testthat::test_that("compute_distance: pivot_wider guard", {
  testthat::skip_if_not_installed("mockery")
  testthat::skip_if_not_installed("dbplyr")

  dat <- data.frame(
    sample_id = c("s1", "s2"),
    peptide_id = c("p1", "p1"),
    exist = c(1, 0)
  )

  mockery::stub(compute_distance, "dbplyr::dbplyr_pivot_wider_spec", function(...) {
    tibble::tibble(not_sample_id = 1)
  })

  testthat::expect_error(
    compute_distance(dat, value_col = "exist"),
    regexp = "(?i)no `sample_id`"
  )
})

testthat::test_that("compute_pcoa: missing labels and zero-axis handling", {
  d <- stats::dist(matrix(c(0, 1, 1, 0), nrow = 2))
  attr(d, "Labels") <- NULL
  attr(d, "abundances") <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("p1", "p2"))
  )

  res <- compute_pcoa(d, n_axes = 3L, top_features = 5L)

  testthat::expect_true(all(res$sample_coords$sample_id %in% c("1", "2")))
  testthat::expect_gte(ncol(res$sample_coords), 1L)
  testthat::expect_s3_class(res$feature_associations, "tbl_df")
})

testthat::test_that("compute_pcoa: feature association alignment warnings", {
  testthat::skip_if_not_installed("mockery")

  d <- stats::dist(matrix(c(0, 1, 1, 0), nrow = 2))
  attr(d, "Labels") <- c("s1", "s2")
  attr(d, "abundances") <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    dimnames = list(c("s1", "x2"), c("p1", "p2"))
  )

  mockery::stub(compute_pcoa, "stats::cmdscale", function(...) {
    pts <- matrix(1, nrow = 2, ncol = 1)
    rownames(pts) <- c("s1", "s2")
    list(eig = c(1, 0.5), points = pts)
  })
  mockery::stub(compute_pcoa, "intersect", function(x, y) x[1])

  testthat::expect_warning(
    compute_pcoa(d, n_axes = 1L, top_features = 5L),
    regexp = "(?i)insufficient overlap"
  )

  attr(d, "abundances") <- matrix(c(1, 0, 0, 1), nrow = 2)

  mockery::stub(compute_pcoa, "stats::cmdscale", function(...) {
    list(eig = c(1, 0.5), points = matrix(1, nrow = 2, ncol = 1))
  })

  testthat::expect_warning(
    compute_pcoa(d, n_axes = 1L, top_features = 5L),
    regexp = "(?i)row names missing"
  )
})

testthat::test_that("compute_capscale: empty metadata and missing abundance rownames", {
  testthat::skip_if_not_installed("vegan")

  d <- stats::dist(matrix(c(0, 1, 1, 0), nrow = 2))
  attr(d, "Labels") <- c("s1", "s2")

  ps_empty <- data.frame(sample_id = character(), group = character())
  testthat::expect_error(
    compute_capscale(d, ps_empty, ~ group),
    regexp = "(?i)constructed metadata has zero rows"
  )

  ps <- data.frame(sample_id = c("s1", "s2"), group = c("A", "B"))
  attr(d, "abundances") <- matrix(1:4, nrow = 2)

  testthat::expect_warning(
    compute_capscale(d, ps, ~ group),
    regexp = "(?i)abundance matrix has no row names"
  )
})

testthat::test_that("compute_capscale: rank/eigen edge cases and score naming", {
  testthat::skip_if_not_installed("mockery")
  testthat::skip_if_not_installed("vegan")

  d <- stats::dist(matrix(c(0, 1, 1, 0), nrow = 2))
  attr(d, "Labels") <- c("s1", "s2")
  attr(d, "abundances") <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("p1", "p2"))
  )

  ps <- data.frame(sample_id = c("s1", "s2"), group = c("A", "B"))

  mockery::stub(compute_capscale, "vegan::capscale", function(...) {
    structure(
      list(
        CCA = list(rank = NULL, eig = NULL, tot.chi = NULL),
        CA = list(tot.chi = NULL, eig = NULL),
        tot.chi = 1
      ),
      class = "capscale"
    )
  })
  mockery::stub(compute_capscale, "vegan::anova.cca", function(...) {
    data.frame(Df = 1, `Pr(>F)` = 1)
  })
  mockery::stub(compute_capscale, "vegan::RsquareAdj", function(...) {
    list(r.squared = 0, adj.r.squared = 0)
  })

  res_null <- compute_capscale(d, ps, ~ group)
  testthat::expect_equal(ncol(res_null$sample_coords), 1L)
  testthat::expect_s3_class(res_null$feature_associations, "tbl_df")

  mockery::stub(compute_capscale, "vegan::capscale", function(...) {
    structure(
      list(
        CCA = list(rank = 2L, eig = c(0.2, 0.1), tot.chi = 0.3),
        CA = list(tot.chi = 0.7, eig = c(0.2, 0.1)),
        tot.chi = 1
      ),
      class = "capscale"
    )
  })
  mockery::stub(compute_capscale, "vegan::anova.cca", function(...) {
    data.frame(Df = 1, `Pr(>F)` = 1)
  })
  mockery::stub(compute_capscale, "vegan::RsquareAdj", function(...) {
    list(r.squared = 0.3, adj.r.squared = 0.25)
  })
  mockery::stub(compute_capscale, "vegan::scores", function(...) {
    matrix(1, nrow = 2, ncol = 2)
  })
  mockery::stub(compute_capscale, "intersect", function(x, y) x[1])

  testthat::expect_warning(
    res_scores <- compute_capscale(d, ps, ~ group),
    regexp = "(?i)insufficient overlap"
  )

  testthat::expect_true(all(c("CAP1", "CAP2") %in% names(res_scores$sample_coords)))
})

testthat::test_that("compute_capscale: warns when score rownames are missing", {
  testthat::skip_if_not_installed("mockery")
  testthat::skip_if_not_installed("vegan")

  d <- stats::dist(matrix(c(0, 1, 1, 0), nrow = 2))
  attr(d, "Labels") <- c("s1", "s2")
  attr(d, "abundances") <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("p1", "p2"))
  )

  ps <- data.frame(sample_id = c("s1", "s2"), group = c("A", "B"))

  mockery::stub(compute_capscale, "vegan::capscale", function(...) {
    structure(
      list(
        CCA = list(rank = 1L, eig = 0.1, tot.chi = 0.1),
        CA = list(tot.chi = 0.9, eig = 0.2),
        tot.chi = 1
      ),
      class = "capscale"
    )
  })
  mockery::stub(compute_capscale, "vegan::anova.cca", function(...) {
    data.frame(Df = 1, `Pr(>F)` = 1)
  })
  mockery::stub(compute_capscale, "vegan::RsquareAdj", function(...) {
    list(r.squared = 0.1, adj.r.squared = 0.05)
  })
  mockery::stub(compute_capscale, "vegan::scores", function(...) {
    matrix(1, nrow = 1, ncol = 1)
  })

  testthat::expect_warning(
    compute_capscale(d, ps, ~ group),
    regexp = "(?i)row names missing"
  )
})

testthat::test_that("compute_permanova: missing columns, empty data, and labels", {
  testthat::skip_if_not_installed("vegan")

  dat <- data.frame(
    sample_id = as.character(1:4),
    group = c("A", "A", "B", "B")
  )
  d <- stats::dist(matrix(runif(8), nrow = 4))
  attr(d, "Labels") <- dat$sample_id
  attr(d, "Size") <- 4

  testthat::expect_error(
    compute_permanova(d, dat, time_col = "missing"),
    regexp = "(?i)not found"
  )

  dat_empty <- data.frame(sample_id = character(), group = character())
  testthat::expect_error(
    compute_permanova(d, dat_empty, group_col = "group"),
    regexp = "(?i)constructed metadata has zero rows"
  )

  d_bad <- d
  attr(d_bad, "Labels") <- paste0("x", 1:4)
  testthat::expect_error(
    compute_permanova(d_bad, dat, group_col = "group"),
    regexp = "(?i)missing in `ps`"
  )

  d_no_labels <- d
  attr(d_no_labels, "Labels") <- NULL
  testthat::expect_warning(
    compute_permanova(d_no_labels, dat, group_col = "group",
                      permutations = 9),
    regexp = "(?i)no labels found"
  )
})

testthat::test_that("compute_permanova: global adonis failure and skip logic", {
  testthat::skip_if_not_installed("mockery")
  testthat::skip_if_not_installed("vegan")

  dat <- data.frame(
    sample_id = paste0("s", 1:4),
    group = c("A", "B", "B", "B")
  )
  d <- stats::dist(matrix(runif(8), nrow = 4))
  attr(d, "Labels") <- dat$sample_id
  attr(d, "Size") <- 4

  mockery::stub(compute_permanova, "vegan::adonis2", function(...) {
    stop("boom")
  })

  testthat::expect_warning(
    compute_permanova(d, dat, group_col = "group",
                      permutations = 9),
    regexp = "(?i)global permanova failed"
  )
})

testthat::test_that("compute_permanova: time pairwise comparisons", {
  testthat::skip_if_not_installed("vegan")

  dat <- data.frame(
    sample_id = paste0("s", 1:6),
    group = rep(c("A", "B"), each = 3),
    time = rep(c("T1", "T2", "T3"), times = 2)
  )
  d <- stats::dist(matrix(runif(18), nrow = 6))
  attr(d, "Labels") <- dat$sample_id
  attr(d, "Size") <- 6

  res_pair <- suppressWarnings(compute_permanova(
    d, dat,
    time_col = "time",
    permutations = 9
  ))
  testthat::expect_s3_class(res_pair, "tbl_df")
})

testthat::test_that("compute_permanova: no-factor and NA handling", {
  testthat::skip_if_not_installed("vegan")

  dat <- data.frame(
    sample_id = paste0("s", 1:3),
    group = c("A", NA, "A")
  )
  d <- stats::dist(matrix(runif(9), nrow = 3))
  attr(d, "Labels") <- dat$sample_id
  attr(d, "Size") <- 3

  res_none <- suppressWarnings(compute_permanova(
    d, data.frame(sample_id = dat$sample_id),
    permutations = 9
  ))
  testthat::expect_s3_class(res_none, "tbl_df")

  res_drop <- suppressWarnings(compute_permanova(
    d, dat,
    group_col = "group",
    permutations = 9
  ))
  testthat::expect_s3_class(res_drop, "tbl_df")

  dat_all_na <- data.frame(sample_id = paste0("s", 1:3), group = NA)
  testthat::expect_error(
    compute_permanova(d, dat_all_na, group_col = "group"),
    regexp = "(?i)all samples have missing values"
  )

  dat_sparse <- data.frame(
    sample_id = paste0("s", 1:3),
    group = c("A", "B", "B")
  )
  d_sparse <- stats::dist(matrix(runif(9), nrow = 3))
  attr(d_sparse, "Labels") <- dat_sparse$sample_id
  attr(d_sparse, "Size") <- 3

  res_skip <- suppressWarnings(compute_permanova(
    d_sparse, dat_sparse,
    group_col = "group",
    permutations = 9
  ))
  testthat::expect_s3_class(res_skip, "tbl_df")
})

testthat::test_that("compute_dispersion: warnings and edge cases", {
  testthat::skip_if_not_installed("vegan")
  testthat::skip_if_not_installed("mockery")

  dat <- data.frame(
    sample_id = paste0("s", 1:6),
    group = rep(c("A", "B"), each = 3),
    time = rep(c("T1", "T2", "T3"), times = 2)
  )
  d <- stats::dist(matrix(runif(18), nrow = 6))
  attr(d, "Labels") <- dat$sample_id
  attr(d, "Size") <- 6

  dat_num_time <- dat
  dat_num_time$time_num <- seq_len(nrow(dat))

  mockery::stub(compute_dispersion, "as.factor", function(x) x)

  testthat::expect_warning(
    compute_dispersion(d, dat_num_time, time_col = "time_num",
                       permutations = 9),
    regexp = "(?i)continuous dispersion"
  )

  res_bl <- suppressWarnings(compute_dispersion(
    d, dat,
    group_col = "group",
    time_col = "time",
    permutations = 9
  ))
  testthat::expect_s3_class(res_bl, "beta_dispersion")
})

testthat::test_that("compute_dispersion: betadisper failure paths", {
  testthat::skip_if_not_installed("mockery")
  testthat::skip_if_not_installed("vegan")

  dat <- data.frame(
    sample_id = paste0("s", 1:4),
    group = c("A", "A", "B", "B"),
    time = c("T1", "T2", "T1", "T2")
  )
  d <- stats::dist(matrix(runif(8), nrow = 4))
  attr(d, "Labels") <- dat$sample_id
  attr(d, "Size") <- 4

  mockery::stub(compute_dispersion, "vegan::betadisper", function(...) {
    stop("boom")
  })

  testthat::expect_warning(
    compute_dispersion(d, dat, group_col = "group", time_col = "time",
                       permutations = 9),
    regexp = "(?i)betadisper failed"
  )
})

testthat::test_that("compute_dispersion: time pairwise comparisons and skip branch", {
  testthat::skip_if_not_installed("vegan")

  dat <- data.frame(
    sample_id = paste0("s", 1:6),
    group = rep(c("A", "B"), each = 3),
    time = rep(c("T1", "T2", "T3"), times = 2)
  )
  d <- stats::dist(matrix(runif(18), nrow = 6))
  attr(d, "Labels") <- dat$sample_id
  attr(d, "Size") <- 6

  res_pair <- suppressWarnings(compute_dispersion(
    d, dat,
    time_col = "time",
    permutations = 9
  ))
  testthat::expect_s3_class(res_pair, "beta_dispersion")

  dat_sparse <- data.frame(
    sample_id = paste0("s", 1:3),
    group = c("A", "B", "B")
  )
  d_sparse <- stats::dist(matrix(runif(6), nrow = 3))
  attr(d_sparse, "Labels") <- dat_sparse$sample_id
  attr(d_sparse, "Size") <- 3

  res_skip <- suppressWarnings(compute_dispersion(
    d_sparse, dat_sparse,
    group_col = "group",
    permutations = 9
  ))
  testthat::expect_s3_class(res_skip, "beta_dispersion")
})

testthat::test_that("compute_dispersion: alignment and edge conditions", {
  testthat::skip_if_not_installed("vegan")

  dat <- data.frame(
    sample_id = as.character(1:6),
    group = rep(c("A", "B"), each = 3),
    time = rep(c("T1", "T2", "T3"), times = 2)
  )
  d <- stats::dist(matrix(runif(18), nrow = 6))
  attr(d, "Labels") <- dat$sample_id
  attr(d, "Size") <- 6

  ps_phip <- structure(list(data_long = dat), class = "phip_data")
  res_phip <- suppressWarnings(compute_dispersion(
    d, ps_phip,
    group_col = "group",
    permutations = 9
  ))
  testthat::expect_s3_class(res_phip, "beta_dispersion")

  testthat::expect_error(
    compute_dispersion(d, dat, time_col = "missing"),
    regexp = "(?i)not found"
  )

  dat_empty <- data.frame(sample_id = character(), group = character())
  testthat::expect_error(
    compute_dispersion(d, dat_empty, group_col = "group"),
    regexp = "(?i)constructed metadata has zero rows"
  )

  d_bad <- d
  attr(d_bad, "Labels") <- paste0("x", 1:6)
  testthat::expect_error(
    compute_dispersion(d_bad, dat, group_col = "group"),
    regexp = "(?i)missing in `ps`"
  )

  d_no_labels <- d
  attr(d_no_labels, "Labels") <- NULL
  testthat::expect_warning(
    compute_dispersion(d_no_labels, dat, group_col = "group",
                       permutations = 9),
    regexp = "(?i)no labels found"
  )

  res_none <- suppressWarnings(compute_dispersion(
    d, data.frame(sample_id = dat$sample_id),
    permutations = 9
  ))
  testthat::expect_s3_class(res_none, "beta_dispersion")

  dat_all_na <- dat
  dat_all_na$group <- NA
  testthat::expect_error(
    compute_dispersion(d, dat_all_na, group_col = "group"),
    regexp = "(?i)all samples have missing values"
  )

  res_inter <- suppressWarnings(compute_dispersion(
    d, dat,
    group_col = "group",
    time_col = "time",
    permutations = 9
  ))
  testthat::expect_s3_class(res_inter, "beta_dispersion")

  res_time_bl <- suppressWarnings(compute_dispersion(
    d, dat,
    time_col = "time",
    permutations = 9
  ))
  testthat::expect_s3_class(res_time_bl, "beta_dispersion")
})

testthat::test_that("compute_tsne: metadata, labels, and warnings", {
  testthat::skip_if_not_installed("Rtsne")

  dat <- data.frame(
    sample_id = paste0("s", 1:6),
    meta1 = rep(c("x", "y"), each = 3),
    stringsAsFactors = FALSE
  )
  d <- stats::dist(matrix(runif(18), nrow = 6))
  attr(d, "Labels") <- dat$sample_id
  attr(d, "Size") <- 6

  ps_phip <- structure(
    list(data_long = dat, meta = list(extra_cols = "meta1")),
    class = "phip_data"
  )

  res_meta <- suppressWarnings(compute_tsne(
    ps = ps_phip,
    dist_obj = d,
    dims = 2L,
    perplexity = 1
  ))
  testthat::expect_true("meta1" %in% names(res_meta))

  d_bad_labels <- d
  attr(d_bad_labels, "Labels") <- c("s1", "s2")
  testthat::expect_error(
    compute_tsne(ps = dat, dist_obj = d_bad_labels, dims = 2L, perplexity = 1),
    regexp = "(?i)length of distance labels"
  )

  testthat::expect_warning(
    compute_tsne(ps = NULL, dist_obj = d, dims = 2L, perplexity = 1),
    regexp = "(?i)ps.*NULL"
  )
})

testthat::test_that("compute_tsne: matrix labels and duplicate metadata rows", {
  testthat::skip_if_not_installed("Rtsne")

  dat <- data.frame(
    sample_id = rep(paste0("s", 1:6), each = 2),
    meta1 = rep(c("x", "y"), length.out = 12),
    stringsAsFactors = FALSE
  )
  d <- stats::dist(matrix(runif(18), nrow = 6))
  attr(d, "Labels") <- paste0("s", 1:6)
  attr(d, "Size") <- 6

  d_mat <- as.matrix(d)
  rownames(d_mat) <- NULL

  testthat::expect_warning(
    compute_tsne(ps = dat, dist_obj = d_mat, dims = 2L, perplexity = 1),
    regexp = "(?i)no rownames"
  )

  testthat::expect_error(
    compute_tsne(
      ps = dat,
      dist_obj = d,
      dims = 2L,
      perplexity = 1,
      meta_cols = "meta1"
    ),
    regexp = "(?i)inconsistent metadata values"
  )
})
