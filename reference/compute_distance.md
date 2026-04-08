# Compute Pairwise Sample Distances

This function builds a sample-by-feature abundance matrix from a
`phip_data` object (using `ps$data_long`), optionally normalizes the
matrix, and then computes pairwise distances between samples.

The normalized abundance matrix used for distance calculation is
attached to the returned `dist` object as attribute `"abundances"`.

Note: this function pivots to a wide matrix in the database (via dbplyr)
and then collects the result into memory. This can be large for big
cohorts and/or large peptide sets.

## Usage

``` r
compute_distance(
  ps,
  value_col = NULL,
  method_normalization = c("auto", "relative", "hellinger", "log", "none"),
  distance = "bray",
  n_threads = 1L
)
```

## Arguments

- ps:

  input data. either:

  - a `phip_data` object, in which case `ps$data_long` is used, or

  - a `data_long`-like table (a data.frame or dplyr `tbl`) containing at
    least `sample_id`, `peptide_id`, and the column given by
    `value_col`.

- value_col:

  character scalar. Name of the abundance column in `ps$data_long`. If
  `NULL`, the function tries (in order) `exist`, `counts_hits`,
  `counts_control`, `fold_change`.

- method_normalization:

  character scalar. Normalization applied to the abundance matrix before
  distance computation. One of:

  - `"auto"`: uses `"none"` for binary (0/1) data, otherwise uses
    `"relative"`.

  - `"relative"`: divide each row by its row sum.

  - `"hellinger"`: `sqrt(relative)`.

  - `"log"`: `log1p` transform.

  - `"none"`: no normalization.

- distance:

  character scalar. Distance method name. The string is lowercased
  internally.

  Supported methods depend on which packages are installed:

  - fast path (if package 'parallelDist' is installed):

    - "bray" (bray-curtis). Computed via threaded l1 distances and
      normalization (equivalent to bray-curtis on the normalized
      matrix).

    - "euclidean"

    - "minkowski"

    - "manhattan"

    - "canberra"

    - "binary"

    - "maximum" (maximum/supremum/chebyshev distance). Note:
      'parallelDist' documents this as method "maximum"; passing
      "chebyshev" may fail unless you map it to "maximum" before calling
      parDist().

    - "cosine"

  - fallback path (requires package 'vegan'). Any method supported by
    vegan::vegdist(), partial match allowed: "manhattan", "euclidean",
    "canberra", "clark", "bray", "kulczynski", "jaccard", "gower",
    "altGower", "morisita", "horn", "mountford", "raup", "binomial",
    "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger",
    "aitchison", "robust.aitchison".

  If 'parallelDist' is installed but the requested method is not in the
  fast list above, the function falls back to vegan::vegdist().

- n_threads:

  integer scalar. Number of cpu threads passed to
  `parallelDist::parDist(threads = ...)`.

## Value

a `dist` object of pairwise sample distances. The attribute
`"abundances"` contains the normalized abundance matrix used for the
calculation (rows are samples, columns are features).

## Examples

``` r
# build an example <phip_data> object from the package example dataset
ps <- load_example_data("small_mixture")

# compute distances (needs either 'parallelDist' or 'vegan')
val_col <- "fold_change"

d <- compute_distance(
  ps,
  value_col = val_col,
  method_normalization = "hellinger",
  distance = "bray",
  n_threads = 2L
)
#> [18:59:41] INFO  building abundance matrix from `ps` using `fold_change`.
#> [18:59:41] INFO  building pivot spec (sample_id x peptide_id).
#> [18:59:41] INFO  Collecting long table (sample_id, peptide_id, value).
#>                  -> compute_distance
#> [18:59:41] INFO  Pivoting to wide abundance matrix in R.
#>                  -> compute_distance
#> [18:59:41] INFO  abundance matrix has 43 samples and 5 features after
#>                  preprocessing.
#> [18:59:41] INFO  computing distance: bray
#> [18:59:41] INFO  distance matrix computation complete.

a <- attr(d, "abundances")
a[1:min(5, nrow(a)), 1:min(5, ncol(a)), drop = FALSE]
#>             16196     16627     18003     24799      5243
#> A_T1_1  0.3191635 0.0000000 0.7566979 0.5705637 0.0000000
#> B_T1_1  0.0000000 0.5056957 0.0000000 0.0000000 0.8627119
#> A_T1_10 0.5965746 0.0000000 0.6709044 0.4404385 0.0000000
#> B_T1_10 0.0000000 0.7612047 0.0000000 0.0000000 0.6485117
#> A_T1_11 0.7955740 0.0000000 0.1061303 0.5964884 0.0000000
```
