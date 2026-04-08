# Principal Components Analysis (PCoA) on a Distance Matrix

Performs PCoA on a distance matrix (typically from
[`compute_distance()`](https://polymerase3.github.io/phiper/reference/compute_distance.md)),
optionally correcting for negative eigenvalues, and returns coordinates,
eigenvalues, variance explained, and feature-axis associations.

## Usage

``` r
compute_pcoa(
  dist_obj,
  neg_correction = c("none", "lingoes", "cailliez"),
  n_axes = 5L
)
```

## Arguments

- dist_obj:

  a `dist` object (for example returned by
  [`compute_distance()`](https://polymerase3.github.io/phiper/reference/compute_distance.md)).
  The normalized abundance matrix used to compute the distances is
  attached as attribute `"abundances"` (numeric matrix with samples in
  rows and features in columns). If missing, feature associations are
  skipped.

- neg_correction:

  character scalar. Method for adjusting negative eigenvalues (if any).
  One of `"none"`, `"lingoes"`, or `"cailliez"`. Default is `"none"`.

- n_axes:

  integer scalar. Number of PCoA axes to return in the sample scores.
  Must be \> 0. Internally, `k = min(n_axes, n_samples - 1)`.

## Value

a list of class `"beta_pcoa"` with elements:

- `sample_coords`: tibble with `sample_id` and columns `PCoA1`, `PCoA2`,
  ... up to `n_axes` (or fewer if `n_samples - 1` is smaller).

- `eigenvalues`: numeric vector of eigenvalues from the PCoA.

- `var_explained`: one-row tibble with percent variance explained by the
  returned axes and `%Other`. percentages are computed from the sum of
  positive eigenvalues.

- `eigen_diagnostics`: one-row tibble with eigenvalue diagnostics:
  `sum_negative`, `sum_positive`, their ratio, the minimum eigenvalue,
  and the number and fraction of negative eigenvalues.

- `correction_infos`: one-row tibble describing any negative eigenvalue
  correction applied, including `corrected`, `correction_method`,
  `correction_add`, and `correction_note`.

## Details

Negative eigenvalues indicate that the distances are not perfectly
euclidean. If `neg_correction` is `"lingoes"` or `"cailliez"`, a
correction is applied via `vegan::wcmdscale(add = ...)`.

## Examples

``` r
# \donttest{
# compute a distance matrix with an attached abundance matrix
# build an example <phip_data> object from the package example dataset
ps <- load_example_data("small_mixture")

# compute distances (needs either 'parallelDist' or 'vegan')
val_col <- "fold_change"

d <- compute_distance(
  ps,
  value_col = val_col,
  distance = "jaccard",
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
#> [18:59:41] INFO  auto normalization selected -> using relative
#> [18:59:41] INFO  computing distance: jaccard
#> [18:59:41] INFO  distance matrix computation complete.

pcoa_res <- compute_pcoa(d, neg_correction = "none", n_axes = 3L)
#> [18:59:41] INFO  performing principal coordinates analysis
#> [18:59:41] INFO  extracting sample coordinates.
#> [18:59:41] INFO  summarizing eigenvalues and variance explained.
#> [18:59:41] INFO  pcoa analysis complete.
pcoa_res$sample_coords
#> # A tibble: 43 × 4
#>    sample_id  PCoA1   PCoA2      PCoA3
#>    <chr>      <dbl>   <dbl>      <dbl>
#>  1 A_T1_1     0.457 -0.0104 -0.287    
#>  2 B_T1_1    -0.389  0.328   0.000742 
#>  3 A_T1_10    0.497 -0.0145 -0.241    
#>  4 B_T1_10   -0.444 -0.126   0.000608 
#>  5 A_T1_11    0.453 -0.0109  0.259    
#>  6 B_T1_11   -0.380  0.353   0.000700 
#>  7 A_T1_12    0.519 -0.0166 -0.109    
#>  8 B_T1_12   -0.404 -0.313   0.0000693
#>  9 A_T1_13    0.500 -0.0146 -0.148    
#> 10 B_T1_13   -0.335 -0.391  -0.000513 
#> # ℹ 33 more rows
pcoa_res$var_explained
#> # A tibble: 1 × 4
#>   `%PCoA1` `%PCoA2` `%PCoA3` `%Other`
#>      <dbl>    <dbl>    <dbl>    <dbl>
#> 1     57.2     15.8     9.45     17.6
# }
```
