# Compute Feature Associations to PCoA Vectors

Calculates feature-axis associations based on given PCoA results (output
of
[`compute_pcoa()`](https://polymerase3.github.io/phiper/reference/compute_pcoa.md))

## Usage

``` r
compute_pcoa_feature_associations(
  dist_obj,
  pcoa_result,
  top_features = 30L,
  association_method = c("weighted_average", "correlation", "regression")
)
```

## Arguments

- dist_obj:

  a `dist` object (for example returned by
  [`compute_distance()`](https://polymerase3.github.io/phiper/reference/compute_distance.md)).
  The normalized abundance matrix used to compute the distances *must*
  be attached as attribute `"abundances"` (numeric matrix with samples
  in rows and features in columns).

- pcoa_result:

  a list of class `"beta_pcoa"`. The result of
  [`compute_pcoa()`](https://polymerase3.github.io/phiper/reference/compute_pcoa.md),
  which contains the resulting PCoA eigen vectors.

- top_features:

  integer scalar. Number of features to keep per axis when reporting
  associations. Features are selected by taking the union of the top
  `top_features` features (by absolute association) for each returned
  axis. Must be \> 0.

- association_method:

  character scalar. Type of feature-axis association to return.
  `"weighted_average"` returns weighted-average feature scores (centroid
  of sample scores weighted by feature abundance). `"correlation"`
  returns feature-axis correlations. `"regression"` returns regression
  slopes for axis scores on feature abundance. `"none"` skips feature
  associations.

## Value

a tibble of feature-axis associations for the returned axes.

## Details

These feature associations are post-hoc summaries of how features relate
to PCoA axes. Weighted-average scores
(`association_method = "weighted_average"`) compute
`t(X) %*% U / colSums(X)`, where `X` is the abundance matrix and `U` are
the sample coordinates. Correlation and regression associations are
computed between feature abundances and axis scores and are not "true"
PCA loadings unless distances are Euclidean and derived compatibly.

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
#> [10:26:57] INFO  building abundance matrix from `ps` using `fold_change`.
#> [10:26:57] INFO  building pivot spec (sample_id x peptide_id).
#> [10:26:57] INFO  Collecting long table (sample_id, peptide_id, value).
#>                  -> compute_distance
#> [10:26:57] INFO  Pivoting to wide abundance matrix in R.
#>                  -> compute_distance
#> [10:26:57] INFO  abundance matrix has 43 samples and 5 features after
#>                  preprocessing.
#> [10:26:57] INFO  auto normalization selected -> using relative
#> [10:26:57] INFO  computing distance: jaccard
#> [10:26:57] INFO  distance matrix computation complete.

# Compute PCoA vectors on these distances
pcoa_res <- compute_pcoa(d, neg_correction = "none", n_axes = 3L)
#> [10:26:57] INFO  performing principal coordinates analysis
#> [10:26:57] INFO  extracting sample coordinates.
#> [10:26:57] INFO  summarizing eigenvalues and variance explained.
#> [10:26:57] INFO  pcoa analysis complete.

feature_associations <- compute_pcoa_feature_associations(d, pcoa_res)
feature_associations
#> # A tibble: 5 × 4
#>   feature  PCoA1   PCoA2     PCoA3
#>   <chr>    <dbl>   <dbl>     <dbl>
#> 1 16196    0.456 -0.0108  0.0469  
#> 2 16627   -0.409 -0.149   0.000312
#> 3 18003    0.451 -0.0100 -0.200   
#> 4 24799    0.456 -0.0106  0.115   
#> 5 5243    -0.383  0.159   0.000490
# }
```
