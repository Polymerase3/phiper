# Compute t-SNE Embeddings for Sample Distances

Performs t-distributed stochastic neighbor embedding (t-SNE) on a sample
distance matrix to create low-dimensional embeddings for visualization.
Returns t-SNE coordinates with optional sample metadata attached.

## Usage

``` r
compute_tsne(
  ps,
  dist_obj,
  dims = 3L,
  perplexity = 30,
  theta = 0.5,
  max_iter = 1000L,
  meta_cols = NULL,
  seed = NULL,
  check_duplicates = FALSE,
  ...
)
```

## Arguments

- ps:

  A `phip_data` object or a table providing sample-level metadata. This
  table must contain `sample_id` and any columns specified in
  `meta_cols`.

- dist_obj:

  A sample distance object. Either:

  - A `dist` object (e.g., from
    [`compute_distance()`](https://polymerase3.github.io/phiper/reference/compute_distance.md)).

  - A numeric, symmetric matrix with row/column names corresponding to
    `sample_id`s.

- dims:

  Integer scalar. Number of t-SNE dimensions to compute (2 or 3).
  Default is 3 to enable both 2D and 3D visualizations.

- perplexity:

  Numeric scalar. Perplexity parameter for t-SNE. Must be smaller than
  the number of samples. If too large, it is automatically reduced with
  a warning.

- theta:

  Numeric scalar. Speed/accuracy tradeoff for Barnes-Hut approximation.
  Default is 0.5.

- max_iter:

  Integer scalar. Maximum number of t-SNE iterations. Default is 1000.

- meta_cols:

  Character vector. Column names from `ps` to attach as metadata. If
  `NULL`, uses `ps$meta$extra_cols` when available.

- seed:

  Integer scalar. Random seed for reproducibility. If provided, the seed
  is temporarily set and restored after computation.

- check_duplicates:

  Logical scalar. Whether to check for duplicate points. Default is
  `FALSE` for distance input.

- ...:

  Additional arguments passed to
  [`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html).

## Value

A tibble of class `c("phip_tsne", "tbl_df", "tbl", "data.frame")` with
columns:

- `sample_id`: Sample identifier from distance labels.

- `tSNE1`, `tSNE2`: First two t-SNE dimensions.

- `tSNE3`: Third dimension if `dims >= 3`, otherwise `NA`.

- Additional metadata columns as specified in `meta_cols`.

Attributes:

- `"distance"`: The original `dist_obj`.

- `"tsne_params"`: List of t-SNE parameters and function call.

- `"meta_cols"`: Character vector of metadata columns used.

## Details

This function runs t-SNE in distance mode (`is_distance = TRUE`), using
the supplied distance matrix directly. The distance computation should
be performed separately using
[`compute_distance()`](https://polymerase3.github.io/phiper/reference/compute_distance.md).

t-SNE is a visualization method: it preserves local neighborhoods rather
than global distances, axes are not directly interpretable, and
embeddings can change with different seeds or perplexity settings.

Samples with missing values in metadata columns are retained in the
t-SNE result but will have `NA` values for those metadata columns.

## Examples

``` r
# \donttest{
# Build example phip_data object
ps <- load_example_data("small_mixture")

# Compute distance matrix
val_col <- "fold_change"

d <- compute_distance(
  ps,
  value_col = val_col,
  method_normalization = "hellinger",
  distance = "bray",
  n_threads = 2L
)
#> [19:43:19] INFO  building abundance matrix from `ps` using `fold_change`.
#> [19:43:19] INFO  building pivot spec (sample_id x peptide_id).
#> [19:43:19] INFO  Collecting long table (sample_id, peptide_id, value).
#>                  -> compute_distance
#> [19:43:19] INFO  Pivoting to wide abundance matrix in R.
#>                  -> compute_distance
#> [19:43:19] INFO  abundance matrix has 43 samples and 5 features after
#>                  preprocessing.
#> [19:43:19] INFO  computing distance: bray
#> [19:43:19] INFO  distance matrix computation complete.

# Compute t-SNE embeddings
tsne_res <- compute_tsne(
  ps = ps,
  dist_obj = d,
  dims = 3L,
  perplexity = 15,
  meta_cols = c("subject_id", "timepoint"),
  seed = 42
)
#> Warning: [19:43:19] WARN  Perplexity (15) is high for n = 43; reducing to 14.
#> [19:43:19] INFO  Running t-SNE with dims = 3, perplexity = 14 on 43 samples
#>                  (distance input).
#> [19:43:19] INFO  Attaching metadata columns to t-SNE result: subject_id,
#>                  timepoint
#> [19:43:20] INFO  t-SNE embedding computation finished.

# View results
head(tsne_res)
#> # A tibble: 6 × 6
#>   sample_id tSNE1 tSNE2 tSNE3 subject_id timepoint
#>   <chr>     <dbl> <dbl> <dbl> <chr>      <chr>    
#> 1 A_T1_1    -6.33  3.08  27.8 1          T1       
#> 2 B_T1_1     5.25 -2.96 -27.0 1          T1       
#> 3 A_T1_10   -6.63  3.05  28.9 10         T1       
#> 4 B_T1_10    4.96 -3.00 -24.6 10         T1       
#> 5 A_T1_11   -5.01  3.78  30.2 11         T1       
#> 6 B_T1_11    4.93 -2.61 -27.1 11         T1       
# }
```
