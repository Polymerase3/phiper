# Test Homogeneity of Dispersion (Beta Dispersion)

Computes distances of samples to group centroids (using
[`vegan::betadisper`](https://vegandevs.github.io/vegan/reference/betadisper.html))
and tests for differences in dispersion among groups or time levels.
Pairwise post-hoc tests are always performed when `group_col` or
`time_col` has \>1 level.

## Usage

``` r
compute_dispersion(
  dist_obj,
  ps,
  group_col = NULL,
  time_col = NULL,
  subject_col = "subject_id",
  permutations = 999,
  p_adjust = "none"
)
```

## Arguments

- dist_obj:

  A `dist` object of sample distances (e.g. from
  [`compute_distance()`](https://polymerase3.github.io/phiper/reference/compute_distance.md)).

- ps:

  A `phip_data` object or a table providing sample-level metadata. This
  table must contain `sample_id` and the columns specified in
  `group_col` and/or `time_col`.

- group_col:

  Name of the group factor column in `ps` (between-subjects). Use `NULL`
  if no group factor.

- time_col:

  Name of the time factor column in `ps` (within-subjects, categorical
  only). Use `NULL` if not applicable.

- subject_col:

  Name of subject identifier column (for reference only; not used
  directly in dispersion test calculations, but kept for API
  consistency). Default `"subject_id"`.

- permutations:

  Number of permutations for significance testing in
  [`vegan::permutest`](https://vegandevs.github.io/vegan/reference/anova.cca.html).
  Default 999.

- p_adjust:

  P-value adjustment method applied within each contrast scope. Use
  `"none"` for raw p-values. Passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A list of class `"beta_dispersion"` with:

- distances:

  Tibble of per-sample distances to centroid. Columns: `sample_id`,
  `distance`, `level` (group/time level for a given scope), `scope`
  (e.g. `"group"`, `"time"`, `"group:time"`), `contrast` (e.g.
  `"<global>"`, `"A vs B"`).

- tests:

  Tibble of dispersion test results. Columns: `scope`, `contrast`,
  `term = "dispersion"`, `p_value`, `p_adjust` (equals `p_value` when
  `p_adjust = "none"`), `n_perm`.

## Examples

``` r
# \donttest{
ps <- load_example_data("small_mixture")

# compute distance matrix
val_col <- "fold_change"

dist_bc <- compute_distance(
  ps,
  value_col = val_col,
  distance = "jaccard",
  n_threads = 2L
)
#> [10:26:56] INFO  building abundance matrix from `ps` using `fold_change`.
#> [10:26:56] INFO  building pivot spec (sample_id x peptide_id).
#> [10:26:56] INFO  Collecting long table (sample_id, peptide_id, value).
#>                  -> compute_distance
#> [10:26:56] INFO  Pivoting to wide abundance matrix in R.
#>                  -> compute_distance
#> [10:26:56] INFO  abundance matrix has 43 samples and 5 features after
#>                  preprocessing.
#> [10:26:56] INFO  auto normalization selected -> using relative
#> [10:26:56] INFO  computing distance: jaccard
#> [10:26:56] INFO  distance matrix computation complete.

dispersion_res <- compute_dispersion(
  dist_bc,
  ps        = ps,
  group_col = "group",
  time_col  = "timepoint",
  p_adjust  = "BH"
)
#> [10:26:57] INFO  preparing distance labels and metadata.
#> [10:26:57] INFO  building metadata from `ps`.
#> [10:26:57] INFO  filtering samples with missing grouping variables.
#> [10:26:57] INFO  computing global dispersion tests.
#> [10:26:57] INFO  running pairwise dispersion contrasts.
dispersion_res$tests
#> # A tibble: 1 × 6
#>   scope contrast term       p_value p_adjust n_perm
#>   <chr> <chr>    <chr>        <dbl>    <dbl>  <dbl>
#> 1 group <global> dispersion   0.098    0.098    999
head(dispersion_res$distances)
#> # A tibble: 6 × 5
#>   sample_id distance level scope contrast
#>   <chr>        <dbl> <chr> <chr> <chr>   
#> 1 A_T1_1       0.397 A     group <global>
#> 2 B_T1_1       0.438 B     group <global>
#> 3 A_T1_10      0.284 A     group <global>
#> 4 B_T1_10      0.115 B     group <global>
#> 5 A_T1_11      0.408 A     group <global>
#> 6 B_T1_11      0.461 B     group <global>
# }
```
