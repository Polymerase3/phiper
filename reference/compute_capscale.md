# Constrained Ordination (db-rda / cap) on Distance Matrix

Performs distance-based redundancy analysis (constrained pcoa, a.k.a.
cap) on a distance matrix using vegan::`capscale`, with optional
negative eigenvalue correction. Returns constrained sample scores,
eigenvalues, variance partitioning, and feature associations.

## Usage

``` r
compute_capscale(
  dist_obj,
  ps,
  formula,
  neg_correction = c("none", "lingoes", "cailliez"),
  top_features = 30L,
  permutations = 999L,
  feature_assoc = c("weighted_average", "correlation", "regression", "none")
)
```

## Arguments

- dist_obj:

  A `dist` object returned by
  [`compute_distance()`](https://polymerase3.github.io/phiper/reference/compute_distance.md).
  The normalized abundance matrix used to compute the distances is
  expected to be attached as an attribute `"abundances"` (samples in
  rows, features in columns).

- ps:

  A `phip_data` object or a table providing sample-level metadata. This
  table must contain `sample_id` and all variables referenced on the
  right-hand side of `formula`. Variable detection uses
  `all.names(terms(formula), functions = FALSE)`, so transformed terms
  like `log(age)` are supported as long as `age` exists.

- formula:

  An R formula specifying the constraints (independent variables) for
  the ordination, e.g. `~ sex + age`. Do not include a response on the
  left-hand side; the distance matrix is provided via `dist_obj`.

- neg_correction:

  One of `"none"`, `"lingoes"`, `"cailliez"`. Method for negative
  eigenvalue correction. Default is `"none"`. Passed to the `add`
  argument of
  [`vegan::capscale()`](https://vegandevs.github.io/vegan/reference/dbrda.html).

- top_features:

  Integer scalar. Number of top features to return in associations
  (selected per constrained axis by absolute association, then unioned).
  Default is 30.

- permutations:

  Integer scalar. Number of permutations for per-term permutation tests
  via `vegan::anova.cca(by = "term")`. Default is 999.

- feature_assoc:

  character scalar. Type of feature-axis association to return.
  `"weighted_average"` returns weighted-average feature scores (centroid
  of sample scores weighted by feature abundance). `"correlation"`
  returns feature-axis correlations. `"regression"` returns regression
  slopes for axis scores on feature abundance. `"none"` skips feature
  associations.

## Value

A list of class `"beta_capscale"` with elements:

- sample_coords:

  Tibble of sample scores on constrained axes (`CAP1`, `CAP2`, ...).
  Contains `sample_id` and coordinates.

- eigenvalues:

  Numeric vector of eigenvalues of the constrained axes.

- variance_partition:

  Tibble with total inertia and inertia partitioned into constrained and
  unconstrained components, with their proportion of total.

- feature_associations:

  Tibble of top feature-axis associations for constrained axes (possibly
  empty if the `"abundances"` attribute is missing or cannot be
  aligned). To limit runtime/memory, associations are computed for at
  most 10 constrained axes.

- r2:

  Numeric scalar. Unadjusted R-squared from
  [`vegan::RsquareAdj()`](https://vegandevs.github.io/vegan/reference/RsquareAdj.html).

- r2_adj:

  Numeric scalar. Adjusted R-squared from
  [`vegan::RsquareAdj()`](https://vegandevs.github.io/vegan/reference/RsquareAdj.html).

- perm_terms:

  Tibble of per-term permutation tests from
  `vegan::anova.cca(by = "term")`.

- cap_model:

  The full
  [`vegan::capscale`](https://vegandevs.github.io/vegan/reference/dbrda.html)
  model object.

## Examples

``` r
# \donttest{
ps <- load_example_data("small_mixture")
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpcfwKE4/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> [11:56:45] INFO  Constructing <phip_data> object
#>                  -> create_data()
#> [11:56:45] INFO  Fetching peptide metadata library via get_peptide_library()
#> [11:56:45] INFO  Retrieving peptide metadata into DuckDB cache
#>                  -> get_peptide_library(force_refresh = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpcfwKE4/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> [11:56:45] INFO  Opened DuckDB connection
#>                    - cache dir:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/phip_cache.duckdb
#>                    - table: peptide_meta
#> [11:56:45] OK    Using cached peptide_meta (fast path)
#> [11:56:45] OK    Retrieving peptide metadata into DuckDB cache - done
#>                  -> elapsed: 0.034s
#> [11:56:45] OK    Peptide metadata acquired
#> [11:56:45] INFO  Validating <phip_data>
#>                  -> validate_phip_data()
#> [11:56:45] INFO  Checking structural requirements (shape & mandatory columns)
#> [11:56:45] INFO  Checking outcome family availability (exist / fold_change /
#>                  raw_counts)
#> [11:56:45] INFO  Checking collisions with reserved names
#>                    - subject_id, sample_id, timepoint, peptide_id, exist,
#>                      fold_change, counts_input, counts_hit
#> [11:56:45] INFO  Ensuring all columns are atomic (no list-cols)
#> [11:56:45] INFO  Checking key uniqueness
#> [11:56:45] INFO  Validating value ranges & types for outcomes
#> [11:56:45] INFO  Assessing sparsity (NA/zero prevalence vs threshold)
#>                    - warn threshold: 50%
#> [11:56:45] INFO  Checking peptide_id coverage against peptide_library
#> [11:56:45] INFO  Checking full grid completeness (peptide * sample)
#> [11:56:45] INFO  Counts table is not a full peptide * sample grid
#>                    - observed rows: 78200
#>                    - expected rows: 156000
#> Warning: [11:56:45] WARN  Grid remains incomplete (auto_expand = FALSE).
#>                  -> grid completeness
#>                    - observed rows: 78200
#>                    - expected rows: 156000.
#> [11:56:45] OK    Validating <phip_data> - done
#>                  -> elapsed: 0.301s
#> [11:56:45] OK    Constructing <phip_data> object - done
#>                  -> elapsed: 0.336s

# compute distance matrix
val_col <- "fold_change"

dist_bc <- compute_distance(
  ps,
  value_col = val_col,
  method_normalization = "hellinger",
  distance = "bray",
  n_threads = 2L
)
#> [11:56:45] INFO  building abundance matrix from `ps` using `fold_change`.
#> [11:56:45] INFO  building pivot spec (sample_id x peptide_id).
#> [11:56:45] INFO  Collecting long table (sample_id, peptide_id, value).
#>                  -> compute_distance
#> [11:56:45] INFO  Pivoting to wide abundance matrix in R.
#>                  -> compute_distance
#> [11:56:45] INFO  abundance matrix has 43 samples and 5 features after
#>                  preprocessing.
#> [11:56:45] INFO  computing distance: bray
#> [11:56:46] INFO  distance matrix computation complete.

# pick a simple constraint that exists in the example data (fallback order)
dat <- ps
cand <- c("group", "big_group", "type_person", "sex", "age")
rhs_var <- cand[cand %in% dplyr::tbl_vars(dat)][1]

cap_res <- compute_capscale(
  dist_bc,
  ps = ps,
  formula = stats::as.formula(paste0("~ ", rhs_var)),
  neg_correction = "none",
  top_features = 30L
)
#> [11:56:46] INFO  building metadata from `ps$data_long`.
#> [11:56:46] INFO  fitting constrained ordination (cap/db-rda)
#>                    - formula: ~group
#> [11:56:46] INFO  extracting constrained sample scores.
#> [11:56:46] INFO  computing variance partitioning and permutation tests.
#> [11:56:46] INFO  computing feature associations: weighted_average.
#> [11:56:46] INFO  cap analysis complete.

cap_res$variance_partition
#> # A tibble: 3 × 3
#>   component     inertia proportion
#>   <chr>           <dbl>      <dbl>
#> 1 Total           11.3       1    
#> 2 Constrained     10.1       0.899
#> 3 Unconstrained    1.27      0.113
cap_res$sample_coords
#> # A tibble: 43 × 2
#>    sample_id   CAP1
#>    <chr>      <dbl>
#>  1 A_T1_1    -0.766
#>  2 B_T1_1     0.669
#>  3 A_T1_10   -0.779
#>  4 B_T1_10    0.678
#>  5 A_T1_11   -0.749
#>  6 B_T1_11    0.666
#>  7 A_T1_12   -0.786
#>  8 B_T1_12    0.664
#>  9 A_T1_13   -0.781
#> 10 B_T1_13    0.634
#> # ℹ 33 more rows
cap_res$feature_associations
#> # A tibble: 5 × 2
#>   feature          CAP1
#>   <chr>           <dbl>
#> 1 agilent_175212 -0.764
#> 2 agilent_196916 -0.763
#> 3 agilent_21155  -0.765
#> 4 agilent_238554  0.662
#> 5 corona2_2925    0.666
# }
```
