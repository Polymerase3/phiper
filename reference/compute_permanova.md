# PERMANOVA with Global and Post-hoc Tests on Beta Diversity

Performs PERMANOVA (adonis2) on a distance matrix for overall group/time
effects, and optionally conducts post-hoc pairwise or contrast tests
(e.g., between each pair of groups, etc.). Supports stratified
permutations for repeated measures.

## Usage

``` r
compute_permanova(
  dist_obj,
  ps,
  group_col = NULL,
  time_col = NULL,
  subject_col = NULL,
  permutations = 999,
  p_adjust = "none"
)
```

## Arguments

- dist_obj:

  A `dist` object of distances between samples (e.g., output of
  [`compute_distance()`](https://polymerase3.github.io/phiper/reference/compute_distance.md)).

- ps:

  A `phip_data` object or a table providing sample-level metadata. This
  table must contain `sample_id` and the columns specified in
  `group_col`, `time_col`, and optionally `subject_col`.

- group_col:

  Name of the grouping column in `ps` (between-subject factor). Use
  `NULL` if no group factor.

- time_col:

  Name of the time factor column in `ps` (within-subject factor for
  longitudinal data). Use `NULL` if not applicable. This should be a
  *categorical* factor for this function (continuous time not
  supported).

- subject_col:

  Name of the subject identifier column in `ps` (for repeated measures).
  Default `NULL`. If provided and this column is present and `time_col`
  is provided, permutations will be stratified by subject. This is a
  simplification and does not implement a full repeated-measures
  permutation design (see Details).

- permutations:

  Number of permutations for significance testing (default 999).

- p_adjust:

  P-value adjustment method applied within each contrast scope. Use
  `"none"` for raw p-values. Passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A tibble with columns:

- scope:

  The scope of the test (e.g., `"global"`, `"group_pairwise"`,
  `"time_pairwise"`).

- contrast:

  Description of the contrast (e.g., `"<global>"` for overall test,
  `"A vs B"` for pairwise group or time comparisons).

- term:

  The term being tested. For global tests, this will be the factor name
  (or interaction). For post-hoc tests, it may be `"group"` or `"time"`
  indicating which factor is being contrasted.

- F_stat:

  F statistic of the PERMANOVA test (for global tests and some contrasts
  where applicable).

- R2:

  R-squared (variance explained) for the term (global tests).

- p_value:

  Permutation p-value for the test.

- p_adjust:

  Adjusted p-value (within scope); equals `p_value` when
  `p_adjust = "none"`.

- n_perm:

  Number of permutations used.

## Details

Global PERMANOVA uses a model with main effects of `group_col` and
`time_col` and their interaction (if both are provided and have \>1
level). Samples with missing values in the constrained variables
(group/time and, when used for stratification, subject) are dropped
*before* fitting the model, and the distance matrix is subset to the
remaining samples so that distances and metadata are always aligned.

Pairwise post-hoc tests are always performed for `group_col` and
`time_col` (when present and with \>1 level), using `adonis2` with
appropriate subsetting and, where applicable, subject stratification.
Stratification via `strata = subject` is a simplified approach and does
not replace a full repeated-measures permutation design (e.g., via
[`permute::how()`](https://rdrr.io/pkg/permute/man/how.html)),
especially for multi-level time with missing timepoints.

P-values are optionally adjusted within each contrast scope (e.g.,
`"group_pairwise"`, `"time_pairwise"`) using
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## Examples

``` r
# \donttest{
ps <- phip_load_example_data("small_mixture")

# compute distance matrix
val_col <- "fold_change"

dist_bc <- compute_distance(
  ps,
  value_col = val_col,
  distance = "jaccard",
  n_threads = 2L
)
#> [10:23:15] INFO  building abundance matrix from `ps` using `fold_change`.
#> [10:23:15] INFO  building pivot spec (sample_id x peptide_id).
#> [10:23:15] INFO  Collecting long table (sample_id, peptide_id, value).
#>                  -> compute_distance
#> [10:23:15] INFO  Pivoting to wide abundance matrix in R.
#>                  -> compute_distance
#> [10:23:15] INFO  abundance matrix has 43 samples and 5 features after
#>                  preprocessing.
#> [10:23:15] INFO  auto normalization selected -> using relative
#> [10:23:15] INFO  computing distance: jaccard
#> [10:23:15] INFO  distance matrix computation complete.

permanova_res <- compute_permanova(
  dist_bc,
  ps        = ps,
  group_col = "group",
  time_col  = "timepoint"
)
#> [10:23:15] INFO  preparing distance labels and metadata.
#> Warning: [10:23:15] WARN  column `subject_id` found in `ps`, but `subject_col` is NULL;
#>                  repeated-measures stratification is disabled.
#> [10:23:15] INFO  building metadata from `ps`.
#> [10:23:15] INFO  filtering samples with missing grouping variables.
#> [10:23:15] INFO  subsetting distance matrix to complete cases.
#> [10:23:15] INFO  preparing global permanova model.
#> [10:23:15] INFO  running global permanova
#>                    - model: d_resp ~ group
#> [10:23:16] INFO  running pairwise permanova contrasts.
#> Warning: number of items to replace is not a multiple of replacement length

permanova_res2 <- compute_permanova(
  dist_bc,
  ps        = ps,
  group_col = "group",
  time_col  = "timepoint",
  p_adjust  = "BH"
)
#> [10:23:16] INFO  preparing distance labels and metadata.
#> Warning: [10:23:16] WARN  column `subject_id` found in `ps`, but `subject_col` is NULL;
#>                  repeated-measures stratification is disabled.
#> [10:23:16] INFO  building metadata from `ps`.
#> [10:23:16] INFO  filtering samples with missing grouping variables.
#> [10:23:16] INFO  subsetting distance matrix to complete cases.
#> [10:23:16] INFO  preparing global permanova model.
#> [10:23:16] INFO  running global permanova
#>                    - model: d_resp ~ group
#> [10:23:16] INFO  running pairwise permanova contrasts.
#> Warning: number of items to replace is not a multiple of replacement length
# }
```
