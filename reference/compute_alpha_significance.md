# Compute statistical significance of alpha diversity between groups

Runs global (Kruskal-Wallis or one-way ANOVA) and pairwise (Wilcoxon or
Tukey HSD) hypothesis tests on each `(rank, metric)` combination present
in the output of
[`compute_alpha()`](https://polymerase3.github.io/phiper/reference/compute_alpha.md).
Pairwise p-values are adjusted with
[`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) and Cohen's d
effect sizes are appended.

## Usage

``` r
compute_alpha_significance(
  x,
  group_col = NULL,
  metric = NULL,
  global_test = c("kruskal", "anova"),
  pairwise_test = c("wilcoxon", "tukey"),
  p_adjust_method = "BH"
)
```

## Arguments

- x:

  A `"phip_alpha_diversity"` list (output of
  [`compute_alpha()`](https://polymerase3.github.io/phiper/reference/compute_alpha.md))
  or a single alpha-diversity data frame.

- group_col:

  Character scalar; name of the grouping column. Inferred from
  `attr(x, "group_cols")` when `NULL` (default).

- metric:

  Character vector; subset of metrics to test. `NULL` (default) uses all
  numeric metric columns present in the data.

- global_test:

  One of `"kruskal"` (Kruskal-Wallis, default) or `"anova"` (one-way
  ANOVA).

- pairwise_test:

  One of `"wilcoxon"` (Wilcoxon rank-sum, default) or `"tukey"` (Tukey
  HSD from [`aov()`](https://rdrr.io/r/stats/aov.html)).

- p_adjust_method:

  Method passed to
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default `"BH"`.

## Value

A named list of class `"phip_alpha_significance"` with two tibbles:

- `$global`:

  One row per `(rank, metric)`: `rank`, `metric`, `statistic`,
  `p_value`, `test`.

- `$pairwise`:

  One row per `(rank, metric, pair)`: `rank`, `metric`, `group1`,
  `group2`, `p_raw`, `p_adj`, `cohens_d`, `stars`, `test`.

Attributes: `group_col`, `global_test`, `pairwise_test`,
`p_adjust_method`, `metrics`, `ranks`.

## Examples

``` r
if (FALSE) { # \dontrun{
alpha <- compute_alpha(phip_obj, group_cols = "group")
sig   <- compute_alpha_significance(alpha)
sig$global
sig$pairwise
} # }
```
