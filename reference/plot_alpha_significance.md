# Visualise alpha diversity significance results

Summarises the pairwise significance table from
[`compute_alpha_significance()`](https://polymerase3.github.io/phiper/reference/compute_alpha_significance.md)
either as a filtered tibble (`type = "table"`) or a Cohen's d heatmap
with significance annotations (`type = "heatmap"`).

## Usage

``` r
plot_alpha_significance(
  x,
  metric = NULL,
  rank = NULL,
  type = c("table", "heatmap"),
  p_threshold = 0.05,
  ...
)
```

## Arguments

- x:

  A `"phip_alpha_significance"` object from
  [`compute_alpha_significance()`](https://polymerase3.github.io/phiper/reference/compute_alpha_significance.md).

- metric:

  Character scalar; metric to display. Defaults to the first metric
  present in `x`.

- rank:

  Character scalar; rank to display. Defaults to the first rank present
  in `x`. Ignored when the data has no `rank` column.

- type:

  One of `"table"` (default) or `"heatmap"`.

- p_threshold:

  Numeric; pairs with `p_adj > p_threshold` are shown in grey in the
  heatmap and excluded from the table. Default `0.05`.

- ...:

  Reserved; ignored.

## Value

- `type = "table"`: a
  [`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  of significant pairwise results.

- `type = "heatmap"`: a `ggplot` heatmap (Cohen's d fill, stars
  overlaid).

## Examples

``` r
if (FALSE) { # \dontrun{
sig <- compute_alpha_significance(alpha_list)
plot_alpha_significance(sig, type = "table")
plot_alpha_significance(sig, type = "heatmap")
} # }
```
