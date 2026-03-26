# Scree Plot for PCoA Eigenvalues

Creates a scree plot from a `"beta_pcoa"` object (output of
[`compute_pcoa()`](https://polymerase3.github.io/phiper/reference/compute_pcoa.md)),
showing the percentage of variance explained for the first `n_axes`
axes.

## Usage

``` r
plot_scree(pcoa_res, n_axes = NULL, type = c("bar", "line"))
```

## Arguments

- pcoa_res:

  A `"beta_pcoa"` object as returned by
  [`compute_pcoa()`](https://polymerase3.github.io/phiper/reference/compute_pcoa.md).
  Must contain an `eigenvalues` component (numeric vector of
  eigenvalues).

- n_axes:

  Integer giving the number of axes to display. If `NULL` (default), the
  function uses `min(10, length(eigenvalues))`.

- type:

  Type of scree plot to draw: `"bar"` (default) for a bar plot of
  variance explained, or `"line"` for a line/point plot.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object representing the scree plot.

## Details

Percentages are computed from the positive part of the eigenvalues, i.e.
`pmax(eigenvalues, 0)`, to handle possible small negative eigenvalues
from non-perfectly Euclidean distance matrices.

The function uses
[`theme_phip()`](https://polymerase3.github.io/phiper/reference/theme_phip.md)
if it exists in the search path; otherwise it falls back to
[`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

## Examples

``` r
if (FALSE) { # \dontrun{
  # pcoa_res <- compute_pcoa(dist_bc)
  plot_scree(pcoa_res)

  # More axes as line plot:
  plot_scree(pcoa_res, n_axes = 15, type = "line")
} # }
```
