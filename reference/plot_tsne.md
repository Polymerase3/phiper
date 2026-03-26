# Plot t-SNE embeddings

Plot t-SNE embeddings computed by
[`compute_tsne()`](https://polymerase3.github.io/phiper/reference/compute_tsne.md).
Supports both 2D (ggplot2) and 3D (plotly) views.

## Usage

``` r
plot_tsne(
  tsne_res,
  view = c("2d", "3d"),
  colour = NULL,
  size = 1.5,
  alpha = 0.8,
  palette = NULL,
  ...
)
```

## Arguments

- tsne_res:

  A `phip_tsne` object returned by
  [`compute_tsne()`](https://polymerase3.github.io/phiper/reference/compute_tsne.md).
  The function also works with a plain tibble/data frame that contains
  at least columns `tSNE1` and `tSNE2` (and `tSNE3` for 3D).

- view:

  Character, either `"2d"` or `"3d"`. `"2d"` returns a `ggplot` object;
  `"3d"` returns a `plotly` HTML widget.

- colour:

  Optional name of a column in `tsne_res` to map to point colour. If
  `NULL`, the function uses the first metadata column stored in
  `attr(tsne_res, "meta_cols")`, if available.

- size:

  Numeric point size for scatter plots. Defaults to `1.5`.

- alpha:

  Numeric transparency for points (0-1). Defaults to `0.8`.

- palette:

  Optional vector of colour values passed to
  [`scale_color_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
  (2D) or `colors` (3D plotly).

- ...:

  Currently ignored; reserved for future extensions.

## Value

For `view = "2d"`, a `ggplot` object. For `view = "3d"`, a `plotly`
object (`htmlwidget`).

## Examples

``` r
if (FALSE) { # \dontrun{
tsne_res <- compute_tsne(ps, compute_distance(ps))

# 2D plot
p2d <- plot_tsne(tsne_res, view = "2d", colour = "type_person")

# 3D interactive plot
p3d <- plot_tsne(tsne_res, view = "3d", colour = "type_person")
} # }
```
