# Plot beta dispersion (distance to centroid)

Plot distances to group/time centroids from a `beta_dispersion` object
(output of
[`compute_dispersion()`](https://polymerase3.github.io/phiper/reference/compute_dispersion.md))
as violins, hollow boxplots and/or jittered points, by level on the
x-axis.

## Usage

``` r
plot_dispersion(
  x,
  scope = "group",
  contrast = "<global>",
  show_violin = TRUE,
  show_box = TRUE,
  show_points = TRUE,
  point_size = 1.5,
  point_alpha = 0.6,
  width_violin = 0.9,
  width_box = 0.5,
  jitter_width = 0.1,
  jitter_height = 0,
  add_pvalue = TRUE
)
```

## Arguments

- x:

  A `beta_dispersion` object returned by
  [`compute_dispersion()`](https://polymerase3.github.io/phiper/reference/compute_dispersion.md).

- scope:

  Character scalar indicating which dispersion scope to plot. Typically
  one of `"group"`, `"time"` or `"group:time"`. Must match values in
  `x$distances$scope`. Default `"group"`.

- contrast:

  Character scalar indicating which contrast within `scope` to plot,
  e.g. `"<global>"` or `"A vs B"`. Must match values in
  `x$distances$contrast`. Default `"<global>"`.

- show_violin:

  Logical; if `TRUE` (default) draw a violin layer.

- show_box:

  Logical; if `TRUE` (default) draw hollow boxplots (`fill = NA`) per
  level.

- show_points:

  Logical; if `TRUE` (default) draw jittered points.

- point_size:

  Numeric point size for jitter layer. Default 1.5.

- point_alpha:

  Numeric alpha for jitter layer. Default 0.6.

- width_violin:

  Width of violin layer. Default 0.9.

- width_box:

  Width of boxplot layer. Default 0.5.

- jitter_width, jitter_height:

  Jitter width/height for points. Defaults: `0.1` and `0`.

- add_pvalue:

  Logical; if `TRUE` (default) and a matching row is found in `x$tests`,
  the p-value is shown in the subtitle.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
  disp_res <- compute_dispersion(dist_bc, ps, group_col = "group_char")

  # Global group dispersion:
  plot_dispersion(disp_res, scope = "group", contrast = "<global>")

  # Pairwise contrast, with only boxplots + points:
  plot_dispersion(
    disp_res,
    scope        = "group",
    contrast     = "dementia vs control",
    show_violin  = FALSE,
    show_box     = TRUE,
    show_points  = TRUE
  )
} # }
```
