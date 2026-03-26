# Plot Principal Coordinates Analysis (PCoA)

Creates a 2D PCoA scatterplot from a `"beta_pcoa"` object (the output of
[`compute_pcoa()`](https://polymerase3.github.io/phiper/reference/compute_pcoa.md)),
with optional grouping by a categorical variable, optional time factor,
centroids, centroid trajectories and ellipses.

## Usage

``` r
plot_pcoa(
  pcoa_res,
  axes = c(1, 2),
  group_col = NULL,
  time_col = NULL,
  show_centroids = TRUE,
  centroid_by = c("auto", "group", "time", "group_time"),
  connect_centroids = c("none", "group", "time"),
  show_ellipses = TRUE,
  ellipse_by = c("group"),
  ellipse_type = c("t", "norm", "euclid"),
  ellipse_level = 0.95,
  point_size = 1.5,
  point_alpha = 0.6,
  centroid_size = 3
)
```

## Arguments

- pcoa_res:

  A `"beta_pcoa"` object as returned by
  [`compute_pcoa()`](https://polymerase3.github.io/phiper/reference/compute_pcoa.md).
  Must contain at least a `sample_coords` component with columns
  `sample_id` and `PCoA1`, `PCoA2`, ..., and a `var_explained` component
  (one-row tibble) with columns such as `"%PCoA1"`, `"%PCoA2"`, etc.

- axes:

  Integer vector of length 2 giving the PCoA axes to plot (e.g.,
  `c(1, 2)` for PCoA1 vs PCoA2). Defaults to `c(1, 2)`.

- group_col:

  Optional name of the grouping column in `pcoa_res$sample_coords`
  (e.g., `"group"`, `"type_person"`). Used for colouring points and for
  group-level centroids / ellipses. If `NULL`, the function attempts to
  auto-detect `"group"` if present.

- time_col:

  Optional name of a **categorical** time column in
  `pcoa_res$sample_coords` (e.g., `"time"` or `"timepoint"`). Used for
  shaping points and for time/interaction centroids / ellipses. If
  `NULL`, the function attempts to auto-detect `"time"` if present.
  Continuous time is **not** supported here.

- show_centroids:

  Logical; if `TRUE` (default), centroids of the chosen grouping are
  shown.

- centroid_by:

  Granularity of centroids. One of:

  - `"auto"` (default): if both `group_col` and `time_col` are present,
    centroids are computed for each group-time combination
    (`"group_time"`); otherwise by `"group"` or `"time"` depending on
    which factor is available.

  - `"group"`: centroids per group (requires `group_col`).

  - `"time"`: centroids per time level (requires `time_col`).

  - `"group_time"`: centroids per group\*time combination (requires both
    `group_col` and `time_col`).

- connect_centroids:

  Character; if centroids are shown and `centroid_by = "group_time"`,
  this controls how they are joined by paths:

  - `"none"` (default): no connecting paths;

  - `"group"`: connect centroids along time within each group;

  - `"time"`: connect centroids across groups within each time level.

  For other `centroid_by` values, connections are ignored.

- show_ellipses:

  Logical; if `TRUE` (default), confidence ellipses (via
  [`ggplot2::stat_ellipse()`](https://ggplot2.tidyverse.org/reference/stat_ellipse.html))
  are drawn for selected factors.

- ellipse_by:

  Character vector specifying which factor(s) to use for ellipses.
  Possible values are `"group"`, `"time"`, `"group_time"`. Default is
  `"group"`. Values that require a missing factor are ignored.

- ellipse_type:

  Type of ellipse passed to
  [`ggplot2::stat_ellipse()`](https://ggplot2.tidyverse.org/reference/stat_ellipse.html),
  one of `"t"`, `"norm"`, `"euclid"`. Defaults to `"t"`.

- ellipse_level:

  Numeric confidence level for ellipses (default `0.95`).

- point_size:

  Numeric size of points. Default `1.5`.

- point_alpha:

  Numeric alpha (opacity) of points in `[0, 1]`. Default `0.6`.

- centroid_size:

  Numeric size of centroid points. Default `3`.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object representing the PCoA scatter plot.

## Details

The function expects that any grouping / time variables you wish to use
(e.g., `"type_person"`, `"age_group"`, `"timepoint"`) have already been
merged into `pcoa_res$sample_coords` (typically by joining with a
metadata table via `sample_id`).

Axis labels are automatically annotated with the percentage of variance
explained, using the `var_explained` component of `pcoa_res` if
available.

The function uses
[`theme_phip()`](https://polymerase3.github.io/phiper/reference/theme_phip.md)
if it exists in the search path; otherwise it falls back to
[`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

## Examples

``` r
if (FALSE) { # \dontrun{
  # pcoa_res <- compute_pcoa(dist_bc)
  # Suppose you have metadata with columns sample_id, type_person, timepoint:
  meta_sc <- dplyr::left_join(
    pcoa_res$sample_coords,
    meta,
    by = "sample_id"
  )
  pcoa_res$sample_coords <- meta_sc

  # Basic PCoA coloured by type_person with group-level ellipses:
  plot_pcoa(
    pcoa_res,
    axes        = c(1, 2),
    group_col   = "type_person",
    time_col    = "timepoint",
    ellipse_by  = "group",
    show_centroids = TRUE,
    centroid_by    = "group_time",
    connect_centroids = "group"
  )
} # }
```
