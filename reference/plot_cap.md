# Plot CAP/db-RDA Results (Constrained Ordination)

Creates a 2D scatter plot from a `"beta_capscale"` object (output of
[`compute_capscale()`](https://polymerase3.github.io/phiper/reference/compute_capscale.md)),
showing sample scores on CAP axes. Points can be coloured by group and
shaped by (categorical) time. Optional ellipses and centroids can be
overlaid for group/time/group\*time summaries.

## Usage

``` r
plot_cap(
  cap_res,
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

- cap_res:

  A `"beta_capscale"` object as returned by
  [`compute_capscale()`](https://polymerase3.github.io/phiper/reference/compute_capscale.md).
  Must contain at least a `sample_coords` tibble with columns
  `sample_id` and `CAP1`, `CAP2`, ..., and a numeric `eigenvalues`
  vector for constrained axes.

- axes:

  Integer vector of length 2 giving the CAP axes to plot (e.g.,
  `c(1, 2)` for CAP1 vs CAP2). Default is `c(1, 2)`.

- group_col:

  Optional name of a grouping column in `cap_res$sample_coords` used to
  colour points (e.g., `"group"`, `"type_person"`). If `NULL`, the
  function auto-detects `"group"` if present.

- time_col:

  Optional name of a **categorical** time column in
  `cap_res$sample_coords` used to shape points (e.g., `"time"` or
  `"timepoint"`). If `NULL`, the function auto-detects `"time"` if
  present. Continuous time is not supported in this plotting function.

- show_centroids:

  Logical; if `TRUE` (default), plots centroids for groups / time levels
  / group\*time combinations depending on `centroid_by`.

- centroid_by:

  How to define centroids. One of `"auto"`, `"group"`, `"time"`,
  `"group_time"`. `"auto"` uses: group\*time if both factors are
  available, otherwise group if present, otherwise time.

- connect_centroids:

  How to connect centroids with lines. One of `"none"` (default),
  `"group"`, `"time"`. The logic is:

  - if centroids are by `group_time` and `connect_centroids = "group"`,
    centroids of each group are connected along time levels;

  - if centroids are by `group_time` and `connect_centroids = "time"`,
    centroids of each time level are connected across groups;

  - for other `centroid_by` values, connections are ignored.

- show_ellipses:

  Logical; if `TRUE` (default), draws ellipses for selected factors.

- ellipse_by:

  Character vector describing which ellipses to draw. Any combination of
  `"group"`, `"time"`, `"group_time"`, or `"none"` (to disable). Default
  is `"group"`.

- ellipse_type:

  Ellipse type passed to
  [`ggplot2::stat_ellipse()`](https://ggplot2.tidyverse.org/reference/stat_ellipse.html),
  one of `"t"`, `"norm"`, `"euclid"`. Default `"t"`.

- ellipse_level:

  Numeric confidence level for ellipses (default `0.95`).

- point_size:

  Numeric size of sample points. Default `1.5`.

- point_alpha:

  Numeric opacity of sample points in `[0, 1]`. Default `0.6`.

- centroid_size:

  Numeric size for centroid points. Default `3`.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object representing the CAP ordination.

## Details

Axis labels are annotated with the percentage of constrained variance
explained by each CAP axis, computed as \\100 \* \lambda_k / \sum_j
\lambda_j\\, where \\\lambda_k\\ are the constrained eigenvalues from
`cap_res$eigenvalues`. Only positive parts of eigenvalues are used when
computing percentages.

The function uses
[`theme_phip()`](https://polymerase3.github.io/phiper/reference/theme_phip.md)
if it exists in the search path; otherwise it falls back to
[`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

Typically you will want to join sample-level metadata into
`cap_res$sample_coords` before plotting, e.g.:

    cap_res$sample_coords <- dplyr::left_join(
      cap_res$sample_coords,
      meta,             # data frame with sample_id, group, time, etc.
      by = "sample_id"
    )

## Examples

``` r
if (FALSE) { # \dontrun{
  # cap_res <- compute_capscale(dist_bc, ps = ps, formula = ~ type_person + age)
  cap_res$sample_coords <- dplyr::left_join(
    cap_res$sample_coords,
    meta,  # contains type_person, timepoint, etc.
    by = "sample_id"
  )

  plot_cap(
    cap_res,
    axes            = c(1, 2),
    group_col       = "type_person",
    time_col        = "timepoint",
    show_centroids  = TRUE,
    centroid_by     = "group_time",
    connect_centroids = "group",
    show_ellipses   = TRUE,
    ellipse_by      = c("group", "group_time")
  )
} # }
```
