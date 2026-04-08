# Delta-prevalence vs Pooled Prevalence

Build a static ggplot showing the per-peptide shift in prevalence
(\\\Delta = group2 - group1\\) as a function of pooled prevalence
(\\(group1 + group2)/2\\). The input should be a tibble/data frame
produced by `ph_compute_prevalence()` or equivalent with columns
`group1`, `group2`, `prop1`, and `prop2`.

## Usage

``` r
deltaplot(
  prev_tbl,
  group_pair_values = NULL,
  group_labels = NULL,
  point_jitter_width = 0.005,
  point_jitter_height = 0.005,
  point_alpha = 0.25,
  point_size = 0.6,
  add_smooth = TRUE,
  smooth_k = 5,
  arrow_color = "red",
  arrow_head_length_mm = 4,
  arrow_x_frac = 0.97,
  plot_title = NULL,
  plot_subtitle = NULL,
  x_label = NULL,
  y_label = NULL
)
```

## Arguments

- prev_tbl:

  Data frame with columns `group1`, `group2`, `prop1`, `prop2`. Optional
  `feature` is used for row identity only.

- group_pair_values:

  Optional length-2 character vector `c(group1, group2)`. Use this when
  `prev_tbl` contains multiple group pairs.

- group_labels:

  Optional length-2 character vector of display labels
  `c(label_group1, label_group2)`. Defaults to `group1`/`group2`.

- point_jitter_width, point_jitter_height:

  Jitter amounts for the points. Defaults 0.005.

- point_alpha:

  Point transparency. Default 0.25.

- point_size:

  Point size. Default 0.6.

- add_smooth:

  Add a GAM smooth curve (`mgcv`). Default `TRUE`.

- smooth_k:

  Basis dimension `k` for the smooth. Default 5.

- arrow_color:

  Color for the directional arrows and labels. Default `"red"`.

- arrow_head_length_mm:

  Arrow head length in mm. Default 4.

- arrow_x_frac:

  Arrow X position as a fraction of the max pooled prevalence. Default
  0.97 (near the right edge).

- plot_title, plot_subtitle:

  Optional plot labels for the title/subtitle.

- x_label, y_label:

  Optional axis labels. Defaults are generated from the group labels.

## Value

A `ggplot` object.

## Details

The plot places each feature (peptide) as a point at:

- x-axis: pooled prevalence `(prop1 + prop2)/2`

- y-axis: prevalence shift `(prop2 - prop1)`

Points are optionally jittered for visibility. A dashed horizontal line
marks \\\Delta = 0\\. Optional arrows and labels indicate the direction
of increased prevalence for `group1` vs `group2`. If
`add_smooth = TRUE`, a GAM smooth is overlaid to summarize the trend
across pooled prevalence.

## Examples

``` r
# \donttest{
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
library(rlang)
ps <- load_example_data("small_mixture")

# pick the grouping column
group_col <- "group"

prev_res <- ph_prevalence_compare(
  ps,
  rank_cols  = "peptide_id",
  group_cols = group_col,
  collect    = TRUE
)
#> [10:26:58] INFO  prevalence_compare (per-rank fdr)
#> Warning: Unknown or uninitialised column: `peptide_library`.
#> Warning: Unknown or uninitialised column: `meta`.
#> [10:26:58] INFO  preparing input data
#>                    - ranks: peptide_id
#>                    - group_cols: group
#>                    - exist_col: exist
#>                    - weight_mode: peptide_count
#>                    - collect: TRUE
#>                    - pop_k_min: 1
#>                    - paired: FALSE
#> [10:26:59] INFO  ranks resolved
#>                    - - available: peptide_id
#> [10:26:59] INFO  grouping universes
#>                    - - per-column only: group
#> [10:26:59] INFO  computing cohort sizes (n) per universe
#> [10:26:59] INFO  computing presence per sample via k-of-n rule
#> [10:26:59] INFO  counting present samples per feature (pop, non-paired)
#> [10:26:59] INFO  fdr accounting
#>                    - pool per rank: peptide_id=5
#>                    - universes: group (k=2, pairs=1)
#>                    - pairs across universes (sum): 1
#>                    - total tests m per rank = pool * pairs: peptide_id=5
#> [10:26:59] INFO  building pairwise comparisons
#> [10:27:00] OK    materialized duckdb table
#>                    - name: ph_prev_20260408_102659
#>                    - computing p-values (fisher-only); then fdr per rank (bh /
#>                      wbh)
#> [10:27:01] OK    prevalence_compare (per-rank fdr) - done
#>                  -> elapsed: 2.142s
prev_tbl <- as.data.frame(prev_res)
pair_tbl <- unique(prev_tbl[, c("group1", "group2")])
group_pair <- c(pair_tbl$group1[1], pair_tbl$group2[1])

p <- deltaplot(
  prev_tbl,
  group_pair_values = group_pair,
  group_labels = group_pair,
  y_label = "Delta prevalence (group2 - group1)"
)
#> [10:27:01] INFO  Preparing delta prevalence plot.

print(p)

# }
```
