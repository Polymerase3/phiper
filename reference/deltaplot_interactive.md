# Interactive Delta-prevalence vs Pooled Prevalence

Build an interactive plotly chart showing the per-peptide shift in
prevalence (\\\Delta = group2 - group1\\) as a function of pooled
prevalence (\\(group1 + group2)/2\\). The input should be a tibble/data
frame produced by `ph_compute_prevalence()` or equivalent with columns
`group1`, `group2`, `prop1`, and `prop2`.

## Usage

``` r
deltaplot_interactive(
  prev_tbl,
  group_pair_values = NULL,
  group_labels = NULL,
  point_alpha = 0.6,
  point_size = 6,
  add_smooth = TRUE,
  smooth_k = 5,
  arrow_color = "red",
  arrow_x_frac = 0.97,
  arrow_length_frac = 0.3,
  label_x_gap_frac = 0.03,
  label_y_gap_frac = 0.02,
  plot_title = NULL,
  plot_subtitle = NULL,
  point_jitter_width = 0.005,
  point_jitter_height = 0.005
)
```

## Arguments

- prev_tbl:

  Data frame with columns `group1`, `group2`, `prop1`, `prop2`. Optional
  `feature` is used for labels.

- group_pair_values:

  Optional length-2 character vector `c(group1, group2)`. Use this when
  `prev_tbl` contains multiple group pairs.

- group_labels:

  Optional length-2 character vector of display labels
  `c(label_group1, label_group2)`. Defaults to `group1`/`group2`.

- point_alpha:

  Point transparency. Default 0.6.

- point_size:

  Point size. Default 6.

- add_smooth:

  Add a GAM smooth curve (`mgcv`). Default `TRUE`.

- smooth_k:

  Basis dimension `k` for the smooth. Default 5.

- arrow_color:

  Color for the directional arrows and labels. Default `"red"`.

- arrow_x_frac:

  Arrow X position as a fraction of the x-range. Default 0.97.

- arrow_length_frac:

  Arrow length as a fraction of the y-range. Default 0.30.

- label_x_gap_frac:

  Horizontal label offset as a fraction of the x-range.

- label_y_gap_frac:

  Vertical label offset as a fraction of the y-range.

- plot_title, plot_subtitle:

  Optional plot labels for the title/subtitle.

- point_jitter_width, point_jitter_height:

  Jitter amounts. Defaults 0.005.

## Value

A plotly object.

## Details

The plot places each feature (peptide) as a point at:

- x-axis: pooled prevalence `(prop1 + prop2)/2`

- y-axis: prevalence shift `(prop2 - prop1)`

Points are optionally jittered for display, and hover text includes the
feature identifier plus prevalence (percent and proportion) and counts
(`n1`, `N1`, `n2`, `N2`) when available in the input table. A dashed
horizontal line marks \\\Delta = 0\\. Optional arrows and labels
indicate the direction of increased prevalence for `group1` vs `group2`.
If `add_smooth = TRUE`, a GAM smooth is overlaid to summarize the trend.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(2)
n <- 40
prev_tbl <- data.frame(
  feature = paste0("pep", seq_len(n)),
  group1  = "A",
  group2  = "B",
  prop1   = runif(n),
  prop2   = runif(n)
)

p <- deltaplot_interactive(
  prev_tbl,
  group_pair_values = c("A", "B"),
  group_labels      = c("Group A", "Group B"),
  add_smooth        = FALSE
)
p
} # }
```
