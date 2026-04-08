# ECDF of Per-peptide Prevalence for Two Groups

Plotly version of `ecdf_prevalence()`, showing ECDF curves for two
groups based on per-peptide prevalence values. The plot can annotate
median shifts and an optional KS test summary.

## Usage

``` r
ecdf_plot_interactive(
  prev_tbl,
  group_pair_values = NULL,
  group_labels = NULL,
  line_width_px = 2,
  line_alpha = 1,
  group1_line_color = "#1f77b4",
  group2_line_color = "#d62728",
  show_median_lines = TRUE,
  show_ks_test = TRUE,
  plot_title = NULL,
  plot_subtitle = NULL
)
```

## Arguments

- prev_tbl:

  Data frame with columns `group1`, `group2`, `prop1`, `prop2`.

- group_pair_values:

  Optional length-2 character vector `c(group1, group2)`. Use this when
  `prev_tbl` contains multiple group pairs.

- group_labels:

  Optional length-2 character vector of display labels
  `c(label_group1, label_group2)`. Defaults to `group1`/`group2`.

- line_width_px:

  Line width for ECDF steps (plotly units). Default 2.0.

- line_alpha:

  Line alpha for ECDF steps. Default 1.0.

- group1_line_color, group2_line_color:

  Line colors for group1 and group2.

- show_median_lines:

  Logical; add median lines. Default `TRUE`.

- show_ks_test:

  Logical; add KS test summary to subtitle. Default `TRUE`.

- plot_title, plot_subtitle:

  Optional plot labels.

## Value

A plotly object.

## Details

Each group is represented by a step curve showing the cumulative
fraction of features with prevalence less than or equal to a given
value. Vertical median lines can be added for each group, and the
subtitle can include the KS statistic and p-value with the median
difference.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(5)
n <- 50
prev_tbl <- data.frame(
  feature = paste0("pep", seq_len(n)),
  group1  = "A",
  group2  = "B",
  prop1   = runif(n),
  prop2   = runif(n)
)

p <- ecdf_plot_interactive(
  prev_tbl,
  group_pair_values = c("A", "B"),
  group_labels      = c("Group A", "Group B"),
  show_ks_test      = TRUE
)
p
} # }
```
