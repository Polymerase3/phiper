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
# \donttest{
library(dplyr)
library(rlang)
ps <- phip_load_example_data("small_mixture")

# pick the grouping column
group_col <- "group"

prev_res <- ph_prevalence_compare(
  ps,
  rank_cols  = "peptide_id",
  group_cols = group_col,
  collect    = TRUE
)
#> [10:23:26] INFO  prevalence_compare (per-rank fdr)
#> Warning: Unknown or uninitialised column: `peptide_library`.
#> Warning: Unknown or uninitialised column: `meta`.
#> [10:23:26] INFO  preparing input data
#>                    - ranks: peptide_id
#>                    - group_cols: group
#>                    - exist_col: exist
#>                    - weight_mode: peptide_count
#>                    - collect: TRUE
#>                    - pop_k_min: 1
#>                    - paired: FALSE
#> [10:23:26] INFO  ranks resolved
#>                    - - available: peptide_id
#> [10:23:26] INFO  grouping universes
#>                    - - per-column only: group
#> [10:23:26] INFO  computing cohort sizes (n) per universe
#> [10:23:27] INFO  computing presence per sample via k-of-n rule
#> [10:23:27] INFO  counting present samples per feature (pop, non-paired)
#> [10:23:27] INFO  fdr accounting
#>                    - pool per rank: peptide_id=5
#>                    - universes: group (k=2, pairs=1)
#>                    - pairs across universes (sum): 1
#>                    - total tests m per rank = pool * pairs: peptide_id=5
#> [10:23:27] INFO  building pairwise comparisons
#> [10:23:28] OK    materialized duckdb table
#>                    - name: ph_prev_20260326_102327
#>                    - computing p-values (fisher-only); then fdr per rank (bh /
#>                      wbh)
#> [10:23:29] OK    prevalence_compare (per-rank fdr) - done
#>                  -> elapsed: 2.236s
prev_tbl <- as.data.frame(prev_res)
pair_tbl <- unique(prev_tbl[, c("group1", "group2")])
group_pair <- c(pair_tbl$group1[1], pair_tbl$group2[1])

p <- ecdf_plot_interactive(
  prev_tbl,
  group_pair_values = group_pair,
  group_labels = group_pair,
  show_ks_test = FALSE
)

p

{"x":{"visdat":{"1dd14c43a8ee":["function () ","plotlyVisDat"]},"cur_data":"1dd14c43a8ee","attrs":{"1dd14c43a8ee":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter","mode":"lines","x":[0,0.050000000000000003,0.10000000000000001,0.69999999999999996],"y":[0.40000000000000002,0.59999999999999998,0.80000000000000004,1],"text":["<b>A<\/b><br>x: 0.0%<br>F(x): 40.0%","<b>A<\/b><br>x: 5.0%<br>F(x): 60.0%","<b>A<\/b><br>x: 10.0%<br>F(x): 80.0%","<b>A<\/b><br>x: 70.0%<br>F(x): 100.0%"],"hoverinfo":"text","line":{"width":2,"color":"#1f77b4","shape":"hv","opacity":1},"name":"A","inherit":true},"1dd14c43a8ee.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter","mode":"lines","x":[0,0.043478260869565216],"y":[0.59999999999999998,1],"text":["<b>B<\/b><br>x: 0.0%<br>F(x): 60.0%","<b>B<\/b><br>x: 4.3%<br>F(x): 100.0%"],"hoverinfo":"text","line":{"width":2,"color":"#d62728","shape":"hv","opacity":1},"name":"B","inherit":true},"1dd14c43a8ee.2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0.050000000000000003,"y":0,"xend":0.050000000000000003,"yend":1,"type":"scatter","mode":"lines","line":{"dash":"dot","width":1,"color":"#1f77b4"},"hoverinfo":"skip","showlegend":false,"inherit":true},"1dd14c43a8ee.3":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0,"y":0,"xend":0,"yend":1,"type":"scatter","mode":"lines","line":{"dash":"dot","width":1,"color":"#d62728"},"hoverinfo":"skip","showlegend":false,"inherit":true}},"layout":{"margin":{"b":60,"l":60,"t":70,"r":40},"title":{"text":"ECDF of per-peptide prevalence ( B vs A )"},"xaxis":{"domain":[0,1],"automargin":true,"title":"Prevalence","tickformat":".0%","range":[0,1]},"yaxis":{"domain":[0,1],"automargin":true,"title":"ECDF","tickformat":".0%","range":[0,1]},"legend":{"orientation":"h","x":0,"y":1.1000000000000001},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"type":"scatter","mode":"lines","x":[0,0.050000000000000003,0.10000000000000001,0.69999999999999996],"y":[0.40000000000000002,0.59999999999999998,0.80000000000000004,1],"text":["<b>A<\/b><br>x: 0.0%<br>F(x): 40.0%","<b>A<\/b><br>x: 5.0%<br>F(x): 60.0%","<b>A<\/b><br>x: 10.0%<br>F(x): 80.0%","<b>A<\/b><br>x: 70.0%<br>F(x): 100.0%"],"hoverinfo":["text","text","text","text"],"line":{"color":"#1f77b4","width":2,"shape":"hv","opacity":1},"name":"A","marker":{"color":"rgba(31,119,180,1)","line":{"color":"rgba(31,119,180,1)"}},"error_y":{"color":"rgba(31,119,180,1)"},"error_x":{"color":"rgba(31,119,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"type":"scatter","mode":"lines","x":[0,0.043478260869565216],"y":[0.59999999999999998,1],"text":["<b>B<\/b><br>x: 0.0%<br>F(x): 60.0%","<b>B<\/b><br>x: 4.3%<br>F(x): 100.0%"],"hoverinfo":["text","text"],"line":{"color":"#d62728","width":2,"shape":"hv","opacity":1},"name":"B","marker":{"color":"rgba(255,127,14,1)","line":{"color":"rgba(255,127,14,1)"}},"error_y":{"color":"rgba(255,127,14,1)"},"error_x":{"color":"rgba(255,127,14,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.050000000000000003,0.050000000000000003],"y":[0,1],"type":"scatter","mode":"lines","line":{"color":"#1f77b4","dash":"dot","width":1},"hoverinfo":["skip","skip"],"showlegend":false,"marker":{"color":"rgba(44,160,44,1)","line":{"color":"rgba(44,160,44,1)"}},"error_y":{"color":"rgba(44,160,44,1)"},"error_x":{"color":"rgba(44,160,44,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,0],"y":[0,1],"type":"scatter","mode":"lines","line":{"color":"#d62728","dash":"dot","width":1},"hoverinfo":["skip","skip"],"showlegend":false,"marker":{"color":"rgba(214,39,40,1)","line":{"color":"rgba(214,39,40,1)"}},"error_y":{"color":"rgba(214,39,40,1)"},"error_x":{"color":"rgba(214,39,40,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}# }
```
