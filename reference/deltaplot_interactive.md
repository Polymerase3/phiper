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
#> [10:23:19] INFO  prevalence_compare (per-rank fdr)
#> Warning: Unknown or uninitialised column: `peptide_library`.
#> Warning: Unknown or uninitialised column: `meta`.
#> [10:23:19] INFO  preparing input data
#>                    - ranks: peptide_id
#>                    - group_cols: group
#>                    - exist_col: exist
#>                    - weight_mode: peptide_count
#>                    - collect: TRUE
#>                    - pop_k_min: 1
#>                    - paired: FALSE
#> [10:23:19] INFO  ranks resolved
#>                    - - available: peptide_id
#> [10:23:20] INFO  grouping universes
#>                    - - per-column only: group
#> [10:23:20] INFO  computing cohort sizes (n) per universe
#> [10:23:20] INFO  computing presence per sample via k-of-n rule
#> [10:23:20] INFO  counting present samples per feature (pop, non-paired)
#> [10:23:20] INFO  fdr accounting
#>                    - pool per rank: peptide_id=5
#>                    - universes: group (k=2, pairs=1)
#>                    - pairs across universes (sum): 1
#>                    - total tests m per rank = pool * pairs: peptide_id=5
#> [10:23:20] INFO  building pairwise comparisons
#> [10:23:22] OK    materialized duckdb table
#>                    - name: ph_prev_20260326_102320
#>                    - computing p-values (fisher-only); then fdr per rank (bh /
#>                      wbh)
#> [10:23:22] OK    prevalence_compare (per-rank fdr) - done
#>                  -> elapsed: 2.368s
prev_tbl <- as.data.frame(prev_res)
pair_tbl <- unique(prev_tbl[, c("group1", "group2")])
group_pair <- c(pair_tbl$group1[1], pair_tbl$group2[1])

p <- deltaplot_interactive(
  prev_tbl,
  group_pair_values = group_pair,
  group_labels = group_pair,
  add_smooth = FALSE
)

p

{"x":{"visdat":{"1dd115472373":["function () ","plotlyVisDat"]},"cur_data":"1dd115472373","attrs":{"1dd115472373":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter","mode":"markers","x":[0.34954589144792403,0.025702491386530354,0.028316040241625162,0.046056873067282143,0.021863734873679829],"y":[-0.70127330887131389,0.04401193219510352,-0.054185718919616196,-0.10216560186585412,0.045090182058452423],"text":["<b>16196<\/b><br>A: 14/20 (70.0%)<br>B: 0/23 (0.0%)<br>p (Fisher, wBH): 2.47e-06<br>pooled: 35.0%<br>Δ: -70.0%","<b>16627<\/b><br>A: 0/20 (0.0%)<br>B: 1/23 (4.3%)<br>p (Fisher, wBH): 1.00e+00<br>pooled: 2.2%<br>Δ: 4.3%","<b>18003<\/b><br>A: 1/20 (5.0%)<br>B: 0/23 (0.0%)<br>p (Fisher, wBH): 7.75e-01<br>pooled: 2.5%<br>Δ: -5.0%","<b>24799<\/b><br>A: 2/20 (10.0%)<br>B: 0/23 (0.0%)<br>p (Fisher, wBH): 5.26e-01<br>pooled: 5.0%<br>Δ: -10.0%","<b>5243<\/b><br>A: 0/20 (0.0%)<br>B: 1/23 (4.3%)<br>p (Fisher, wBH): 1.00e+00<br>pooled: 2.2%<br>Δ: 4.3%"],"hoverinfo":"text","marker":{"size":6,"opacity":0.59999999999999998},"inherit":true},"1dd115472373.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter","mode":"lines","x":[0,1],"y":[0,0],"hoverinfo":"skip","line":{"dash":"dash"},"showlegend":false,"inherit":true}},"layout":{"margin":{"b":60,"l":60,"t":60,"r":80},"title":{"text":"Per-peptide shift vs pooled prevalence<br />                                ( B − A )"},"xaxis":{"domain":[0,1],"automargin":true,"title":"Pooled prevalence (A & B)","tickformat":".0%%","range":[0,1]},"yaxis":{"domain":[0,1],"automargin":true,"title":"Δ prevalence (B − A)","tickformat":".1%"},"annotations":[{"x":0.96999999999999997,"y":0.22304347826086954,"xref":"x","yref":"y","ax":0.96999999999999997,"ay":0,"axref":"x","ayref":"y","text":"","showarrow":true,"arrowhead":2,"arrowsize":0.59999999999999998,"arrowwidth":2,"arrowcolor":"red"},{"x":0.96999999999999997,"y":-0.22304347826086954,"xref":"x","yref":"y","ax":0.96999999999999997,"ay":0,"axref":"x","ayref":"y","text":"","showarrow":true,"arrowhead":2,"arrowsize":0.59999999999999998,"arrowwidth":2,"arrowcolor":"red"},{"x":0.93999999999999995,"y":0.028608695652173916,"xref":"x","yref":"y","text":"More in B","showarrow":false,"xanchor":"right","font":{"color":"red","size":12},"bgcolor":"rgba(255,255,255,0.65)","bordercolor":"red","borderwidth":1},{"x":0.93999999999999995,"y":-0.13382608695652171,"xref":"x","yref":"y","text":"More in A","showarrow":false,"xanchor":"right","font":{"color":"red","size":12},"bgcolor":"rgba(255,255,255,0.65)","bordercolor":"red","borderwidth":1}],"hovermode":"closest","showlegend":false},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"type":"scatter","mode":"markers","x":[0.34954589144792403,0.025702491386530354,0.028316040241625162,0.046056873067282143,0.021863734873679829],"y":[-0.70127330887131389,0.04401193219510352,-0.054185718919616196,-0.10216560186585412,0.045090182058452423],"text":["<b>16196<\/b><br>A: 14/20 (70.0%)<br>B: 0/23 (0.0%)<br>p (Fisher, wBH): 2.47e-06<br>pooled: 35.0%<br>Δ: -70.0%","<b>16627<\/b><br>A: 0/20 (0.0%)<br>B: 1/23 (4.3%)<br>p (Fisher, wBH): 1.00e+00<br>pooled: 2.2%<br>Δ: 4.3%","<b>18003<\/b><br>A: 1/20 (5.0%)<br>B: 0/23 (0.0%)<br>p (Fisher, wBH): 7.75e-01<br>pooled: 2.5%<br>Δ: -5.0%","<b>24799<\/b><br>A: 2/20 (10.0%)<br>B: 0/23 (0.0%)<br>p (Fisher, wBH): 5.26e-01<br>pooled: 5.0%<br>Δ: -10.0%","<b>5243<\/b><br>A: 0/20 (0.0%)<br>B: 1/23 (4.3%)<br>p (Fisher, wBH): 1.00e+00<br>pooled: 2.2%<br>Δ: 4.3%"],"hoverinfo":["text","text","text","text","text"],"marker":{"color":"rgba(31,119,180,1)","size":6,"opacity":0.59999999999999998,"line":{"color":"rgba(31,119,180,1)"}},"error_y":{"color":"rgba(31,119,180,1)"},"error_x":{"color":"rgba(31,119,180,1)"},"line":{"color":"rgba(31,119,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"type":"scatter","mode":"lines","x":[0,1],"y":[0,0],"hoverinfo":["skip","skip"],"line":{"color":"rgba(255,127,14,1)","dash":"dash"},"showlegend":false,"marker":{"color":"rgba(255,127,14,1)","line":{"color":"rgba(255,127,14,1)"}},"error_y":{"color":"rgba(255,127,14,1)"},"error_x":{"color":"rgba(255,127,14,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}# }
```
