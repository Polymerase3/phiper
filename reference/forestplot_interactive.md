# Interactive Forest Plot of Top/Bottom DELTA/Stouffer Statistics

Plotly version of
[`forestplot()`](https://polymerase3.github.io/phiper/reference/forestplot.md)
that highlights the most extreme positive and negative features within a
chosen rank. The plot shows the top `n_pos_each` features on the
positive side and the top `n_neg_each` features on the negative side,
ordered by the selected statistic (`statistic_to_plot`).

## Usage

``` r
forestplot_interactive(
  results_tbl,
  rank_of_interest,
  statistic_to_plot = c("T", "T_stand", "Z_from_p"),
  n_neg_each = 15,
  n_pos_each = 15,
  filter_significant = "none",
  sig_level = 0.05,
  left_label = "More in group1",
  right_label = "More in group2",
  arrow_length_frac = 0.35,
  label_x_gap_frac = 0.06,
  label_y_offset = 0,
  arrow_color = "red",
  arrow_linewidth = 0.6,
  arrow_head_length_mm = 3,
  use_diverging_colors = FALSE,
  show_grid = FALSE,
  base_text_pt = 12,
  font_family = "Montserrat",
  seg_width = 1.6,
  point_size = 11
)
```

## Arguments

- results_tbl:

  Data frame/tibble with at least:
  `rank, feature, group1, group2, design, T_obs, p_perm`.

- rank_of_interest:

  Character scalar specifying the rank to plot (e.g., `"species"`).

- statistic_to_plot:

  Which statistic to rank/plot: `"T"` (raw `T_obs`), `"T_stand"`
  (permutation-standardized), or `"Z_from_p"` (signed Z from permutation
  p).

- n_neg_each:

  Number of most negative features to show. Default 15.

- n_pos_each:

  Number of most positive features to show. Default 15.

- filter_significant:

  Column name to filter on, or `"none"` to disable filtering. If the
  column is numeric, keep rows where `col <= sig_level`.

- sig_level:

  Significance threshold used when `filter_significant` is numeric.
  Default 0.05.

- left_label:

  Text for the left arrow/side label.

- right_label:

  Text for the right arrow/side label.

- arrow_length_frac:

  Fraction of max \\\|T\|\\ used as half-length of arrows.

- label_x_gap_frac:

  Horizontal label gap for arrow labels beyond arrow tips (fraction of
  max \\\|T\|\\).

- label_y_offset:

  Additional vertical offset for arrow-end labels (y-axis units).

- arrow_color:

  Arrow/label color.

- arrow_linewidth:

  Arrow line width (plotly units).

- arrow_head_length_mm:

  Arrow head size (approximate, plotly units).

- use_diverging_colors:

  Logical; if `TRUE`, lines/points are shaded blue (negative) to red
  (positive) with higher contrast; otherwise monochrome.

- show_grid:

  Logical; show grid lines.

- base_text_pt:

  Base text size.

- font_family:

  Font family name.

- seg_width:

  Segment line width.

- point_size:

  Point size.

## Value

A list with `data` and `plot` (plotly object).

## Details

The y-axis lists feature names (sorted by the chosen statistic), and the
x-axis shows the signed effect size for each feature. Each feature is
drawn as a horizontal segment from zero to its statistic value, with a
point at the end of the segment. A dashed vertical line marks zero to
separate negative from positive shifts. The title and subtitle report
the contrast and how many features are shown.

If `use_diverging_colors = TRUE`, segments/points are colored by
magnitude on a blue-to-red scale (negative to positive), otherwise a
single color is used. Arrow annotations label the direction of
enrichment for each group.

`statistic_to_plot` controls the statistic used for both ranking and
plotting: raw `T_obs`, permutation-standardized `T_obs_stand`, or
`Z_from_p` (signed Z from permutation p-values).

## Examples

``` r
# \donttest{
  set.seed(1)
  n <- 20
  results_tbl <- data.frame(
    rank = rep("species", n),
    feature = paste0("feat_", seq_len(n)),
    group1 = "control",
    group2 = "treated",
    design = "case-control",
    T_obs = rnorm(n, sd = 2),
    p_perm = runif(n),
    T_obs_stand = rnorm(n),
    Z_from_p = qnorm(1 - runif(n) / 2) * sign(rnorm(n))
  )

  out <- forestplot_interactive(
    results_tbl,
    rank_of_interest = "species",
    statistic_to_plot = "T",
    n_neg_each = 5,
    n_pos_each = 5,
    left_label = "More in control",
    right_label = "More in treated",
    use_diverging_colors = TRUE,
    label_x_gap_frac = 0,
    label_y_offset = -0.05
  )

  out$plot

{"x":{"visdat":{"1c3ac2cf66c":["function () ","plotlyVisDat"]},"cur_data":"1c3ac2cf66c","attrs":{"1c3ac2cf66c":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0,"y":"feat_14","xend":-4.4293997743549998,"yend":"feat_14","type":"scatter","mode":"lines","showlegend":false,"hoverinfo":"text","text":"Feature: feat_14<br>n_peptides_used: NA<br>Stouffer T:<br />    -4.4294<br>Interpretation: More in control","line":{"color":"#1F4E79","width":1.6000000000000001},"opacity":0.94999999999999996,"inherit":true},"1c3ac2cf66c.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0,"y":"feat_3","xend":-1.6712572248200943,"yend":"feat_3","type":"scatter","mode":"lines","showlegend":false,"hoverinfo":"text","text":"Feature: feat_3<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.6713<br>Interpretation: More in control","line":{"color":"#7F9FC7","width":1.6000000000000001},"opacity":0.94999999999999996,"inherit":true},"1c3ac2cf66c.2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0,"y":"feat_6","xend":-1.6409367682360305,"yend":"feat_6","type":"scatter","mode":"lines","showlegend":false,"hoverinfo":"text","text":"Feature: feat_6<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.6409<br>Interpretation: More in control","line":{"color":"#81A1C8","width":1.6000000000000001},"opacity":0.94999999999999996,"inherit":true},"1c3ac2cf66c.3":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0,"y":"feat_1","xend":-1.2529076214846648,"yend":"feat_1","type":"scatter","mode":"lines","showlegend":false,"hoverinfo":"text","text":"Feature: feat_1<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.2529<br>Interpretation: More in control","line":{"color":"#97B0CF","width":1.6000000000000001},"opacity":0.94999999999999996,"inherit":true},"1c3ac2cf66c.4":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0,"y":"feat_13","xend":-1.2424811610836075,"yend":"feat_13","type":"scatter","mode":"lines","showlegend":false,"hoverinfo":"text","text":"Feature: feat_13<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.2425<br>Interpretation: More in control","line":{"color":"#97B0CF","width":1.6000000000000001},"opacity":0.94999999999999996,"inherit":true},"1c3ac2cf66c.5":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0,"y":"feat_4","xend":3.1905616042755831,"yend":"feat_4","type":"scatter","mode":"lines","showlegend":false,"hoverinfo":"text","text":"Feature: feat_4<br>n_peptides_used: NA<br>Stouffer T:<br />     3.1906<br>Interpretation: More in treated","line":{"color":"#8E1B10","width":1.6000000000000001},"opacity":0.94999999999999996,"inherit":true},"1c3ac2cf66c.6":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0,"y":"feat_11","xend":3.023562336901696,"yend":"feat_11","type":"scatter","mode":"lines","showlegend":false,"hoverinfo":"text","text":"Feature: feat_11<br>n_peptides_used: NA<br>Stouffer T:<br />     3.0236<br>Interpretation: More in treated","line":{"color":"#96271D","width":1.6000000000000001},"opacity":0.94999999999999996,"inherit":true},"1c3ac2cf66c.7":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0,"y":"feat_15","xend":2.2498618362862164,"yend":"feat_15","type":"scatter","mode":"lines","showlegend":false,"hoverinfo":"text","text":"Feature: feat_15<br>n_peptides_used: NA<br>Stouffer T:<br />     2.2499<br>Interpretation: More in treated","line":{"color":"#BC6158","width":1.6000000000000001},"opacity":0.94999999999999996,"inherit":true},"1c3ac2cf66c.8":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0,"y":"feat_18","xend":1.8876724213705984,"yend":"feat_18","type":"scatter","mode":"lines","showlegend":false,"hoverinfo":"text","text":"Feature: feat_18<br>n_peptides_used: NA<br>Stouffer T:<br />     1.8877<br>Interpretation: More in treated","line":{"color":"#CE7D75","width":1.6000000000000001},"opacity":0.94999999999999996,"inherit":true},"1c3ac2cf66c.9":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":0,"y":"feat_19","xend":1.6424423901961771,"yend":"feat_19","type":"scatter","mode":"lines","showlegend":false,"hoverinfo":"text","text":"Feature: feat_19<br>n_peptides_used: NA<br>Stouffer T:<br />     1.6424<br>Interpretation: More in treated","line":{"color":"#DA9088","width":1.6000000000000001},"opacity":0.94999999999999996,"inherit":true},"1c3ac2cf66c.10":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":[-4.4293997743549998,-1.6712572248200943,-1.6409367682360305,-1.2529076214846648,-1.2424811610836075,3.1905616042755831,3.023562336901696,2.2498618362862164,1.8876724213705984,1.6424423901961771],"y":["feat_14","feat_3","feat_6","feat_1","feat_13","feat_4","feat_11","feat_15","feat_18","feat_19"],"type":"scatter","mode":"markers","showlegend":false,"hoverinfo":"text","text":["Feature: feat_14<br>n_peptides_used: NA<br>Stouffer T:<br />    -4.4294<br>Interpretation: More in control","Feature: feat_3<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.6713<br>Interpretation: More in control","Feature: feat_6<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.6409<br>Interpretation: More in control","Feature: feat_1<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.2529<br>Interpretation: More in control","Feature: feat_13<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.2425<br>Interpretation: More in control","Feature: feat_4<br>n_peptides_used: NA<br>Stouffer T:<br />     3.1906<br>Interpretation: More in treated","Feature: feat_11<br>n_peptides_used: NA<br>Stouffer T:<br />     3.0236<br>Interpretation: More in treated","Feature: feat_15<br>n_peptides_used: NA<br>Stouffer T:<br />     2.2499<br>Interpretation: More in treated","Feature: feat_18<br>n_peptides_used: NA<br>Stouffer T:<br />     1.8877<br>Interpretation: More in treated","Feature: feat_19<br>n_peptides_used: NA<br>Stouffer T:<br />     1.6424<br>Interpretation: More in treated"],"marker":{"size":11,"color":["#1F4E79","#7F9FC7","#81A1C8","#97B0CF","#97B0CF","#8E1B10","#96271D","#BC6158","#CE7D75","#DA9088"],"opacity":1},"inherit":true}},"layout":{"margin":{"b":40,"l":45,"t":90,"r":10},"shapes":[{"type":"line","x0":0,"x1":0,"xref":"x","y0":0,"y1":1,"yref":"paper","line":{"dash":"dash","width":1,"color":"rgba(0,0,0,0.5)"}}],"title":{"text":"Top/Bottom Stouffer T - rank: species<br><sub>Contrast: More in control vs More in treated; shown: 5 neg + 5 pos<\/sub>","font":{"size":14,"family":"Montserrat"}},"xaxis":{"domain":[0,1],"automargin":true,"title":"Stouffer T (raw)","zeroline":false,"showgrid":false,"titlefont":{"size":12,"family":"Montserrat"},"tickfont":{"size":11,"family":"Montserrat"}},"yaxis":{"domain":[0,1],"automargin":true,"title":[],"categoryorder":"array","categoryarray":["feat_14","feat_3","feat_6","feat_1","feat_13","feat_19","feat_18","feat_15","feat_11","feat_4"],"showgrid":false,"tickfont":{"size":11,"family":"Montserrat"},"type":"category"},"font":{"family":"Montserrat","size":12},"annotations":[{"x":-1.5502899210242498,"y":1.02,"ax":0,"ay":1.02,"xref":"x","yref":"paper","axref":"x","ayref":"paper","showarrow":true,"arrowhead":3,"arrowsize":1.6000000000000001,"arrowwidth":2.6999999999999997,"arrowcolor":"red"},{"x":1.5502899210242498,"y":1.02,"ax":0,"ay":1.02,"xref":"x","yref":"paper","axref":"x","ayref":"paper","showarrow":true,"arrowhead":3,"arrowsize":1.6000000000000001,"arrowwidth":2.6999999999999997,"arrowcolor":"red"},{"x":-1.5502899210242498,"y":0.98999999999999999,"xref":"x","yref":"paper","text":"More in control","showarrow":false,"font":{"size":13,"color":"red","family":"Montserrat"},"xanchor":"right","yanchor":"bottom"},{"x":1.5502899210242498,"y":0.98999999999999999,"xref":"x","yref":"paper","text":"More in treated","showarrow":false,"font":{"size":13,"color":"red","family":"Montserrat"},"xanchor":"left","yanchor":"bottom"}],"hovermode":"closest","showlegend":false},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[0,-4.4293997743549998],"y":["feat_14","feat_14"],"type":"scatter","mode":"lines","showlegend":false,"hoverinfo":["text","text"],"text":["Feature: feat_14<br>n_peptides_used: NA<br>Stouffer T:<br />    -4.4294<br>Interpretation: More in control","Feature: feat_14<br>n_peptides_used: NA<br>Stouffer T:<br />    -4.4294<br>Interpretation: More in control"],"line":{"color":"#1F4E79","width":1.6000000000000001},"opacity":0.94999999999999996,"marker":{"color":"rgba(31,119,180,1)","line":{"color":"rgba(31,119,180,1)"}},"error_y":{"color":"rgba(31,119,180,1)"},"error_x":{"color":"rgba(31,119,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,-1.6712572248200943],"y":["feat_3","feat_3"],"type":"scatter","mode":"lines","showlegend":false,"hoverinfo":["text","text"],"text":["Feature: feat_3<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.6713<br>Interpretation: More in control","Feature: feat_3<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.6713<br>Interpretation: More in control"],"line":{"color":"#7F9FC7","width":1.6000000000000001},"opacity":0.94999999999999996,"marker":{"color":"rgba(255,127,14,1)","line":{"color":"rgba(255,127,14,1)"}},"error_y":{"color":"rgba(255,127,14,1)"},"error_x":{"color":"rgba(255,127,14,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,-1.6409367682360305],"y":["feat_6","feat_6"],"type":"scatter","mode":"lines","showlegend":false,"hoverinfo":["text","text"],"text":["Feature: feat_6<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.6409<br>Interpretation: More in control","Feature: feat_6<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.6409<br>Interpretation: More in control"],"line":{"color":"#81A1C8","width":1.6000000000000001},"opacity":0.94999999999999996,"marker":{"color":"rgba(44,160,44,1)","line":{"color":"rgba(44,160,44,1)"}},"error_y":{"color":"rgba(44,160,44,1)"},"error_x":{"color":"rgba(44,160,44,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,-1.2529076214846648],"y":["feat_1","feat_1"],"type":"scatter","mode":"lines","showlegend":false,"hoverinfo":["text","text"],"text":["Feature: feat_1<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.2529<br>Interpretation: More in control","Feature: feat_1<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.2529<br>Interpretation: More in control"],"line":{"color":"#97B0CF","width":1.6000000000000001},"opacity":0.94999999999999996,"marker":{"color":"rgba(214,39,40,1)","line":{"color":"rgba(214,39,40,1)"}},"error_y":{"color":"rgba(214,39,40,1)"},"error_x":{"color":"rgba(214,39,40,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,-1.2424811610836075],"y":["feat_13","feat_13"],"type":"scatter","mode":"lines","showlegend":false,"hoverinfo":["text","text"],"text":["Feature: feat_13<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.2425<br>Interpretation: More in control","Feature: feat_13<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.2425<br>Interpretation: More in control"],"line":{"color":"#97B0CF","width":1.6000000000000001},"opacity":0.94999999999999996,"marker":{"color":"rgba(148,103,189,1)","line":{"color":"rgba(148,103,189,1)"}},"error_y":{"color":"rgba(148,103,189,1)"},"error_x":{"color":"rgba(148,103,189,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,3.1905616042755831],"y":["feat_4","feat_4"],"type":"scatter","mode":"lines","showlegend":false,"hoverinfo":["text","text"],"text":["Feature: feat_4<br>n_peptides_used: NA<br>Stouffer T:<br />     3.1906<br>Interpretation: More in treated","Feature: feat_4<br>n_peptides_used: NA<br>Stouffer T:<br />     3.1906<br>Interpretation: More in treated"],"line":{"color":"#8E1B10","width":1.6000000000000001},"opacity":0.94999999999999996,"marker":{"color":"rgba(140,86,75,1)","line":{"color":"rgba(140,86,75,1)"}},"error_y":{"color":"rgba(140,86,75,1)"},"error_x":{"color":"rgba(140,86,75,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,3.023562336901696],"y":["feat_11","feat_11"],"type":"scatter","mode":"lines","showlegend":false,"hoverinfo":["text","text"],"text":["Feature: feat_11<br>n_peptides_used: NA<br>Stouffer T:<br />     3.0236<br>Interpretation: More in treated","Feature: feat_11<br>n_peptides_used: NA<br>Stouffer T:<br />     3.0236<br>Interpretation: More in treated"],"line":{"color":"#96271D","width":1.6000000000000001},"opacity":0.94999999999999996,"marker":{"color":"rgba(227,119,194,1)","line":{"color":"rgba(227,119,194,1)"}},"error_y":{"color":"rgba(227,119,194,1)"},"error_x":{"color":"rgba(227,119,194,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,2.2498618362862164],"y":["feat_15","feat_15"],"type":"scatter","mode":"lines","showlegend":false,"hoverinfo":["text","text"],"text":["Feature: feat_15<br>n_peptides_used: NA<br>Stouffer T:<br />     2.2499<br>Interpretation: More in treated","Feature: feat_15<br>n_peptides_used: NA<br>Stouffer T:<br />     2.2499<br>Interpretation: More in treated"],"line":{"color":"#BC6158","width":1.6000000000000001},"opacity":0.94999999999999996,"marker":{"color":"rgba(127,127,127,1)","line":{"color":"rgba(127,127,127,1)"}},"error_y":{"color":"rgba(127,127,127,1)"},"error_x":{"color":"rgba(127,127,127,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,1.8876724213705984],"y":["feat_18","feat_18"],"type":"scatter","mode":"lines","showlegend":false,"hoverinfo":["text","text"],"text":["Feature: feat_18<br>n_peptides_used: NA<br>Stouffer T:<br />     1.8877<br>Interpretation: More in treated","Feature: feat_18<br>n_peptides_used: NA<br>Stouffer T:<br />     1.8877<br>Interpretation: More in treated"],"line":{"color":"#CE7D75","width":1.6000000000000001},"opacity":0.94999999999999996,"marker":{"color":"rgba(188,189,34,1)","line":{"color":"rgba(188,189,34,1)"}},"error_y":{"color":"rgba(188,189,34,1)"},"error_x":{"color":"rgba(188,189,34,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,1.6424423901961771],"y":["feat_19","feat_19"],"type":"scatter","mode":"lines","showlegend":false,"hoverinfo":["text","text"],"text":["Feature: feat_19<br>n_peptides_used: NA<br>Stouffer T:<br />     1.6424<br>Interpretation: More in treated","Feature: feat_19<br>n_peptides_used: NA<br>Stouffer T:<br />     1.6424<br>Interpretation: More in treated"],"line":{"color":"#DA9088","width":1.6000000000000001},"opacity":0.94999999999999996,"marker":{"color":"rgba(23,190,207,1)","line":{"color":"rgba(23,190,207,1)"}},"error_y":{"color":"rgba(23,190,207,1)"},"error_x":{"color":"rgba(23,190,207,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-4.4293997743549998,-1.6712572248200943,-1.6409367682360305,-1.2529076214846648,-1.2424811610836075,3.1905616042755831,3.023562336901696,2.2498618362862164,1.8876724213705984,1.6424423901961771],"y":["feat_14","feat_3","feat_6","feat_1","feat_13","feat_4","feat_11","feat_15","feat_18","feat_19"],"type":"scatter","mode":"markers","showlegend":false,"hoverinfo":["text","text","text","text","text","text","text","text","text","text"],"text":["Feature: feat_14<br>n_peptides_used: NA<br>Stouffer T:<br />    -4.4294<br>Interpretation: More in control","Feature: feat_3<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.6713<br>Interpretation: More in control","Feature: feat_6<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.6409<br>Interpretation: More in control","Feature: feat_1<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.2529<br>Interpretation: More in control","Feature: feat_13<br>n_peptides_used: NA<br>Stouffer T:<br />    -1.2425<br>Interpretation: More in control","Feature: feat_4<br>n_peptides_used: NA<br>Stouffer T:<br />     3.1906<br>Interpretation: More in treated","Feature: feat_11<br>n_peptides_used: NA<br>Stouffer T:<br />     3.0236<br>Interpretation: More in treated","Feature: feat_15<br>n_peptides_used: NA<br>Stouffer T:<br />     2.2499<br>Interpretation: More in treated","Feature: feat_18<br>n_peptides_used: NA<br>Stouffer T:<br />     1.8877<br>Interpretation: More in treated","Feature: feat_19<br>n_peptides_used: NA<br>Stouffer T:<br />     1.6424<br>Interpretation: More in treated"],"marker":{"color":["#1F4E79","#7F9FC7","#81A1C8","#97B0CF","#97B0CF","#8E1B10","#96271D","#BC6158","#CE7D75","#DA9088"],"size":11,"opacity":1,"line":{"color":"rgba(31,119,180,1)"}},"error_y":{"color":"rgba(31,119,180,1)"},"error_x":{"color":"rgba(31,119,180,1)"},"line":{"color":"rgba(31,119,180,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}# }
```
