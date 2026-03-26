# Static scatterplot of percent1 vs percent2 from `ph_prevalence_compare()`

A ggplot2 clone of
[`scatter_interactive()`](https://polymerase3.github.io/phiper/reference/scatter_interactive.md)
with the same filtering behavior and default coloring rules:

- If `color_by = NULL`, points are colored by significance category. If
  `category_rank_wbh` exists, it is used directly. Otherwise, categories
  are derived with priority: p_adj_rank_wbh -\> p_adj_rank -\> p_raw
  (threshold = `alpha`). If only one category remains, it auto-falls
  back to p-value bins for informative coloring.

- If `color_by` is provided, it joins peptide metadata (like
  interactive).

## Usage

``` r
scatter_static(
  df,
  pair = NULL,
  rank = NULL,
  universe = NULL,
  features = NULL,
  features_regex = FALSE,
  universe_regex = FALSE,
  xlab = NULL,
  ylab = NULL,
  alpha = 0.05,
  prefer_flags = TRUE,
  color_by = NULL,
  color_title = NULL,
  point_size = 2,
  point_alpha = 0.85,
  jitter_width_pp = 0,
  jitter_height_pp = 0,
  font_family = NULL,
  font_size = 12
)
```

## Arguments

- df:

  A `ph_prev_result` object or a data frame with prevalence results.

- pair:

  optional group pair (character length-2).

- rank:

  optional single rank (character) to keep.

- universe:

  optional `group_col` value or regex (if `universe_regex = TRUE`).

- features:

  optional character vector or regex patterns (if
  `features_regex = TRUE`).

- features_regex:

  logical; treat `features` as regex patterns.

- universe_regex:

  logical; treat `universe` as regex pattern(s).

- xlab, ylab:

  axis labels; if missing and `pair` is provided, they default to
  `pair[1]` and `pair[2]`.

- alpha:

  numeric in (0,1\]; used only for nominal labels.

- prefer_flags:

  logical; reserved for future use (kept for back-compat).

- color_by:

  optional peptide-level meta column name to color by.

- color_title:

  optional legend title for `color_by`.

- point_size:

  Numeric; marker size for points.

- point_alpha:

  Numeric in (0,1); marker opacity.

- jitter_width_pp:

  Numeric; jitter width in percentage points.

- jitter_height_pp:

  Numeric; jitter height in percentage points.

- font_family:

  Character; font family for plot text.

- font_size:

  Numeric; font size for plot text.

## Value

A ggplot object.
