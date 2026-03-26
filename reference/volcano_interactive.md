# Interactive volcano plot (log2 ratio vs -log10 p)

Interactive volcano plot (log2 ratio vs -log10 p)

## Usage

``` r
volcano_interactive(
  df,
  pair = NULL,
  rank = NULL,
  universe = NULL,
  features = NULL,
  features_regex = FALSE,
  universe_regex = FALSE,
  color_by = NULL,
  color_title = NULL,
  fc_cut = 1,
  p_cut = 0.05,
  p_mode = c("raw", "bh", "wbh"),
  significant_colors = c(`not significant` = "#386cb0", `significant prior correction` =
    "#1b9e77", `significant post fdr correction` = "#e31a1c")
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

- color_by:

  optional peptide-level meta column name to color by.

- color_title:

  optional legend title for `color_by`.

- fc_cut:

  Numeric; absolute log2 fold-change cutoff.

- p_cut:

  Numeric; p-value cutoff.

- p_mode:

  One of `c("raw","bh","wbh")` controlling which p-values to use.

- significant_colors:

  Named vector of colors for significance categories.

## Value

A `plotly` htmlwidget.
