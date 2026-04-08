# Interactive volcano plot (log2 ratio vs -log10 p)

Interactive volcano plot (log2 ratio vs -log10 p)

## Usage

``` r
volcano_interactive(
  df,
  pair = NULL,
  rank = NULL,
  color_by = NULL,
  color_title = NULL,
  fc_cut = 1,
  p_cut = 0.05,
  p_mode = c("raw", "bh"),
  significant_colors = c(`not significant` = "#386cb0", `significant prior correction` =
    "#1b9e77", `significant post fdr correction` = "#e31a1c")
)
```

## Arguments

- df:

  A data frame with prevalence results.

- pair:

  optional group pair (character length-2).

- rank:

  optional single rank (character) to keep.

- color_by:

  optional named vector identifying peptide-library values to highlight,
  e.g. `c("is_flagellum" = TRUE)`.

- color_title:

  optional legend title when `color_by` is used.

- fc_cut:

  Numeric; absolute log2 fold-change cutoff.

- p_cut:

  Numeric; p-value cutoff.

- p_mode:

  One of `c("raw","bh")`; `"bh"` applies BH correction per-plot from
  `p_raw`.

- significant_colors:

  Named vector of colors for significance categories.

## Value

A `plotly` htmlwidget.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(3)
prev <- data.frame(
  rank       = "peptide_id",
  feature    = paste0("pep", 1:40),
  group1     = "A",
  group2     = "B",
  prop1      = runif(40),
  prop2      = runif(40),
  percent1   = runif(40, 0, 100),
  percent2   = runif(40, 0, 100),
  ratio      = runif(40, 0.1, 10),
  p_raw      = c(runif(10, 0, 0.01), runif(30, 0.1, 1)),
  n_peptides = 1L
)

# interactive volcano — hover to inspect each peptide
volcano_interactive(prev)

# BH correction
volcano_interactive(prev, p_mode = "bh", fc_cut = 1.5, p_cut = 0.01)
} # }
```
