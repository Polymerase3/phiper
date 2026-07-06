# Interactive prevalence scatter for prevalence results

Creates an interactive scatter (plotly) comparing prevalence in **group
a** vs **group b** for per-feature results produced by
[`compute_pop()`](https://polymerase3.github.io/phiper/reference/compute_pop.md).
Accepts a plain data.frame with columns `percent1`, `percent2`,
`feature`, `group1`, `group2`, `p_raw`, etc.

When `pair` is given, subsetting uses
[`.ph_filter_pairs()`](https://polymerase3.github.io/phiper/reference/dot-ph_filter_pairs.md).

Color mapping:

- Default (`color_by = NULL`): BH-corrected p-values computed per-plot
  from `p_raw` ("significant (BH)", "nominal only", "not significant").

- `color_by` as a named vector highlights points matching the specified
  peptide-library values:


        color_by = c("is_flagellum" = TRUE)
        color_by = c("species" = "Staphylococcus aureus")
        color_by = c("is_flagellum" = TRUE, "species" = "Staphylococcus aureus")

## Usage

``` r
scatter_interactive(
  df,
  pair = NULL,
  rank = NULL,
  xlab = NULL,
  ylab = NULL,
  alpha = 0.05,
  color_by = NULL,
  color_title = NULL,
  peplib = NULL,
  background_df = NULL,
  ...
)
```

## Arguments

- df:

  A `data.frame` with columns `percent1`, `percent2`, `feature`,
  `group1`, `group2`, `p_raw`, optionally `n_peptides`, `rank`,
  `peptide_id`.

- pair:

  optional length-2 character, e.g.
  `c("kid_serum::T2","kid_serum::T8")`.

- rank:

  optional single rank (character) to keep.

- xlab, ylab:

  axis labels; defaults to `pair[1]`/`pair[2]` when `pair` is given.

- alpha:

  numeric in (0,1\]; significance threshold; not the plotly alpha.

- color_by:

  optional named vector identifying peptide-library values to highlight,
  e.g. `c("is_flagellum" = TRUE)` or
  `c("species" = "Staphylococcus aureus")`.

- color_title:

  optional legend title when `color_by` is used.

- peplib:

  Optional peptide metadata table used to resolve `color_by` when not
  available via the global library.

- background_df:

  Optional data frame of background points.

- ...:

  graphical parameters: `category_colors`, `show_background`,
  `background_name`, `background_color`, `background_size`,
  `background_alpha`, `background_max_n`, `background_seed`,
  `point_line_width`, `point_line_color`, `point_size`, `point_alpha`,
  `jitter_width_pp`, `jitter_height_pp`, `font_family`, `font_size`.

## Value

a `plotly` object.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- scatter_interactive(scatters,
  pair        = c("kid_serum::T2", "kid_serum::T8"),
  rank        = "peptide_id",
  color_by    = c("is_flagellum" = TRUE),
  color_title = "Flagellum"
)
p
} # }
```
