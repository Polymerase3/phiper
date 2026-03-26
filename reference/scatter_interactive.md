# Interactive prevalence scatter for `ph_prev_result`

creates an interactive scatter (plotly) comparing prevalence in **group
a** vs **group b** for per-feature results produced by
[`ph_prevalence_compare()`](https://polymerase3.github.io/phiper/reference/ph_prevalence_compare.md).
accepts either a `ph_prev_result` (recommended) or a plain data.frame
with the same columns (`percent1`, `percent2`, `feature`, `group1`,
`group2`, etc.).

when a `ph_prev_result` is passed, you may specify a concrete pair of
group levels via `pair = c("A","B")`, and optionally restrict to a
`rank`, a `universe` (`group_col`) and/or `features`. internally,
subsetting uses
[`prev_filter_pairs()`](https://polymerase3.github.io/phiper/reference/prev_filter_pairs.md),
preserving metadata.

color mapping:

- by default (when `color_by = NULL`), points are colored by
  `category_rank_wbh` with three levels: "significant (wBH, per rank)",
  "nominal only", "not significant".

- when `color_by` is a peptide-level meta column name (e.g.
  `"is_flagellum"`), the function tries to join peptide metadata using
  `peptide_id` (or `feature` if `rank == "peptide_id"`), and color by
  that variable.

## Usage

``` r
scatter_interactive(
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
  peplib = NULL,
  category_colors = c(`significant (wBH, per rank)` = "#009E73", `nominal only` =
    "#E69F00", `not significant` = "#0072B2", `other/NA` = "#999999"),
  show_background = FALSE,
  background_df = NULL,
  background_name = "background peptides",
  background_color = "#808080",
  background_size = 6,
  background_alpha = 0.3,
  background_max_n = 3000,
  background_seed = 1L,
  point_size = 7,
  point_alpha = 0.85,
  point_line_width = 0.7,
  point_line_color = "rgba(0,0,0,0.6)",
  jitter_width_pp = 0,
  jitter_height_pp = 0,
  font_family = NULL,
  font_size = 12
)
```

## Arguments

- df:

  a `ph_prev_result` or a `data.frame` with columns: `percent1`,
  `percent2`, `feature`, `group1`, `group2`, `p_raw`, `p_adj_rank`,
  `p_adj_rank_wbh`, `category_rank_wbh`, optionally `n_peptides`,
  `rank`, `peptide_id`, `group_col`.

- pair:

  optional length-2 character, e.g.
  `c("kid_serum::T2","kid_serum::T8")`. only used when `df` is a
  `ph_prev_result`; filters rows to that exact pair (order-agnostic).

- rank:

  optional single rank (character) to keep (only used with
  `ph_prev_result`).

- universe:

  optional `group_col` value or regex (if `universe_regex = TRUE`); only
  used with `ph_prev_result`.

- features:

  optional character vector or regex patterns (if
  `features_regex = TRUE`); only used with `ph_prev_result`.

- features_regex:

  logical; treat `features` as regex patterns (or).

- universe_regex:

  logical; treat `universe` as regex pattern(s) (or).

- xlab, ylab:

  axis labels; if missing and `pair` is provided, they default to
  `pair[1]` and `pair[2]`.

- alpha:

  numeric in (0,1\]; used only for nominal labels; not the plotly alpha.

- prefer_flags:

  logical; reserved for future use (kept for back-compat).

- color_by:

  optional peptide-level meta column name to color by.

- color_title:

  optional legend title for `color_by`.

- peplib:

  Optional peptide metadata table used to color points by `color_by`
  when not present in `df`.

- category_colors:

  Optional named vector mapping categories to colors.

- show_background:

  Logical; plot background peptides if provided.

- background_df:

  Optional data frame of background points.

- background_name:

  Character; legend name for background points.

- background_color:

  Character; color for background points.

- background_size:

  Numeric; size for background points.

- background_alpha:

  Numeric in (0,1); opacity for background points.

- background_max_n:

  Integer; maximum number of background points to plot.

- background_seed:

  Integer; RNG seed for background downsampling.

- point_size:

  Numeric; marker size for points.

- point_alpha:

  Numeric in (0,1); marker opacity.

- point_line_width:

  Numeric; outline width for points.

- point_line_color:

  Character; outline color for points.

- jitter_width_pp:

  Numeric; jitter width in percentage points.

- jitter_height_pp:

  Numeric; jitter height in percentage points.

- font_family:

  Character; font family for plot text.

- font_size:

  Numeric; font size for plot text.

## Value

a `plotly` object.

## Examples

``` r
if (FALSE) { # \dontrun{
# typical usage with ph_prev_result:
p <- scatter_interactive(scatters,
  pair     = c("kid_serum::T2","kid_serum::T8"),
  rank     = "peptide_id",
  color_by = "is_flagellum",
  color_title = "Flagellum"
)
p
} # }
```
