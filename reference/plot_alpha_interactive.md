# Plot alpha diversity (precomputed) — interactive (plotly)

native plotly version of
[`plot_alpha()`](https://polymerase3.github.io/phiper/reference/plot_alpha.md).
mirrors the ggplot look: per-group boxplots with jittered points and
optional faceting by rank.

## Usage

``` r
plot_alpha_interactive(
  x,
  metric = c("richness", "shannon_diversity", "simpson_diversity", "pielou_evenness",
    "berger_parker_dominance"),
  group_col = "group",
  rank_col = "rank",
  filter_groups = NULL,
  filter_ranks = NULL,
  custom_colors = NULL,
  facet_by_rank = TRUE,
  ncol = 2,
  facet_scales = "fixed",
  interaction_only = FALSE,
  interaction_sep = NULL,
  x_order = NULL,
  x_labels = NULL,
  y_range = NULL,
  x_tickangle = 0,
  quartile_method = c("exclusive", "inclusive", "linear"),
  jitter_width = 0.25,
  point_size = 6,
  point_alpha = 0.85,
  text_size = 12,
  font_family = "Montserrat",
  show_grids = TRUE,
  ...
)
```

## Arguments

- x:

  a `"phip_alpha_diversity"` list (output of
  [`compute_alpha()`](https://polymerase3.github.io/phiper/reference/compute_alpha.md))
  or a single alpha-diversity data frame.

- metric:

  one of `"richness"`, `"shannon_diversity"`, `"simpson_diversity"`,
  `"pielou_evenness"`, or `"berger_parker_dominance"`.

- group_col:

  name of the grouping column in the alpha table (default `"group"` when
  `group_cols = NULL` in the computation step).

- rank_col:

  name of the rank column (default `"rank"`).

- filter_groups:

  optional character vector; keep only these group levels.

- filter_ranks:

  optional character vector; keep only these ranks.

- custom_colors:

  named character vector of hex colors for groups (like ggplot
  scale_fill_manual()).

- facet_by_rank:

  logical; facet by `rank_col` if multiple ranks present.

- ncol:

  integer; number of columns in facet wrap (default `2`).

- facet_scales:

  `"fixed"`, `"free_x"`, `"free_y"`, or `"free"`.

- interaction_only:

  logical; if `TRUE`, plot only the interaction table (when available in
  `x`). useful when
  [`compute_alpha()`](https://polymerase3.github.io/phiper/reference/compute_alpha.md)
  was run with `group_interaction = TRUE`.

- interaction_sep:

  character; separator used to join interaction labels. default is taken
  from `attr(x, "interaction_sep")` if present, otherwise `" * "`.

- x_order:

  optional character vector with desired group order (levels).

- x_labels:

  optional named character vector mapping group -\> label (used on x
  axis).

- y_range:

  optional numeric length-2; y axis range (e.g., c(0, 2300)).

- x_tickangle:

  numeric; tick label rotation in degrees (default 0; e.g., 25).

- quartile_method:

  one of c("exclusive","inclusive","linear"); passed to plotly box
  (default "exclusive" ~ ggplot).

- jitter_width:

  Numeric; horizontal jitter width for points.

- point_size:

  Numeric; size of jittered points.

- point_alpha:

  Numeric in (0,1); alpha for jittered points.

- text_size:

  Numeric; base text size for plot labels.

- font_family:

  Character; font family for plot text.

- show_grids:

  Logical; whether to show panel grid lines.

- ...:

  Reserved for future extensions; ignored.

## Value

A plotly htmlwidget.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_alpha_interactive(df, metric = "richness", group_col = "group")
} # }
```
