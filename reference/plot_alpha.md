# Plot alpha diversity (richness/Shannon/Simpson) from precomputed results

Plot alpha diversity metrics from the precomputed output of
[`compute_alpha()`](https://polymerase3.github.io/phiper/reference/compute_alpha.md).
Supports filtering groups/ranks and optional faceting by rank. If the
input contains an interaction table, you can request only that via
`interaction_only = TRUE`.

## Usage

``` r
plot_alpha(
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
  jitter_width = 0.25,
  point_size = 1.8,
  point_alpha = 0.7,
  text_size = 12,
  font_family = "Montserrat",
  show_grids = TRUE,
  x_order = NULL,
  x_labels = NULL,
  y_range = NULL,
  x_tickangle = 0,
  significance = NULL,
  show_significance = FALSE,
  sig_p_threshold = 0.05,
  sig_step_increase = 0.05,
  sig_tip_length = 0.01,
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

  optional named vector of colors for groups.

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

- x_order:

  Optional character vector specifying the x-axis order.

- x_labels:

  Optional named character vector mapping x-axis labels.

- y_range:

  Optional numeric length-2 vector for y-axis limits.

- x_tickangle:

  Numeric; x-axis label angle in degrees.

- significance:

  Optional `"phip_alpha_significance"` object from
  [`compute_alpha_significance()`](https://polymerase3.github.io/phiper/reference/compute_alpha_significance.md);
  used to add significance brackets when `show_significance = TRUE`.

- show_significance:

  Logical; if `TRUE` and `significance` is supplied, add pairwise
  significance brackets via `ggsignif` (package must be installed).
  Default `FALSE`.

- sig_p_threshold:

  Numeric; only pairs with `p_adj <= sig_p_threshold` receive a bracket.
  Default `0.05`.

- sig_step_increase:

  Numeric; vertical step between stacked brackets. Passed to
  [`ggsignif::geom_signif()`](https://const-ae.github.io/ggsignif/reference/stat_signif.html).
  Default `0.05`.

- sig_tip_length:

  Numeric; length of bracket tips. Default `0.01`.

- ...:

  Reserved for future extensions; ignored.

## Value

a `ggplot` object.

## Details

- `x` can be:

  - the named list returned by
    [`compute_alpha()`](https://polymerase3.github.io/phiper/reference/compute_alpha.md)
    (class `"phip_alpha_diversity"`), or

  - a single data frame taken from that list.

- when `interaction_only = TRUE`, the function tries to select the
  element named by
  `paste(attr(x, "group_cols"), collapse = interaction_sep)`. if not
  found, it falls back to the first element whose name contains the
  separator.

## Examples

``` r
if (FALSE) { # \dontrun{
# precomputed alpha (list) -> boxplot per group
p <- plot_alpha(alpha_list, metric = "richness", group_col = "Cohort")

# only the interaction table (if available)
p_int <- plot_alpha(
  alpha_list,
  metric = "shannon_diversity",
  group_col = "Cohort * timepoint",
  interaction_only = TRUE
)
} # }
```
