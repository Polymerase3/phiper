# Internal helper: .ph_filter_pairs

Filter prevalence result rows by group pair(s), rank, features,
universe, and optional significance thresholds. Preserves
`ph_prev_result` class and `prev_meta` attribute when the input carries
them.

## Usage

``` r
.ph_filter_pairs(
  df,
  gA,
  gB,
  ranks = NULL,
  features = NULL,
  features_regex = FALSE,
  group_universe = NULL,
  universe_regex = FALSE,
  col_rank = "rank",
  col_feature = "feature",
  col_g1 = "group1",
  col_g2 = "group2",
  col_groupcol = "group_col",
  p_raw_max = NULL,
  q_bh_max = NULL,
  q_wbh_max = NULL,
  passed_only = FALSE,
  drop_na = TRUE,
  keep_cols = NULL
)
```
