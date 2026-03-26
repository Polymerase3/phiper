# Filter pairwise results by groups/ranks/features with optional q-value gates

Order-agnostic filter for rows comparing two group levels (e.g. "B" vs
"M12"), z opcjonalnym zawężeniem do ranków/feature'ów/uniwersów oraz
progów istotności. Działa zarówno na obiekcie klasy `ph_prev_result`,
jak i zwykłym `data.frame`.

## Usage

``` r
prev_filter_pairs(
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

## Arguments

- df:

  A `ph_prev_result` or `data.frame` with columns `group1`, `group2` and
  optionally `group_col`, `rank`, `feature`, p/q columns.

- gA, gB:

  character vector or length-1 character (pojedyncza para). Można też
  przekazać listę par (np. `list(c("B","M12"), c("T0","T2"))`) albo
  dwu-kolumnowy `data.frame`/`tibble` z nazwami kolumn `gA`, `gB`.

- ranks:

  Optional character vector of ranks to keep.

- features:

  Optional character vector or regex pattern(s).

- features_regex:

  Logical; treat `features` as regex pattern(s) (OR).

- group_universe:

  Optional character vector or regex pattern(s) for `group_col`.

- universe_regex:

  Logical; treat `group_universe` as regex pattern(s) (OR).

- col_rank, col_feature, col_g1, col_g2, col_groupcol:

  Column names in `df`.

- p_raw_max, q_bh_max, q_wbh_max:

  Optional numeric thresholds (e.g. 0.05).

- passed_only:

  Logical; if TRUE, keep rows with any `passed_* == TRUE` found.

- drop_na:

  Logical; drop rows with NA in group1/group2 before pair match.

- keep_cols:

  Optional character vector of columns to retain (NULL = all).

## Value

Filtered object of the same base type as `df`; if input is a
`ph_prev_result`, output keeps the class and augmented metadata.

## Details

**pary grup (symetrycznie).** Dla każdej zadanej pary {gA, gB} dobierane
są wiersze, w których `(min(group1,group2), max(group1,group2))` zgadza
się z `sort(c(gA, gB))`, niezależnie od kolejności w danych.

**filtrowanie.**

- `ranks`: wektor rang do zatrzymania.

- `features`: wektor lub regex; jeśli `features_regex=TRUE`, wszystkie
  wzorce są łączone logicznym OR (dopasowanie nieczułe na wielkość
  liter).

- `group_universe`: dokładne wartości `group_col` albo regex (gdy
  `universe_regex=TRUE`).

- progi istotności: `p_raw_max`, `q_bh_max`, `q_wbh_max` oraz
  `passed_only` (korzysta z kolumn `passed_rank_bh`/`passed_rank_wbh`,
  jeśli są dostępne).

**klasa.** Jeżeli wejście ma klasę `ph_prev_result`, wyjście zachowuje
ją oraz atrybut `prev_meta` (dodaje `subsetted = TRUE`, oraz ewentualną
notkę o filtrze).

## Examples

``` r
if (FALSE) { # \dontrun{
# single pair, regex on features:
# prev_filter_pairs(res, "B", "M12", ranks = "species",
#                   features = "flagellin|fliC", features_regex = TRUE,
#                   q_bh_max = 0.1)

# multiple pairs:
# pairs <- list(c("B","M12"), c("T0","T2"))
# prev_filter_pairs(res, pairs[[1]][1], pairs[[1]][2], ranks="peptide_id")
} # }
```
