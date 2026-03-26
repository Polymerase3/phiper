# Prevalence by group with pairwise tests (POP; per-rank FDR; BH & weighted BH)

Computes prevalence (counts, proportions, percentages) for features
defined by the requested `rank_cols` across one or more grouping
"universes" (`group_cols`, optional interaction), and performs
**pairwise** statistical tests between all group levels. Presence is
computed with a k-of-n POP rule within each sample (`pop_k_min`, default
1). Two testing modes are supported:

- **Unpaired:** Fisher’s exact test (2×2) per (rank, feature, group
  pair).

- **Paired:** McNemar’s exact (binomial) per subject (`paired = TRUE`
  requires `subject_id`).

P-values are adjusted **per rank** (single FDR family for each rank
across all requested universes/pairs/features) using BH and optional
weighted BH (`weight_mode = "peptide_count"`).

## Usage

``` r
ph_prevalence_compare(
  x,
  rank_cols,
  group_cols,
  exist_col = "exist",
  weight_mode = c("peptide_count", "none"),
  parallel = NULL,
  compute_ratios_db = TRUE,
  interaction = FALSE,
  combine_cols = NULL,
  interaction_sep = "::",
  collect = TRUE,
  register_name = NULL,
  pop_k_min = 1L,
  paired = FALSE,
  peptide_library = NULL
)
```

## Arguments

- x:

  A `phip_data` object with DuckDB backend **or** a data.frame/tibble
  with at least: `sample_id`, `peptide_id`, `exist_col`, all
  `group_cols`, and if `paired=TRUE` also `subject_id`. If `rank_cols`
  include non-peptide taxa, `x` must provide `peptide_library` with
  those columns.

- rank_cols:

  Character vector of rank columns, e.g. `c("peptide_id","species")`.

- group_cols:

  Character vector of grouping columns defining universes.

- exist_col:

  Name of the binary presence column (default `"exist"`).

- weight_mode:

  `"peptide_count"` (default) or `"none"`.

- parallel:

  Logical; compute Fisher p-values in parallel if possible (default
  `NULL` = auto).

- compute_ratios_db:

  Logical; compute simple ratios/delta ratios in SQL (unpaired only).

- interaction:

  Logical; also create an interaction universe of the first two
  `group_cols`.

- combine_cols:

  Optional length-2 character vector to build only that interaction
  universe.

- interaction_sep:

  Separator for interaction labels (default `"::"`).

- collect:

  Logical; if `TRUE` (default) return a collected tibble; otherwise a
  lazy table.

- register_name:

  Optional DuckDB table name for materialization (unpaired path).

- pop_k_min:

  Integer \\\ge\\ 1; k-of-n POP threshold per sample (default 1).

- paired:

  Logical; use paired design (McNemar exact) with `subject_id` (default
  `FALSE`). NOTE: can also be a character scalar naming the column that
  links related samples (e.g. "subject_id" or "dyade"). If so, only
  samples present in both groups for that identifier will be used for
  paired McNemar tests.

- peptide_library:

  Optional peptide metadata table used to map non-peptide rank columns.
  If `NULL`, the function will use `x$peptide_library`.

## Value

An object of class `ph_prev_result`, i.e., a tibble (or lazy table if
`collect = FALSE` on the unpaired path) with attributes:

- `prev_meta$m_by_rank`: named integer vector of `m_r` per rank,

- `prev_meta$pairs_by_universe`: tibble with `group_col`, `k_levels`,
  `n_pairs`,

- `prev_meta$pool_by_rank`: tibble with `rank`, `POOL`,

- other bookkeeping: `paired`, `weight_mode`, `pop_k_min`,
  `fdr_scope = "per-rank"`,

- `register_name` (unpaired path) and `view` (if available).

Columns include (subset may differ between paths): `view`, `rank`,
`feature`, `n_peptides`, `group_col`, `group1`, `group2`, `n1`, `N1`,
`prop1`, `percent1`, `n2`, `N2`, `prop2`, `percent2`, optional `ratio`,
`delta_ratio`, `p_raw`, `p_adj_rank`, `passed_rank_bh`,
`category_rank_bh`, `p_adj_rank_wbh`, `passed_rank_wbh`,
`category_rank_wbh`.

## Details

**universe construction.** For each `group_col` we create a universe of
levels present in the data (non-missing). If `interaction = TRUE` (or
`combine_cols` is provided), we also build a combined universe where the
group value is `<col1>::<col2>`. Denominators `N` are distinct sample
counts per `(group_col, group_value)`.

**presence (pop).** For each sample and each `(rank, feature)`, we count
the number of positive peptides contributing to that feature; a feature
is marked present if `k >= pop_k_min` (default 1). This yields
`n_present` per `(group_col, group_value, rank, feature)` and prevalence
`prop = n_present / N`.

**pairwise tests.** For each universe with `K` levels we form all
unordered pairs `K*(K-1)/2`. For each `(rank, feature, pair)`:

- unpaired: build the 2×2 table from (`n1`, `N1 - n1`, `n2`, `N2 - n2`)
  and run Fisher’s exact test (two-sided).

- paired: compute discordant counts `(n01, n10)` per subject and run
  McNemar’s exact binomial test (`binom.test(n01, n01+n10)`).

**fdr families and the number of comparisons.** Let `POOL_r` be the
number of unique features for a given **rank** `r` **across all
requested universes** (after presence aggregation and any missing value
filtering). Let `PAIRS = sum_u K_u*(K_u-1)/2` be the total number of
unordered level-pairs summed over all universes `u` (each `u` is one
`group_col` or an interaction universe). Then the size of the **single
FDR family** for rank `r` is:

\$\$ m_r \\=\\ POOL_r \times PAIRS. \$\$

All p-values produced for rank `r` (covering **all** universes and
**all** level pairs) are adjusted together as one family of size `m_r`.
This makes the per-rank FDR **stricter** when (i) there are many
features for that rank, and/or (ii) many universes or levels generate
many pairwise comparisons.

**bh (benjamini–hochberg).** Within each rank `r` (and `view` if
present), BH q-values are computed via `p.adjust(method="BH")` on the
vector of p-values of length `m_r` (excluding NAs). Reported columns:
`p_adj_rank`, `passed_rank_bh`, `category_rank_bh`.

**weighted bh.** If `weight_mode="peptide_count"`, each `(rank,feature)`
gets a base weight equal to the number of distinct peptides mapping to
that feature. Within each rank `r` (and `view` if present), **weights
are scaled to sum to** `m_r`:

\$\$ w_i^{\*} \\=\\ w_i \cdot \frac{m_r}{\sum_j w_j}. \$\$

We adjust using the standard weighted step-up rule on \\p_i /
w_i^{\*}\\. The resulting q-values are reported in `p_adj_rank_wbh`,
with flags `passed_rank_wbh` and labels `category_rank_wbh`. If
`weight_mode="none"`, all weights are 1 and wBH reduces to BH.

**logging of comparisons.** The function logs:

- `POOL_r` per rank,

- number of levels `K_u` and pair counts per universe,

- `PAIRS = sum_u K_u*(K_u-1)/2`,

- `m_r = POOL_r * PAIRS` for each rank. These values are also returned
  in the object metadata (see Value).

## Examples

``` r
if (FALSE) { # \dontrun{
res <- ph_prevalence_compare(pd, rank_cols=c("species"), group_cols=c("big_group"))
print(res)
} # }
```
