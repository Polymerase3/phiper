# Prevalence comparison by group (POP framework)

Computes feature-level prevalence across binary group columns and
performs pairwise statistical tests:

- **Unpaired:** Fisher's exact test (2x2) per (rank, feature, group
  pair).

- **Paired:** McNemar's exact binomial test per subject.

Presence per sample is determined by a k-of-n rule: a feature is
considered present in a sample if at least `pop_k_min` of its
contributing peptides are positive. Each `group_col` must have exactly 2
non-missing levels in the data.

## Usage

``` r
compute_pop(
  x,
  rank_cols,
  group_cols,
  exist_col = "exist",
  pop_k_min = 1L,
  paired = FALSE,
  peptide_library = NULL
)
```

## Arguments

- x:

  A `phip_data` object or a data.frame/tibble with at least columns
  `sample_id`, `peptide_id`, `exist_col`, and all `group_cols`. If
  `paired` is set, the paired linking column must also be present. If
  `rank_cols` include non-peptide taxa, `x` must provide
  `peptide_library`.

- rank_cols:

  Character vector of rank columns, e.g. `c("peptide_id", "species")`.

- group_cols:

  Character vector of binary grouping columns. Each column must have
  exactly 2 non-missing levels in the data.

- exist_col:

  Name of the binary presence column (default `"exist"`).

- pop_k_min:

  Integer \>= 1; k-of-n POP threshold per sample (default 1).

- paired:

  `FALSE` (default) or a single string naming the column that links
  related samples (e.g. `"subject_id"`). When set, paired McNemar tests
  are used instead of Fisher.

- peptide_library:

  Optional peptide metadata table with `peptide_id` and any non-peptide
  rank columns. If `NULL`, taken from `x$peptide_library`.

## Value

A data.frame with one row per (rank, feature, group_col) comparison:
`rank`, `feature`, `group_col`, `group1`, `group2`, `n1`, `N1`, `prop1`,
`percent1`, `n2`, `N2`, `prop2`, `percent2`, `ratio`, `delta_ratio`
(unpaired only), `p_raw`, `n_peptides`. A `view` column is prepended
when the input carries a view attribute.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- compute_pop(pd, rank_cols = "species", group_cols = "group")
res
} # }
```
