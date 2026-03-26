# Summarize FDR accounting and design for `ph_prev_result`

Provides a compact overview of the testing design and FDR accounting:

- per-rank pool sizes (`POOL_r`), per-universe level counts and pair
  counts,

- the number of pairwise comparisons aggregated across universes
  (`PAIRS`),

- the per-rank family sizes `m_r = POOL_r * PAIRS`,

- design flags (paired/unpaired), weight mode, `pop_k_min`, `fdr_scope`,
  and `view`.

## Usage

``` r
# S3 method for class 'ph_prev_result'
summary(object, ...)
```

## Arguments

- object:

  A `ph_prev_result` object returned by
  [`ph_prevalence_compare()`](https://polymerase3.github.io/phiper/reference/ph_prevalence_compare.md).

- ...:

  Unused, kept for S3 compatibility.

## Value

An object of class `summary.ph_prev_result` containing:

- `overview`: tibble with one row per rank (`rank`, `POOL`, `PAIRS`,
  `m_r`).

- `pairs_by_universe`: tibble with `group_col`, `k_levels`, `n_pairs`.

- `pool_by_rank`: tibble with `rank`, `POOL`.

- `totals`: named list with `PAIRS_total` and `m_total`.

- `design`: named list with `paired`, `weight_mode`, `pop_k_min`,
  `fdr_scope`, `view`.

The object has a custom [`print()`](https://rdrr.io/r/base/print.html)
for nice console output.

## Examples

``` r
if (FALSE) { # \dontrun{
s <- summary(res)   # where res is a ph_prev_result
s
} # }
```
