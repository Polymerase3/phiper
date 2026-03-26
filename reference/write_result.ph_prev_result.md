# Write method for ph_prev_result (ph_prevalence_compare output)

Write method for ph_prev_result (ph_prevalence_compare output)

## Usage

``` r
# S3 method for class 'ph_prev_result'
write_result(
  x,
  path,
  filter = c("none", "only_significant", "only_significant_fdr"),
  format = NULL,
  sheet_by_rank = TRUE,
  overwrite = FALSE,
  ...
)
```

## Arguments

- x:

  A `ph_prev_result` object (tibble/data.frame with attribute
  `prev_meta`).

- path:

  File path. Extension determines format (.xlsx/.csv/.parquet) unless
  format is provided.

- filter:

  One of c("none","only_significant","only_significant_fdr").

- format:

  Optional override for format: "xlsx","csv","parquet".

- sheet_by_rank:

  Logical; if TRUE create one sheet/file per rank when multiple ranks
  present.

- overwrite:

  Logical; allow overwriting existing files (default FALSE).

- ...:

  Reserved for future use.

## Value

Invisible `NULL`.
