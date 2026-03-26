# Prepare sample metadata and (optional) comparisons for legacy import

Reads the `samples_file`, optionally reads `comparisons_file` and
`timepoints_file`, verifies that any dummy-coded grouping columns are
mutually exclusive (row sum == 1), and collapses them into a single
`group` column. When `timepoints_file` present, it merges tha `samples`
with `timepoints` and adds an person identifier and timepoint variable.

## Usage

``` r
.ph_legacy_prepare_metadata(
  samples_file,
  comparisons_file = NULL,
  timepoints_file = NULL,
  extra_cols = character()
)
```

## Arguments

- samples_file:

  Absolute path to the samples CSV/Parquet.

- comparisons_file:

  Absolute path to comparisons CSV/Parquet, or `NULL`.

- timepoints_file:

  Absolute path to timepoints CSV/Parquet, or `NULL`.

- extra_cols:

  Character vector of extra metadata columns to keep.

## Value

A list with elements `samples`, `comparisons` (or `NULL`), and
`extra_cols` (possibly augmented).
