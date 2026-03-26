# Retrieve the comparisons definition table

Returns the two-way comparison specifications stored inside a
**phip_data** object.

## Usage

``` r
get_comparisons(x)
```

## Arguments

- x:

  A valid `phip_data` object.

## Value

A tibble with columns: `comparison`, `group1`, `group2`, and `variable`.
