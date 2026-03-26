# Merge or join a `phip_data` object

Merge or join a `phip_data` object

## Usage

``` r
# S3 method for class 'phip_data'
merge(x, y, ..., confirm = interactive())
```

## Arguments

- x:

  A `phip_data` object.

- y:

  A data-frame-like object *or* another `phip_data`.

- ...:

  Arguments forwarded to either
  [`base::merge()`](https://rdrr.io/r/base/merge.html) or the chosen
  **dplyr** join (e.g. `by =`, `suffix =`, etc.).

- confirm:

  Logical. When `TRUE` *and* `type = "base"` *and* the session is
  interactive, the user is asked to confirm. Set to `FALSE` to skip the
  prompt (use sparingly-OOM risk remains).

## Value

A new `phip_data` whose `data_long` contains the merged / joined tibble.
