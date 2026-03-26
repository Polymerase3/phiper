# Construct a **phip_data** object

Creates a fully-validated S3 object that bundles the tidy PhIP-Seq
counts (`data_long`), optional comparison definitions, a peptide-library
annotation table, and other metadata. The function performs a minimal
sanity check on *comparisons* before returning the object (validation of
the data itself happens via
[`validate_phip_data()`](https://polymerase3.github.io/phiper/reference/validate_phip_data.md)
helper).

## Usage

``` r
new_phip_data(
  data_long,
  comparisons,
  peptide_library = TRUE,
  auto_expand = TRUE,
  materialise_table = TRUE,
  meta = list()
)
```

## Arguments

- data_long:

  A tidy data frame (or `tbl_lazy`) with one row per `peptide_id` x
  `sample_id` combination. **Required.**

- comparisons:

  A data frame describing two-way contrasts (`comparison`, `group1`,
  `group2`, `variable`); defaults to an empty tibble if `NULL`.

- peptide_library:

  A data frame with one row per `peptide_id` and its annotations. If
  `NULL`, the package’s current default library is used.

- auto_expand:

  Logical. If `TRUE` and the input is **not** already the full Cartesian
  product of `sample_id` x `peptide_id`, the function fills in the
  missing combinations.

  - Columns that are constant within a `sample_id` (metadata) are
    duplicated to the newly created rows.

  - Measurement columns such as `fold_change`, `exist`, raw counts, or
    any other non-recyclable fields are initialised to 0. The expanded
    table replaces `data_long` in place.

- materialise_table:

  Logical. If `FALSE` (default) the result is registered as a **view**.
  If `TRUE` the result is fully **materialised** and stored as a
  physical table, which speeds up repeated queries at the cost of extra
  memory/disk.

- meta:

  Optional named list of metadata flags to pre-populate the `meta` slot
  (rarely needed by users).

## Value

An object of class `"phip_data"`.

## Examples

``` r
if (FALSE) { # \dontrun{
## minimal constructor call
pd <- new_phip_data(
  data_long = tidy_counts,
  comparisons = NULL,
  peptide_library = TRUE
)
} # }
```
