# Convert raw PhIP-Seq output into a `phip_data` object

`phip_convert()` ingests a "long" table of PhIPsSeq read counts /
enrichment statistics, optionally expands it to the full
`sample_id x peptide_id` grid, and registers the result in DuckDB. The
function returns a fully initialised **`phip_data`** object that can be
queried with the tidy API used throughout the package.

## Usage

``` r
phip_convert(
  data_long_path,
  sample_id = NULL,
  peptide_id = NULL,
  subject_id = NULL,
  timepoint = NULL,
  exist = NULL,
  fold_change = NULL,
  counts_input = NULL,
  counts_hit = NULL,
  n_cores = 8,
  materialise_table = TRUE,
  auto_expand = FALSE,
  peptide_library = TRUE
)
```

## Arguments

- data_long_path:

  Character scalar. File or directory containing the *long-format*
  PhIP-Seq data. Allowed extensions are **`.csv`** and **`.parquet`**.
  Directories are treated as partitions of a parquet set.

- sample_id, peptide_id, subject_id, timepoint, exist, fold_change,
  counts_input, counts_hit:

  Optional character strings. Supply these only if your column names
  differ from the defaults (`"sample_id"`, `"peptide_id"`,
  `"subject_id"`, `"timepoint"`, `"exist"`, `"fold_change"`,
  `"counts_input"`, `"counts_hit"`). Each argument should contain the
  *name* of the column in the incoming data; `NULL` lets the default
  stand.

- n_cores:

  Integer \>= 1. Number of CPU threads DuckDB may use while reading and
  writing files.

- materialise_table:

  Logical. If `FALSE` the result is registered as a **view**; if `TRUE`
  the table is fully **materialised** and stored on disk, trading higher
  load time and storage for faster repeated queries.

- auto_expand:

  Logical. If `TRUE` and the incoming data are **not** a complete
  Cartesian product of `sample_id x peptide_id`, missing combinations
  are generated:

  - Columns that are constant within each `sample_id` (metadata) are
    copied to the new rows.

  - Non-recyclable measurement columns (`fold_change`, `exist`,
    `counts_input`, `counts_hit`, etc.) are initialised to 0. The
    expanded table replaces the original *in place*.

- peptide_library:

  Logical. If `TRUE` (default) `phip_convert()` will attempt to locate
  and attach the matching peptide-library metadata for downstream
  annotation. Set to `FALSE` to skip this step.

## Value

An S3 object of class **`phip_data`** containing:

- `data_long`:

  The (possibly expanded) long-format table.

- `comparisons`:

  A tibble of pre-computed group comparisons or `NULL` if none were
  supplied.

- `peptide_library`:

  Loaded peptide-library metadata (if `peptide_library = TRUE`).

- `meta`:

  List with DuckDB connection handles.

## Details

*Paths are resolved to absolute form* before any work begins, and
explicit checks confirm existence as well as extension validity.

## See also

- [`new_phip_data()`](https://polymerase3.github.io/phiper/reference/new_phip_data.md)
  for the object constructor.

- [`dplyr::tbl()`](https://dplyr.tidyverse.org/reference/tbl.html) to
  query DuckDB tables lazily.

## Examples

``` r
# Basic import, auto-detecting default column names
phip_obj <- phip_convert(
  data_long_path = phip_example_path("phip_mixture"),
  n_cores = 4,
  materialise_table = TRUE
)
#> [10:23:32] INFO  Constructing <phip_data> object
#>                  -> new_phip_data()
#> [10:23:32] INFO  Fetching peptide metadata library via get_peptide_meta()
#> [10:23:32] INFO  Retrieving peptide metadata into DuckDB cache
#>                  -> get_peptide_meta(force_refresh = FALSE)
#> [10:23:32] INFO  Opened DuckDB connection
#>                    - cache dir:
#>                      /tmp/RtmptPywE9/phiper_cache1dd17993b5e3/phip_cache.duckdb
#>                    - table: peptide_meta
#> [10:23:32] INFO  Starting download
#>                    - dest: /tmp/RtmptPywE9/file1dd14fbcbfab.rds
#> [10:23:32] OK    Download succeeded (method = <getOption()>)
#> [10:23:32] OK    Checksum verified (SHA-256 match)
#> [10:23:35] OK    Download complete and loaded into R
#> [10:23:40] INFO  Importing sanitized metadata into DuckDB cache...
#> [10:23:42] OK    peptide_meta table created in DuckDB cache
#> [10:23:42] OK    Retrieving peptide metadata into DuckDB cache - done
#>                  -> elapsed: 9.476s
#> [10:23:42] OK    Peptide metadata acquired
#> [10:23:42] INFO  Validating <phip_data>
#>                  -> validate_phip_data()
#> [10:23:42] INFO  Checking structural requirements (shape & mandatory columns)
#> [10:23:42] INFO  Checking outcome family availability (exist / fold_change /
#>                  raw_counts)
#> [10:23:42] INFO  Checking collisions with reserved names
#>                    - subject_id, sample_id, timepoint, peptide_id, exist,
#>                      fold_change, counts_input, counts_hit
#> [10:23:42] INFO  Ensuring all columns are atomic (no list-cols)
#> [10:23:42] INFO  Checking key uniqueness
#> [10:23:42] INFO  Validating value ranges & types for outcomes
#> [10:23:42] INFO  Assessing sparsity (NA/zero prevalence vs threshold)
#>                    - warn threshold: 50%
#> [10:23:42] INFO  Checking peptide_id coverage against peptide_library
#> Warning: [10:23:42] WARN  peptide_id not found in peptide_library (e.g. 10003)
#>                  -> peptide library coverage.
#> [10:23:42] INFO  Validating comparisons table (if provided)
#> [10:23:42] INFO  Checking full grid completeness (peptide * sample)
#> Warning: [10:23:43] WARN  Counts table is not a full peptide * sample grid.
#>                  -> grid completeness
#>                    - observed rows: 78200
#>                    - expected rows: 156000.
#> Warning: [10:23:43] WARN  Grid remains incomplete (auto_expand = FALSE).
#>                  -> grid completeness
#>                    - observed rows: 78200
#>                    - expected rows: 156000.
#> [10:23:43] OK    Validating <phip_data> - done
#>                  -> elapsed: 0.64s
#> [10:23:43] OK    Constructing <phip_data> object - done
#>                  -> elapsed: 10.207s

if (FALSE) { # \dontrun{
# Import a CSV and rename columns
phip_mem <- phip_convert(
  data_long_path = "data/phip_long.csv",
  sample_id      = "sample",
  peptide_id     = "pep"
)
} # }
```
