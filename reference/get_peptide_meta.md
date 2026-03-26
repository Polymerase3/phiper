# Retrieve the peptide metadata table into DuckDB, forcing atomic types

This function uses the phiper logging utilities for consistent,
ASCII-only progress messages and timing. Long-running steps are
bracketed with `.ph_with_timing()`, and informational/warning/error
messages are emitted via `.ph_log_info()`, `.ph_log_ok()`, `.ph_warn()`,
and `.ph_abort()`.

- Downloads the RDS once, sanitizes types (logical, character, numeric),
  and writes into a DuckDB cache on disk.

- Subsequent calls return a lazy `tbl_dbi` without loading into R
  memory.

## Usage

``` r
get_peptide_meta(force_refresh = FALSE)
```

## Arguments

- force_refresh:

  Logical. If `TRUE`, re-downloads and rebuilds the cache.

## Value

A `dplyr::tbl_dbi` pointing to the `peptide_meta` table. The returned
object carries an attribute `"duckdb_con"` with the open `DBI`
connection.

## Details

**Caching:** A temporary DuckDB database is created within a temp
directory (via
[`withr::local_tempdir()`](https://withr.r-lib.org/reference/with_tempfile.html)),
so each R session gets an isolated cache. The `force_refresh` argument
bypasses the fast path and rebuilds the cache.

**Sanitization:** Columns are stripped of attributes, list-columns are
flattened, textual `"NaN"` and numeric `NaN` are coerced to `NA`. Binary
0/1 fields are converted to `logical`, `"TRUE"/"FALSE"`
(case-insensitive) are converted to `logical`, and numeric-looking
character columns (beyond trivial 0/1) are converted to `numeric`. All
other atomic types are preserved.

**Integrity check:** If a SHA-256 checksum is provided, a warning is
logged when the downloaded file’s checksum does not match the expected
value.

## See also

[`dplyr::tbl()`](https://dplyr.tidyverse.org/reference/tbl.html),
[`DBI::dbConnect()`](https://dbi.r-dbi.org/reference/dbConnect.html),
[`duckdb::duckdb()`](https://r.duckdb.org/reference/duckdb.html)
