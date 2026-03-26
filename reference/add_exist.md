# Ensure an existence flag (all ones) on `data_long`

Appends/overwrites a column (default: "exist") filled with 1L on the
lazy `data_long` table. Preserves laziness; no collection is forced.

## Usage

``` r
add_exist(phip_data, exist_col = "exist", overwrite = FALSE)
```

## Arguments

- phip_data:

  A \<phip_data\> object.

- exist_col:

  Name of the existence column to append/overwrite.

- overwrite:

  If FALSE and the column exists, abort with a phiper-style error.

## Value

Modified \<phip_data\> with updated `data_long`.

## Examples

``` r
pd <- phip_load_example_data()
#> [10:22:51] INFO  Constructing <phip_data> object
#>                  -> new_phip_data()
#> [10:22:51] INFO  Fetching peptide metadata library via get_peptide_meta()
#> [10:22:51] INFO  Retrieving peptide metadata into DuckDB cache
#>                  -> get_peptide_meta(force_refresh = FALSE)
#> [10:22:51] INFO  Opened DuckDB connection
#>                    - cache dir:
#>                      /tmp/RtmptPywE9/phiper_cache1dd150858922/phip_cache.duckdb
#>                    - table: peptide_meta
#> [10:22:51] INFO  Starting download
#>                    - dest: /tmp/RtmptPywE9/file1dd17017956f.rds
#> [10:22:51] OK    Download succeeded (method = <getOption()>)
#> [10:22:51] OK    Checksum verified (SHA-256 match)
#> [10:22:54] OK    Download complete and loaded into R
#> [10:22:59] INFO  Importing sanitized metadata into DuckDB cache...
#> [10:23:01] OK    peptide_meta table created in DuckDB cache
#> [10:23:01] OK    Retrieving peptide metadata into DuckDB cache - done
#>                  -> elapsed: 10.259s
#> [10:23:01] OK    Peptide metadata acquired
#> [10:23:01] INFO  Validating <phip_data>
#>                  -> validate_phip_data()
#> [10:23:01] INFO  Checking structural requirements (shape & mandatory columns)
#> [10:23:01] INFO  Checking outcome family availability (exist / fold_change /
#>                  raw_counts)
#> [10:23:01] INFO  Checking collisions with reserved names
#>                    - subject_id, sample_id, timepoint, peptide_id, exist,
#>                      fold_change, counts_input, counts_hit
#> [10:23:01] INFO  Ensuring all columns are atomic (no list-cols)
#> [10:23:01] INFO  Checking key uniqueness
#> [10:23:01] INFO  Validating value ranges & types for outcomes
#> Warning: Missing values are always removed in SQL aggregation functions.
#> Use `na.rm = TRUE` to silence this warning
#> This warning is displayed once every 8 hours.
#> [10:23:02] INFO  Assessing sparsity (NA/zero prevalence vs threshold)
#>                    - warn threshold: 50%
#> [10:23:02] INFO  Checking peptide_id coverage against peptide_library
#> Warning: [10:23:02] WARN  peptide_id not found in peptide_library (e.g. 10003)
#>                  -> peptide library coverage.
#> [10:23:02] INFO  Validating comparisons table (if provided)
#> [10:23:02] INFO  Checking full grid completeness (peptide * sample)
#> Warning: [10:23:02] WARN  Counts table is not a full peptide * sample grid.
#>                  -> grid completeness
#>                    - observed rows: 78200
#>                    - expected rows: 156000.
#> Warning: [10:23:02] WARN  Grid remains incomplete (auto_expand = FALSE).
#>                  -> grid completeness
#>                    - observed rows: 78200
#>                    - expected rows: 156000.
#> [10:23:02] OK    Validating <phip_data> - done
#>                  -> elapsed: 0.613s
#> [10:23:02] OK    Constructing <phip_data> object - done
#>                  -> elapsed: 10.978s
pd <- add_exist(pd, overwrite = TRUE) # overwrites if present
#> [10:23:02] INFO  Ensuring existence flag on data_long
#>                  -> column: 'exist'; overwrite: TRUE
#> Warning: [10:23:02] WARN  Overwriting existing existence flag.
#>                  -> adding existence indicator
#>                    - column: "exist".
#> [10:23:02] OK    Ensuring existence flag on data_long - done
#>                  -> elapsed: 0.044s
```
