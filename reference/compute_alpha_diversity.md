# Compute alpha diversity per sample / group across ranks

Computes **richness**, **Shannon**, and **Simpson** diversity per sample
and per grouping variable at one or more **ranks** (columns describing
peptides).

## Usage

``` r
compute_alpha_diversity(
  x,
  group_cols = NULL,
  ranks = "peptide_id",
  fc_threshold = NULL,
  shannon_log = c("ln", "log2", "log10"),
  carry_cols = NULL,
  group_interaction = FALSE,
  interaction_only = FALSE,
  interaction_sep = " * "
)
```

## Arguments

- x:

  A `<phip_data>` object or a long `data.frame`.

- group_cols:

  Character vector of grouping columns, or `NULL` for a single aggregate
  (non-facetted) table. All columns must be present on the input.

- ranks:

  Character vector of **exact column names** to aggregate by. Typical
  values: `"peptide_id"` or taxonomy/lineage columns.

- fc_threshold:

  Numeric or `NULL`. If `NULL` (default), presence is `exist > 0`. If
  numeric, presence is `fold_change > fc_threshold`.

- shannon_log:

  One of `"ln"`, `"log2"`, `"log10"`; reporting base for the Shannon
  index (via base change from natural log).

- carry_cols:

  Optional character vector of extra columns to carry forward into the
  output if present (e.g. sample metadata in your table).

- group_interaction:

  Logical; also compute the interaction of all `group_cols` (default
  `FALSE`).

- interaction_only:

  Logical; if `TRUE`, return only the interaction table (requires
  `group_interaction = TRUE` and at least two `group_cols`).

- interaction_sep:

  Separator used for the interaction label (default `" * "`).

## Value

A **named list** of data frames with S3 class `"phip_alpha_diversity"`.
Each element (per `group_col`, plus optional interaction or
`"all_samples"`) contains: `rank`, `sample_id`, the grouping column (or
`group` when `group_cols = NULL`), any `carry_cols`, and the metrics:
`richness`, `shannon_diversity`, `simpson_diversity`.

## Details

### Ranks

Ranks are the peptide identities or characteristics you aggregate by.
They must be **exact column names**:

- For `<phip_data>`: columns in the PHIPER peptide library (e.g.
  `peptide_id`, lineage/taxa fields).

- For `data.frame`: columns present on your long count table.

### Presence rule

- Default: `exist > 0`.

- If `fc_threshold` is numeric, presence is
  `fold_change > fc_threshold`.

### Grouping, interactions, and interaction-only mode

- `group_cols` can be a character vector; the return value is a **named
  list** of data frames, one per `group_col`.

- If `group_cols = NULL`, a single non-facetted table is returned under
  the name `"all_samples"`.

- If `group_interaction = TRUE` (and you supplied \>= 2 `group_cols`),
  an additional element is computed for the interaction of all group
  columns, with labels joined by `interaction_sep`.

- If `interaction_only = TRUE`, you get **only** that interaction
  element (requires `group_interaction = TRUE` and at least 2
  `group_cols`).

## Examples

``` r
pd <- load_example_data()
#> [11:24:20] INFO  Constructing <phip_data> object
#>                  -> create_data()
#> [11:24:20] INFO  Fetching peptide metadata library via get_peptide_library()
#> [11:24:20] INFO  Retrieving peptide metadata into DuckDB cache
#>                  -> get_peptide_library(force_refresh = FALSE)
#> [11:24:20] INFO  Opened DuckDB connection
#>                    - cache dir:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/phip_cache.duckdb
#>                    - table: peptide_meta
#> [11:24:20] INFO  Starting download
#>                    - dest:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/combined_library_15.01.26.rds
#> [11:24:20] OK    Download succeeded (method = <getOption()>)
#> [11:24:20] OK    Checksum verified (SHA-256 match)
#> [11:24:23] OK    Download complete and loaded into R
#> [11:24:28] INFO  Importing sanitized metadata into DuckDB cache...
#> [11:24:30] OK    peptide_meta table created in DuckDB cache
#> [11:24:30] OK    Retrieving peptide metadata into DuckDB cache - done
#>                  -> elapsed: 10.256s
#> [11:24:30] OK    Peptide metadata acquired
#> [11:24:30] INFO  Validating <phip_data>
#>                  -> validate_phip_data()
#> [11:24:30] INFO  Checking structural requirements (shape & mandatory columns)
#> [11:24:30] INFO  Checking outcome family availability (exist / fold_change /
#>                  raw_counts)
#> [11:24:30] INFO  Checking collisions with reserved names
#>                    - subject_id, sample_id, timepoint, peptide_id, exist,
#>                      fold_change, counts_input, counts_hit
#> [11:24:30] INFO  Ensuring all columns are atomic (no list-cols)
#> [11:24:30] INFO  Checking key uniqueness
#> [11:24:30] INFO  Validating value ranges & types for outcomes
#> Warning: Missing values are always removed in SQL aggregation functions.
#> Use `na.rm = TRUE` to silence this warning
#> This warning is displayed once every 8 hours.
#> [11:24:30] INFO  Assessing sparsity (NA/zero prevalence vs threshold)
#>                    - warn threshold: 50%
#> [11:24:30] INFO  Checking peptide_id coverage against peptide_library
#> Warning: [11:24:31] WARN  peptide_id not found in peptide_library (e.g. 10003)
#>                  -> peptide library coverage.
#> [11:24:31] INFO  Checking full grid completeness (peptide * sample)
#> Warning: [11:24:31] WARN  Counts table is not a full peptide * sample grid.
#>                  -> grid completeness
#>                    - observed rows: 78200
#>                    - expected rows: 156000.
#> Warning: [11:24:31] WARN  Grid remains incomplete (auto_expand = FALSE).
#>                  -> grid completeness
#>                    - observed rows: 78200
#>                    - expected rows: 156000.
#> [11:24:31] OK    Validating <phip_data> - done
#>                  -> elapsed: 0.637s
#> [11:24:31] OK    Constructing <phip_data> object - done
#>                  -> elapsed: 10.896s
# phip_data input: peptide-level diversity by group
out <- compute_alpha_diversity(
  pd, group_cols = "group", ranks = "peptide_id"
)
#> [11:24:31] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group'; ranks: 'peptide_id'
#> [11:24:31] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 0.213s

# include interaction of multiple grouping variables
out2 <- compute_alpha_diversity(
  pd,
  group_cols = c("group", "timepoint"),
  ranks = c("peptide_id", "family", "genus"),
  group_interaction = TRUE
)
#> [11:24:31] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group', 'timepoint'; ranks: 'peptide_id',
#>                     'family', 'genus'
#> [11:24:33] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 1.768s

# interaction only (returns a single element named "group * timepoint")
out3 <- compute_alpha_diversity(
  pd,
  group_cols = c("group", "timepoint"),
  ranks = "peptide_id",
  group_interaction = TRUE,
  interaction_only = TRUE
)
#> [11:24:33] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group', 'timepoint'; ranks: 'peptide_id'
#> [11:24:33] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 0.208s

if (FALSE) { # \dontrun{
# data.frame input: ranks must be columns in the data
out_df <- compute_alpha_diversity(
  df_long, group_cols = NULL, ranks = "peptide_id"
)
} # }

# presence via fold-change
out_fc <- compute_alpha_diversity(
  pd, group_cols = "group", ranks = "peptide_id", fc_threshold = 1.5
)
#> [11:24:33] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group'; ranks: 'peptide_id'
#> [11:24:33] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 0.207s
```
