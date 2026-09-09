# Compute alpha diversity per sample / group across ranks

Computes **richness**, **Shannon**, **Simpson**, **Pielou's evenness**,
and **Berger-Parker dominance** per sample and per grouping variable at
one or more **ranks** (columns describing peptides).

## Usage

``` r
compute_alpha(
  x,
  group_cols = NULL,
  ranks = "peptide_id",
  mode = c("binary", "threshold", "abundance"),
  abundance_col = NULL,
  threshold = NULL,
  abundance_agg = c("mean", "sum", "max"),
  metrics = .alpha_metric_names,
  shannon_base = c("ln", "log2", "log10"),
  carry_cols = NULL,
  group_interaction = FALSE,
  interaction_only = FALSE,
  interaction_sep = " * ",
  shannon_log = NULL
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

- mode:

  One of `"binary"` (default), `"threshold"`, or `"abundance"`. Controls
  how per-category counts are computed; see **Modes** in Details.

- abundance_col:

  Character scalar naming the abundance column. Required for
  `mode = "abundance"`; optional for `mode = "threshold"` (defaults to
  `"fold_change"`); ignored for `mode = "binary"`.

- threshold:

  Numeric scalar. Required for `mode = "threshold"`. Ignored for other
  modes.

- abundance_agg:

  One of `"mean"` (default), `"sum"`, or `"max"`. Used only in
  `mode = "abundance"` when `rank != "peptide_id"` to aggregate
  peptide-level abundances within each rank category.

- metrics:

  Character vector of diversity metrics to compute. Any subset of:
  `"richness"`, `"shannon"`, `"simpson"`, `"pielou_evenness"`,
  `"berger_parker"`. Defaults to all five. Output columns follow the
  canonical names: `richness`, `shannon_diversity`, `simpson_diversity`,
  `pielou_evenness`, `berger_parker_dominance`.

- shannon_base:

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

  Separator used for the interaction label (default `" * "`). The
  resulting string (e.g. `"GroupA * T1"`) becomes both a list name and a
  value in the output column, so choose a separator that does not appear
  in your group levels.

- shannon_log:

  Deprecated. Use `shannon_base` instead.

## Value

A **named list** of data frames with S3 class `"phip_alpha_diversity"`.
Each element (per `group_col`, plus optional interaction or
`"all_samples"`) contains: `rank`, `sample_id`, the grouping column (or
`group` when `group_cols = NULL`), any `carry_cols`, and one column per
requested metric (see `metrics`). `pielou_evenness` is `NA` when
richness \<= 1; `berger_parker_dominance` is `NA` when richness == 0.

## Details

### Ranks

Ranks are the peptide identities or characteristics you aggregate by.
They must be **exact column names**:

- For `<phip_data>`: columns in the PHIPER peptide library (e.g.
  `peptide_id`, lineage/taxa fields).

- For `data.frame`: columns present on your long count table.

### Modes

The `mode` parameter controls what `n` represents in every diversity
formula:

- `"binary"` (default): presence = `exist > 0`. `n` = count of distinct
  enriched peptides per `(sample, rank_val)`.

- `"threshold"`: presence = `abundance_col > threshold`. `n` = count of
  peptides passing the threshold. `threshold` is required;
  `abundance_col` defaults to `"fold_change"` when `NULL`.

- `"abundance"`: no binarisation. `abundance_col` is required. At
  `rank = "peptide_id"`, `n` = the raw abundance value per peptide. At
  higher ranks, peptides are aggregated within each `(sample, rank_val)`
  using `abundance_agg` (`"mean"`, `"sum"`, or `"max"`).

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
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpcfwKE4/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> [11:56:34] INFO  Constructing <phip_data> object
#>                  -> create_data()
#> [11:56:34] INFO  Fetching peptide metadata library via get_peptide_library()
#> [11:56:34] INFO  Retrieving peptide metadata into DuckDB cache
#>                  -> get_peptide_library(force_refresh = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpcfwKE4/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> [11:56:34] INFO  Opened DuckDB connection
#>                    - cache dir:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/phip_cache.duckdb
#>                    - table: peptide_meta
#> [11:56:34] INFO  Starting download
#>                    - dest:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/combined_library_06.07.26.rds
#> [11:56:34] OK    Download succeeded (method = <getOption()>)
#> [11:56:34] OK    Checksum verified (SHA-256 match)
#> [11:56:37] OK    Download complete and loaded into R
#> [11:56:41] INFO  Importing sanitized metadata into DuckDB cache...
#> [11:56:42] OK    peptide_meta table created in DuckDB cache
#> [11:56:42] OK    Retrieving peptide metadata into DuckDB cache - done
#>                  -> elapsed: 8.768s
#> [11:56:42] OK    Peptide metadata acquired
#> [11:56:42] INFO  Validating <phip_data>
#>                  -> validate_phip_data()
#> [11:56:42] INFO  Checking structural requirements (shape & mandatory columns)
#> [11:56:42] INFO  Checking outcome family availability (exist / fold_change /
#>                  raw_counts)
#> [11:56:42] INFO  Checking collisions with reserved names
#>                    - subject_id, sample_id, timepoint, peptide_id, exist,
#>                      fold_change, counts_input, counts_hit
#> [11:56:42] INFO  Ensuring all columns are atomic (no list-cols)
#> [11:56:42] INFO  Checking key uniqueness
#> [11:56:42] INFO  Validating value ranges & types for outcomes
#> Warning: Missing values are always removed in SQL aggregation functions.
#> Use `na.rm = TRUE` to silence this warning
#> This warning is displayed once every 8 hours.
#> [11:56:43] INFO  Assessing sparsity (NA/zero prevalence vs threshold)
#>                    - warn threshold: 50%
#> [11:56:43] INFO  Checking peptide_id coverage against peptide_library
#> [11:56:43] INFO  Checking full grid completeness (peptide * sample)
#> [11:56:43] INFO  Counts table is not a full peptide * sample grid
#>                    - observed rows: 78200
#>                    - expected rows: 156000
#> Warning: [11:56:43] WARN  Grid remains incomplete (auto_expand = FALSE).
#>                  -> grid completeness
#>                    - observed rows: 78200
#>                    - expected rows: 156000.
#> [11:56:43] OK    Validating <phip_data> - done
#>                  -> elapsed: 0.336s
#> [11:56:43] OK    Constructing <phip_data> object - done
#>                  -> elapsed: 9.107s
# phip_data input: peptide-level diversity by group
out <- compute_alpha(
  pd, group_cols = "group", ranks = "peptide_id"
)
#> [11:56:43] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group'; ranks: 'peptide_id'
#> [11:56:43] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 0.182s

# include interaction of multiple grouping variables
out2 <- compute_alpha(
  pd,
  group_cols = c("group", "timepoint"),
  ranks = c("peptide_id", "family", "genus"),
  group_interaction = TRUE
)
#> [11:56:43] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group', 'timepoint'; ranks: 'peptide_id',
#>                     'family', 'genus'
#> [11:56:44] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 1.206s

# interaction only (returns a single element named "group * timepoint")
out3 <- compute_alpha(
  pd,
  group_cols = c("group", "timepoint"),
  ranks = "peptide_id",
  group_interaction = TRUE,
  interaction_only = TRUE
)
#> [11:56:44] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group', 'timepoint'; ranks: 'peptide_id'
#> [11:56:44] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 0.172s

if (FALSE) { # \dontrun{
# data.frame input: ranks must be columns in the data
out_df <- compute_alpha(
  df_long, group_cols = NULL, ranks = "peptide_id"
)
} # }

# threshold mode: presence = fold_change > 1.5
out_thr <- compute_alpha(
  pd, group_cols = "group", ranks = "peptide_id",
  mode = "threshold", threshold = 1.5
)
#> [11:56:44] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group'; ranks: 'peptide_id'
#> [11:56:45] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 0.17s
```
