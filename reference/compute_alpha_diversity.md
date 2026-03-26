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
pd <- phip_load_example_data()
# phip_data input: peptide-level diversity by group
out <- compute_alpha_diversity(
  pd, group_cols = "group", ranks = "peptide_id"
)
#> [10:23:04] INFO  Peptide library attached on main connection
#>                    - available columns: peptide_id, aa_seq, pos, len_seq,
#>                      Fullname, Description, is_IEDB_or_cntrl, is_auto ...(+36)
#> [10:23:04] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group'; ranks: 'peptide_id'
#> [10:23:04] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 0.213s

# include interaction of multiple grouping variables
out2 <- compute_alpha_diversity(
  pd,
  group_cols = c("group", "timepoint"),
  ranks = c("peptide_id", "family", "genus"),
  group_interaction = TRUE
)
#> [10:23:06] INFO  Peptide library attached on main connection
#>                    - available columns: peptide_id, aa_seq, pos, len_seq,
#>                      Fullname, Description, is_IEDB_or_cntrl, is_auto ...(+36)
#> [10:23:06] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group', 'timepoint'; ranks: 'peptide_id',
#>                     'family', 'genus'
#> [10:23:08] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 2.131s

# interaction only (returns a single element named "group * timepoint")
out3 <- compute_alpha_diversity(
  pd,
  group_cols = c("group", "timepoint"),
  ranks = "peptide_id",
  group_interaction = TRUE,
  interaction_only = TRUE
)
#> [10:23:10] INFO  Peptide library attached on main connection
#>                    - available columns: peptide_id, aa_seq, pos, len_seq,
#>                      Fullname, Description, is_IEDB_or_cntrl, is_auto ...(+36)
#> [10:23:10] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group', 'timepoint'; ranks: 'peptide_id'
#> [10:23:10] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 0.17s

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
#> [10:23:12] INFO  Peptide library attached on main connection
#>                    - available columns: peptide_id, aa_seq, pos, len_seq,
#>                      Fullname, Description, is_IEDB_or_cntrl, is_auto ...(+36)
#> [10:23:12] INFO  Computing alpha diversity (<phip_data>)
#>                  -> group_cols: 'group'; ranks: 'peptide_id'
#> [10:23:12] OK    Computing alpha diversity (<phip_data>) - done
#>                  -> elapsed: 0.174s
```
