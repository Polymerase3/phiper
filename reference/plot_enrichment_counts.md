# Plot enrichment counts per group (and optional interaction)

Visualizes per-sample peptide enrichment counts across groups in a
`<phip_data>` object, with optional interaction of multiple grouping
variables. Presence is defined by `exist > 0`. The plot(s) are produced
by an internal helper `.plot_enrichment_counts_one()`.

## Usage

``` r
plot_enrichment_counts(
  phip_data,
  group_cols = NULL,
  prevalence_threshold = 0.05,
  custom_colors = NULL,
  binwidth = 1,
  group_interaction = FALSE,
  interaction_only = FALSE,
  interaction_sep = " * ",
  annotation_size = 4,
  ...
)
```

## Arguments

- phip_data:

  A `<phip_data>` object with `data_long` containing at least
  `sample_id`, `peptide_id`, and `exist`.

- group_cols:

  Character vector of grouping columns in `data_long`, or `NULL` to plot
  all samples together.

- prevalence_threshold:

  Numeric in `[0,1]`; minimum prevalence used by
  `.plot_enrichment_counts_one()` to filter/annotate bins (default
  `0.05`).

- custom_colors:

  Optional named vector for group colors passed through to
  `.plot_enrichment_counts_one()` (default `NULL`).

- binwidth:

  Numeric bin width for histograms (default `1`).

- group_interaction:

  Logical; also compute a plot for the interaction of all `group_cols`
  (default `FALSE`).

- interaction_only:

  Logical; if `TRUE`, return only the interaction plot (requires
  `group_interaction = TRUE` and at least two `group_cols`).

- interaction_sep:

  Character separator for interaction labels (default `" * "`).

- annotation_size:

  Numeric; size of the in-plot threshold annotations (passed to
  `geom_text(size = ...)`). Typical range 3–6. Default `4`.

- ...:

  Reserved for future extensions; ignored.

## Value

- A single plot object (when `group_cols = NULL`, or when only one plot
  is produced), or

- A named list of plot objects (when multiple plots are produced).

## Details

- If `group_cols = NULL`, a single plot is returned for all samples.

- If `group_cols` is a character vector, a list of plots is returned
  (one per grouping column) unless `interaction_only = TRUE`.

- If `group_interaction = TRUE` and at least two `group_cols` are
  supplied, an additional interaction plot is created whose label joins
  groups using `interaction_sep`.

- If `interaction_only = TRUE`, only the interaction plot is returned
  (this requires `group_interaction = TRUE` and at least two
  `group_cols`).

## Examples

``` r
# per-group plots
pd <- load_example_data()
p <- plot_enrichment_counts(pd, group_cols = c("group","timepoint"))
#> [13:05:50] INFO  Plotting enrichment counts (<phip_data>)
#>                  -> group_cols: 'group', 'timepoint'
#> [13:05:50] INFO  building enrichment count plot
#>                  -> grouping variable: 'group'
#> [13:05:51] OK    plot built
#> [13:05:51] OK    building enrichment count plot - done
#>                  -> elapsed: 0.253s
#> [13:05:51] INFO  building enrichment count plot
#>                  -> grouping variable: 'timepoint'
#> [13:05:51] OK    plot built
#> [13:05:51] OK    building enrichment count plot - done
#>                  -> elapsed: 0.25s
#> [13:05:51] OK    Plotting enrichment counts (<phip_data>) - done
#>                  -> elapsed: 0.504s

# add interaction plot
p2 <- plot_enrichment_counts(pd,
  group_cols = c("group","timepoint"),
  group_interaction = TRUE
)
#> [13:05:51] INFO  Plotting enrichment counts (<phip_data>)
#>                  -> group_cols: 'group', 'timepoint'
#> [13:05:51] INFO  building enrichment count plot
#>                  -> grouping variable: 'group'
#> [13:05:51] OK    plot built
#> [13:05:51] OK    building enrichment count plot - done
#>                  -> elapsed: 0.249s
#> [13:05:51] INFO  building enrichment count plot
#>                  -> grouping variable: 'timepoint'
#> [13:05:51] OK    plot built
#> [13:05:51] OK    building enrichment count plot - done
#>                  -> elapsed: 0.229s
#> [13:05:51] INFO  building enrichment count plot
#>                  -> grouping variable: '..phip_interaction..'
#> [13:05:52] OK    plot built
#> [13:05:52] OK    building enrichment count plot - done
#>                  -> elapsed: 0.247s
#> [13:05:52] OK    Plotting enrichment counts (<phip_data>) - done
#>                  -> elapsed: 0.764s

# interaction only
p3 <- plot_enrichment_counts(pd,
  group_cols = c("group","timepoint"),
  group_interaction = TRUE,
  interaction_only = TRUE
)
#> [13:05:52] INFO  Plotting enrichment counts (<phip_data>)
#>                  -> group_cols: 'group', 'timepoint'
#> [13:05:52] INFO  building enrichment count plot
#>                  -> grouping variable: '..phip_interaction..'
#> [13:05:52] OK    plot built
#> [13:05:52] OK    building enrichment count plot - done
#>                  -> elapsed: 0.265s
#> [13:05:52] OK    Plotting enrichment counts (<phip_data>) - done
#>                  -> elapsed: 0.273s
```
