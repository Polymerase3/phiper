# Load Example PhIP-Seq Dataset as \<phip_data\>

Convenience helper to quickly load a shipped example dataset
("phip_mixture") into a `<phip_data>` object, suitable for downstream
analysis and visualization. This function wraps
[`phip_convert`](https://polymerase3.github.io/phiper/reference/phip_convert.md),
automatically supplying the correct parameters for the included example
data.

## Usage

``` r
phip_load_example_data(name = c("phip_mixture", "small_mixture"))
```

## Arguments

- name:

  Character scalar. Name of the shipped example dataset. Currently
  supported: `"phip_mixture"`, `"small_mixture"`.

## Value

A `<phip_data>` object created from the chosen example dataset.

## Examples

``` r
# Load the example data shipped with the package:
ex <- phip_load_example_data()
# ex is now a <phip_data> object ready for analysis

# Specify the dataset name explicitly
ex2 <- phip_load_example_data("small_mixture")

# Use with plotting functions
p = plot_enrichment_counts(ex, group_cols = "timepoint")
#> [10:23:43] INFO  Plotting enrichment counts (<phip_data>)
#>                  -> group_cols: 'timepoint'
#> [10:23:43] INFO  building enrichment count plot
#>                  -> grouping variable: 'timepoint'
#> [10:23:43] OK    plot built
#> [10:23:43] OK    building enrichment count plot - done
#>                  -> elapsed: 0.224s
#> [10:23:43] OK    Plotting enrichment counts (<phip_data>) - done
#>                  -> elapsed: 0.225s
```
