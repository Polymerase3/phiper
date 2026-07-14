# phiper

`phiper` is an end-to-end R toolkit for Phage ImmunoPrecipitation
sequencing (PhIP-seq) data analysis. It covers everything from quality
control and enrichment statistics to diversity metrics and
publication-ready figures — supporting both cross-sectional and
longitudinal study designs.

## Installation

You can install the development version of `phiper` from GitHub with
either `pak` or `devtools`:

``` r

# install.packages("pak")
pak::pak("Polymerase3/phiper")

# or, using devtools:
# install.packages("devtools")
devtools::install_github("Polymerase3/phiper")
```

> **Note:** `phiper` depends on
> [`phiperio`](https://github.com/Polymerase3/phiperio) for data import
> and class construction. It is installed automatically by the commands
> above.

## What can phiper do?

`phiper` is organised into four analysis domains, each with dedicated
computation and visualisation functions:

- **Enrichment counts** —
  [`plot_enrichment_counts()`](https://polymerase3.github.io/phiper/reference/plot_enrichment_counts.md)
  gives a quick overview of per-group peptide hit counts.
- **Alpha diversity** — richness, Shannon, Simpson, Pielou evenness, and
  Berger-Parker dominance, with optional significance testing and
  bracket annotations
  ([`compute_alpha()`](https://polymerase3.github.io/phiper/reference/compute_alpha.md),
  [`compute_alpha_significance()`](https://polymerase3.github.io/phiper/reference/compute_alpha_significance.md)).
- **Beta diversity** — pairwise distances, PCoA, CAP (constrained
  ordination), PERMANOVA, dispersion tests, and t-SNE
  ([`compute_distance()`](https://polymerase3.github.io/phiper/reference/compute_distance.md),
  [`compute_pcoa()`](https://polymerase3.github.io/phiper/reference/compute_pcoa.md),
  [`compute_capscale()`](https://polymerase3.github.io/phiper/reference/compute_capscale.md),
  [`compute_permanova()`](https://polymerase3.github.io/phiper/reference/compute_permanova.md),
  [`compute_dispersion()`](https://polymerase3.github.io/phiper/reference/compute_dispersion.md),
  [`compute_tsne()`](https://polymerase3.github.io/phiper/reference/compute_tsne.md)).
- **Global prevalence shift (delta)** — subject-level permutation tests
  for antigen-wide prevalence shifts between groups, with forest plots
  and ECDF visualisations
  ([`compute_delta()`](https://polymerase3.github.io/phiper/reference/compute_delta.md),
  [`forestplot()`](https://polymerase3.github.io/phiper/reference/forestplot.md),
  [`ecdf_plot()`](https://polymerase3.github.io/phiper/reference/ecdf_plot.md)).
- **Per-feature prevalence comparison (POP)** — Fisher’s exact and
  McNemar tests for every peptide or taxonomic rank, visualised as
  scatter plots and volcano plots
  ([`compute_pop()`](https://polymerase3.github.io/phiper/reference/compute_pop.md),
  [`scatter_static()`](https://polymerase3.github.io/phiper/reference/scatter_static.md),
  [`volcano_static()`](https://polymerase3.github.io/phiper/reference/volcano_static.md)).

### The `compute_` / `plot_` pattern

Every domain follows the same two-step workflow: a `compute_*()`
function returns a tidy data structure, which is then passed to one or
more `plot_*()` functions. This keeps computation and visualisation
cleanly separated and makes it easy to inspect, filter, or export
results before plotting.

### Static and interactive outputs

All visualisation functions come in two flavours: a static version built
on **ggplot2** (composable with the full tidyverse graphics ecosystem)
and an interactive version built on **plotly** (useful for exploratory
data analysis and HTML reports). Both use the built-in
[`theme_phip()`](https://polymerase3.github.io/phiper/reference/theme_phip.md)
and phiper colour scales for a consistent visual style.

## Documentation, examples, and vignettes

Full documentation, worked examples, and vignettes are available on the
package website:

**<https://polymerase3.github.io/phiper/>**

A hands-on [**Typical workflow
tutorial**](https://polymerase3.github.io/phiper/articles/typical-workflow.html)
walks through a full cross-sectional analysis — loading data, alpha/beta
diversity, POP, and DELTA — using a bundled dummy dataset.

The [**Peptide Library
Metadata**](https://polymerase3.github.io/phiper/articles/peptide-library.html)
vignette describes the Agilent, Twist, and Corona2 peptide library
annotation columns.

## Contributing

Bug reports, feature requests, and questions are welcome. Please open an
issue on the GitHub tracker:

- <https://github.com/Polymerase3/phiper/issues>

If you’d like to contribute, please read
[CONTRIBUTING.md](https://polymerase3.github.io/phiper/CONTRIBUTING.md).
