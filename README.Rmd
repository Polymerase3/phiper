---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![Codecov test coverage](https://codecov.io/gh/Polymerase3/phiper/graph/badge.svg)](https://app.codecov.io/gh/Polymerase3/phiper)
[![R-CMD-check](https://github.com/Polymerase3/phiper/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Polymerase3/phiper/actions/workflows/R-CMD-check.yaml)
[![pkgcheck](https://github.com/Polymerase3/phiper/workflows/pkgcheck/badge.svg)](https://github.com/Polymerase3/phiper/actions?query=workflow%3Apkgcheck)
[![version](https://img.shields.io/github/r-package/v/Polymerase3/phiper?label=version)](https://github.com/Polymerase3/phiper)
[![Project Status: Active/Unstable](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

# phiper

`phiper` is an end-to-end R toolkit for Phage ImmunoPrecipitation sequencing
(PhIP-seq) data analysis. It covers everything from quality control and
enrichment statistics to diversity metrics and publication-ready figures —
supporting both cross-sectional and longitudinal study designs.

## Installation

You can install the development version of `phiper` from GitHub with either `pak` or `devtools`:

``` r
# install.packages("pak")
pak::pak("Polymerase3/phiper")

# or, using devtools:
# install.packages("devtools")
devtools::install_github("Polymerase3/phiper")
```

> **Note:** `phiper` depends on
> [`phiperio`](https://github.com/Polymerase3/phiperio) for data import and
> class construction. It is installed automatically by the commands above.

## What can phiper do?

`phiper` is organised into four analysis domains, each with dedicated
computation and visualisation functions:

- **Enrichment counts** — `plot_enrichment_counts()` gives a quick overview of
  per-group peptide hit counts.
- **Alpha diversity** — richness, Shannon, Simpson, Pielou evenness, and
  Berger-Parker dominance, with optional significance testing and bracket
  annotations (`compute_alpha()`, `compute_alpha_significance()`).
- **Beta diversity** — pairwise distances, PCoA, CAP (constrained ordination),
  PERMANOVA, dispersion tests, and t-SNE (`compute_distance()`,
  `compute_pcoa()`, `compute_capscale()`, `compute_permanova()`,
  `compute_dispersion()`, `compute_tsne()`).
- **Global prevalence shift (delta)** — subject-level permutation tests for
  antigen-wide prevalence shifts between groups, with forest plots and ECDF
  visualisations (`compute_delta()`, `forestplot()`, `ecdf_plot()`).
- **Per-feature prevalence comparison (POP)** — Fisher's exact and McNemar
  tests for every peptide or taxonomic rank, visualised as scatter plots and
  volcano plots (`compute_pop()`, `scatter_static()`, `volcano_static()`).

### The `compute_` / `plot_` pattern

Every domain follows the same two-step workflow: a `compute_*()` function
returns a tidy data structure, which is then passed to one or more `plot_*()`
functions. This keeps computation and visualisation cleanly separated and makes
it easy to inspect, filter, or export results before plotting.

### Static and interactive outputs

All visualisation functions come in two flavours: a static version built on
**ggplot2** (composable with the full tidyverse graphics ecosystem) and an
interactive version built on **plotly** (useful for exploratory data analysis
and HTML reports). Both use the built-in `theme_phip()` and phiper colour
scales for a consistent visual style.

## Documentation, examples, and vignettes

Full documentation, worked examples, and vignettes are available on the
package website:

**<https://polymerase3.github.io/phiper/>**

A hands-on **[Typical workflow tutorial](https://polymerase3.github.io/phiper/articles/typical-workflow.html)**
walks through a full cross-sectional analysis — loading data, alpha/beta diversity,
POP, and DELTA — using a bundled dummy dataset.

## Contributing

Bug reports, feature requests, and questions are welcome. Please open an issue
on the GitHub tracker:

- <https://github.com/Polymerase3/phiper/issues>

If you'd like to contribute, please read [CONTRIBUTING.md](CONTRIBUTING.md).
