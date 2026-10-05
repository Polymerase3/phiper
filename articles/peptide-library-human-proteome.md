# Human Proteome Library Metadata

## Overview

This vignette describes the **human proteome** PhIP-seq peptide library
(300,000 `humanProteome_*` peptides). The library can be retrieved via
`phiperio`’s
[`get_peptide_library()`](https://polymerase3.github.io/phiperio/reference/get_peptide_library.html):

``` r

library(phiperio)
library(phiper)
library(dplyr)

peplib <- get_peptide_library("human_proteome") %>%
  collect()
#> [13:04:43] INFO  Retrieving peptide metadata into DuckDB cache
#>                  -> get_peptide_library(library = human_proteome, force_refresh
#>                     = FALSE)
#> [13:04:43] INFO  Opened DuckDB connection
#>                    - cache dir:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/phip_cache.duckdb
#>                    - tables: peptide_meta_human_proteome
#> [13:04:43] INFO  Starting download
#>                    - dest:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/human_proteome_library_16.09.26.rds
#> [13:04:43] OK    Download succeeded (method = <getOption()>)
#> [13:04:43] OK    Checksum verified (SHA-256 match)
#> [13:04:48] OK    Download complete and loaded into R
#> [13:04:53] INFO  Importing sanitized metadata into DuckDB cache...
#> [13:04:55] OK    peptide_meta_human_proteome table created in DuckDB cache
#> [13:04:55] OK    Retrieving peptide metadata into DuckDB cache - done
#>                  -> elapsed: 12.6s
```

## Metadata fields

Detailed documentation of the metadata fields is coming soon.

## Credits

The human proteome library was curated by **Nicolai Hörstke**.
