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
#> [10:31:20] INFO  Retrieving peptide metadata into DuckDB cache
#>                  -> get_peptide_library(library = human_proteome, force_refresh
#>                     = FALSE)
#> [10:31:21] INFO  Opened DuckDB connection
#>                    - cache dir:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/phip_cache.duckdb
#>                    - tables: peptide_meta_human_proteome
#> [10:31:21] INFO  Starting download
#>                    - dest:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/human_proteome_library_16.09.26.rds
#> [10:31:21] OK    Download succeeded (method = <getOption()>)
#> [10:31:21] OK    Checksum verified (SHA-256 match)
#> [10:31:26] OK    Download complete and loaded into R
#> [10:31:31] INFO  Importing sanitized metadata into DuckDB cache...
#> [10:31:33] OK    peptide_meta_human_proteome table created in DuckDB cache
#> [10:31:33] OK    Retrieving peptide metadata into DuckDB cache - done
#>                  -> elapsed: 12.501s
```

## Metadata fields

Detailed documentation of the metadata fields is coming soon.

## Credits

The human proteome library was curated by **Nicolai Hörstke**.
