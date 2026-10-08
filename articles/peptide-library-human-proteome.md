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
#> [08:52:57] INFO  Retrieving peptide metadata into DuckDB cache
#>                  -> get_peptide_library(library = human_proteome, force_refresh
#>                     = FALSE)
#> [08:52:57] INFO  Opened DuckDB connection
#>                    - cache dir:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/phip_cache.duckdb
#>                    - tables: peptide_meta_human_proteome
#> [08:52:57] INFO  Starting download
#>                    - dest:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/human_proteome_library_16.09.26.rds
#> [08:52:57] OK    Download succeeded (method = <getOption()>)
#> [08:52:57] OK    Checksum verified (SHA-256 match)
#> [08:53:02] OK    Download complete and loaded into R
#> [08:53:06] INFO  Importing sanitized metadata into DuckDB cache...
#> [08:53:08] OK    peptide_meta_human_proteome table created in DuckDB cache
#> [08:53:08] OK    Retrieving peptide metadata into DuckDB cache - done
#>                  -> elapsed: 11.733s
```

All columns beginning with `is_` are logical flags: `TRUE` indicates
that the corresponding annotation applies to the peptide, `FALSE`
indicates that it does not.

## General library information

The following columns describe peptide-level and source-protein
information.

| Column | Description |
|----|----|
| `peptide_id` | Unique peptide identifier in the human proteome library. Final IDs use the format `humanProteome_<number>`. |
| `aa_seq` | Amino-acid sequence of the encoded peptide/oligo. |
| `barcode_0` | Nucleotide barcode sequence of 100 nt length. Used to identify the enriched peptide after sequencing. |
| `pos` | 0-based start position of the selected peptide mapping within the full original source-protein sequence. When multiple valid mappings existed, the first retained position is used in the final table. |
| `origin` | Original source protein and position of the peptide in the original protein. |
| `mapped` | All proteins and positions the peptide was mapped to in the complete library. The positions are filtered to only include positions that are in-frame of the peptide tiling. |
| `len_seq` | Length, in amino acids, of the full original source-protein sequence. |
| `Fullname` | Name and description of the original source protein. This includes additional information based on the source protein category, for example UniProt IDs, HLA nomenclature or mutation description for neoantigens. |
| `full_aa_seq` | Full amino-acid sequence of the original source protein. Not included in the library served by `phiperio` because of large file size. |
| `Description` | UniRef annotation of the source protein. |

## Taxonomic information

The taxonomic information is mainly of importance for the control
sequences, as the majority of the library is comprised of human
sequences.

| Column | Description |
|----|----|
| `domain` | Highest-level taxonomic assignment used in the metadata, e.g. Bacteria, Archaea, Eukaryota, or Viruses. |
| `kingdom` | Taxonomic kingdom, when available in the NCBI lineage. |
| `phylum` | Taxonomic phylum. |
| `class` | Taxonomic class. |
| `order` | Taxonomic order. |
| `family` | Taxonomic family. |
| `genus` | Taxonomic genus. |
| `species` | Taxonomic species. |
| `common` | Common name of the source organism. |

## Library annotation flags

The following fields derive from annotation labels carried through the
creation of the human proteome library. They can be used to filter
peptides by category of the source protein.

| Column | Description |
|----|----|
| `is_proteome` | One of the source sequences this peptide was mapped to is part of the UniProt-sourced proteome, including isoforms. |
| `is_mitochond_protein` | One of the source sequences this peptide was mapped to was annotated as a mitochondrial protein using the MitoProteome Database build from 18.01.2022. |
| `is_HLA` | One of the source sequences this peptide was mapped to is an HLA allele sequence sourced from the IPD-IMGT/HLA database version 3.56. |
| `is_HLA_eplet` | One of the source sequences this peptide was mapped to is an HLA eplet sequence received in April 2025 from Konstantin Doberer and Sebastian Kapps of the Vienna Transplant and Complement Lab. |
| `is_control` | One of the source sequences this peptide was mapped to is a control sequence. The controls include all control protein sequences of the PhIP-Seq library described in Vogl et al., 2021 (<https://doi.org/10.1038/s41591-021-01409-3>) and all peptides that were 5-100% prevalent in previous PhIP-Seq studies (<https://doi.org/10.1038/s41591-021-01409-3>, <https://doi.org/10.1016/j.immuni.2023.04.003>, <https://doi.org/10.1126/sciimmunol.abe9950>). |
| `is_neoantigen` | One of the source sequences this peptide was mapped to is a cancer neoantigen sourced from the TCGA in November 2024. |
| `is_surface_neoantigen` | One of the source sequences this peptide was mapped to is one of the 6,500 most frequent cancer surface neoantigens sourced from the TCGA in November 2024. |
| `is_common_neoantigen_in_selected_cancers` | One of the source sequences this peptide was mapped to is a cancer neoantigen prevalent in more than 1% of bladder, breast, colorectal, kidney, liver and lung cancer as well as melanoma cases. Sourced from the TCGA in November 2024. |
| `is_cryptic_peptide` | One of the source sequences this peptide was mapped to is a cryptic peptide sequence received in January 2025 from Andreas Schlosser of the Rudolf-Virchow-Zentrum - Center for Integrative and Translational Bioimaging in Würzburg. |
| `is_stable_transpos_ORF` | One of the source sequences this peptide was mapped to is encoded by a stable transposable element ORF. The sequences were sourced from Arribas et al. (<https://doi.org/10.1016/j.cell.2024.11.011>). |
| `is_therapeutic` | One of the source sequences this peptide was mapped to is the variable fragment of a therapeutic antibody or the TNF alpha inhibitor Etanercept. The sequences were sourced in October 2024 from the Thera-SAbDab (<https://doi.org/10.1093/nar/gkz827>). |

## Credits

The human proteome library was curated by **Nicolai Hörstke**, who also
wrote this documentation of its metadata fields.
