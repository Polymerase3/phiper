# ICAM Library Metadata

## Overview

This vignette describes the metadata fields provided with the **ICAM**
PhIP-seq peptide library. The metadata combine peptide-level
information, source-protein sequences, taxonomic annotations,
library-construction labels, and information about fused protein
fragments. The library can be retrieved via `phiperio`’s
[`get_peptide_library()`](https://polymerase3.github.io/phiperio/reference/get_peptide_library.html):

``` r

library(phiperio)
library(phiper)
library(dplyr)

peplib <- get_peptide_library("icam") %>%
  collect()
#> [14:27:39] INFO  Retrieving peptide metadata into DuckDB cache
#>                  -> get_peptide_library(library = icam, force_refresh = FALSE)
#> [14:27:39] INFO  Opened DuckDB connection
#>                    - cache dir:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/phip_cache.duckdb
#>                    - tables: peptide_meta_icam
#> [14:27:40] INFO  Starting download
#>                    - dest:
#>                      /home/runner/.cache/R/phiperio/peptide_meta/icam_library_01.10.26.rds
#> [14:27:41] OK    Download succeeded (method = <getOption()>)
#> [14:27:41] OK    Checksum verified (SHA-256 match)
#> [14:27:45] OK    Download complete and loaded into R
#> [14:27:50] INFO  Importing sanitized metadata into DuckDB cache...
#> [14:27:52] OK    peptide_meta_icam table created in DuckDB cache
#> [14:27:52] OK    Retrieving peptide metadata into DuckDB cache - done
#>                  -> elapsed: 13.071s
```

All columns beginning with `is_` are binary dummy variables: `1`
indicates that the corresponding annotation applies to the
peptide/protein, `0` indicates that it does not.

## General library information

The following columns describe peptide-level and source-protein
information.

| Column | Description |
|----|----|
| `peptide_id` | Unique peptide identifier in the ICAM library. Final IDs use the format `icam_<number>`. |
| `aa_seq` | Amino-acid sequence of the encoded peptide/oligo. |
| `pos` | 0-based start position of the selected peptide mapping within the full original source-protein sequence. When multiple valid mappings existed, the first retained position is used in the final table. |
| `len_seq` | Length, in amino acids, of the full original source-protein sequence. |
| `full_aa_seq` | Full amino-acid sequence of the original source protein. Not included in the library served by `phiperio` because of large file size. |
| `Description` | Protein annotation/description associated with the peptide, derived from the original PGAP annotation within each proteome downloaded from GenBank or RefSeq. This corresponds to the protein description previously stored as `Fullname`. |

## Taxonomic information

Taxonomic lineage fields were harmonized using NCBI taxonomy where
possible. The normalized `species` field and the original
`species_exact` field are intentionally retained separately.

| Column | Description |
|----|----|
| `tax_id` | NCBI Taxonomy identifier assigned during taxonomic harmonization, when available. |
| `domain` | Highest-level taxonomic assignment used in the metadata, e.g. Bacteria, Archaea, Eukaryota, or Viruses. |
| `kingdom` | Taxonomic kingdom, when available in the NCBI lineage. |
| `phylum` | Taxonomic phylum. |
| `class` | Taxonomic class. |
| `order` | Taxonomic order. |
| `family` | Taxonomic family. |
| `genus` | Taxonomic genus. |
| `species` | Taxonomic species. |
| `species_exact` | Original organism name retained from the source annotation (the strain name found in the GenBank file of the downloaded strain). This may preserve strain, subspecies, isolate, or historical naming information not present in the normalized `species` field. |

## Library annotation flags

The following fields derive from annotation labels carried through the
ICAM library-construction workflow. They can be used to filter peptides
by control set, disease-associated mimotope set, protein function,
predicted function, localization, or curated library component.

| Column | Description |
|----|----|
| `is_Controls_Agilent` | Peptides belonging to the Agilent control set. |
| `is_Controls_Israeli_Common_epitopes` | Peptides belonging to the Israeli common-epitope control set (epitopes found at 5-95% prevalence in an Israeli cohort). |
| `is_Controls_Random` | Random control peptides. |
| `is_EBV_variants` | Epstein-Barr virus (EBV) variant peptides. |
| `is_MS_mimotope` | Peptides annotated as multiple sclerosis (MS) mimotopes. |
| `is_RA_mimotope` | Peptides annotated as rheumatoid arthritis (RA) mimotopes. |
| `is_SARS_CoV_variants` | SARS-CoV variant peptides. |
| `is_SLE_mimotope` | Peptides annotated as systemic lupus erythematosus (SLE) mimotopes. |
| `is_SSc_mimotope_100` | Peptides in the systemic sclerosis (SSc) mimotope set labelled `100`. |
| `is_SSc_mimotope_90` | Peptides in the systemic sclerosis (SSc) mimotope set labelled `90`. |
| `is_adhesin` | Peptides derived from proteins identified as adhesins. |
| `is_aeromonas_aerolysin` | Peptides belonging to the Aeromonas aerolysin annotation set. |
| `is_agilent_flagellins` | Peptides belonging to the Agilent flagellin annotation set. |
| `is_antibiotic_resistance` | Proteins predicted as antibiotic-resistance related using DIAMOND against the curated CARD database. |
| `is_axSpA_mimotope` | Peptides annotated as axial spondyloarthritis (axSpA) mimotopes. |
| `is_crc_proteins` | Proteins included in the colorectal cancer (CRC)-associated protein set used during library construction. |
| `is_flagellin` | Proteins identified as flagellins within the initial PGAP annotations. |
| `is_fusobacterium_orthologues` | Peptides/proteins belonging to the Fusobacterium subspecies orthologue set. |
| `is_immunogenic_proteins` | Proteins included in the immunogenic-protein annotation set. |
| `is_invasion_proteins` | Proteins annotated as invasion-related proteins. |
| `is_ley_flagellins` | Peptides belonging to the Ley flagellin annotation set. |
| `is_membrane` | Proteins annotated as membrane-associated. |
| `is_predicted_antibiotic_resistance` | Proteins computationally predicted to be antibiotic-resistance related. |
| `is_predicted_flagellin` | Proteins predicted as flagellins using DIAMOND against a bacterial flagella database downloaded from UniProt. |
| `is_predicted_toxin` | Proteins computationally predicted to be toxins. |
| `is_secreted` | Proteins annotated or predicted as secreted. |
| `is_ureases` | Proteins belonging to the urease annotation set. |

## Fused-protein information

Some ICAM library entries were constructed by joining multiple protein
fragments with a linker. The following field preserves the original
location of each component fragment.

| Column | Description |
|----|----|
| `origin_fused_chunks` | For peptides derived from fused library proteins, the source chunk origins as a list of `(original_protein_id, start_position)` tuples. Start positions are 0-based coordinates in the full original protein. The field may contain two, three, or more tuples depending on the number of fragments used in the fusion. For non-fused proteins it is not applicable and may be empty/`NA`. |

## Credits

This documentation of the ICAM library metadata fields was written by
**Gabriel Innocenti**.
