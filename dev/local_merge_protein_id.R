library(readr)

input_lib  <- "~/Documents/R/phiper/library-metadata/combined_library_15.01.26.rds"
input_map  <- "~/Documents/R/phiper/library-metadata/peptide_to_protein_map.csv"
output_lib <- "~/Documents/R/phiper/library-metadata/combined_library_06.07.26.rds"

peplib <- readRDS(input_lib)
prot_map <- read_delim(input_map, delim = ";", col_types = cols(.default = col_character()))

peplib2 <- merge(peplib, prot_map, by = "peptide_id", all.x = TRUE)

saveRDS(peplib2, output_lib)
