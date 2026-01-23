library(phiper)
library(readxl)
library(dplyr)
library(purrr)
library(readr)

# Local file paths for the Sasha annotation export.
input_path <- "~/Downloads/combined_antigen_data_LLNEXT_with_extra_food_bacteriophages_SZ.xlsx - combined_antigen_data_LLNEXT_wi.csv"
output_path <- "~/Documents/R/phiper/library-metadata/combined_library_15.01.26.rds"

peplib <- get_peptide_meta() %>% collect()

cols_keep <- c(
  "peptide_id","aa_seq","pos","len_seq","Fullname","Description",
  "is_IEDB_or_cntrl","is_auto","is_infect","is_EBV","is_toxin","is_PNP",
  "is_EM","is_MPA","is_patho","is_probio","is_IgA","is_flagellum",
  "signalp6_slow","is_topgraph_new","is_allergens","domain","kingdom",
  "phylum","class","order","family","genus","species","common"
)

peplib <- peplib %>%
  dplyr::select(dplyr::any_of(cols_keep))

sasha_peplib <- read_csv(
  input_path,
  col_types = cols(
    SZ_fly_worm_bee_cocroach_mite_mosquito = col_character(),
    SZ_fungi = col_character(),
    SZ_fungi_type = col_character(),
    SZ_food = col_character(),
    SZ_food_subcaterogy = col_character(),
    SZ_food_item = col_character(),
    SZ_HomoSapi = col_character(),
    SZ_Lacto_phage = col_character(),
    SZ_comments = col_character()
  ),
  skip = 1
)

# Keep the annotation block, normalize column names, and clean ad-hoc values.
sasha_peplib <- sasha_peplib %>%
  select(antigen, SZ_other:SZ_comments) %>%
  rename_with(~ sub("^SZ", "anno", .x), starts_with("SZ")) %>%
  mutate(anno_other = na_if(anno_other, "other")) %>%
  mutate(anno_other = na_if(anno_other, "others"))

# Convert known label columns to logical flags.
label_map <- c(
  anno_HomoSapi = "yes",
  anno_fungi = "yes",
  anno_food = "food",
  anno_Lacto_phage = "Lacto_phage"
)

sasha_peplib <- sasha_peplib %>%
  mutate(across(all_of(names(label_map)), ~ .x == label_map[cur_column()])) %>%
  select(-anno_comments) %>%
  rename(
    anno_food_subcategory = anno_food_subcaterogy,
    anno_is_homo_sapiens = anno_HomoSapi,
    anno_viruses_bacteriophage = anno_Viruses_bacteriophage,
    anno_is_fungi = anno_fungi,
    anno_is_food = anno_food,
    anno_is_lacto_phage = anno_Lacto_phage
  )

# quick profiling
cols <- sasha_peplib %>%
  select(antigen, starts_with("anno_")) %>%
  names()

counts_list <- set_names(cols) %>%
  map(~ sasha_peplib %>% count(.data[[.x]], sort = TRUE, name = "n"))
counts_list

peplib2 <- peplib %>%
  left_join(
    sasha_peplib,
    by = c("peptide_id" = "antigen")
  )

saveRDS(as.data.frame(peplib2), output_path)
