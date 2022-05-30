# Creation of a strains specific table

data_cyto <- readxl::read_xlsx("data/data_robot/data_cyto.xlsx")
bak <- readxl::read_xlsx("data/_cerevisiae_phenotypage_bak.xlsx")

# > integrate the metadata file data to the strains.csv table (before see with Delphine its structure to see what is unique)
metadata <- readxl::read_xlsx("data/_metadata_saccharomyces_cerevisiae.xlsx") %>%
  select(c("strain_name", "strain", "ploidy", "phenotyping", "genome_seq", 
           "plaque", "puit", "genus", "DB_name", "score", "identity")) %>%
  unique()

strains <- data_cyto %>%
  select(strain_name, species, isolation, plaque, puit, baker_id, backslopping, genome_seq, flour_id) %>%
  unique

strains2 <- bak %>%
  mutate(strain_name = ifelse(is.na(strain_name...5), strain_name...17, strain_name...5)) %>%
  select(strain_name, species, isolation, plaque, puit, baker_id, collection_date, backslopping, to_sequence, séquencage, project, flour_id) %>%
  unique

strains_tot <- full_join(strains, strains2, by = "strain_name") %>% relocate(., names(.)[order(names(.))]) %>%
  mutate(backslopping = ifelse(is.na(backslopping.x), backslopping.y, backslopping.x),
         baker_id = ifelse(is.na(baker_id.x), baker_id.y, baker_id.x),
         isolation = ifelse(is.na(isolation.x), isolation.y, isolation.x),
         plaque = ifelse(is.na(plaque.x), plaque.y, plaque.x),
         puit = ifelse(is.na(puit.x), puit.y, puit.x),
         species = "Saccharomyces cerevisiae",
         flour_id = flour_id.x,
         collection_date = if_else(nchar(collection_date) == 4, 
                                   as.Date(as.character(collection_date), format = "%Y"),
                                   as.Date(collection_date, origin = "1899-12-30"))) %>%
  select(-ends_with(".y"), -ends_with(".x")) %>%
  relocate(c("strain_name","species","baker_id","flour_id","backslopping",
             "isolation","collection_date","plaque","puit", "project")) %>%
  arrange(strain_name)

strains3 <- full_join(strains_tot, metadata, by = "strain_name")

write.csv(file = "data/strains.csv", row.names = F, na = "")
