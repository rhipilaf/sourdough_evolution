
# These SraRunInfo-*.csv files were downloaded from the NCBI after searching for 
# SRA study ids ERP119191, ERP014555, SRP109074, SRP010663 in the SRA browser 
# and clicking on Send to > File > RunInfo.

sra <- list()
sra$bigey <- readr::read_csv("data/other_strains/SraRunInfo-Bigey.csv")
sra$coreColl <- readr::read_csv("data/other_strains/SraRunInfo-CoreCollection.csv")
sra$gallone <- readr::read_csv("data/other_strains/SraRunInfo-Gallone.csv")
sra$dunn <- readr::read_csv("data/other_strains/SraRunInfo-Dunn.csv")
sra$zhu1 <- readr::read_csv("data/other_strains/SraRunInfo-Zhu1.csv")
sra$zhu2 <- readr::read_csv("data/other_strains/SraRunInfo-Zhu2.csv")


# Just below , filters for each dataset if needed. For Bigey et al. data, we include 
# everything, so, no filter.

fltr <- list()
fltr$coreColl <- readr::read_lines("data/other_strains/coreCollection_list.txt")
fltr$dunn <- c("MoroccoBreadG17","SardSourdoughS11")

col_names <- c("strain_name", "habitat", "to_dl", "ncbi_bioproject", 
               "ncbi_biosample", "ncbi_sra")

other_strains <- bind_rows(sra, .id = "project") %>%
  select(where(function(x) any(!is.na(x)))) %>% # Columns that are not empty
  bind_cols(as_tibble(matrix(ncol = length(col_names), 
                             nrow = nrow(.), 
                             dimnames = list(NULL, col_names))), .) %>%
  mutate(strain_name = case_when(project == "coreColl" ~ str_extract(LibraryName, 
                                        '(?<=^BCM_)[A-Z]{3}'),
                                 TRUE ~ LibraryName),
    to_dl = case_when(!is.na(download_path) ~ 1,
                      TRUE ~ 0)) %>%
  filter(project == "bigey" 
         | (project == "coreColl" & strain_name %in% fltr$coreColl)
         | (project == "gallone" & grepl("bread", strain_name))
         | (project == "dunn" & strain_name %in% fltr$dunn)
         | project == "zhu1"
         | project == "zhu2")


write.table(other_strains, file='data/other_strains/strains_RunInfo.tsv', 
            quote=FALSE, sep='\t', row.names = FALSE)
xlsx::write.xlsx(x = as.data.frame(other_strains), file='data/other_strains/strains_RunInfo.xlsx', showNA = F, row.names = F)
