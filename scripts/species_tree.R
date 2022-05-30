library(taxize)
library(ggtree)
library(ape)

sps <- unique(divspe$species)

sps_namesupd <- gnr_resolve(names = sps, canonical = T, fields = "all", 
                            data_source_ids = "4", best_match_only = T)
sps_classif <- classification(sps_namesupd$taxon_id %>% unique(), db = 'ncbi')
fungi
sps_order <- class2tree(classif)$phylo %>%
  fortify() %>%
  filter(isTip) %>%
  arrange(y) %>%
  select(label) %>%
  flatten_chr()


sps_phylo <- ggtree(sps_tree$phylo)
sps_phylo$labels <- sps_

names_upd %>%
  select(user_supplied_name, submitted_name, matched_name, score) %>%
  unique()
