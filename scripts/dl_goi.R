# Script to download info concerning genes of interest

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("InterMineR")

library(InterMineR)

gene_list <- readxl::read_xlsx("data/genes.xlsx")


im <- initInterMine(mine=listMines()["YeastMine"])
im

template = getTemplates(im)
template
template[grep("gene", template$name, ignore.case=TRUE),]


# Query for gene pathways
queryGenePath = getTemplateQuery(
  im = im, 
  name = "Gene_Phenotype_New"
)

# Specify a gene name (here PDR18)
queryGenePath$where = setConstraints(
  modifyQueryConstraints = queryGenePath,
  m.index = 2,
  values = list("PDR18")
)

# Add a new column to the output dataframe
queryGenePath.InterMineR = setQuery(
  inheritQuery = queryGenePath,
  select = c(queryGenePath$select, 
             "Gene.diseases.name")
)
getSelect(queryGenePath.InterMineR)

# Run the query
resGenePath.InterMineR <- runQuery(im, queryGenePath.InterMineR)
resGenePath <- runQuery(im, queryGenePath)
head(resGenePath)




