## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load, warning=FALSE, message=FALSE---------------------------------------
library(InterMineR)
listMines()

## ----humanmine, warning=FALSE, message=FALSE----------------------------------

# load HumaMine
im <- initInterMine(mine=listMines()["HumanMine"])
im


## ----get_queries, warning=FALSE, message=FALSE--------------------------------
# Get template (collection of pre-defined queries)
template = getTemplates(im)
head(template)

## ----gene_templates, warning=FALSE, message=FALSE-----------------------------

# Get gene-related templates
template[grep("gene", template$name, ignore.case=TRUE),]


## ----gene_pathway, warning=FALSE, message=FALSE-------------------------------
# Query for gene pathways
queryGenePath = getTemplateQuery(
  im = im, 
  name = "Gene_Pathway"
)
queryGenePath


## ----run_genepath, warning=FALSE, message=FALSE-------------------------------

resGenePath <- runQuery(im, queryGenePath)
head(resGenePath)


## ----modify_query, warning=FALSE, message=FALSE-------------------------------

# modify directly the value of the first constraint from the list query
queryGenePath$where[[1]][["value"]] <- "ABO"

# or modify the value of the first constraint from the list query with setConstraints
queryGenePath$where = setConstraints(
  modifyQueryConstraints = queryGenePath,
  m.index = 1,
  values = list("ABO")
)

queryGenePath$where


## ----modify_query2, warning=FALSE, message=FALSE------------------------------

resGenePath <- runQuery(im, queryGenePath)
head(resGenePath)


## ----new_constraint, warning=FALSE, message=FALSE-----------------------------
newConstraint <- list(
  path=c("Gene.pathways.name"),
  op=c("="), 
  value=c("ABO blood group biosynthesis"), 
  code=c("B")
)

queryGenePath$where[[2]] <- newConstraint
queryGenePath$where


## ----new_constraint2, warning=FALSE, message=FALSE----------------------------
resGenePath <- runQuery(im, queryGenePath)
resGenePath

## ----add_column3, warning=FALSE, message=FALSE--------------------------------

# use setQuery function which will create an InterMineR-class query
queryGenePath.InterMineR = setQuery(
  inheritQuery = queryGenePath,
  select = c(queryGenePath$select, 
             "Gene.diseases.name")
)

getSelect(queryGenePath.InterMineR)
#queryGenePath.InterMineR@select

# or assign new column directly to the existing list query
queryGenePath$select[[8]] <- "Gene.diseases.name"
queryGenePath$select

# run queries
resGenePath.InterMineR <- runQuery(im, queryGenePath.InterMineR)
resGenePath <- runQuery(im, queryGenePath)

all(resGenePath == resGenePath.InterMineR)

head(resGenePath, 3)


## ----constrant_logic, warning=FALSE, message=FALSE----------------------------

queryGenePath$constraintLogic <- "A and B"
queryGenePath$constraintLogic


## ----constrant_logic2, warning=FALSE, message=FALSE---------------------------

resGenePath <- runQuery(im, queryGenePath)
resGenePath


## ----gene_go, warning=FALSE, message=FALSE------------------------------------

queryGeneGO <- getTemplateQuery(im, "Gene_GO")
queryGeneGO


## ----gene_go2, warning=FALSE, message=FALSE-----------------------------------

queryGeneGO$select <- queryGeneGO$select[2:5]
queryGeneGO$select


## ----gene_go3, warning=FALSE, message=FALSE-----------------------------------
queryGeneGO$where[[1]][["value"]] <- "ABO"
queryGeneGO$where

## ----gene_go4, warning=FALSE, message=FALSE-----------------------------------
resGeneGO <- runQuery(im, queryGeneGO )
head(resGeneGO)

## ----metal_ion_binding, warning=FALSE, message=FALSE--------------------------
queryGOGene <- getTemplateQuery(im, "GOterm_Gene")
queryGOGene

## ----metal_ion_binding2, warning=FALSE, message=FALSE-------------------------
queryGOGene$select <- queryGOGene$select[2:5]
queryGOGene$select

## ----metal_ion_binding3, warning=FALSE, message=FALSE-------------------------

queryGOGene$where[[1]]$value = "metal ion binding"
queryGOGene$where


## ----metal_ion_binding4, warning=FALSE, message=FALSE-------------------------

resGOGene <- runQuery(im, queryGOGene )
head(resGOGene)


## ----neighbor_genes, warning=FALSE, message=FALSE-----------------------------

queryGeneLoc = getTemplateQuery(im, "Gene_Location")
queryGeneLoc$where[[2]][["value"]] = "ABCA6"
resGeneLoc= runQuery(im, queryGeneLoc)

resGeneLoc


## ----neighbor_genes2, warning=FALSE, message=FALSE----------------------------

# set constraints
constraints = setConstraints(
  paths = c(
    "Gene.chromosome.primaryIdentifier",
    "Gene.locations.start",
    "Gene.locations.end",
    "Gene.organism.name"
  ),
  operators = c(
    "=",
    ">=",
    "<=",
    "="
  ),
  values = list(
    resGeneLoc[1, "Gene.chromosome.primaryIdentifier"],
    as.character(as.numeric(resGeneLoc[1, "Gene.locations.start"])-50000),
    as.character(as.numeric(resGeneLoc[1, "Gene.locations.end"])+50000),
    "Homo sapiens"
  )
)

# set InterMineR-class query
queryNeighborGene = setQuery(
  select = c("Gene.primaryIdentifier", 
             "Gene.symbol",
             "Gene.chromosome.primaryIdentifier",
             "Gene.locations.start", 
             "Gene.locations.end", 
             "Gene.locations.strand"),
  where = constraints
)

summary(queryNeighborGene)


## ----neighbor_genes6, warning=FALSE, message=FALSE----------------------------

resNeighborGene <- runQuery(im, queryNeighborGene)
resNeighborGene


## ----neighbor_genes7, warning=FALSE, message=FALSE----------------------------

resNeighborGene$Gene.locations.strand[which(resNeighborGene$Gene.locations.strand==1)]="+"

resNeighborGene$Gene.locations.strand[which(resNeighborGene$Gene.locations.strand==-1)]="-"

gene.idx = which(nchar(resNeighborGene$Gene.symbol)==0)

resNeighborGene$Gene.symbol[gene.idx]=resNeighborGene$Gene.primaryIdentifier[gene.idx]

## ----load_gviz, warning=FALSE, message=FALSE----------------------------------
require(Gviz)

## ----plottracks, warning=FALSE, message=FALSE---------------------------------
annTrack = AnnotationTrack(
  start=resNeighborGene$Gene.locations.start,
  end=resNeighborGene$Gene.locations.end,
  strand=resNeighborGene$Gene.locations.strand,
  chromosome=resNeighborGene$Gene.chromosome.primaryIdentifier[1],
  genome="GRCh38", 
  name="around ABCA6",
  id=resNeighborGene$Gene.symbol)

gtr <- GenomeAxisTrack()
itr <- IdeogramTrack(genome="hg38", chromosome="chr17")

plotTracks(list(gtr, itr, annTrack), shape="box", showFeatureId=TRUE, fontcolor="black")


## ----sessioInfo---------------------------------------------------------------

sessionInfo()


