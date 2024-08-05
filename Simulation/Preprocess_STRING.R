## Preprocess
## Mapped to gene
## Remove dupicated and loop pairs

## STRING data (Homo Sapiens) is downloaded on 8/30/2023 from: https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz
## Introduction of scores in STRING: http://version10.string-db.org/help/getting_started/
##                                   http://version10.string-db.org/help/faq/
## * explore basic statistics of the STRING database
## * check the relationship between confidence in STRING and the gene pair expression level
## * compare the connectivity of confident pairs with gene pair expression level
## p-values for gene coexpression network

setwd("/gpfs/gibbs/pi/zhao/xs282/validation/STRING/")
set.seed(8162022)
library(tidyverse)
library(magrittr)
library(rtracklayer)
library(biomaRt)
library(ggvenn)
library(ggplot2)
library(igraph)

hs <- read.table("9606.protein.links.detailed.v12.0.txt.gz")
colnames(hs) <- hs[1,]
hs <- hs[-1,]

pre_explore_fun <- function(x){
  hist(x)
  print(paste0("Number of zero: ", sum(x==0)," Prop of zero: ", sum(x==0)/length(x)))
  summary(x)
}

hs$combined_score <- as.numeric(hs$combined_score)
pre_explore_fun(hs$combined_score)

hs$neighborhood <- as.numeric(hs$neighborhood)
pre_explore_fun(hs$neighborhood)

hs$fusion <- as.numeric(hs$fusion)
pre_explore_fun(hs$fusion)

hs$cooccurence <- as.numeric(hs$cooccurence)
pre_explore_fun(hs$cooccurence)

hs$coexpression <- as.numeric(hs$coexpression)
pre_explore_fun(hs$coexpression)

hs$experimental <- as.numeric(hs$experimental)
pre_explore_fun(hs$experimental)

hs$database <- as.numeric(hs$database)
pre_explore_fun(hs$database)

hs$textmining <- as.numeric(hs$textmining)
pre_explore_fun(hs$textmining)

hs$protein1 <- gsub("9606\\.", "", hs$protein1)
hs$protein2 <- gsub("9606\\.", "", hs$protein2)

#### string contains id1 id2 and id2 id1, exactly half
# map Ensembl protein to genes (codes from: https://github.com/skinnider/SCT-MoA/blob/master/R/networks/preprocess-string.R)
gtf <- as.data.frame(import("Homo_sapiens.GRCh38.110.gtf.gz"))
g2p <- gtf %>% dplyr::filter(type == 'CDS') %>%
  dplyr::select(gene_id, protein_id) %>%
  group_by(gene_id, protein_id) %>%
  dplyr::slice(1) %>%
  ungroup()
map <- with(g2p, setNames(gene_id, protein_id))
hs$protein1 = map[hs$protein1]
hs$protein2 = map[hs$protein2]

hs_filter <- hs[!(is.na(hs$protein1)|is.na(hs$protein2)),]
# identify and remove self-connected pairs
idx <- which(hs_filter$protein1==hs_filter$protein2)
if (length(idx)>0){hs_filter <- hs_filter[-idx,]}
# identify and remove repeated pairs but with different order
hs_filter <- hs_filter %>%
  group_by(grp = paste(pmax(protein1, protein2), pmin(protein1, protein2), sep = "_")) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(-grp)

saveRDS(hs_filter, "hs_filter_10_4_2023.rds")


