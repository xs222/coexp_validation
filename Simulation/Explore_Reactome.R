
# explore the bias in Reactome
# Reactome: downloaded on 10/6/2023, from https://reactome.org/download/current/Ensembl2Reactome.txt
# generate correlated pairs based on the pathway: if they are in the same pathway, they are correlated
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/Reactome/")
set.seed(8162022)
library(tidyverse)
library(magrittr)
library(rtracklayer)
library(biomaRt)
library(ggvenn)
library(ggplot2)
library(igraph)
library(venn)

# https://reactome.org/download-data
reactome <- read_tsv("Ensembl2Reactome.txt", col_names = F)
colnames(reactome) <- c("EnsemblID", "PathwayID", "URL", "PathwayName", "EvidenceCode", "Species")

table(reactome$Species)
table(reactome$EvidenceCode)
# IEA (Inferred from Electronic Annotation): This evidence code is used when information
# is inferred through automated computational methods or algorithms.
# TAS (Traceable Author Statement): TAS is used when information is provided by an
# author or expert and can be traced back to a specific literature source.

reactome_human <- reactome[reactome$Species=="Homo sapiens",]

# current data contains ENSP, ENST
human_enst <- reactome_human[substr(reactome_human$EnsemblID, 1,4)=="ENST",]
human_ensp <- reactome_human[substr(reactome_human$EnsemblID, 1,4)=="ENSP",]
human_ensg <- reactome_human[substr(reactome_human$EnsemblID, 1,4)=="ENSG",]

gtf <- as.data.frame(import("/gpfs/gibbs/pi/zhao/xs282/validation/STRING/Homo_sapiens.GRCh38.110.gtf.gz"))
g2p <- gtf %>% dplyr::filter(type == 'CDS') %>%
  dplyr::select(gene_id, protein_id) %>%
  group_by(gene_id, protein_id) %>%
  dplyr::slice(1) %>%
  ungroup()
map_g2p <- with(g2p, setNames(gene_id, protein_id))
human_ensp$EnsemblID <- map_g2p[human_ensp$EnsemblID]

g2t <- gtf %>% dplyr::filter(type == 'CDS') %>%
  dplyr::select(gene_id, transcript_id) %>%
  group_by(gene_id, transcript_id) %>%
  dplyr::slice(1) %>%
  ungroup()
map_g2t <- with(g2t, setNames(gene_id, transcript_id))
human_enst$EnsemblID <- map_g2t[human_enst$EnsemblID]

human <- rbind(human_enst, human_ensp, human_ensg)
human <- human[1:2]
human$EnsemblID <- gsub("\\..*", "", human$EnsemblID)
human <- human %>% group_by(EnsemblID, PathwayID) %>% dplyr::slice(1) %>% ungroup()
human <- human[!is.na(human$EnsemblID),]

pathwayID <- unique(human$PathwayID)
pathway_list <- list()
for (i in 1:length(pathwayID)){
  if(i%%100==0){print(i)}
  sub_data <- human[human$PathwayID==pathwayID[i],]
  pathway_list[[pathwayID[i]]] <- sub_data$EnsemblID
}
idx <- sapply(pathway_list, function(x){length(x)<2})
sum(idx)
which(idx)

pathway_list_filter <- pathway_list[!idx]
all_pair <- lapply(pathway_list_filter, function(x){t(combn(x,2))})
all_pair_df <- do.call(rbind, all_pair)
all_pair_df <- as.data.frame(all_pair_df) %>% group_by(V1,V2) %>% dplyr::slice(1) %>% ungroup()
sum(grepl("\\.", all_pair_df$V1, perl = TRUE))
saveRDS(all_pair_df, "Processed_Reactome_10_07_2023.rds")
