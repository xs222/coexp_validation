library(ggplot2)
library(dplyr)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")

# save count matrix
# sc_obj <- readRDS("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/MIT_ROSMAP_Multiomics/Gene_Expression/snRNAseq/snRNAseq-10x/processed/Oligodendrocytes.rds")
# clic = read.csv("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/ROSMAP/Metadata/ROSMAP_clinical.csv")
# sc_obj_meta = sc_obj@meta.data
# sc_obj_meta = left_join(sc_obj_meta, clic)
# sc_obj$braaksc <- sc_obj_meta$braaksc
# count_ROSMAP_oli <- as.matrix(sc_obj[["RNA"]]@counts[, which(sc_obj$braaksc==0)])
# saveRDS(count_ROSMAP_oli, "marginal_fit/ROSMAP_NC_Oli_count.rds")
count_ROSMAP_oli <- readRDS("marginal_fit/ROSMAP_NC_Oli_count.rds")

# sc_obj <- readRDS("/gpfs/gibbs/pi/zhao/cs2629/AD_Nancy/seurat_obj_cell_type_labelled.rds")
# table(sc_obj@meta.data$cell_type)
# sc_obj$disease = sapply(sc_obj$orig.ident, function(or) {substr(strsplit(or, '_')[[1]][2], 1, 2)})
# table(sc_obj$disease)
# count_PNAS_oli <- as.matrix(sc_obj[["RNA"]]@counts[, which(sc_obj$cell_type=="Oli" & sc_obj$disease=="NC")])
# saveRDS(count_PNAS_oli, "marginal_fit/PNAS_NC_Oli_count.rds")
count_PNAS_oli <- readRDS("marginal_fit/PNAS_NC_Oli_count.rds")

# estimate correlation structure
marginal_fit_ROSMAP = readRDS('marginal_fit/ROSMAP_NC_Oli_marginal_fit.rds')
marginal_fit_PNAS = readRDS('marginal_fit/PNAS_NC_Oli_marginal_fit.rds')
marginal_join = inner_join(marginal_fit_ROSMAP, marginal_fit_PNAS, by="gene")

## 1000 gene cor mat
ngene <- 1000
marginal_join <- marginal_join[order(marginal_join$mu.x, decreasing = T),]
gene_name <- marginal_join$gene[1:ngene]

set.seed(11122023)
dim(count_ROSMAP_oli)

source("/gpfs/gibbs/pi/zhao/xs282/validation/coexp_function.R")
sc_obj <- CreateSeuratObject(counts = count_ROSMAP_oli)
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc.sel <- subset(sc_obj, features = gene_name)
ROSMAP_sct_prn <- sct_cor(sc_obj, sc.sel, gene_name)
saveRDS(ROSMAP_sct_prn, "mean_cor/semi_PD/simu/ROSMAP_NC_Oli_sct1000.rds")

sc_obj <- CreateSeuratObject(counts = count_PNAS_oli)
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc.sel <- subset(sc_obj, features = gene_name)
PNAS_sct_prn <- sct_cor(sc_obj, sc.sel, gene_name)
saveRDS(PNAS_sct_prn, "mean_cor/semi_PD/simu/PNAS_NC_Oli_sct1000.rds")
