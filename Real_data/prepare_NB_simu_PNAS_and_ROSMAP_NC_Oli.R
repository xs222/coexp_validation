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

source("compare_simulation/NB_copula/NB_copula_function.R")
source("/gpfs/gibbs/pi/zhao/xs282/coexp-sc/IRLS_CSCORE/CscoreSimplifiedIRLS.R")
source("/gpfs/gibbs/pi/zhao/xs282/validation/cscore_real_data_function.R")
set.seed(11122023)
dim(count_ROSMAP_oli)
ROSMAP_ests <- CscoreSimplifiedIRLS(count_ROSMAP_oli[gene_name,] %>% as.matrix %>% t,
                                    colSums(count_ROSMAP_oli), covar_weight="regularized")
ROSMAP_ests$est <- post_process_est(ROSMAP_ests$est)
saveRDS(ROSMAP_ests, "marginal_fit/ROSMAP_NC_Oli_cscore_cor1000.rds")
filtered_ROSMAP_ests <- ROSMAP_ests$est
filtered_ROSMAP_ests[MatrixBH(ROSMAP_ests$p_value) > 0.05] <- 0

dim(count_PNAS_oli)
PNAS_ests <- CscoreSimplifiedIRLS(count_PNAS_oli[gene_name,] %>% as.matrix %>% t,
                                  colSums(count_PNAS_oli), covar_weight="regularized")
PNAS_ests$est <- post_process_est(PNAS_ests$est)
saveRDS(PNAS_ests, "marginal_fit/PNAS_NC_Oli_cscore_cor1000.rds")
filtered_PNAS_ests <- PNAS_ests$est
filtered_PNAS_ests[MatrixBH(PNAS_ests$p_value) > 0.05] <- 0



