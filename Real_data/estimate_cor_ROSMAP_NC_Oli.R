library(ggpubr)
library(SingleCellExperiment)
library(jsonlite)
library(optparse)
library(DESeq2)
library(tidyr)
library(parallel)
set.seed(10222023)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("compare_simulation/NB_copula/NB_copula_function.R")
source("/gpfs/gibbs/pi/zhao/xs282/coexp-sc/IRLS_CSCORE/CscoreSimplifiedIRLS.R")
source("/gpfs/gibbs/pi/zhao/xs282/validation/cscore_real_data_function.R")
source("/gpfs/gibbs/pi/zhao/xs282/validation/coexp_function.R")
seed <- 10222023

# ROSMAP-------------------------------------------------------------------
ROSMAP_oli_ct <- readRDS("marginal_fit/ROSMAP_NC_Oli_count.rds")
ori_ests <- readRDS("marginal_fit/ROSMAP_NC_Oli_cscore_cor1000.rds")
cor_gene_name <- colnames(ori_ests$est)

extract_upp_tri <- function(data, gene_name){
  data <- data[gene_name, gene_name]
  return(data[upper.tri(data, diag = FALSE)])
}

est_mat_ROSMAP <- as.data.frame(matrix(NA,ncol=9,nrow=length(cor_gene_name)*(length(cor_gene_name)-1)/2))

# sct
sc_obj <- CreateSeuratObject(counts = ROSMAP_oli_ct)
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc.sel <- subset(sc_obj, features = cor_gene_name)
ROSMAP_sct_prn <- sct_cor(sc_obj, sc.sel, cor_gene_name)
est_mat_ROSMAP$sct <- extract_upp_tri(ROSMAP_sct_prn, cor_gene_name)

# noise regularization
path2 <- paste0("real/oli/noise/")
if(!file.exists(path2)){
  dir.create(path2,recursive = T)
}
noise_cor <- noise_fun(sc_obj, sc.sel, sel.gene=cor_gene_name,
                       seed=seed, path2 = path2)
est_mat_ROSMAP$noise <- extract_upp_tri(noise_cor, cor_gene_name)

# cscore
ROSMAP_cscore <- CscoreSimplifiedIRLS(ROSMAP_oli_ct[cor_gene_name, ] %>% as.matrix %>% t,
                                      colSums(ROSMAP_oli_ct), covar_weight="regularized")
ROSMAP_cscore$est <- post_process_est(ROSMAP_cscore$est)
est_mat_ROSMAP$cscore_p <- extract_upp_tri(ROSMAP_cscore$p_value, cor_gene_name)
est_mat_ROSMAP$cscore_est <- extract_upp_tri(ROSMAP_cscore$est, cor_gene_name)
est_mat_ROSMAP$cscore_stat <- extract_upp_tri(ROSMAP_cscore$test_stat, cor_gene_name)


## analytic pearson
ROSMAP_ana_prn <- ana_prn(ROSMAP_oli_ct, cor_gene_name, colSums(ROSMAP_oli_ct))
est_mat_ROSMAP$ana_prn <- extract_upp_tri(ROSMAP_ana_prn, cor_gene_name)


## propr
pr <- propr(counts = t(ROSMAP_oli_ct), metric = "rho", select =cor_gene_name, alpha = NA,p = 100)
est_mat_ROSMAP$propr <- extract_upp_tri(pr@matrix, cor_gene_name)


## pearson
norm.data <- as.matrix(GetAssayData(sc.sel, assay = "RNA", slot = "data"))
cor_m_pearson <- cor(t(norm.data),method = "pearson")
est_mat_ROSMAP$prn <- extract_upp_tri(cor_m_pearson, cor_gene_name)
## spearman
cor_m_spr <- cor(t(norm.data),method = "spearman")
est_mat_ROSMAP$spr <- extract_upp_tri(cor_m_spr, cor_gene_name)

saveRDS(est_mat_ROSMAP, paste0("real/oli/", "/est_cor", ".rds"))
