library(ggpubr)
library(SingleCellExperiment)
library(jsonlite)
library(optparse)
library(DESeq2)
library(tidyr)
library(propr)
library(parallel)
set.seed(10222023)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("/gpfs/gibbs/pi/zhao/xs282/coexp-sc/IRLS_CSCORE/CscoreSimplifiedIRLS.R")
source("/gpfs/gibbs/pi/zhao/xs282/validation/cscore_real_data_function.R")


# ROSMAP-------------------------------------------------------------------
ROSMAP_oli_ct <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_sct_cor_NB_simu1000_abs_thresh.rds")
ROSMAP_ori_ests <- readRDS("mean_cor/semi_PD/simu/ROSMAP_NC_Oli_sct1000.rds")
gene_name <- rownames(ROSMAP_ori_ests)

# cscore
ROSMAP_cscore <- CscoreSimplifiedIRLS(ROSMAP_oli_ct[gene_name, ] %>% as.matrix %>% t,
                                 colSums(ROSMAP_oli_ct), covar_weight="regularized")
ROSMAP_cscore$est <- post_process_est(ROSMAP_cscore$est)
saveRDS(ROSMAP_cscore, "mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_cscore1000_abs_thresh.rds")


# sctransfrom
source("/gpfs/gibbs/pi/zhao/xs282/validation/coexp_function.R")
sc_obj <- CreateSeuratObject(counts = ROSMAP_oli_ct)
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc.sel <- subset(sc_obj, features = gene_name)
ROSMAP_sct_prn <- sct_cor(sc_obj, sc.sel, gene_name)
saveRDS(ROSMAP_sct_prn, "mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_sct1000_abs_thresh.rds")


## analytic pearson
ROSMAP_ana_prn <- ana_prn(ROSMAP_oli_ct, gene_name, colSums(ROSMAP_oli_ct))
saveRDS(ROSMAP_ana_prn, "mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_ana_prn1000_abs_thresh.rds")

## propr
ROSMAP_pr <- propr(counts = t(ROSMAP_oli_ct), metric = "rho", select =gene_name, alpha = NA,p = 100)
saveRDS(ROSMAP_pr, "mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_propr1000_abs_thresh.rds")

## pearson
norm.data <- as.matrix(sc.sel@assays$RNA@data)
ROSMAP_cor_m_pearson <- cor(t(norm.data),method = "pearson")
saveRDS(ROSMAP_cor_m_pearson, "mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_prn1000_abs_thresh.rds")
## spearman
ROSMAP_cor_m_spr <- cor(t(norm.data),method = "spearman")
saveRDS(ROSMAP_cor_m_spr, "mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_spr1000_abs_thresh.rds")

# noise regularization
path2 <- paste0("mean_cor/semi_PD_sparse/simu/noise/ROSMAP/")
if(!file.exists(path2)){
  dir.create(path2,recursive = T)
}
ROSMAP_noise_cor <- noise_fun(sc_obj, sc.sel, sel.gene=gene_name,
                            seed=12052023, path2 = path2)
saveRDS(ROSMAP_noise_cor, "mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_noise1000_abs_thresh.rds")



# PNAS-------------------------------------------------------------------
PNAS_oli_ct <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_sct_cor_NB_simu1000_abs_thresh.rds")
ROSMAP_ori_ests <- readRDS("mean_cor/semi_PD/simu/ROSMAP_NC_Oli_sct1000.rds")
gene_name <- rownames(ROSMAP_ori_ests)

# cscore
PNAS_cscore <- CscoreSimplifiedIRLS(PNAS_oli_ct[gene_name, ] %>% as.matrix %>% t,
                                      colSums(PNAS_oli_ct), covar_weight="regularized")
PNAS_cscore$est <- post_process_est(PNAS_cscore$est)
saveRDS(PNAS_cscore, "mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_cscore1000_abs_thresh.rds")


# sctransfrom
source("/gpfs/gibbs/pi/zhao/xs282/validation/coexp_function.R")
sc_obj <- CreateSeuratObject(counts = PNAS_oli_ct)
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc.sel <- subset(sc_obj, features = gene_name)
PNAS_sct_prn <- sct_cor(sc_obj, sc.sel, gene_name)
saveRDS(PNAS_sct_prn, "mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_sct1000_abs_thresh.rds")

## analytic pearson
PNAS_ana_prn <- ana_prn(PNAS_oli_ct, gene_name, colSums(PNAS_oli_ct))
saveRDS(PNAS_ana_prn, "mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_ana_prn1000_abs_thresh.rds")

## propr
PNAS_pr <- propr(counts = t(PNAS_oli_ct), metric = "rho", select =gene_name, alpha = NA,p = 100)
saveRDS(PNAS_pr, "mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_propr1000_abs_thresh.rds")

## pearson
norm.data <- as.matrix(sc.sel@assays$RNA@data)
PNAS_cor_m_pearson <- cor(t(norm.data),method = "pearson")
saveRDS(PNAS_cor_m_pearson, "mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_prn1000_abs_thresh.rds")
## spearman
PNAS_cor_m_spr <- cor(t(norm.data),method = "spearman")
saveRDS(PNAS_cor_m_spr, "mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_spr1000_abs_thresh.rds")

# noise regularization
path2 <- paste0("mean_cor/semi_PD_sparse/simu/noise/PNAS/")
if(!file.exists(path2)){
  dir.create(path2,recursive = T)
}
PNAS_noise_cor <- noise_fun(sc_obj, sc.sel, sel.gene=gene_name,
                       seed=12052023, path2 = path2)
saveRDS(PNAS_noise_cor, "mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_noise1000_abs_thresh.rds")

