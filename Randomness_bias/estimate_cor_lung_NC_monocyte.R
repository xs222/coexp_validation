library(ggpubr)
library(SingleCellExperiment)
library(jsonlite)
library(optparse)
library(DESeq2)
library(tidyr)
library(parallel)
library(propr)
library(Seurat)
set.seed(10222023)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("AFinal/NB_copula_function.R")
source("AFinal/CscoreSimplifiedIRLS.R")
source("AFinal/cscore_real_data_function.R")
source("AFinal/coexp_function.R")
seed <- 10222023

# lung-------------------------------------------------------------------
sc_obj = readRDS("/gpfs/gibbs/pi/zhao/xs282/lung_data/lung_nc_myeloid.rds")
sc_obj <- subset(sc_obj, subset = Subclass_Cell_Identity=="Monocyte")
lung_mye_ct <- as.matrix(sc_obj[["RNA"]]@counts)
vanilla_mye <- readRDS(paste0("revision/estimate_cor_lung/lung_Monocyte_marginal_fit.rds"))
vanilla_mye <- vanilla_mye[order(vanilla_mye$mu, decreasing = T),]
cor_gene_name <- vanilla_mye$gene[1:1000]
saveRDS(cor_gene_name, "revision/estimate_cor_lung/lung_NC_mono_cor_gene.rds")


extract_upp_tri <- function(data, gene_name){
  data <- data[gene_name, gene_name]
  return(data[upper.tri(data, diag = FALSE)])
}

est_mat_lung <- as.data.frame(matrix(NA,ncol=9,nrow=length(cor_gene_name)*(length(cor_gene_name)-1)/2))

# sct
sc_obj <- CreateSeuratObject(counts = lung_mye_ct)
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc.sel <- subset(sc_obj, features = cor_gene_name)
lung_sct_prn <- sct_cor(sc_obj, sc.sel, cor_gene_name)
est_mat_lung$sct <- extract_upp_tri(lung_sct_prn, cor_gene_name)

# noise regularization
path2 <- paste0("revision/estimate_cor_lung/mono/noise/")
if(!file.exists(path2)){
  dir.create(path2,recursive = T)
}
noise_cor <- noise_fun(sc_obj, sc.sel, sel.gene=cor_gene_name,
                       seed=seed, path2 = path2)
est_mat_lung$noise <- extract_upp_tri(noise_cor, cor_gene_name)

# cscore
lung_cscore <- CscoreSimplifiedIRLS(lung_mye_ct[cor_gene_name, ] %>% as.matrix %>% t,
                                      colSums(lung_mye_ct), covar_weight="regularized")
lung_cscore$est <- post_process_est(lung_cscore$est)
est_mat_lung$cscore_p <- extract_upp_tri(lung_cscore$p_value, cor_gene_name)
est_mat_lung$cscore_est <- extract_upp_tri(lung_cscore$est, cor_gene_name)
est_mat_lung$cscore_stat <- extract_upp_tri(lung_cscore$test_stat, cor_gene_name)


## analytic pearson
lung_ana_prn <- ana_prn(lung_mye_ct, cor_gene_name, colSums(lung_mye_ct))
est_mat_lung$ana_prn <- extract_upp_tri(lung_ana_prn, cor_gene_name)


## propr
pr <- propr(counts = t(lung_mye_ct), metric = "rho", select =cor_gene_name, alpha = NA,p = 100)
est_mat_lung$propr <- extract_upp_tri(pr@matrix, cor_gene_name)


## pearson
norm.data <- as.matrix(GetAssayData(sc.sel, assay = "RNA", slot = "data"))
cor_m_pearson <- cor(t(norm.data),method = "pearson")
est_mat_lung$prn <- extract_upp_tri(cor_m_pearson, cor_gene_name)
## spearman
cor_m_spr <- cor(t(norm.data),method = "spearman")
est_mat_lung$spr <- extract_upp_tri(cor_m_spr, cor_gene_name)

saveRDS(est_mat_lung, paste0("revision/estimate_cor_lung/mono/est_cor.rds"))