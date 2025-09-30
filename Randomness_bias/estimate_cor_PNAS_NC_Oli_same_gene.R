library(ggpubr)
library(SingleCellExperiment)
library(jsonlite)
library(optparse)
library(DESeq2)
library(tidyr)
library(parallel)
library(Seurat)
library(propr)
set.seed(10222023)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("AFinal/NB_copula_function.R")
source("AFinal/CscoreSimplifiedIRLS.R")
source("AFinal/cscore_real_data_function.R")
source("AFinal/coexp_function.R")
seed <- 10222023

# PNAS-------------------------------------------------------------------
PNAS_oli_ct <- readRDS("marginal_fit/PNAS_NC_Oli_count.rds")
# vanilla <- readRDS(paste0("marginal_fit/PNAS_NC_Oli_marginal_fit.rds"))
# vanilla <- vanilla[order(vanilla$mu, decreasing = T),]
# cor_gene_name <- vanilla$gene[1:1000]
# saveRDS(cor_gene_name, "revision/estimate_cor_PNAS/PNAS_NC_oli_cor_gene.rds")
rosmap_ests <- readRDS("marginal_fit/ROSMAP_NC_Oli_cscore_cor1000.rds")
cor_gene_name <- colnames(rosmap_ests$est)

extract_upp_tri <- function(data, gene_name){
  data <- data[gene_name, gene_name]
  return(data[upper.tri(data, diag = FALSE)])
}

est_mat_PNAS <- as.data.frame(matrix(NA,ncol=9,nrow=length(cor_gene_name)*(length(cor_gene_name)-1)/2))

# sct
sc_obj <- CreateSeuratObject(counts = PNAS_oli_ct)
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc.sel <- subset(sc_obj, features = cor_gene_name)
PNAS_sct_prn <- sct_cor(sc_obj, sc.sel, cor_gene_name)
est_mat_PNAS$sct <- extract_upp_tri(PNAS_sct_prn, cor_gene_name)

# noise regularization
path2 <- paste0("revision/estimate_cor_PNAS_same_gene/oli/noise/")
if(!file.exists(path2)){
  dir.create(path2,recursive = T)
}
noise_cor <- noise_fun(sc_obj, sc.sel, sel.gene=cor_gene_name,
                       seed=seed, path2 = path2)
est_mat_PNAS$noise <- extract_upp_tri(noise_cor, cor_gene_name)

# cscore
PNAS_cscore <- CscoreSimplifiedIRLS(PNAS_oli_ct[cor_gene_name, ] %>% as.matrix %>% t,
                                      colSums(PNAS_oli_ct), covar_weight="regularized")
PNAS_cscore$est <- post_process_est(PNAS_cscore$est)
est_mat_PNAS$cscore_p <- extract_upp_tri(PNAS_cscore$p_value, cor_gene_name)
est_mat_PNAS$cscore_est <- extract_upp_tri(PNAS_cscore$est, cor_gene_name)
est_mat_PNAS$cscore_stat <- extract_upp_tri(PNAS_cscore$test_stat, cor_gene_name)


## analytic pearson
PNAS_ana_prn <- ana_prn(PNAS_oli_ct, cor_gene_name, colSums(PNAS_oli_ct))
est_mat_PNAS$ana_prn <- extract_upp_tri(PNAS_ana_prn, cor_gene_name)


## propr
pr <- propr(counts = t(PNAS_oli_ct), metric = "rho", select =cor_gene_name, alpha = NA,p = 100)
est_mat_PNAS$propr <- extract_upp_tri(pr@matrix, cor_gene_name)


## pearson
norm.data <- as.matrix(GetAssayData(sc.sel, assay = "RNA", slot = "data"))
cor_m_pearson <- cor(t(norm.data),method = "pearson")
est_mat_PNAS$prn <- extract_upp_tri(cor_m_pearson, cor_gene_name)
## spearman
cor_m_spr <- cor(t(norm.data),method = "spearman")
est_mat_PNAS$spr <- extract_upp_tri(cor_m_spr, cor_gene_name)

saveRDS(est_mat_PNAS, paste0("revision/estimate_cor_PNAS_same_gene/oli/est_cor.rds"))