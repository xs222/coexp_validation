
# using NB copula to simulate PNAS and ROSMAP NC Oli

setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("compare_simulation/NB_copula/NB_copula_function.R")
source("/gpfs/gibbs/pi/zhao/xs282/coexp-sc/IRLS_CSCORE/CscoreSimplifiedIRLS.R")
source("/gpfs/gibbs/pi/zhao/xs282/validation/cscore_real_data_function.R")

library(pheatmap)
library(Seurat)
library(dplyr)


################################ ROSMAP ###################################
set.seed(10212023)
count_ex <- readRDS("marginal_fit/ROSMAP_NC_Oli_count.rds")
marginal_fit_fn <- 'marginal_fit/ROSMAP_NC_Oli_marginal_fit.rds'
vanilla_ex <- readRDS(marginal_fit_fn)
log10mu <- vanilla_ex$mu
gene_name <- vanilla_ex$gene
names(log10mu) <- gene_name
mu <- 10^log10mu
sum(is.na(mu))

cell_idx <- 1:ncol(count_ex)
cell_name <- colnames(count_ex)[cell_idx]
seq_depth <- colSums(count_ex)[cell_idx]

# use the sampled mean and the trend between mean and alpha to get corresponding alpha
km_Ex5 <- readRDS('marginal_fit/ROSMAP_NC_Oli_ks_fit_5.rds')
fitted_trend <- data.frame(mu=km_Ex5$x, alpha=km_Ex5$y)
log10alpha <- rep(NA,nrow(vanilla_ex))
names(log10alpha) <- vanilla_ex$gene
for (i in 1:nrow(vanilla_ex)){
  idx <- which.min(abs(log10mu[i]-fitted_trend$mu))
  log10alpha[i] <- fitted_trend$alpha[idx]
}
alpha <- 10^log10alpha

# use correlation from SCT
ori_ests <- readRDS("mean_cor/semi_PD/simu/ROSMAP_NC_Oli_sct1000.rds")
cor_mat <- ori_ests
cor_mat[abs(cor_mat)<0.015] <- 0
is.positive.semi.definite(cor_mat)
mean(cor_mat==0)

simu_cscore_cor <- NB_copula(mu, gene_name, seq_depth, cell_name, alpha,
                             cor_mat, ind=F, seed=11132023)
saveRDS(simu_cscore_cor, "mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_sct_cor_NB_simu1000_abs_thresh.rds")

################################ PNAS ###################################
set.seed(10212023)
count_ex <- readRDS("marginal_fit/PNAS_NC_Oli_count.rds")
marginal_fit_fn <- 'marginal_fit/PNAS_NC_Oli_marginal_fit.rds'
vanilla_ex <- readRDS(marginal_fit_fn)
log10mu <- vanilla_ex$mu
gene_name <- vanilla_ex$gene
names(log10mu) <- gene_name
mu <- 10^log10mu
sum(is.na(mu))

cell_idx <- 1:ncol(count_ex)
cell_name <- colnames(count_ex)[cell_idx]
seq_depth <- colSums(count_ex)[cell_idx]

# use the sampled mean and the trend between mean and alpha to get corresponding alpha
km_Ex5 <- readRDS('marginal_fit/PNAS_NC_Oli_ks_fit_5.rds')
fitted_trend <- data.frame(mu=km_Ex5$x, alpha=km_Ex5$y)
log10alpha <- rep(NA,nrow(vanilla_ex))
names(log10alpha) <- vanilla_ex$gene
for (i in 1:nrow(vanilla_ex)){
  idx <- which.min(abs(log10mu[i]-fitted_trend$mu))
  log10alpha[i] <- fitted_trend$alpha[idx]
}
alpha <- 10^log10alpha

ori_ests <- readRDS("mean_cor/semi_PD/simu/PNAS_NC_Oli_sct1000.rds")
cor_mat <- ori_ests
cor_mat[abs(cor_mat)<0.017] <- 0
is.positive.semi.definite(cor_mat)
mean(cor_mat==0)

simu_cscore_cor <- NB_copula(mu, gene_name, seq_depth, cell_name, alpha,
                             cor_mat, ind=F, seed=11132023)
saveRDS(simu_cscore_cor, "mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_sct_cor_NB_simu1000_abs_thresh.rds")



