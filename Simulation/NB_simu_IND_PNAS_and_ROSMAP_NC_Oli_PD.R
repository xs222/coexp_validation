
# using NB to simulate IND PNAS and ROSMAP NC Oli

setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("compare_simulation/NB_copula/NB_copula_function.R")
source("/gpfs/gibbs/pi/zhao/xs282/coexp-sc/IRLS_CSCORE/CscoreSimplifiedIRLS.R")
source("/gpfs/gibbs/pi/zhao/xs282/validation/cscore_real_data_function.R")
source("/gpfs/gibbs/pi/zhao/xs282/validation/coexp_function.R")

library(pheatmap)
library(Seurat)
library(dplyr)
library(optparse)
library(propr)
library(SingleCellExperiment)
library(optparse)
library(DESeq2)
library(tidyr)
library(parallel)

option_list <- list(
  make_option(c("--prefix"), type="character",
              default='s1',
              help="prefix of the experiment result to look at", metavar="character"),
  make_option(c("--data"), type="character",
              default='ROSMAP',
              help="data name", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

prefix <- opt$prefix
dname <- opt$data
i <- as.numeric(gsub("s","",prefix))
seed <- 10212023+i

set.seed(seed)
count_ex <- readRDS(paste0("marginal_fit/",dname,"_NC_Oli_count.rds"))
marginal_fit_fn <- paste0("marginal_fit/",dname,"_NC_Oli_PD_sparse_marginal_fit_simulated.rds")
vanilla_ex <- readRDS(marginal_fit_fn)
log10mu <- vanilla_ex$mu
gene_name <- vanilla_ex$gene
names(log10mu) <- gene_name
mu <- 10^log10mu

cell_idx <- 1:ncol(count_ex)
cell_name <- colnames(count_ex)[cell_idx]
seq_depth <- colSums(count_ex)[cell_idx]

# use the sampled mean and the trend between mean and alpha to get corresponding alpha
km_Ex5 <- readRDS(paste0('marginal_fit/',dname,'_NC_Oli_PD_sparse_ks_fit_5_simulated.rds'))
fitted_trend <- data.frame(mu=km_Ex5$x, alpha=km_Ex5$y)
log10alpha <- rep(NA,nrow(vanilla_ex))
names(log10alpha) <- vanilla_ex$gene
for (i in 1:nrow(vanilla_ex)){
  idx <- which.min(abs(log10mu[i]-fitted_trend$mu))
  log10alpha[i] <- fitted_trend$alpha[idx]
}
alpha <- 10^log10alpha

simu_ROSMAP <- NB_copula(mu, gene_name, seq_depth, cell_name, alpha,
                         ind=T, seed=seed)
rm(gene_name)
if(!file.exists(paste0("mean_cor/semi_PD_sparse/simu/IND/", dname,"/"))){
  dir.create(paste0("mean_cor/semi_PD_sparse/simu/IND/", dname,"/"),recursive = T)
}
saveRDS(simu_ROSMAP, paste0("mean_cor/semi_PD_sparse/simu/IND/", dname,"/simu_", prefix, ".rds"))

### estimate correlation --------------------------------------------------------
ori_ests <- readRDS("mean_cor/semi_PD/simu/ROSMAP_NC_Oli_sct1000.rds")
cor_gene_name <- colnames(ori_ests)

extract_upp_tri <- function(data, gene_name){
  data <- data[gene_name, gene_name]
  return(data[upper.tri(data, diag = FALSE)])
}

est_mat_ROSMAP <- as.data.frame(matrix(NA,ncol=9,nrow=length(cor_gene_name)*(length(cor_gene_name)-1)/2))

# sct
sc_obj <- CreateSeuratObject(counts = simu_ROSMAP)
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc.sel <- subset(sc_obj, features = cor_gene_name)
ROSMAP_sct_prn <- sct_cor(sc_obj, sc.sel, cor_gene_name)
est_mat_ROSMAP$sct <- extract_upp_tri(ROSMAP_sct_prn, cor_gene_name)

# noise regularization
path2 <- paste0("mean_cor/semi_PD_sparse/simu/IND/",dname,"/",prefix,"/")
if(!file.exists(path2)){
  dir.create(path2,recursive = T)
}
noise_cor <- noise_fun(sc_obj, sc.sel, sel.gene=cor_gene_name,
                       seed=seed, path2 = path2)
est_mat_ROSMAP$noise <- extract_upp_tri(noise_cor, cor_gene_name)

# cscore
ROSMAP_cscore <- CscoreSimplifiedIRLS(simu_ROSMAP[cor_gene_name, ] %>% as.matrix %>% t,
                                      colSums(simu_ROSMAP), covar_weight="regularized")
ROSMAP_cscore$est <- post_process_est(ROSMAP_cscore$est)
est_mat_ROSMAP$cscore_p <- extract_upp_tri(ROSMAP_cscore$p_value, cor_gene_name)
est_mat_ROSMAP$cscore_est <- extract_upp_tri(ROSMAP_cscore$est, cor_gene_name)
est_mat_ROSMAP$cscore_stat <- extract_upp_tri(ROSMAP_cscore$test_stat, cor_gene_name)


## analytic pearson
ROSMAP_ana_prn <- ana_prn(simu_ROSMAP, cor_gene_name, colSums(simu_ROSMAP))
est_mat_ROSMAP$ana_prn <- extract_upp_tri(ROSMAP_ana_prn, cor_gene_name)


## propr
pr <- propr(counts = t(simu_ROSMAP), metric = "rho", select =cor_gene_name, alpha = NA,p = 100)
est_mat_ROSMAP$propr <- extract_upp_tri(pr@matrix, cor_gene_name)


## pearson
norm.data <- as.matrix(GetAssayData(sc.sel, assay = "RNA", slot = "data"))
cor_m_pearson <- cor(t(norm.data),method = "pearson")
est_mat_ROSMAP$prn <- extract_upp_tri(cor_m_pearson, cor_gene_name)
## spearman
cor_m_spr <- cor(t(norm.data),method = "spearman")
est_mat_ROSMAP$spr <- extract_upp_tri(cor_m_spr, cor_gene_name)

saveRDS(est_mat_ROSMAP, paste0("mean_cor/semi_PD_sparse/simu/IND/", dname,"/est_cor_", prefix, ".rds"))




