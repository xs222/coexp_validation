library(optparse)
library(glmGamPoi)
library(ggplot2)
library(peakRAM)
library(Seurat)
library(dplyr)
library(DESeq2)
library(SingleCellExperiment)
library(tidyr)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("compare_simulation/NB_copula/NB_copula_function.R")
source("AFinal/CscoreSimplifiedIRLS.R")
source("AFinal/cscore_real_data_function.R")
source("AFinal/coexp_function.R")

option_list <- list(
  make_option(c("--ncell"), type="integer"),
  make_option(c("--ngene"), type="integer")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

ncell <- opt$ncell
ngene <- opt$ngene

# data from: https://www.pnas.org/doi/10.1073/pnas.2008762117
sc_obj <- readRDS("/gpfs/gibbs/pi/zhao/cs2629/AD_Nancy/seurat_obj_cell_type_labelled.rds")
table(sc_obj@meta.data$cell_type)

sc_obj$disease = sapply(sc_obj$orig.ident, function(or) {substr(strsplit(or, '_')[[1]][2], 1, 2)})
table(sc_obj$disease)

length(unique(sc_obj$orig.ident[sc_obj$disease=="NC"]))
sc_obj = subset(sc_obj, subset = disease == "NC")

count_ex <- as.matrix(sc_obj[["RNA"]]@counts[, which(sc_obj$cell_type=="Ex" & sc_obj$disease=="NC")])
print(dim(count_ex))

# sample genes and cells
set.seed(2242026)
cell_idx = sample(1:ncol(count_ex), ncell)
gene_idx = sample(1:nrow(count_ex), ngene)
ori_ct = count_ex[gene_idx, cell_idx]

# para est --------------------------------------------------------------------
bench_para <- peakRAM({
  # estimate marginal parameter
  size_factors_sel <- colSums(ori_ct)
  gp_ex <- glm_gp(ori_ct, size_factors = size_factors_sel, verbose = T, overdispersion_shrinkage = T, do_cox_reid_adjustment = T)
  vanilla <- data.frame(mu = log10(exp(gp_ex$Beta[,1])),
                        alpha = -log10(gp_ex$overdispersions),
                        deviances = gp_ex$deviances,
                        gene=rownames(ori_ct))

  # fit the kernel regression
  vanilla$up <- ifelse(vanilla$alpha>2.5, "upper","lower")
  ## fit line (only use the lower cluster)
  marginal_fit_PNAS_sel <- vanilla[vanilla$up=="lower",]
  ## kernel smooth
  km5 <- ksmooth(marginal_fit_PNAS_sel$mu, marginal_fit_PNAS_sel$alpha,
                 kernel="normal", bandwidth = bw.SJ(marginal_fit_PNAS_sel$mu)*5)
})

cpu_para <- system.time({
  # estimate marginal parameter
  size_factors_sel <- colSums(ori_ct)
  gp_ex <- glm_gp(ori_ct, size_factors = size_factors_sel, verbose = T, overdispersion_shrinkage = T, do_cox_reid_adjustment = T)
  vanilla <- data.frame(mu = log10(exp(gp_ex$Beta[,1])),
                        alpha = -log10(gp_ex$overdispersions),
                        deviances = gp_ex$deviances,
                        gene=rownames(ori_ct))

  # fit the kernel regression
  vanilla$up <- ifelse(vanilla$alpha>2.5, "upper","lower")
  ## fit line (only use the lower cluster)
  marginal_fit_PNAS_sel <- vanilla[vanilla$up=="lower",]
  ## kernel smooth
  km5 <- ksmooth(marginal_fit_PNAS_sel$mu, marginal_fit_PNAS_sel$alpha,
                 kernel="normal", bandwidth = bw.SJ(marginal_fit_PNAS_sel$mu)*5)
})


# estimate marginal parameter
size_factors_sel <- colSums(ori_ct)
gp_ex <- glm_gp(ori_ct, size_factors = size_factors_sel, verbose = T, overdispersion_shrinkage = T, do_cox_reid_adjustment = T)
vanilla <- data.frame(mu = log10(exp(gp_ex$Beta[,1])),
                      alpha = -log10(gp_ex$overdispersions),
                      deviances = gp_ex$deviances,
                      gene=rownames(ori_ct))

# fit the kernel regression
vanilla$up <- ifelse(vanilla$alpha>2.5, "upper","lower")
## fit line (only use the lower cluster)
marginal_fit_PNAS_sel <- vanilla[vanilla$up=="lower",]
## kernel smooth
km5 <- ksmooth(marginal_fit_PNAS_sel$mu, marginal_fit_PNAS_sel$alpha,
               kernel="normal", bandwidth = bw.SJ(marginal_fit_PNAS_sel$mu)*5)


# simulation -------------------------------------------------------------------
bench_sim <- peakRAM({
  # simulate
  source("AFinal/NB_copula_function.R")
  log10mu <- vanilla$mu
  gene_name <- vanilla$gene
  names(log10mu) <- gene_name
  mu <- 10^log10mu

  cell_name <- colnames(ori_ct)
  seq_depth <- colSums(ori_ct)

  # use the sampled mean and the trend between mean and alpha to get the corresponding alpha
  fitted_trend <- data.frame(mu=km5$x, alpha=km5$y)
  log10alpha <- rep(NA,nrow(vanilla))
  names(log10alpha) <- vanilla$gene
  for (i in 1:nrow(vanilla)){
    idx <- which.min(abs(log10mu[i]-fitted_trend$mu))
    log10alpha[i] <- fitted_trend$alpha[idx]
  }
  alpha <- 10^log10alpha

  simu_nb <- NB_copula(mu, gene_name, seq_depth, cell_name, alpha,
                       cor_mat=NULL, ind=T, seed=11132023)
})

cpu_sim <- system.time({
  # simulate
  source("AFinal/NB_copula_function.R")
  log10mu <- vanilla$mu
  gene_name <- vanilla$gene
  names(log10mu) <- gene_name
  mu <- 10^log10mu

  cell_name <- colnames(ori_ct)
  seq_depth <- colSums(ori_ct)

  # use the sampled mean and the trend between mean and alpha to get the corresponding alpha
  fitted_trend <- data.frame(mu=km5$x, alpha=km5$y)
  log10alpha <- rep(NA,nrow(vanilla))
  names(log10alpha) <- vanilla$gene
  for (i in 1:nrow(vanilla)){
    idx <- which.min(abs(log10mu[i]-fitted_trend$mu))
    log10alpha[i] <- fitted_trend$alpha[idx]
  }
  alpha <- 10^log10alpha

  simu_nb <- NB_copula(mu, gene_name, seq_depth, cell_name, alpha,
                       cor_mat=NULL, ind=T, seed=11132023)
})

# simulate
source("AFinal/NB_copula_function.R")
log10mu <- vanilla$mu
gene_name <- vanilla$gene
names(log10mu) <- gene_name
mu <- 10^log10mu

cell_name <- colnames(ori_ct)
seq_depth <- colSums(ori_ct)

# use the sampled mean and the trend between mean and alpha to get the corresponding alpha
fitted_trend <- data.frame(mu=km5$x, alpha=km5$y)
log10alpha <- rep(NA,nrow(vanilla))
names(log10alpha) <- vanilla$gene
for (i in 1:nrow(vanilla)){
  idx <- which.min(abs(log10mu[i]-fitted_trend$mu))
  log10alpha[i] <- fitted_trend$alpha[idx]
}
alpha <- 10^log10alpha

simu_nb <- NB_copula(mu, gene_name, seq_depth, cell_name, alpha,
                     cor_mat=NULL, ind=T, seed=11132023)


# cor est -------------------------------------------------------------------
ncor_gene = c(500, 1000, 1500)
bench_noise_ls = list()
bench_cscore_ls = list()
bench_ana_ls = list()
cpu_noise_ls = list()
cpu_cscore_ls = list()
cpu_ana_ls = list()

for (i in ncor_gene){
  set.seed(2242026)
  cor_gene_name = sample(rownames(simu_nb), i)

  # noise
  path2 <- paste0("BIB_R1/scalability/res/simu/ncell",ncell,"_ngene",ngene,"/")
  if(!file.exists(path2)){
    dir.create(path2,recursive = T)
  }
  bench_noise <- peakRAM({
    sc_obj <- CreateSeuratObject(counts = simu_nb)
    sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    sc.sel <- subset(sc_obj, features = cor_gene_name)
    noise_cor <- noise_fun(sc_obj, sc.sel, sel.gene=cor_gene_name,
                           seed=2242026, path2 = path2)
  })
  bench_noise_ls[[paste0("ncor", i)]] = bench_noise

  cpu_noise <- system.time({
    sc_obj <- CreateSeuratObject(counts = simu_nb)
    sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    sc.sel <- subset(sc_obj, features = cor_gene_name)
    noise_cor <- noise_fun(sc_obj, sc.sel, sel.gene=cor_gene_name,
                           seed=2242026, path2 = path2)
  })
  cpu_noise_ls[[paste0("ncor", i)]] = cpu_noise

  # cscore
  bench_cscore <- peakRAM({
    ROSMAP_cscore <- CscoreSimplifiedIRLS(simu_nb[cor_gene_name, ] %>% as.matrix %>% t,
                                          colSums(simu_nb), covar_weight="regularized")
    ROSMAP_cscore$est <- post_process_est(ROSMAP_cscore$est)
  })
  bench_cscore_ls[[paste0("ncor", i)]] = bench_cscore

  cpu_cscore <- system.time({
    ROSMAP_cscore <- CscoreSimplifiedIRLS(simu_nb[cor_gene_name, ] %>% as.matrix %>% t,
                                          colSums(simu_nb), covar_weight="regularized")
    ROSMAP_cscore$est <- post_process_est(ROSMAP_cscore$est)
  })
  cpu_cscore_ls[[paste0("ncor", i)]] = cpu_cscore

  # analytic pearson
  bench_ana <- peakRAM({
    ROSMAP_ana_prn <- ana_prn(simu_nb, cor_gene_name, colSums(simu_nb))
  })
  bench_ana_ls[[paste0("ncor", i)]] = bench_ana

  cpu_ana <- system.time({
    ROSMAP_ana_prn <- ana_prn(simu_nb, cor_gene_name, colSums(simu_nb))
  })
  cpu_ana_ls[[paste0("ncor", i)]] = cpu_ana
}

res = list(bench_para=bench_para, cpu_para=cpu_para, bench_sim=bench_sim,
           cpu_sim=cpu_sim, bench_noise_ls=bench_noise_ls,
           bench_cscore_ls=bench_cscore_ls, bench_ana_ls=bench_ana_ls,
           cpu_noise_ls=cpu_noise_ls, cpu_cscore_ls=cpu_cscore_ls,
           cpu_ana_ls=cpu_ana_ls)

saveRDS(res, paste0("BIB_R1/scalability/res/ncell",ncell, "_ngene", ngene, ".rds"))

##################################################################33
# path <- "ml miniconda; conda activate validation; Rscript /gpfs/gibbs/pi/zhao/xs282/validation/BIB_R1/scalability/scalability.r --ncell "
# ncells = c(8000, 16000, 24000)
# ngenes = c(9000, 18000, 27000)
#
# job <- c()
# for (ncell in ncells){
#   for (ngene in ngenes){
#     job <- c(job, paste0(path, ncell, " --ngene ", ngene))
#   }
# }
# job_mx <- matrix(job, nrow=length(job))
# write.table(job_mx, file = "/gpfs/gibbs/pi/zhao/xs282/validation/BIB_R1/scalability/joblist_scalability.txt", sep = "\t", row.names = F, col.names = F, quote = F)
# dsq --job-file /gpfs/gibbs/pi/zhao/xs282/validation/BIB_R1/scalability/joblist_scalability.txt --mem-per-cpu 200g -t 2:00:00 --mail-type ALL --partition scavenge
#
#
