
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(glmGamPoi)
library(optparse)


option_list <- list(
  make_option(c("--path"), type="character",
              help="real data path", metavar="character",
              default="compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_ex.rds")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

data_path <- opt$path
data_name <- basename(data_path)
prefix <- gsub(data_name, "", data_path)
data_name <- gsub(".rds", "", data_name)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")


# read real data ----------------------------------------------------------------
ori_ct <- readRDS(data_path)
sce <- SingleCellExperiment(assays = list(counts = ori_ct))
sce$celltype <- "TypeI"

time_table <- matrix(NA, nrow=3, ncol=2)
path <- paste0(prefix, data_name,"_IND/")
if(!file.exists(path)){
  dir.create(path,recursive = T)
}

# NB ----------------------------------------------------------------------------
print("NB ----------------------------------------------------------------------------")
set.seed(1162024)
NB_start_time <- Sys.time()
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

# simulate
source("compare_simulation/NB_copula/NB_copula_function.R")
log10mu <- vanilla$mu
gene_name <- vanilla$gene
names(log10mu) <- gene_name
mu <- 10^log10mu

cell_name <- colnames(ori_ct)
seq_depth <- colSums(ori_ct)

# use the sampled mean and the trend between mean and alpha to get corresponding alpha
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
NB_end_time <- Sys.time()
saveRDS(simu_nb, paste0(path,data_name,"_NB.rds"))
time_table[1,] <- c("NB", difftime(NB_end_time,NB_start_time, units="secs"))
rm(simu_nb)


# scDesign2 ---------------------------------------------------------------------
print("scdesign2 ----------------------------------------------------------------------------")
library(scDesign2)

set.seed(10212023)
scDesign2_start_time <- Sys.time()
ori_ct_scdesign <- ori_ct
colnames(ori_ct_scdesign) <- rep("TpyeA", ncol(ori_ct_scdesign))
copula_result <- fit_model_scDesign2(ori_ct_scdesign, 'TpyeA', sim_method = 'ind', marginal="nb")
sim_count_copula <- simulate_count_scDesign2(copula_result, ncol(ori_ct_scdesign), sim_method = 'ind')
scDesign2_end_time <- Sys.time()

saveRDS(sim_count_copula, paste0(path,data_name,"_scDesign2.rds"))
time_table[2,] <- c("scDesign2", difftime(scDesign2_end_time,scDesign2_start_time, units="secs"))
rm(sim_count_copula)
rm(copula_result)

# SPsimSeq ---------------------------------------------------------------------
# https://github.com/CenterForStatistics-UGent/SPsimSeq
print("SPsimSeq ----------------------------------------------------------------------------")
library(SPsimSeq)

set.seed(10212023)
SPsimSeq_start_time <- Sys.time()
ori_ct_sub <- ori_ct[rowSums(ori_ct)>0, ]
# need to modify the function a little bit
# otherwise:Error in FUN(X[[i]], ...) : object 'corMats.batch' not found
# sim.data.sc <- SPsimSeq(n.sim = 1, s.data = ori_ct_sub,n.genes = nrow(ori_ct_sub),
#                         batch.config = 1, model.zero.prob = TRUE,
#                         result.format = "SCE", genewiseCor=FALSE)
sim.data.sc <- tryCatch({
  SPsimSeq(n.sim = 1, s.data = ori_ct_sub,n.genes = nrow(ori_ct_sub),
           batch.config = 1, model.zero.prob = TRUE,
           result.format = "SCE", genewiseCor=FALSE)
}, error = function(e) {
  # Print the error message
  cat("An error occurred:", e$message, "\n")
  NULL  # Return NULL if there's an error
})

SPsimSeq_end_time <- Sys.time()

saveRDS(sim.data.sc, paste0(path,data_name,"_SPsimSeq.rds"))
time_table[3,] <- c("SPsimSeq", difftime(SPsimSeq_end_time,SPsimSeq_start_time, units="secs"))
rm(sim.data.sc)


# ESCO -------------------------------------------------------------------------
# https://github.com/HelenaLC/simulation-comparison/blob/master/code/03-est_pars-ESCO.R
# https://github.com/HelenaLC/simulation-comparison/blob/master/code/04-sim_data-ESCO.R
# randomly selected genes to be correlated
print("ESCO ----------------------------------------------------------------------------")
library(ESCO)

set.seed(10212023)
ESCO_start_time <- Sys.time()
z <- escoEstimate(ori_ct)
y <- escoSimulate(z, "single",withcorr=F)
ESCO_end_time <- Sys.time()
saveRDS(y, paste0(path,data_name,"_ESCO.rds"))
time_table[3,] <- c("ESCO", difftime(ESCO_end_time,ESCO_start_time, units="secs"))


saveRDS(time_table, paste0(path,"time.rds"))
