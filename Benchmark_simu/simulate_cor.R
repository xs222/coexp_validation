
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(glmGamPoi)
# library(SPARSim)
library(optparse)
library(splatter)

option_list <- list(
  make_option(c("--path"), type="character",
              help="real data path", metavar="character",
              default="compare_simulation/compare_simu_multi_data/PNAS_NC/PNAS_NC_EX_sel5000.rds")
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

time_table <- matrix(NA, nrow=14, ncol=2)
path <- paste0(prefix, data_name,"/")
if(!file.exists(path)){
  dir.create(path,recursive = T)
}

# muscat -----------------------------------------------------------------------
# https://bioconductor.org/packages/release/bioc/vignettes/muscat/inst/doc/simulation.html#simdata-simulating-complex-designs
print("muscat ----------------------------------------------------------------------------")
library(muscat)

set.seed(10212023)
muscat_start_time <- Sys.time()

sce$cluster_id <- "1"
sce$sample_id <- sampleID[colnames(sce)]
sce$group_id <- "control"
ref <- prepSim(sce,min_size = NULL)
sim <- simData(ref, dd = FALSE)

muscat_end_time <- Sys.time()
saveRDS(sim, paste0(path,data_name,"_muscat.rds"))
time_table[13,] <- c("muscat", difftime(muscat_end_time,muscat_start_time, units="secs"))
rm(sim)


# SCRIP ------------------------------------------------------------------------
# https://github.com/thecailab/SCRIP/blob/main/vignettes/SCRIPsimu.pdf
print("SCRIP ----------------------------------------------------------------------------")
library(SCRIP)
set.seed(10212023)
SCRIP_start_time <- Sys.time()
params <- splatEstimate(sce)
sim_trend <-  SCRIPsimu(data=ori_ct, params=params, mode="GP-trendedBCV")
SCRIP_end_time <- Sys.time()

saveRDS(sim_trend, paste0(path,data_name,"_SCRIP.rds"))
time_table[14,] <- c("SCRIP", difftime(SCRIP_end_time,SCRIP_start_time, units="secs"))
rm(sim_trend)

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

## 1000 gene cor mat
ncor_gene <- 1000
vanilla <- vanilla[order(vanilla$mu, decreasing = T),]
cor_gene_name <- vanilla$gene[1:ncor_gene]

source("/gpfs/gibbs/pi/zhao/xs282/coexp-sc/IRLS_CSCORE/CscoreSimplifiedIRLS.R")
source("/gpfs/gibbs/pi/zhao/xs282/validation/cscore_real_data_function.R")
cor_ests <- CscoreSimplifiedIRLS(ori_ct[cor_gene_name,] %>% as.matrix %>% t,
                                 colSums(ori_ct), covar_weight="regularized")
cor_ests$est <- post_process_est(cor_ests$est)
filtered_cor_ests <- cor_ests$est
filtered_cor_ests[MatrixBH(cor_ests$p_value) > 0.05] <- 0
filtered_cor_ests[abs(filtered_cor_ests)<0.1] <- 0

cor_mat <- (filtered_cor_ests+t(filtered_cor_ests))/2
simu_nb <- NB_copula(mu, gene_name, seq_depth, cell_name, alpha,
                     cor_mat, ind=F, seed=11132023)
NB_end_time <- Sys.time()
saveRDS(simu_nb, paste0(path,data_name,"_NB.rds"))
time_table[1,] <- c("NB", difftime(NB_end_time,NB_start_time, units="secs"))
rm(simu_nb)


# Splat ------------------------------------------------------------------------
# https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html
# x: SingleCellExperiment
print("Splat ----------------------------------------------------------------------------")

set.seed(10212023)
Splat_start_time <- Sys.time()
params <- splatEstimate(sce)
sim <- splatSimulate(params)
Splat_end_time <- Sys.time()

saveRDS(sim, paste0(path,data_name,"_Splat.rds"))
time_table[3,] <- c("Splat", difftime(Splat_end_time,Splat_start_time, units="secs"))
rm(sim)

# powsimR ----------------------------------------------------------------------
# https://bvieth.github.io/powsimR/articles/powsimR.html#introduction-1
print("powsimR ----------------------------------------------------------------------------")
library(powsimR)

set.seed(10212023)
powsimR_start_time <- Sys.time()

z <- estimateParam(countData = ori_ct, RNAseq = "singlecell",
  Protocol = "UMI", Distribution = "NB",
  Normalisation = "scran", GeneFilter = 0,
  SampleFilter = Inf)
setupres <- Setup(ngenes = z$totalG, estParamRes = z,
  n1 = z$totalS, n2 = 2, # has to be at least 2
  p.DE = 0, pLFC = 0, nsims = 1, setup.seed = 1222024)
y <- simulateDE(SetupRes = setupres,
  Normalisation = "scran",DEmethod = "DESeq2", Counts = TRUE)
powsimR_end_time <- Sys.time()

saveRDS(y$Counts[[1]][[1]][, -c(z$totalS+1, z$totalS+2)],
        paste0(path,data_name,"_powsimR.rds"))
time_table[4,] <- c("powsimR", difftime(powsimR_end_time,powsimR_start_time, units="secs"))
rm(y)


# ZINB-WaVE ---------------------------------------------------------------------
print("ZINB-WaVE ----------------------------------------------------------------------------")
library(zinbwave)

set.seed(10212023)
ZINB_WaVE_start_time <- Sys.time()
params = splatter::zinbEstimate(sce)
newcount = try(splatter::zinbSimulate(params, verbose = FALSE),
               silent = TRUE)
ZINB_WaVE_end_time <- Sys.time()

saveRDS(newcount, paste0(path,data_name,"_ZINB_WaVE.rds"))
time_table[6,] <- c("ZINB_WaVE", difftime(ZINB_WaVE_end_time,ZINB_WaVE_start_time, units="secs"))
rm(newcount)

# SymSim -----------------------------------------------------------------------
# https://github.com/HelenaLC/simulation-comparison/blob/master/code/03-est_pars-SymSim.R
# https://github.com/HelenaLC/simulation-comparison/blob/master/code/04-sim_data-SymSim.R
# https://rpubs.com/zhenmiao/658664
print("Symsim ----------------------------------------------------------------------------")
library(SymSim)

set.seed(10212023)
SymSim_start_time <- Sys.time()

ori_ct_sub <- ori_ct[rowSums(ori_ct)>0, ]
z <- BestMatchParams(tech = "UMI", counts = ori_ct_sub, plotfilename = "foo", n_optimal = 1)
file.remove("foo.pdf")
z$ngenes <- nrow(ori_ct_sub)
z$ncells_total <- ncol(ori_ct_sub)

args <- intersect(names(z),names(formals(SimulateTrueCounts)))
args <- c(z[args], list(randseed = 10122023))
y <- do.call(SimulateTrueCounts, args)

data("gene_len_pool")
gene_len <- sample(gene_len_pool, args$ngenes, TRUE)

args <- intersect(names(z), names(formals(True2ObservedCounts)))
args <- c(z[args], list(true_counts = y$counts,
  meta_cell = y$cell_meta, gene_len = gene_len))
z <- do.call(True2ObservedCounts, args)

SymSim_end_time <- Sys.time()

saveRDS(z, paste0(path,data_name,"_SymSim.rds"))
time_table[7,] <- c("SymSim", difftime(SymSim_end_time,SymSim_start_time, units="secs"))
rm(z)

# scDesign2 ---------------------------------------------------------------------
print("scdesign2 ----------------------------------------------------------------------------")
library(scDesign2)

set.seed(10212023)
scDesign2_start_time <- Sys.time()
ori_ct_scdesign <- ori_ct
colnames(ori_ct_scdesign) <- rep("TpyeA", ncol(ori_ct_scdesign))
copula_result <- fit_model_scDesign2(ori_ct_scdesign, 'TpyeA', sim_method = 'copula', marginal="nb")
sim_count_copula <- simulate_count_scDesign2(copula_result, ncol(ori_ct_scdesign), sim_method = 'copula')
scDesign2_end_time <- Sys.time()

saveRDS(sim_count_copula, paste0(path,data_name,"_scDesign2.rds"))
time_table[8,] <- c("scDesign2", difftime(scDesign2_end_time,scDesign2_start_time, units="secs"))
rm(sim_count_copula)
rm(copula_result)

# # SPsimSeq ---------------------------------------------------------------------
# # https://github.com/CenterForStatistics-UGent/SPsimSeq
# print("SPsimSeq ----------------------------------------------------------------------------")
# library(SPsimSeq)
#
# set.seed(10212023)
# SPsimSeq_start_time <- Sys.time()
# ori_ct_sub <- ori_ct[rowSums(ori_ct)>0, ]
# sim.data.sc <- SPsimSeq(n.sim = 1, s.data = ori_ct_sub,n.genes = nrow(ori_ct_sub),
#                         batch.config = 1, model.zero.prob = TRUE,
#                         result.format = "SCE")
# SPsimSeq_end_time <- Sys.time()
#
# saveRDS(sim.data.sc, paste0(path,data_name,"_SPsimSeq.rds"))
# time_table[9,] <- c("SPsimSeq", difftime(SPsimSeq_end_time,SPsimSeq_start_time, units="secs"))
# rm(sim.data.sc)

# POWSC ------------------------------------------------------------------------
# https://github.com/HelenaLC/simulation-comparison/blob/master/code/04-sim_data-POWSC.R
# https://github.com/suke18/POWSC
print("POWSC ----------------------------------------------------------------------------")
library(POWSC)

set.seed(10212023)
POWSC_start_time <- Sys.time()
y <- Est2Phase(ori_ct)
simu <- Simulate2SCE(n = ncol(ori_ct), perDE = 0, estParas1 = y, estParas2 = y)
POWSC_end_time <- Sys.time()

saveRDS(simu, paste0(path,data_name,"_POWSC.rds"))
time_table[10,] <- c("POWSC", difftime(POWSC_end_time,POWSC_start_time, units="secs"))
rm(simu)

# # ESCO -------------------------------------------------------------------------
# # https://github.com/HelenaLC/simulation-comparison/blob/master/code/03-est_pars-ESCO.R
# # https://github.com/HelenaLC/simulation-comparison/blob/master/code/04-sim_data-ESCO.R
# print("ESCO ----------------------------------------------------------------------------")
# library(ESCO)
#
# set.seed(10212023)
# ESCO_start_time <- Sys.time()
# z <- escoEstimate(ori_ct)
# y <- escoSimulate(z, "single")
# ESCO_end_time <- Sys.time()
#
# saveRDS(y, paste0(path,data_name,"_ESCO.rds"))
# time_table[11,] <- c("ESCO", difftime(ESCO_end_time,ESCO_start_time, units="secs"))
# rm(y)


saveRDS(time_table, paste0(path,"time.rds"))
