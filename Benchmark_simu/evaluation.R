library(RANN)
library(scater)
library(scran)
library(SingleCellExperiment)
library(Seurat)
library(optparse)
library(dplyr)

option_list <- list(
  make_option(c("--path"), type="character",
              help="real data path", metavar="character",
              default="compare_simulation/compare_simu_multi_data/PNAS_NC/PNAS_NC_END.rds")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

data_path <- opt$path
data_name <- basename(data_path)
prefix <- gsub(data_name, "", data_path)
data_name <- gsub(".rds", "", data_name)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
path <- paste0(prefix, data_name,"/")

eval_fun <- function(ct){
  eval_metrics <- list()

  sce <- SingleCellExperiment(assays = list(counts = ct))
  sc.obj <- CreateSeuratObject(counts = ct)
  sc.obj <- NormalizeData(sc.obj, scale = 1e+06)
  normalized_data <- GetAssayData(sc.obj, layer = "data")

  # average of logCPM
  eval_metrics[["gene_mean"]] <- rowMeans(normalized_data)

  # variance of logCPM
  eval_metrics[["gene_var"]] <- rowVars(normalized_data)

  # coefficient of variation
  eval_metrics[["gene_cv"]] <- sqrt(eval_metrics[["gene_var"]])/eval_metrics[["gene_mean"]]

  # fraction of genes with zero counts
  eval_metrics[["gene_frq_zero"]] <- rowMeans(ct==0)

  # fraction of cells with zero counts
  eval_metrics[["cell_frq_zero"]] <- colMeans(ct==0)

  # library size (total counts)
  eval_metrics[["lib_size"]] <- colSums(ct)

  # cell-to-cell correlation
  print("cell_cor")
  if (max(sample_idx)>ncol(normalized_data)){
    new_sample_idx <- sample_idx[sample_idx<=ncol(normalized_data)]
    cell_cor <- cor(as.matrix(normalized_data)[,new_sample_idx], method = "spearman", use = "pairwise.complete.obs")
  } else{
    cell_cor <- cor(as.matrix(normalized_data)[,sample_idx], method = "spearman", use = "pairwise.complete.obs")
  }
  eval_metrics[["cell_cor"]] <- cell_cor[upper.tri(cell_cor)]

  # gene-gene correlation
  print("gene_cor")
  source("/gpfs/gibbs/pi/zhao/xs282/coexp-sc/IRLS_CSCORE/CscoreSimplifiedIRLS.R")
  source("/gpfs/gibbs/pi/zhao/xs282/validation/cscore_real_data_function.R")
  genes_selected <- names(sort(eval_metrics[["gene_mean"]], decreasing = T))[1:1000]
  NB_ests <- CscoreSimplifiedIRLS(GetAssayData(sc.obj, layer = "counts")[genes_selected,] %>% as.matrix %>% t,
                                  colSums(GetAssayData(sc.obj, layer = "counts")), covar_weight="regularized")
  NB_ests$est <- post_process_est(NB_ests$est)
  filtered_NB_ests <- NB_ests$est
  filtered_NB_ests[MatrixBH(NB_ests$p_value) > 0.05] <- 0
  eval_metrics[["gene_cor"]] <- filtered_NB_ests[upper.tri(filtered_NB_ests)]


  # pcd: cell-to-cell distance (in PCA space)
  # https://github.com/HelenaLC/simulation-comparison/blob/master/code/05-calc_qc-cell_pcd.R
  print("pcd")
  sce <- logNormCounts(sce)
  stats <- modelGeneVar(sce)
  hvgs <- getTopHVGs(stats, n = 500)
  sce <- runPCA(sce, subset_row = hvgs)
  pca <- reducedDim(sce, "PCA")
  eval_metrics[["cell_distance"]] <- c(dist(pca, upper = TRUE))

  # knn: number of KNN occurrences
  # https://github.com/HelenaLC/simulation-comparison/blob/master/code/05-calc_qc-cell_knn.R
  print("knn")
  k <- round(0.05*ncol(sce)) # (where k = 5% of cells)
  knn <- nn2(pca, k = k+1)
  idx <- knn$nn.idx[, seq(2, k+1)]
  eval_metrics[["cell_knn"]] <- vapply(seq(ncol(sce)), function(i) sum(idx == i), numeric(1))

  return(eval_metrics)
}

# read data -------------------------------------------------------------------
simu_method <- c("NB", "SPsimSeq", "POWSC", "scDesign2", "SymSim", "ZINB_WaVE",
                 "powsimR", "Splat", "SCRIP", "muscat", "ESCO", "SPARSim", "hierarchicell")
eval_result <- list()

# ori data
ori_ct <- readRDS(data_path)
set.seed(1272024)
sample_idx <- sample(1:ncol(ori_ct), min(500, ncol(ori_ct)))
eval_result[["ori"]] <- eval_fun(ori_ct)

# NB--------------------------------------------------------------------
print("NB --------------------------------------------------------------------")
simu_ct_NB <- readRDS(paste0(path,data_name,"_NB.rds"))
eval_result[["NB"]] <- eval_fun(simu_ct_NB)

# # SPsimSeq--------------------------------------------------------------------
# simu_ct_SPsimSeq <- readRDS(paste0(path,data_name,"_SPsimSeq.rds"))
# simu_ct_SPsimSeq <- counts(simu_ct_SPsimSeq[[1]])

# POWSC--------------------------------------------------------------------
print("POWSC --------------------------------------------------------------------")
simu_ct_POWSC <- readRDS(paste0(path,data_name,"_POWSC.rds"))
simu_ct_POWSC <- counts(simu_ct_POWSC$sce)
colnames(simu_ct_POWSC) <- paste0("cell", 1:ncol(simu_ct_POWSC))
eval_result[["POWSC"]] <- eval_fun(simu_ct_POWSC)

# scDesign2 --------------------------------------------------------------------
print("scDesign2 --------------------------------------------------------------------")
simu_ct_scDesign2 <- readRDS(paste0(path,data_name,"_scDesign2.rds"))
colnames(simu_ct_scDesign2) <- paste0("cell", 1:ncol(simu_ct_scDesign2))
rownames(simu_ct_scDesign2) <- paste0("gene",1:nrow(simu_ct_scDesign2))
eval_result[["scDesign2"]] <- eval_fun(simu_ct_scDesign2)

# SymSim--------------------------------------------------------------------
print("SymSim --------------------------------------------------------------------")
simu_ct_SymSim <- readRDS(paste0(path,data_name,"_SymSim.rds"))
simu_ct_SymSim <- simu_ct_SymSim$counts
colnames(simu_ct_SymSim) <- paste0("cell", 1:ncol(simu_ct_SymSim))
rownames(simu_ct_SymSim) <- paste0("gene",1:nrow(simu_ct_SymSim))
eval_result[["SymSim"]] <- eval_fun(simu_ct_SymSim)

# ZINB_WaVE --------------------------------------------------------------------
print("zinb-wave --------------------------------------------------------------------")
simu_ct_ZINB_WaVE <- readRDS(paste0(path,data_name,"_ZINB_WaVE.rds"))
simu_ct_ZINB_WaVE <- as.matrix(counts(simu_ct_ZINB_WaVE))
eval_result[["ZINB_WaVE"]] <- eval_fun(simu_ct_ZINB_WaVE)

# powsimR --------------------------------------------------------------------
print("powsimr --------------------------------------------------------------------")
simu_ct_powsimR <- readRDS(paste0(path,data_name,"_powsimR.rds"))
eval_result[["powsimR"]] <- eval_fun(simu_ct_powsimR)

# Splat --------------------------------------------------------------------
print("splat --------------------------------------------------------------------")
simu_ct_Splat <- readRDS(paste0(path,data_name,"_Splat.rds"))
simu_ct_Splat <- as.matrix(counts(simu_ct_Splat))
eval_result[["Splat"]] <- eval_fun(simu_ct_Splat)

# SCRIP --------------------------------------------------------------------
print("scrip --------------------------------------------------------------------")
simu_ct_SCRIP <- readRDS(paste0(path,data_name,"_SCRIP.rds"))
simu_ct_SCRIP <- counts(simu_ct_SCRIP)
eval_result[["SCRIP"]] <- eval_fun(simu_ct_SCRIP)

# muscat --------------------------------------------------------------------
print("muscat --------------------------------------------------------------------")
simu_ct_muscat <- readRDS(paste0(path,data_name,"_muscat.rds"))
simu_ct_muscat <- counts(simu_ct_muscat)
eval_result[["muscat"]] <- eval_fun(simu_ct_muscat)

# # ESCO --------------------------------------------------------------------
# print("esco --------------------------------------------------------------------")
# simu_ct_ESCO <- readRDS(paste0(path,data_name,"_ESCO.rds"))
# simu_ct_ESCO <- counts(simu_ct_ESCO)
# eval_result[["ESCO"]] <- eval_fun(simu_ct_ESCO)

# SPARSim --------------------------------------------------------------------
print("sparsim --------------------------------------------------------------------")
simu_ct_SPARSim <- readRDS(paste0(path,data_name,"_SPARSim.rds"))
simu_ct_SPARSim <- simu_ct_SPARSim$count_matrix
eval_result[["SPARSim"]] <- eval_fun(simu_ct_SPARSim)

# hierarchicell --------------------------------------------------------------------
print("hierarchicell --------------------------------------------------------------------")
simu_ct_hierarchicell <- readRDS(paste0(path,data_name,"_hierarchicell.rds"))
eval_result[["hierarchicell"]] <- eval_fun(simu_ct_hierarchicell)

saveRDS(eval_result, paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval.rds"))




