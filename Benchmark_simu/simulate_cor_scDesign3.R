
library(scDesign3)
library(SingleCellExperiment)
library(scran)
library(tidyverse)
library(optparse)

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
colData(sce)$library = colSums(counts(sce))

time_table <- matrix(NA, nrow=1, ncol=2)

# scDesign3 ---------------------------------------------------------------------
set.seed(10222023)

scDesign3_start_time <- Sys.time()
example_simu <- scdesign3(sce = sce, celltype = "celltype",
  pseudotime = NULL, spatial = NULL,
  other_covariates = "library",
  mu_formula = "offset(log(library))",
  corr_formula = "1",
  parallelization = "pbmcmapply")

scDesign3_end_time <- Sys.time()

path <- paste0(prefix, data_name,"/")
if(!file.exists(path)){
  dir.create(path,recursive = T)
}
saveRDS(example_simu, paste0(path,data_name,"_scDesign3.rds"))
time_table[1,] <- c("scDesign3", scDesign3_end_time-scDesign3_start_time)
saveRDS(time_table, paste0(path,"time_scDesign3.rds"))

