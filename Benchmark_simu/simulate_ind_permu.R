
library(dplyr)
library(Seurat)
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

time_table <- matrix(NA, nrow=1, ncol=2)
path <- paste0(prefix, data_name,"_IND/")
if(!file.exists(path)){
  dir.create(path,recursive = T)
}

# permutation -------------------------------------------------------------------
print("permutation ----------------------------------------------------------------------------")
source("/gpfs/gibbs/pi/zhao/xs282/validation/permutation_our/permute_our_function.R")
set.seed(1162024)
permu_start_time <- Sys.time()
perm_count <- permute_our(t(ori_ct), colSums(ori_ct), seed=10232023)
perm_count <- t(perm_count)
permu_end_time <- Sys.time()
saveRDS(perm_count, paste0(path,data_name,"_permu.rds"))
time_table[1,] <- c("permu", difftime(permu_end_time,permu_start_time, units="secs"))

saveRDS(time_table, paste0(path,"time_permu.rds"))

