.libPaths("/gpfs/gibbs/project/zhao/xs282/R/4.3")

library(dplyr)
library(SPARSim)
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
path <- paste0(prefix, data_name,"/")
if(!file.exists(path)){
  dir.create(path,recursive = T)
}

# SPARSim ----------------------------------------------------------------------
# https://gitlab.com/sysbiobig/sparsim/-/blob/master/vignettes/sparsim.Rmd?ref_type=heads
set.seed(10212023)
SPARSim_start_time <- Sys.time()
Example_count_matrix_norm <- scran_normalization(ori_ct)
cond_A_column_index <- c(1:ncol(ori_ct))
Example_count_matrix_conditions <- list(cond_A = cond_A_column_index)
SPARSim_sim_param <- SPARSim_estimate_parameter_from_data(raw_data = ori_ct,
                                                          norm_data = Example_count_matrix_norm,
                                                          conditions = Example_count_matrix_conditions)
sim_result <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param)
SPARSim_end_time <- Sys.time()

saveRDS(sim_result, paste0(path,data_name,"_SPARSim.rds"))
time_table[1,] <- c("SPARSim", difftime(SPARSim_end_time,SPARSim_start_time, units="secs"))
saveRDS(time_table, paste0(path,"time_SPARSim.rds"))

