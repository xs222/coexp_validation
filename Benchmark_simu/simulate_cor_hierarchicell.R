
library(dplyr)
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

# hierarchicell ----------------------------------------------------------------
# https://github.com/kdzimm/hierarchicell
# https://github.com/HelenaLC/simulation-comparison/blob/master/code/03-est_pars-hierarchicell.R
print("hierarchicell ----------------------------------------------------------------------------")
library(hierarchicell)

set.seed(10212023)
if (grepl("PNAS", data_name)){
  meta <- read.csv("compare_simulation/compare_simu_multi_data/PNAS_NC/PNAS_NC_meta.csv")
  sampleID <- meta$orig.ident
  names(sampleID) <- meta$X
} else if (grepl("ROSMAP", data_name)){
  meta <- read.csv(paste0("compare_simulation/compare_simu_multi_data/ROSMAP/", data_name, "_meta.csv"),row.names = 1)
  sampleID <- meta$projid
  names(sampleID) <- rownames(meta)
}

hierarchicell_start_time <- Sys.time()
df <- data.frame(colnames(ori_ct), sampleID[colnames(ori_ct)],t(ori_ct))
clean_expr_data <- filter_counts(df)
data_summaries <- compute_data_summaries(clean_expr_data, type="Raw")
simulated_counts <- simulate_hierarchicell(data_summaries, n_genes = nrow(ori_ct),
                                           n_cases = 1,cells_per_case = 1,
                                           n_controls = length(unique(sampleID[colnames(ori_ct)])),
                                           cells_per_control = tabulate(as.factor(sampleID[colnames(ori_ct)])),
                                           ncells_variation_type = "Fixed",)
simulated_counts <- simulated_counts[simulated_counts$Status == "Control", ]
simulated_counts <- t(as.matrix(simulated_counts[, grep("^Gene", names(simulated_counts))]))

hierarchicell_end_time <- Sys.time()
saveRDS(simulated_counts, paste0(path,data_name,"_hierarchicell.rds"))
time_table[1,] <- c("hierarchicell", difftime(hierarchicell_end_time,hierarchicell_start_time, units="secs"))
rm(simulated_counts)
saveRDS(time_table, paste0(path,"time_hierarchicell.rds"))

