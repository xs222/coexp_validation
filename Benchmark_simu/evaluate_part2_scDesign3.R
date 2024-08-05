library(optparse)
library(dplyr)
library(ks)

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

eval_result <- readRDS(paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval_scDesign3.rds"))
ori_dist <- readRDS(paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval.rds"))[["ori"]]


simu_method <- c("scDesign3")
eval_metrics <- c("gene_mean", "gene_var", "gene_cv", "gene_frq_zero", "cell_frq_zero",
                  "lib_size", "cell_cor", "gene_cor", "cell_distance", "cell_knn")

kde_stat <- matrix(nrow=length(simu_method), ncol=length(eval_metrics))
rownames(kde_stat) <- simu_method
colnames(kde_stat) <- c(eval_metrics)

ks_stat <- matrix(nrow=length(simu_method), ncol=length(eval_metrics))
rownames(ks_stat) <- simu_method
colnames(ks_stat) <- c(eval_metrics)


for (i in simu_method){
  print(paste(i, "----------------------------------------------"))
  simu_dist <- eval_result[[i]]
  for (j in eval_metrics){
    print(j)
    x <- ori_dist[[j]]
    x <- ifelse(is.na(x), 0, x)
    y <- simu_dist[[j]]
    y <- ifelse(is.na(y), 0, y)

    # ks test
    z <- ks.test(x, y)
    ks_stat[i,j] <- z$statistic

    # kde test
    set.seed(1272024)
    if (length(x)>10000 | length(y)>10000){
      if (length(x)==length(y)){
        idx <- sample(1:length(x), 10000)
        x_new <- x[idx]
        y_new <- y[idx]
      } else{
        if (length(x)>10000){
          x_new <- x[sample(1:length(x), 10000)]
        } else{
          x_new <- x
        }

        if (length(y)>10000){
          y_new <- y[sample(1:length(y), 10000)]
        } else{
          y_new <- y
        }
      }
      if (length(unique(x_new))==1){
        x_new[1] <- x_new[1]+0.0001
      }

      if (length(unique(y_new))==1){
        y_new[1] <- y_new[1]+0.0001
      }
    } else{
      x_new <- x
      y_new <- y

      if (length(unique(x_new))==1){
        x_new[1] <- x_new[1]+0.0001
      }

      if (length(unique(y_new))==1){
        y_new[1] <- y_new[1]+0.0001
      }
    }
    kde_stat[i,j] <- kde.test(x_new, y_new)$zstat
  }

}

ks_stat <- as.data.frame(ks_stat)
ks_stat$data <- data_name
saveRDS(ks_stat, paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval_ks_scDesign3.rds"))

kde_stat <- as.data.frame(kde_stat)
kde_stat$data <- data_name
saveRDS(kde_stat, paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval_kde_scDesign3.rds"))
