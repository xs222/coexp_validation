library(matrixStats)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")

seed <- 12052023
set.seed(seed)

ori_ests <- readRDS("marginal_fit/ROSMAP_NC_Oli_cscore_cor1000.rds")
gene_name <- colnames(ori_ests$est)

############# ROSMAP
# cor estimations
ROSMAP_est <- readRDS("real/oli/est_cor.rds")

n_pair <- 1000*(1000-1)/2
n_permu <- 500
nulll <- list()
nulll[["sct"]] <-  matrix(nrow=n_pair, ncol=n_permu)
nulll[["cscore_est"]] <- matrix(nrow=n_pair, ncol=n_permu)
nulll[["propr"]] <- matrix(nrow=n_pair, ncol=n_permu)
nulll[["spr"]] <- matrix(nrow=n_pair, ncol=n_permu)
nulll[["prn"]] <- matrix(nrow=n_pair, ncol=n_permu)
nulll[["ana_prn"]] <- matrix(nrow=n_pair, ncol=n_permu)
nulll[["noise"]] <- matrix(nrow=n_pair, ncol=n_permu)

for (i in 1:n_permu){
  print(i)
  path <- paste0("real/oli/simu_IND/ROSMAP/est_cor_s", i, ".rds")
  est_mat <- readRDS(path)

  nulll[["sct"]][,i] <- est_mat[,"sct"]
  nulll[["cscore_est"]][,i] <- est_mat[,"cscore_est"]
  nulll[["prn"]][,i] <- est_mat[,"prn"]
  nulll[["propr"]][,i] <- est_mat[,"propr"]
  nulll[["spr"]][,i] <- est_mat[,"spr"]
  nulll[["ana_prn"]][,i] <- est_mat[,"ana_prn"]
  nulll[["noise"]][,i] <- est_mat[,"noise"]
}

saveRDS(nulll, paste0("real/oli/ROSMAP_simu_null.rds"))

nulll <- readRDS("real/oli/ROSMAP_simu_null.rds")
mean_mat <- data.frame(sct=rowMeans(nulll$sct), prn=rowMeans(nulll$prn),
                       spr=rowMeans(nulll$spr), propr=rowMeans(nulll$propr),
                       cscore_est=rowMeans(nulll$cscore_est),
                       ana_prn=rowMeans(nulll$ana_prn), noise=rowMeans(nulll$noise))

sd_mat <- data.frame(sct=rowSds(nulll$sct), prn=rowSds(nulll$prn),
                     spr=rowSds(nulll$spr), propr=rowSds(nulll$propr),
                     cscore_est=rowSds(nulll$cscore_est),
                     ana_prn=rowSds(nulll$ana_prn), noise=rowSds(nulll$noise))

est_mat <- ROSMAP_est
p_val <- data.frame(sct=pnorm(abs(est_mat$sct), mean = mean_mat$sct, sd=sd_mat$sct, lower.tail = F)+
                      pnorm(-abs(est_mat$sct), mean = mean_mat$sct, sd=sd_mat$sct, lower.tail = T),
                   prn=pnorm(abs(est_mat$prn), mean = mean_mat$prn, sd=sd_mat$prn, lower.tail = F)+
                     pnorm(-abs(est_mat$prn), mean = mean_mat$prn, sd=sd_mat$prn, lower.tail = T),
                   spr=pnorm(abs(est_mat$spr), mean = mean_mat$spr, sd=sd_mat$spr, lower.tail = F)+
                     pnorm(-abs(est_mat$spr), mean = mean_mat$spr, sd=sd_mat$spr, lower.tail = T),
                   propr=pnorm(abs(est_mat$propr), mean = mean_mat$propr, sd=sd_mat$propr, lower.tail = F)+
                     pnorm(-abs(est_mat$propr), mean = mean_mat$propr, sd=sd_mat$propr, lower.tail = T),
                   cscore_est=pnorm(abs(est_mat$cscore_est), mean = mean_mat$cscore_est, sd=sd_mat$cscore_est, lower.tail = F)+
                     pnorm(-abs(est_mat$cscore_est), mean = mean_mat$cscore_est, sd=sd_mat$cscore_est, lower.tail = T),
                   ana_prn=pnorm(abs(est_mat$ana_prn), mean = mean_mat$ana_prn, sd=sd_mat$ana_prn, lower.tail = F)+
                     pnorm(-abs(est_mat$ana_prn), mean = mean_mat$ana_prn, sd=sd_mat$ana_prn, lower.tail = T),
                   noise=pnorm(abs(est_mat$noise), mean = mean_mat$noise, sd=sd_mat$noise, lower.tail = F)+
                     pnorm(-abs(est_mat$noise), mean = mean_mat$noise, sd=sd_mat$noise, lower.tail = T))
saveRDS(p_val, paste0("real/oli/ROSMAP_oli_norm_p.rds"))

