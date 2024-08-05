# null distribution from /gpfs/gibbs/pi/zhao/xs282/validation/mean_cor/NB_simu_IND_PNAS_and_ROSMAP_NC_Oli_11_14_2023.R
# cor est from /gpfs/gibbs/pi/zhao/xs282/validation/mean_cor/NB_simu_PNAS_and_ROSMAP_NC_Oli_11_13_2023.R
#              /gpfs/gibbs/pi/zhao/xs282/validation/mean_cor/estimate_cor_1000_abs_thresh_PNAS_and_ROSMAP_NC_Oli_11_15_2023.R
library(matrixStats)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")

seed <- 12052023
set.seed(seed)

ROSMAP_ori_ests <- readRDS("mean_cor/semi_PD/simu/ROSMAP_NC_Oli_sct1000.rds")
gene_name <- rownames(ROSMAP_ori_ests)

############# ROSMAP
# cor estimations
ROSMAP_cscore <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_cscore1000_abs_thresh.rds")
ROSMAP_cscore_est <- ROSMAP_cscore$est[gene_name, gene_name]
PNAS_cscore <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_cscore1000_abs_thresh.rds")
PNAS_cscore_est <- PNAS_cscore$est[gene_name, gene_name]

ROSMAP_sct <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_sct1000_abs_thresh.rds")
ROSMAP_sct_est <- ROSMAP_sct[gene_name, gene_name]
PNAS_sct <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_sct1000_abs_thresh.rds")
PNAS_sct_est <- PNAS_sct[gene_name, gene_name]

ROSMAP_ana_prn <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_ana_prn1000_abs_thresh.rds")
ROSMAP_ana_prn_est <- ROSMAP_ana_prn[gene_name, gene_name]
PNAS_ana_prn <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_ana_prn1000_abs_thresh.rds")
PNAS_ana_prn_est <- PNAS_ana_prn[gene_name, gene_name]

ROSMAP_noise <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_noise1000_abs_thresh.rds")
ROSMAP_noise_est <- ROSMAP_noise[gene_name, gene_name]
PNAS_noise <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_noise1000_abs_thresh.rds")
PNAS_noise_est <- PNAS_noise[gene_name, gene_name]

ROSMAP_prn <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_prn1000_abs_thresh.rds")
ROSMAP_prn_est <- ROSMAP_prn[gene_name, gene_name]
PNAS_prn <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_prn1000_abs_thresh.rds")
PNAS_prn_est <- PNAS_prn[gene_name, gene_name]

ROSMAP_spr <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_spr1000_abs_thresh.rds")
ROSMAP_spr_est <- ROSMAP_spr[gene_name, gene_name]
PNAS_spr <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_spr1000_abs_thresh.rds")
PNAS_spr_est <- PNAS_spr[gene_name, gene_name]

ROSMAP_propr <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_propr1000_abs_thresh.rds")
ROSMAP_propr_est <- ROSMAP_propr@matrix[gene_name, gene_name]
PNAS_propr <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_propr1000_abs_thresh.rds")
PNAS_propr_est <- PNAS_propr@matrix[gene_name, gene_name]

ROSMAP_est <- data.frame(sct=ROSMAP_sct_est[upper.tri(ROSMAP_sct_est)],
                         cscore=ROSMAP_cscore_est[upper.tri(ROSMAP_cscore_est)],
                         ana_prn=ROSMAP_ana_prn_est[upper.tri(ROSMAP_ana_prn_est)],
                         noise=ROSMAP_noise_est[upper.tri(ROSMAP_noise_est)],
                         prn=ROSMAP_prn_est[upper.tri(ROSMAP_prn_est)],
                         spr=ROSMAP_spr_est[upper.tri(ROSMAP_spr_est)],
                         propr=ROSMAP_propr_est[upper.tri(ROSMAP_propr_est)])

PNAS_est <- data.frame(sct=PNAS_sct_est[upper.tri(PNAS_sct_est)],
                         cscore=PNAS_cscore_est[upper.tri(PNAS_cscore_est)],
                         ana_prn=PNAS_ana_prn_est[upper.tri(PNAS_ana_prn_est)],
                         noise=PNAS_noise_est[upper.tri(PNAS_noise_est)],
                         prn=PNAS_prn_est[upper.tri(PNAS_prn_est)],
                         spr=PNAS_spr_est[upper.tri(PNAS_spr_est)],
                         propr=PNAS_propr_est[upper.tri(PNAS_propr_est)])

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
  path <- paste0("mean_cor/semi_PD_sparse/simu/IND/ROSMAP/est_cor_s", i, ".rds")
  est_mat <- readRDS(path)

  nulll[["sct"]][,i] <- est_mat[,"sct"]
  nulll[["cscore_est"]][,i] <- est_mat[,"cscore_est"]
  nulll[["prn"]][,i] <- est_mat[,"prn"]
  nulll[["propr"]][,i] <- est_mat[,"propr"]
  nulll[["spr"]][,i] <- est_mat[,"spr"]
  nulll[["ana_prn"]][,i] <- est_mat[,"ana_prn"]
  nulll[["noise"]][,i] <- est_mat[,"noise"]
}

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
                   cscore_est=pnorm(abs(est_mat$cscore), mean = mean_mat$cscore_est, sd=sd_mat$cscore_est, lower.tail = F)+
                     pnorm(-abs(est_mat$cscore), mean = mean_mat$cscore_est, sd=sd_mat$cscore_est, lower.tail = T),
                   ana_prn=pnorm(abs(est_mat$ana_prn), mean = mean_mat$ana_prn, sd=sd_mat$ana_prn, lower.tail = F)+
                     pnorm(-abs(est_mat$ana_prn), mean = mean_mat$ana_prn, sd=sd_mat$ana_prn, lower.tail = T),
                   noise=pnorm(abs(est_mat$noise), mean = mean_mat$noise, sd=sd_mat$noise, lower.tail = F)+
                     pnorm(-abs(est_mat$noise), mean = mean_mat$noise, sd=sd_mat$noise, lower.tail = T))
saveRDS(p_val, paste0("mean_cor/semi_PD_sparse/simu/ROSMAP_simu_norm_p.rds"))

plot_dat <- data.frame(pair2 = nulll$sct[17770,],
                       pair3 = nulll$sct[163976,])


p2 <- ggplot(plot_dat, aes(x=pair2))+geom_histogram(alpha=0.6)+
  geom_vline(aes(xintercept = est_mat$sct[17770]), color = "darkred", linetype="dashed", size=1)+
  theme_bw()+labs(x="Gene pair 1", y="Counts")+
  annotate("text", x = -0.027, y = 20, label = paste0("Est = ", round(est_mat$sct[17770],4)),
           color = "darkred", angle = 90, size=5)+
  annotate("text", x = -0.0015, y = 40.5, label = expression(italic("Adjusted t-dist p = 0.045")),
           color = "darkblue", size=4.5)+
  annotate("text", x = 0.0015, y = 35, label = expression(italic("Adjusted empirical p = 0.064")),
           color = "darkblue", size=4.5)+
  theme(text = element_text(size = 14))
p2

p3 <- ggplot(plot_dat, aes(x=pair3))+geom_histogram(alpha=0.6)+
  geom_vline(aes(xintercept = est_mat$sct[163976]), color = "darkred", linetype="dashed", size=1)+
  theme_bw()+labs(x="Gene pair 2", y="Counts")+
  annotate("text", x = -0.027, y = 25, label = paste0("Est = ", round(est_mat$sct[163976],4)),
           color = "darkred", angle = 90, size=5)+
  annotate("text", x = -0.002, y = 50, label = expression(italic("Adjusted t-dist p = 0.045")),
           color = "darkblue", size=4.5)+
  annotate("text", x = 0.001, y = 43, label = expression(italic("Adjusted empirical p = 0.031")),
           color = "darkblue", size=4.5)+
  theme(text = element_text(size = 14))
p3

pdf('mean_cor/semi_PD_sparse/figures/t_dist_improper_v2.pdf', width = 12, height = 2.5, onefile = T)
ggarrange(NULL,p2,p3, NULL, ncol=4, nrow=1, labels = c("","A", "B",""), widths = c(0.8,1,1,0.8))
dev.off()


ROSMAP_ori_ests[abs(ROSMAP_ori_ests)<0.015] <- 0
p_val$true <- ROSMAP_ori_ests[upper.tri(ROSMAP_ori_ests)]
p_val_true <- p_val[p_val$true==0,]
prepare_data <- function(x, idx, name){
  return(data.frame(emp=sort(x),
                    unif=(1:length(idx))/(length(idx)),
                    cor_method=name))
}

set.seed(5112024)
idx <- sample(1:nrow(p_val_true), 10000)
p_dat <- rbind(prepare_data(p_val_true$sct[idx], idx, "sctransform"),
               prepare_data(p_val_true$prn[idx], idx, "Pearson"),
               prepare_data(p_val_true$spr[idx], idx, "Spearman"),
               prepare_data(p_val_true$ana_prn[idx], idx, "Analytic PR"),
               prepare_data(p_val_true$propr[idx], idx, "propr"),
               prepare_data(p_val_true$cscore_est[idx], idx, "CS-CORE"),
               prepare_data(p_val_true$noise[idx], idx, "Noise \nRegularization"))

pdf('mean_cor/semi_PD_sparse/figures/qq_our_rosmap.pdf', width = 12, height = 2.5, onefile = T)
ggplot(p_dat, aes(x=unif, y=emp))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, col="red")+
  facet_wrap(.~cor_method, nrow=1)+theme_bw()+
  labs(x="Theoretical Quantiles", y="Sample Quantiles")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14))
dev.off()



# PNAS
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
  path <- paste0("mean_cor/semi_PD_sparse/simu/IND/PNAS/est_cor_s", i, ".rds")
  est_mat <- readRDS(path)

  nulll[["sct"]][,i] <- est_mat[,"sct"]
  nulll[["cscore_est"]][,i] <- est_mat[,"cscore_est"]
  nulll[["prn"]][,i] <- est_mat[,"prn"]
  nulll[["propr"]][,i] <- est_mat[,"propr"]
  nulll[["spr"]][,i] <- est_mat[,"spr"]
  nulll[["ana_prn"]][,i] <- est_mat[,"ana_prn"]
  nulll[["noise"]][,i] <- est_mat[,"noise"]
}

mean_mat <- data.frame(sct=rowMeans(nulll$sct), prn=rowMeans(nulll$prn),
                       spr=rowMeans(nulll$spr), propr=rowMeans(nulll$propr),
                       cscore_est=rowMeans(nulll$cscore_est),
                       ana_prn=rowMeans(nulll$ana_prn), noise=rowMeans(nulll$noise))
sd_mat <- data.frame(sct=rowSds(nulll$sct), prn=rowSds(nulll$prn),
                     spr=rowSds(nulll$spr), propr=rowSds(nulll$propr),
                     cscore_est=rowSds(nulll$cscore_est),
                     ana_prn=rowSds(nulll$ana_prn), noise=rowSds(nulll$noise))

est_mat <- PNAS_est
p_val <- data.frame(sct=pnorm(abs(est_mat$sct), mean = mean_mat$sct, sd=sd_mat$sct, lower.tail = F)+
                      pnorm(-abs(est_mat$sct), mean = mean_mat$sct, sd=sd_mat$sct, lower.tail = T),
                    prn=pnorm(abs(est_mat$prn), mean = mean_mat$prn, sd=sd_mat$prn, lower.tail = F)+
                      pnorm(-abs(est_mat$prn), mean = mean_mat$prn, sd=sd_mat$prn, lower.tail = T),
                    spr=pnorm(abs(est_mat$spr), mean = mean_mat$spr, sd=sd_mat$spr, lower.tail = F)+
                      pnorm(-abs(est_mat$spr), mean = mean_mat$spr, sd=sd_mat$spr, lower.tail = T),
                    propr=pnorm(abs(est_mat$propr), mean = mean_mat$propr, sd=sd_mat$propr, lower.tail = F)+
                      pnorm(-abs(est_mat$propr), mean = mean_mat$propr, sd=sd_mat$propr, lower.tail = T),
                    cscore_est=pnorm(abs(est_mat$cscore), mean = mean_mat$cscore_est, sd=sd_mat$cscore_est, lower.tail = F)+
                      pnorm(-abs(est_mat$cscore), mean = mean_mat$cscore_est, sd=sd_mat$cscore_est, lower.tail = T),
                    ana_prn=pnorm(abs(est_mat$ana_prn), mean = mean_mat$ana_prn, sd=sd_mat$ana_prn, lower.tail = F)+
                      pnorm(-abs(est_mat$ana_prn), mean = mean_mat$ana_prn, sd=sd_mat$ana_prn, lower.tail = T),
                    noise=pnorm(abs(est_mat$noise), mean = mean_mat$noise, sd=sd_mat$noise, lower.tail = F)+
                      pnorm(-abs(est_mat$noise), mean = mean_mat$noise, sd=sd_mat$noise, lower.tail = T))
saveRDS(p_val, paste0("mean_cor/semi_PD_sparse/simu/PNAS_simu_norm_p.rds"))

PNAS_ori_ests <- readRDS("mean_cor/semi_PD/simu/PNAS_NC_Oli_sct1000.rds")[gene_name, gene_name]
PNAS_ori_ests[abs(PNAS_ori_ests)<0.017] <- 0

p_val$true <- PNAS_ori_ests[upper.tri(PNAS_ori_ests)]

p_val_true <- p_val[p_val$true==0,]
prepare_data <- function(x, idx, name){
  return(data.frame(emp=sort(x),
                    unif=(1:length(idx))/(length(idx)),
                    cor_method=name))
}

set.seed(5112024)
idx <- sample(1:nrow(p_val_true), 10000)
p_dat <- rbind(prepare_data(p_val_true$sct[idx], idx, "sctransform"),
               prepare_data(p_val_true$prn[idx], idx, "Pearson"),
               prepare_data(p_val_true$spr[idx], idx, "Spearman"),
               prepare_data(p_val_true$ana_prn[idx], idx, "Analytic PR"),
               prepare_data(p_val_true$propr[idx], idx, "propr"),
               prepare_data(p_val_true$cscore_est[idx], idx, "CS-CORE"),
               prepare_data(p_val_true$noise[idx], idx, "Noise \nRegularization"))

pdf('mean_cor/semi_PD_sparse/figures/qq_our_pnas.pdf', width = 12, height = 2.5, onefile = T)
ggplot(p_dat, aes(x=unif, y=emp))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, col="red")+
  facet_wrap(.~cor_method, nrow=1)+theme_bw()+
  labs(x="Theoretical Quantiles", y="Sample Quantiles")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14))
dev.off()

