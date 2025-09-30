# reproducibility based if the gene pair is correlated
# /gpfs/gibbs/pi/zhao/xs282/validation/mean_cor/p_value_PNAS_and_ROSMAP_NC_Oli_12_5_2023.R

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(magrittr)
library(rtracklayer)
library(biomaRt)
library(data.table)
library(matrixcalc)
set.seed(11272023)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("AFinal/cscore_real_data_function.R")

ROSMAP_oli_ct <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_sct_cor_NB_simu1000_abs_thresh.rds")
PNAS_oli_ct <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_sct_cor_NB_simu1000_abs_thresh.rds")


ROSMAP_ori_ests <- readRDS("mean_cor/semi_PD/simu/ROSMAP_NC_Oli_sct1000.rds")
ROSMAP_ori_ests[abs(ROSMAP_ori_ests)<0.015] <- 0
gene_name <- rownames(ROSMAP_ori_ests)
mean(ROSMAP_ori_ests!=0)

PNAS_ori_ests <- readRDS("mean_cor/semi_PD/simu/PNAS_NC_Oli_sct1000.rds")
PNAS_ori_ests[abs(PNAS_ori_ests)<0.017] <- 0
PNAS_ori_ests <- PNAS_ori_ests[gene_name, gene_name]
mean(PNAS_ori_ests!=0)

# original reproduce
ori_mat <- data.frame(ROSMAP=ROSMAP_ori_ests[upper.tri(ROSMAP_ori_ests)],
                      PNAS=PNAS_ori_ests[upper.tri(PNAS_ori_ests)])
ori_mat <- as.data.frame(apply(ori_mat, 2, function(x){ifelse(x!=0,1,0)}))
true_reproduce <- sum(ori_mat$PNAS==1 & ori_mat$ROSMAP==1)
sum(ori_mat$ROSMAP>0)
sum(ori_mat$PNAS>0)
true_reproduce_pair <- ori_mat$PNAS==1 & ori_mat$ROSMAP==1

# cor estimations
ROSMAP_cscore <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_cscore1000_abs_thresh.rds")
ROSMAP_cscore_p <- MatrixBH(ROSMAP_cscore$p_value)
ROSMAP_cscore_p <- ROSMAP_cscore_p[gene_name, gene_name]
ROSMAP_cscore_est <- ROSMAP_cscore$est[gene_name, gene_name]

PNAS_cscore <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_cscore1000_abs_thresh.rds")
PNAS_cscore_p <- MatrixBH(PNAS_cscore$p_value)
PNAS_cscore_p <- PNAS_cscore_p[gene_name, gene_name]
PNAS_cscore_est <- PNAS_cscore$est[gene_name, gene_name]

ROSMAP_p <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_simu_norm_p.rds")
ROSMAP_p_adj <- as.data.frame(apply(ROSMAP_p, 2, function(x){p.adjust(x, method = "BH")}))
ROSMAP_p_adj$cscore_p <- ROSMAP_cscore_p[upper.tri(ROSMAP_cscore_p)]
PNAS_p <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_simu_norm_p.rds")
PNAS_p_adj <- as.data.frame(apply(PNAS_p, 2, function(x){p.adjust(x, method = "BH")}))
PNAS_p_adj$cscore_p <- PNAS_cscore_p[upper.tri(PNAS_cscore_p)]
ROSMAP_p$cscore_p <- ROSMAP_cscore$p_value[gene_name, gene_name][upper.tri(ROSMAP_cscore$p_value[gene_name, gene_name])]
PNAS_p$cscore_p <- PNAS_cscore$p_value[gene_name, gene_name][upper.tri(PNAS_cscore$p_value[gene_name, gene_name])]

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

ROSMAP_propr <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_propr1000_abs_thresh.rds")
ROSMAP_propr_est <- ROSMAP_propr@matrix[gene_name, gene_name]
PNAS_propr <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_propr1000_abs_thresh.rds")
PNAS_propr_est <- PNAS_propr@matrix[gene_name, gene_name]

ROSMAP_prn <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_prn1000_abs_thresh.rds")
ROSMAP_prn_est <- ROSMAP_prn[gene_name, gene_name]
PNAS_prn <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_prn1000_abs_thresh.rds")
PNAS_prn_est <- PNAS_prn[gene_name, gene_name]

ROSMAP_spr <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_spr1000_abs_thresh.rds")
ROSMAP_spr_est <- ROSMAP_spr[gene_name, gene_name]
PNAS_spr <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_spr1000_abs_thresh.rds")
PNAS_spr_est <- PNAS_spr[gene_name, gene_name]

ROSMAP_est <- data.frame(sct=ROSMAP_sct_est[upper.tri(ROSMAP_sct_est)],
                         prn=ROSMAP_prn_est[upper.tri(ROSMAP_prn_est)],
                         spr=ROSMAP_spr_est[upper.tri(ROSMAP_spr_est)],
                         propr=ROSMAP_propr_est[upper.tri(ROSMAP_propr_est)],
                         cscore_est=ROSMAP_cscore_est[upper.tri(ROSMAP_cscore_est)],
                         ana_prn=ROSMAP_ana_prn_est[upper.tri(ROSMAP_ana_prn_est)],
                         noise=ROSMAP_noise_est[upper.tri(ROSMAP_noise_est)])

PNAS_est <- data.frame(sct=PNAS_sct_est[upper.tri(PNAS_sct_est)],
                       prn=PNAS_prn_est[upper.tri(PNAS_prn_est)],
                       spr=PNAS_spr_est[upper.tri(PNAS_spr_est)],
                       propr=PNAS_propr_est[upper.tri(PNAS_propr_est)],
                       cscore_est=PNAS_cscore_est[upper.tri(PNAS_cscore_est)],
                       ana_prn=PNAS_ana_prn_est[upper.tri(PNAS_ana_prn_est)],
                       noise=PNAS_noise_est[upper.tri(PNAS_noise_est)])
colnames(ROSMAP_p)


# based on p-value--------------------------------------------------------------
p_cutoff <- c(0.001, 0.005, 0.01, 0.05, 0.1)
total_cor_PNAS <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(total_cor_PNAS) <- colnames(ROSMAP_p_adj)
colnames(total_cor_PNAS) <- p_cutoff

total_cor_ROSMAP <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(total_cor_ROSMAP) <- colnames(ROSMAP_p_adj)
colnames(total_cor_ROSMAP) <- p_cutoff

repduc_p <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(repduc_p) <- colnames(ROSMAP_p_adj)
colnames(repduc_p) <- p_cutoff

repduc_truth_p <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(repduc_truth_p) <- colnames(ROSMAP_p_adj)
colnames(repduc_truth_p) <- p_cutoff

prec_ROSMAP_p <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(prec_ROSMAP_p) <- colnames(ROSMAP_p_adj)
colnames(prec_ROSMAP_p) <- p_cutoff

prec_PNAS_p <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(prec_PNAS_p) <- colnames(ROSMAP_p_adj)
colnames(prec_PNAS_p) <- p_cutoff

for (i in 1:length(p_cutoff)){
  thresh <- p_cutoff[i]
  print(thresh)

  for (j in 1:ncol(ROSMAP_p_adj)){
    ROSMAP_deci <- as.numeric(ROSMAP_p_adj[,j]<thresh)
    PNAS_deci <- as.numeric(PNAS_p_adj[,j]<thresh)
    repduc_p[j,i] <- sum(ROSMAP_deci==1 & PNAS_deci==1)
    repduc_truth_p[j,i] <- sum(ROSMAP_deci==1 & PNAS_deci==1 & ori_mat$PNAS==1 & ori_mat$ROSMAP==1)
    prec_PNAS_p[j,i] <- sum(PNAS_deci==1 & ori_mat$PNAS==1)/sum(PNAS_deci==1)
    prec_ROSMAP_p[j,i] <- sum(ROSMAP_deci==1 & ori_mat$ROSMAP==1)/sum(ROSMAP_deci==1)
    total_cor_PNAS[j,i] <- sum(PNAS_deci==1)
    total_cor_ROSMAP[j,i] <- sum(ROSMAP_deci==1)
  }
}


prec_PNAS_p_long <- melt(prec_PNAS_p)
colnames(prec_PNAS_p_long) <- c("Method", "Top", "PNAS")
prec_ROSMAP_p_long <- melt(prec_ROSMAP_p)
colnames(prec_ROSMAP_p_long) <- c("Method", "Top", "ROSMAP")
repduc_p_long <- melt(repduc_p)
colnames(repduc_p_long) <- c("Method", "Top", "Reproduce")

repduc_prec_overlap_p <- left_join(prec_PNAS_p_long, prec_ROSMAP_p_long, by=c("Method", "Top"))
repduc_prec_overlap_p <- left_join(repduc_prec_overlap_p, repduc_p_long, by=c("Method", "Top"))
repduc_prec_overlap_p$avg <- (repduc_prec_overlap_p$PNAS+repduc_prec_overlap_p$ROSMAP)/2
repduc_prec_overlap_p$Method <- recode(repduc_prec_overlap_p$Method,
                                       sct="sctransform", prn="Pearson", spr="Spearman",
                                       propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                                       noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

color_setting <- c("CS-CORE \n(Empirical)"="brown", "Noise \nRegularization"="#AF58BA",
                   "CS-CORE"="#339933", "sctransform"="#ff6699",
                   "Pearson"="#F28522", "Spearman"="#ffff66","Analytic PR"="#99ccff",
                   "propr"="#3366cc")


inflation_p <- as.data.frame(as.table((repduc_p-repduc_truth_p)/repduc_p))
inflation_p$Var1 <- recode(inflation_p$Var1,
                           sct="sctransform", prn="Pearson", spr="Spearman",
                           propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                           noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")


inflation_p2 <- as.data.frame(as.table(repduc_truth_p))
inflation_p2$group <- "True reproducible pairs"
inflation_p2$Var1 <- recode(inflation_p2$Var1,
                            sct="sctransform", prn="Pearson", spr="Spearman",
                            propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                            noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

format_supp <- theme(text = element_text(size = 14),
                     legend.position="none")
repduc_prec_overlap_p_sub = repduc_prec_overlap_p[repduc_prec_overlap_p$Method=="CS-CORE",]
repduc_scaling_factor_overlap_p <- max(repduc_prec_overlap_p_sub$Reproduce) 
p_repduc_prec_overlap_p_cscore = ggplot(repduc_prec_overlap_p_sub, aes(x = as.factor(Top), group = 1)) +
  geom_col(aes(y = Reproduce), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = ROSMAP * repduc_scaling_factor_overlap_p), color = "darkred", size = 1.5) +
  geom_point(aes(y = ROSMAP * repduc_scaling_factor_overlap_p), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of reproducible pairs",
    sec.axis = sec_axis(trans = ~ . / repduc_scaling_factor_overlap_p, name = "Precision")) +
  labs(title = "P-value",x = "P-value cutoffs") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_repduc_prec_overlap_p_cscore

p_repduc_prec_overlap_p_cscore_PNAS = ggplot(repduc_prec_overlap_p_sub, aes(x = as.factor(Top), group = 1)) +
  geom_col(aes(y = Reproduce), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = PNAS * repduc_scaling_factor_overlap_p), color = "darkred", size = 1.5) +
  geom_point(aes(y = PNAS * repduc_scaling_factor_overlap_p), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of reproducible pairs",
    sec.axis = sec_axis(trans = ~ . / repduc_scaling_factor_overlap_p, name = "Precision (PNAS)")) +
  labs(title = "P-value",x = "P-value cutoffs") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_repduc_prec_overlap_p_cscore_PNAS

repduc_inflation_p_cscore = inflation_p
colnames(repduc_inflation_p_cscore)[3] = "infla"
repduc_inflation_p_cscore = left_join(repduc_inflation_p_cscore, inflation_p2, by=c("Var1", "Var2"))
repduc_inflation_p_cscore = repduc_inflation_p_cscore[repduc_inflation_p_cscore$Var1=="CS-CORE",]
repduc_scaling_factor_infla_p <- max(repduc_inflation_p_cscore$Freq) 

p_repduc_infla_p_cscore = ggplot(repduc_inflation_p_cscore, aes(x = Var2, group = 1)) +
  geom_col(aes(y = Freq), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = infla * repduc_scaling_factor_infla_p), color = "darkred", size = 1.5) +
  geom_point(aes(y = infla * repduc_scaling_factor_infla_p), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of true reproducible pairs",
    sec.axis = sec_axis(trans = ~ . / repduc_scaling_factor_infla_p, name = "Prop of misidentified pairs")) +
  labs(title = "P-value",x = "P-value cutoffs") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_repduc_infla_p_cscore

# based on the unfiltered order instead of the gene cor-------------------------
top_cutoff <- sort(c(1000, 5000, 10000, 50000, 100000, sum(ori_mat$ROSMAP==1), sum(ori_mat$PNAS==1)))
repduc_top <- matrix(NA, nrow=7, ncol=length(top_cutoff))
rownames(repduc_top) <- colnames(ROSMAP_est)
colnames(repduc_top) <- top_cutoff

repduc_truth_top <- matrix(NA, nrow=7, ncol=length(top_cutoff))
rownames(repduc_truth_top) <- colnames(ROSMAP_est)
colnames(repduc_truth_top) <- top_cutoff

total_pair <- ncol(ROSMAP_sct_est)*(ncol(ROSMAP_sct_est)-1)/2

prec_PNAS <- matrix(NA, nrow=7, ncol=length(top_cutoff))
rownames(prec_PNAS) <- colnames(ROSMAP_est)
colnames(prec_PNAS) <- top_cutoff

prec_ROSMAP <- matrix(NA, nrow=7, ncol=length(top_cutoff))
rownames(prec_ROSMAP) <- colnames(ROSMAP_est)
colnames(prec_ROSMAP) <- top_cutoff

for (i in 1:length(top_cutoff)){
  thresh <- top_cutoff[i]
  print(thresh)
  for (j in 1:ncol(ROSMAP_est)){
    cor_ROSMAP <- abs(ROSMAP_est[,j])
    deci_ROSMAP <- as.numeric(cor_ROSMAP>quantile(cor_ROSMAP, 1-thresh/total_pair))
    cor_PNAS <- abs(PNAS_est[,j])
    deci_PNAS <- as.numeric(cor_PNAS>quantile(cor_PNAS, 1-thresh/total_pair))

    repduc_top[j,i] <- sum(deci_ROSMAP==1 & deci_PNAS==1)
    repduc_truth_top[j,i] <- sum(deci_ROSMAP==1 & deci_PNAS==1 & ori_mat$PNAS==1 & ori_mat$ROSMAP==1)
    prec_PNAS[j,i] <- sum(deci_PNAS==1 & ori_mat$PNAS==1)/sum(deci_PNAS==1)
    prec_ROSMAP[j,i] <- sum(deci_ROSMAP==1 & ori_mat$ROSMAP==1)/sum(deci_ROSMAP==1)
  }
}


prec_PNAS_long <- melt(prec_PNAS)
colnames(prec_PNAS_long) <- c("Method", "Top", "PNAS")
prec_ROSMAP_long <- melt(prec_ROSMAP)
colnames(prec_ROSMAP_long) <- c("Method", "Top", "ROSMAP")
repduc_top_long <- melt(repduc_top)
colnames(repduc_top_long) <- c("Method", "Top", "Reproduce")

repduc_prec_overlap <- left_join(prec_PNAS_long, prec_ROSMAP_long, by=c("Method", "Top"))
repduc_prec_overlap <- left_join(repduc_prec_overlap, repduc_top_long, by=c("Method", "Top"))
repduc_prec_overlap$Method <- recode(repduc_prec_overlap$Method,
                            sct="sctransform", prn="Pearson", spr="Spearman",
                            propr="propr",ana_prn="Analytic PR",
                            noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")


inflation_unfil_top <- as.data.frame(as.table((repduc_top-repduc_truth_top)/repduc_top))
inflation_unfil_top$Var1 <- recode(inflation_unfil_top$Var1,
                           sct="sctransform", prn="Pearson", spr="Spearman",
                           propr="propr",ana_prn="Analytic PR",
                           noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

inflation_unfil_top2 <- as.data.frame(as.table(repduc_truth_top))
inflation_unfil_top2$group <- "True reproducible pairs"
inflation_unfil_top2$Var1 <- recode(inflation_unfil_top2$Var1,
                            sct="sctransform", prn="Pearson", spr="Spearman",
                            propr="propr",ana_prn="Analytic PR",
                            noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")


repduc_prec_overlap_sub = repduc_prec_overlap[repduc_prec_overlap$Method=="CS-CORE \n(Empirical)",]
repduc_scaling_factor_overlap <- max(repduc_prec_overlap_sub$Reproduce) 
p_repduc_prec_overlap_cscore = ggplot(repduc_prec_overlap_sub, aes(x = as.factor(Top), group = 1)) +
  geom_col(aes(y = Reproduce), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = ROSMAP * repduc_scaling_factor_overlap), color = "darkred", size = 1.5) +
  geom_point(aes(y = ROSMAP * repduc_scaling_factor_overlap), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of reproducible pairs",
    sec.axis = sec_axis(trans = ~ . / repduc_scaling_factor_overlap, name = "Precision")) +
  labs(title = "Correlation strength",x = "Top") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_repduc_prec_overlap_cscore

p_repduc_prec_overlap_cscore_PNAS = ggplot(repduc_prec_overlap_sub, aes(x = as.factor(Top), group = 1)) +
  geom_col(aes(y = Reproduce), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = PNAS * repduc_scaling_factor_overlap), color = "darkred", size = 1.5) +
  geom_point(aes(y = PNAS * repduc_scaling_factor_overlap), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of reproducible pairs",
    sec.axis = sec_axis(trans = ~ . / repduc_scaling_factor_overlap, name = "Precision (PNAS)")) +
  labs(title = "Correlation strength",x = "Top") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_repduc_prec_overlap_cscore_PNAS    

repduc_inflation_cscore = inflation_unfil_top
colnames(repduc_inflation_cscore)[3] = "infla"
repduc_inflation_cscore = left_join(repduc_inflation_cscore, inflation_unfil_top2, by=c("Var1", "Var2"))
repduc_inflation_cscore = repduc_inflation_cscore[repduc_inflation_cscore$Var1=="CS-CORE \n(Empirical)",]
repduc_scaling_factor_infla <- max(repduc_inflation_cscore$Freq) 

p_repduc_infla_cscore = ggplot(repduc_inflation_cscore, aes(x = Var2, group = 1)) +
  geom_col(aes(y = Freq), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = infla * repduc_scaling_factor_infla), color = "darkred", size = 1.5) +
  geom_point(aes(y = infla * repduc_scaling_factor_infla), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of true reproducible pairs",
    sec.axis = sec_axis(trans = ~ . / repduc_scaling_factor_infla, name = "Prop of misidentified pairs")) +
  labs(title = "Correlation strength",x = "Top") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_repduc_infla_cscore




# fixed mis prop --------------------------------------------------------------
# based on p-value--------------------------------------------------------------
p_cutoff <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
total_cor_PNAS <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(total_cor_PNAS) <- colnames(ROSMAP_p_adj)
colnames(total_cor_PNAS) <- p_cutoff

total_cor_ROSMAP <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(total_cor_ROSMAP) <- colnames(ROSMAP_p_adj)
colnames(total_cor_ROSMAP) <- p_cutoff

repduc_p <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(repduc_p) <- colnames(ROSMAP_p_adj)
colnames(repduc_p) <- p_cutoff

repduc_truth_p <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(repduc_truth_p) <- colnames(ROSMAP_p_adj)
colnames(repduc_truth_p) <- p_cutoff

prec_ROSMAP_p <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(prec_ROSMAP_p) <- colnames(ROSMAP_p_adj)
colnames(prec_ROSMAP_p) <- p_cutoff

prec_PNAS_p <- matrix(NA, nrow=8, ncol=length(p_cutoff))
rownames(prec_PNAS_p) <- colnames(ROSMAP_p_adj)
colnames(prec_PNAS_p) <- p_cutoff

for (i in 1:length(p_cutoff)){
  thresh <- p_cutoff[i]
  print(thresh)

  for (j in 1:ncol(ROSMAP_p_adj)){
    ROSMAP_deci <- as.numeric(ROSMAP_p_adj[,j]<thresh)
    PNAS_deci <- as.numeric(PNAS_p_adj[,j]<thresh)
    repduc_p[j,i] <- sum(ROSMAP_deci==1 & PNAS_deci==1)
    repduc_truth_p[j,i] <- sum(ROSMAP_deci==1 & PNAS_deci==1 & ori_mat$PNAS==1 & ori_mat$ROSMAP==1)
    prec_PNAS_p[j,i] <- sum(PNAS_deci==1 & ori_mat$PNAS==1)/sum(PNAS_deci==1)
    prec_ROSMAP_p[j,i] <- sum(ROSMAP_deci==1 & ori_mat$ROSMAP==1)/sum(ROSMAP_deci==1)
    total_cor_PNAS[j,i] <- sum(PNAS_deci==1)
    total_cor_ROSMAP[j,i] <- sum(ROSMAP_deci==1)
  }
}


repduc_p_long <- reshape2::melt(repduc_truth_p)
colnames(repduc_p_long) <- c("Method", "Top", "Reproduce")

inflation_p <- as.data.frame(as.table((repduc_p-repduc_truth_p)/repduc_p))
colnames(inflation_p) <- c("Method", "cutoff", "Mis")
inflation_p$true <- repduc_p_long$Reproduce


# based on the unfiltered order instead of the gene cor-------------------------
top_cutoff <- sort(c(seq(1000, 10000, by=1000), 20000,25000,30000,35000,40000,45000,50000))
repduc_top <- matrix(NA, nrow=7, ncol=length(top_cutoff))
rownames(repduc_top) <- colnames(ROSMAP_est)
colnames(repduc_top) <- top_cutoff

repduc_truth_top <- matrix(NA, nrow=7, ncol=length(top_cutoff))
rownames(repduc_truth_top) <- colnames(ROSMAP_est)
colnames(repduc_truth_top) <- top_cutoff

total_pair <- ncol(ROSMAP_sct_est)*(ncol(ROSMAP_sct_est)-1)/2

prec_PNAS <- matrix(NA, nrow=7, ncol=length(top_cutoff))
rownames(prec_PNAS) <- colnames(ROSMAP_est)
colnames(prec_PNAS) <- top_cutoff

prec_ROSMAP <- matrix(NA, nrow=7, ncol=length(top_cutoff))
rownames(prec_ROSMAP) <- colnames(ROSMAP_est)
colnames(prec_ROSMAP) <- top_cutoff

for (i in 1:length(top_cutoff)){
  thresh <- top_cutoff[i]
  print(thresh)
  for (j in 1:ncol(ROSMAP_est)){
    cor_ROSMAP <- abs(ROSMAP_est[,j])
    deci_ROSMAP <- as.numeric(cor_ROSMAP>quantile(cor_ROSMAP, 1-thresh/total_pair))
    cor_PNAS <- abs(PNAS_est[,j])
    deci_PNAS <- as.numeric(cor_PNAS>quantile(cor_PNAS, 1-thresh/total_pair))

    repduc_top[j,i] <- sum(deci_ROSMAP==1 & deci_PNAS==1)
    repduc_truth_top[j,i] <- sum(deci_ROSMAP==1 & deci_PNAS==1 & ori_mat$PNAS==1 & ori_mat$ROSMAP==1)
    prec_PNAS[j,i] <- sum(deci_PNAS==1 & ori_mat$PNAS==1)/sum(deci_PNAS==1)
    prec_ROSMAP[j,i] <- sum(deci_ROSMAP==1 & ori_mat$ROSMAP==1)/sum(deci_ROSMAP==1)
  }
}


repduc_top_long <- melt(repduc_truth_top)
colnames(repduc_top_long) <- c("Method", "Top", "Reproduce")

inflation_unfil_top <- as.data.frame(as.table((repduc_top-repduc_truth_top)/repduc_top))
colnames(inflation_unfil_top) <- c("Method", "cutoff", "Mis")
inflation_unfil_top$true <- repduc_top_long$Reproduce
inflation_unfil_top_sub <- inflation_unfil_top[inflation_unfil_top$Method=="cscore_est",]
inflation_unfil_top_sub$Method <- "cscore_p"
inflation_unfil_top <- rbind(inflation_unfil_top_sub, inflation_unfil_top)

inflation_unfil_top$Method <- recode(inflation_unfil_top$Method,
                            sct="sctransform", prn="Pearson", spr="Spearman",
                            propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                            noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")
inflation_p$Method <- recode(inflation_p$Method,
                             sct="sctransform", prn="Pearson", spr="Spearman",
                             propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                             noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")


i = "CS-CORE"
plot_dat1 <- inflation_p[inflation_p$Method==i,]
plot_dat1$group <- "P-value"
plot_dat2 <- inflation_unfil_top[inflation_unfil_top$Method==i,]
plot_dat2$group <- "Cor-strength"
plot_dat <- rbind(plot_dat1, plot_dat2)
fix_mis_cscore <- ggplot(plot_dat, aes(x=Mis, y=true, color=group))+
    geom_point(size=2)+geom_line(size=1)+labs( x="Prop of \nmisidentified pairs", y="# of true reproducible pairs", color="")+
    theme_bw()+guides(color = guide_legend(nrow = 2)) +
    # scale_colour_manual(values = c("P-value"="darkblue", "Cor-strength"="darkred"))+
    theme(text = element_text(size = 14),
         plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1),
        #  legend.position = c(0.98,1.1),
        # legend.justification = c("right", "top"),
        # legend.background = element_blank()
         legend.position = "top", legend.box.margin = margin(b = -10))
fix_mis_cscore

plot_rep_ori <- ggarrange(p_repduc_prec_overlap_cscore, p_repduc_prec_overlap_p_cscore,
                          p_repduc_infla_cscore, p_repduc_infla_p_cscore, fix_mis_cscore, ncol=5,nrow=1,
                          widths = c(1, 0.9, 1, 0.9,0.8),
          labels = c("A", "B", "C", "D","E"))
p_reproduce = annotate_figure(ggarrange(plot_rep_ori, ncol=1, nrow=1),
  top = text_grob("Reproducibility", color = "black", face = "bold", size = 16))

##############################################################################################################
# overlap with string ----------------------------------------------------------------------------------------
##############################################################################################################
# prepare biological network----------------------------------------------------
# STRING
hs_filter <- readRDS("STRING/hs_filter_10_4_2023.rds")
hs_filter <- hs_filter %>%
  mutate(grp = paste(pmax(protein1, protein2), pmin(protein1, protein2), sep = "_"))
ROSMAP_p <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_simu_norm_p.rds")

# cor estimations
ROSMAP_cscore_est_filter <- ROSMAP_cscore_est
ROSMAP_cscore_est_filter[ROSMAP_cscore_p >= 0.05] <- 0

upper2matrix <- function(est_mat, col_name){
  p_mat <- est_mat-est_mat
  p_mat[upper.tri(p_mat)] <- ROSMAP_p_adj[,col_name]
  p_mat <- p_mat + t(p_mat)
  return(p_mat)
}

ROSMAP_sct_p <- upper2matrix(ROSMAP_sct_est, "sct")
ROSMAP_sct_est_filter <- ROSMAP_sct_est
ROSMAP_sct_est_filter[ROSMAP_sct_p >= 0.05] <- 0

ROSMAP_ana_prn_p <- upper2matrix(ROSMAP_ana_prn_est, "ana_prn")
ROSMAP_ana_prn_est_filter <- ROSMAP_ana_prn_est
ROSMAP_ana_prn_est_filter[ROSMAP_ana_prn_p >= 0.05] <- 0

ROSMAP_noise_est <- apply(ROSMAP_noise_est, c(1,2), as.numeric)
ROSMAP_noise_p <- upper2matrix(ROSMAP_noise_est, "noise")
ROSMAP_noise_est_filter <- ROSMAP_noise_est
ROSMAP_noise_est_filter[ROSMAP_noise_p >= 0.05] <- 0

ROSMAP_propr_p <- upper2matrix(ROSMAP_propr_est, "propr")
ROSMAP_propr_est_filter <- ROSMAP_propr_est
ROSMAP_propr_est_filter[ROSMAP_propr_p >= 0.05] <- 0

ROSMAP_prn_p <- upper2matrix(ROSMAP_prn_est, "prn")
ROSMAP_prn_est_filter <- ROSMAP_prn_est
ROSMAP_prn_est_filter[ROSMAP_prn_p >= 0.05] <- 0

ROSMAP_spr_p <- upper2matrix(ROSMAP_spr_est, "spr")
ROSMAP_spr_est_filter <- ROSMAP_spr_est
ROSMAP_spr_est_filter[ROSMAP_spr_p >= 0.05] <- 0

ROSMAP_cscore_simu_p <- upper2matrix(ROSMAP_cscore_est, "cscore_est")
ROSMAP_cscore_est_filter_simu <- ROSMAP_cscore_est
ROSMAP_cscore_est_filter_simu[ROSMAP_cscore_simu_p >= 0.05] <- 0

#-----------------------------------------------------------
marginal_fit_ROSMAP = readRDS('marginal_fit/ROSMAP_NC_Oli_marginal_fit.rds')
mu_ROSMAP <- marginal_fit_ROSMAP[gene_name,]$mu

ncor_gene <- length(gene_name)
mu_col_ROSMAP <- matrix(mu_ROSMAP,ncor_gene,ncor_gene,byrow = T)
mu_row_ROSMAP <- matrix(mu_ROSMAP,ncor_gene,ncor_gene)

# filtered
tri = upper.tri(ROSMAP_ori_ests, diag = FALSE)
idxs = which(tri, arr.ind = T)
estimate_ROSMAP <- data.frame(id1=rownames(ROSMAP_ori_ests)[idxs[,1]],
                              id2=colnames(ROSMAP_ori_ests)[idxs[,2]],
                              ROSMAP_ori=ROSMAP_ori_ests[tri],
                              mu_col_ROSMAP=mu_col_ROSMAP[tri],
                              mu_row_ROSMAP=mu_row_ROSMAP[tri],
                              ROSMAP_ori_ests=ROSMAP_ori_ests[tri])
estimate_ROSMAP$log10mean_mu_ROSMAP <- log10(sqrt(10^estimate_ROSMAP$mu_col_ROSMAP*10^estimate_ROSMAP$mu_row_ROSMAP))
estimate_ROSMAP$true_cor <- ifelse(abs(estimate_ROSMAP$ROSMAP_ori_ests)!=0, 1, 0)


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=rownames(ROSMAP_ori_ests),mart= mart)
G_list <- G_list %>% group_by(hgnc_symbol) %>%
  dplyr::slice(1) %>% ungroup()
unmapped <- rownames(ROSMAP_ori_ests)[!rownames(ROSMAP_ori_ests) %in% G_list$hgnc_symbol]
ensembl <- c(unmapped, G_list$ensembl_gene_id)
names(ensembl) <- c(unmapped, G_list$hgnc_symbol)
estimate_ROSMAP$id1 <- ensembl[estimate_ROSMAP$id1]
estimate_ROSMAP$id2 <- ensembl[estimate_ROSMAP$id2]
estimate_ROSMAP <- estimate_ROSMAP %>%
  mutate(grp = paste(pmax(id1, id2), pmin(id1, id2), sep = "_"))

estimate_ROSMAP_filter <- estimate_ROSMAP

estimate_ROSMAP$cscore_est <- ROSMAP_cscore_est[tri]
estimate_ROSMAP$sct <- ROSMAP_sct_est[tri]
estimate_ROSMAP$ana_prn <- ROSMAP_ana_prn_est[tri]
estimate_ROSMAP$noise <- ROSMAP_noise_est[tri]
estimate_ROSMAP$propr <- ROSMAP_propr_est[tri]
estimate_ROSMAP$prn <- ROSMAP_prn_est[tri]
estimate_ROSMAP$spr <- ROSMAP_spr_est[tri]
estimate_ROSMAP$cscore_p <- estimate_ROSMAP$cscore_est

estimate_ROSMAP_filter$cscore_p <- ROSMAP_cscore_est_filter[tri]
estimate_ROSMAP_filter$sct <- ROSMAP_sct_est_filter[tri]
estimate_ROSMAP_filter$ana_prn <- ROSMAP_ana_prn_est_filter[tri]
estimate_ROSMAP_filter$noise <- ROSMAP_noise_est_filter[tri]
estimate_ROSMAP_filter$propr <- ROSMAP_propr_est_filter[tri]
estimate_ROSMAP_filter$prn <- ROSMAP_prn_est_filter[tri]
estimate_ROSMAP_filter$spr <- ROSMAP_spr_est_filter[tri]
estimate_ROSMAP_filter$cscore_est <- ROSMAP_cscore_est_filter_simu[tri]

estimate_ROSMAP_p_adj <- estimate_ROSMAP_filter
estimate_ROSMAP_p_adj[,colnames(ROSMAP_p_adj)] <- ROSMAP_p_adj
estimate_ROSMAP_p_adj$cscore_p <- ROSMAP_cscore_p[tri]


# number of total overlap ------------------------------------------------------
biological_net <- hs_filter[hs_filter$combined_score>500,]

## overlap with biological network with different top threshold -----------------
top_cutoff <- c(1000,5000, 10000, 50000,100000, sum(estimate_ROSMAP$true_cor))
overlap_string <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_string) <- top_cutoff
rownames(overlap_string) <- colnames(ROSMAP_p)

overlap_true <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_true) <- top_cutoff
rownames(overlap_true) <- colnames(ROSMAP_p)

overlap_prec <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_prec) <- top_cutoff
rownames(overlap_prec) <- colnames(ROSMAP_p)


total_pair <- nrow(estimate_ROSMAP)
for (i in 1:length(top_cutoff)){
  thresh <- top_cutoff[i]
  print(thresh)

  selected_genes_thresh <- list()
  for (j in rownames(overlap_string)){
    cor_ROSMAP <- abs(estimate_ROSMAP[,j])
    deci_ROSMAP <- estimate_ROSMAP[cor_ROSMAP>quantile(cor_ROSMAP, 1-thresh/total_pair),]
    overlap_string[j,i] <- sum(deci_ROSMAP$grp %in% biological_net$grp)
    overlap_true[j,i] <- sum((deci_ROSMAP$grp %in% biological_net$grp)&deci_ROSMAP$true_cor==1)
    overlap_prec[j,i] <- sum(deci_ROSMAP$true_cor==1)/nrow(deci_ROSMAP)
  }
}

overlap_prec_long <- reshape2::melt(overlap_prec)
colnames(overlap_prec_long) <- c("Method", "Top", "Precision")
overlap_string_long <- reshape2::melt(overlap_string)
colnames(overlap_string_long) <- c("Method", "Top", "Overlap")

overlap_prec_overlap <- left_join(overlap_prec_long, overlap_string_long, by=c("Method", "Top"))
overlap_prec_overlap$Top <- as.factor(overlap_prec_overlap$Top)
overlap_prec_overlap$Group <- "Original"
overlap_prec_overlap$Method <- recode(overlap_prec_overlap$Method,
                                           sct="sctransform", prn="Pearson", spr="Spearman",
                                           propr="propr",ana_prn="Analytic PR",
                                           noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")


inflation_unfil_top <- as.data.frame(as.table((overlap_string-overlap_true)/overlap_string))
inflation_unfil_top$Var1 <- recode(inflation_unfil_top$Var1,
                                      sct="sctransform", prn="Pearson", spr="Spearman",
                                      propr="propr",ana_prn="Analytic PR",
                                      noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

inflation_unfil_top2 <- as.data.frame(as.table(overlap_true))
inflation_unfil_top2$group <- "True overlaps"
inflation_unfil_top2$Var1 <- recode(inflation_unfil_top2$Var1,
                                    sct="sctransform", prn="Pearson", spr="Spearman",
                                    propr="propr",ana_prn="Analytic PR",
                                    noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")



## overlap with biological network with different p-value threshold --------------------
p_cutoff <- c(0.001, 0.005, 0.01, 0.05, 0.1)
overlap_string_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_string_p) <- p_cutoff
rownames(overlap_string_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_true_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_true_p) <- p_cutoff
rownames(overlap_true_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_prec_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_prec_p) <- p_cutoff
rownames(overlap_prec_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_cor_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_cor_p) <- p_cutoff
rownames(overlap_cor_p) <- c(colnames(ROSMAP_p), "cscore_p")


estimate_p <- estimate_ROSMAP[,1:9]
estimate_p <- cbind(estimate_p, ROSMAP_p_adj)
estimate_p$cscore_p <- ROSMAP_cscore_p[upper.tri(ROSMAP_cscore_p)]
total_pair <- nrow(estimate_p)

for (i in 1:length(p_cutoff)){
  thresh <- p_cutoff[i]
  print(thresh)

  for (j in rownames(overlap_true_p)){
    ROSMAP_deci <- estimate_p[estimate_p[,j]<thresh,]
    overlap_string_p[j,i] <- sum(ROSMAP_deci$grp %in% biological_net$grp)
    overlap_true_p[j,i] <- sum((ROSMAP_deci$grp %in% biological_net$grp)&ROSMAP_deci$true_cor==1)
    overlap_prec_p[j,i] <- sum(ROSMAP_deci$true_cor==1)/nrow(ROSMAP_deci)
    overlap_cor_p[j,i] <- nrow(ROSMAP_deci)
  }
}


overlap_prec_p_long <- melt(overlap_prec_p)
colnames(overlap_prec_p_long) <- c("Method", "Top", "Precision")
overlap_string_p_long <- melt(overlap_string_p)
colnames(overlap_string_p_long) <- c("Method", "Top", "Overlap")

overlap_prec_overlap_p <- left_join(overlap_prec_p_long, overlap_string_p_long, by=c("Method", "Top"))
overlap_prec_overlap_p$Top <- as.factor(overlap_prec_overlap_p$Top)
overlap_prec_overlap_p$Method <- recode(overlap_prec_overlap_p$Method,
                                       sct="sctransform", prn="Pearson", spr="Spearman",
                                       propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                                       noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

inflation_p <- as.data.frame(as.table((overlap_string_p-overlap_true_p)/overlap_string_p))
inflation_p$Var1 <- recode(inflation_p$Var1,
                                 sct="sctransform", prn="Pearson", spr="Spearman",
                                 propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                                 noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

inflation_p2 <- as.data.frame(as.table(overlap_true_p))
inflation_p2$group <- "True overlaps"
inflation_p2$Var1 <- recode(inflation_p2$Var1,
                                  sct="sctransform", prn="Pearson", spr="Spearman",
                                  propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                                  noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

format_supp <- theme(text = element_text(size = 14),
                     legend.position="none")
string_prec_overlap_p_sub = overlap_prec_overlap_p[overlap_prec_overlap_p$Method=="CS-CORE",]
string_scaling_factor_overlap_p <- max(string_prec_overlap_p_sub$Overlap) 
p_string_prec_overlap_p_cscore = ggplot(string_prec_overlap_p_sub, aes(x = as.factor(Top), group = 1)) +
  geom_col(aes(y = Overlap), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = Precision * string_scaling_factor_overlap_p), color = "darkred", size = 1.5) +
  geom_point(aes(y = Precision * string_scaling_factor_overlap_p), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of overlaps with STRING",
    sec.axis = sec_axis(trans = ~ . / string_scaling_factor_overlap_p, name = "Precision")) +
  labs(title = "P-value",x = "P-value cutoffs") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_string_prec_overlap_p_cscore

string_inflation_p_cscore = inflation_p
colnames(string_inflation_p_cscore)[3] = "infla"
string_inflation_p_cscore = left_join(string_inflation_p_cscore, inflation_p2, by=c("Var1", "Var2"))
string_inflation_p_cscore = string_inflation_p_cscore[string_inflation_p_cscore$Var1=="CS-CORE",]
string_scaling_factor_infla_p <- max(string_inflation_p_cscore$Freq) 

p_string_infla_p_cscore = ggplot(string_inflation_p_cscore, aes(x = Var2, group = 1)) +
  geom_col(aes(y = Freq), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = infla * string_scaling_factor_infla_p), color = "darkred", size = 1.5) +
  geom_point(aes(y = infla * string_scaling_factor_infla_p), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of true overlaps with STRING",
    sec.axis = sec_axis(trans = ~ . / string_scaling_factor_infla_p, name = "Prop of misidentified overlaps")) +
  labs(title = "P-value",x = "P-value cutoffs") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_string_infla_p_cscore


string_prec_overlap_sub = overlap_prec_overlap[overlap_prec_overlap$Method=="CS-CORE \n(Empirical)",]
string_scaling_factor_overlap <- max(string_prec_overlap_sub$Overlap) 
p_string_prec_overlap_cscore = ggplot(string_prec_overlap_sub, aes(x = as.factor(Top), group = 1)) +
  geom_col(aes(y = Overlap), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = Precision * string_scaling_factor_overlap), color = "darkred", size = 1.5) +
  geom_point(aes(y = Precision * string_scaling_factor_overlap), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of overlaps with STRING",
    sec.axis = sec_axis(trans = ~ . / string_scaling_factor_overlap, name = "Precision")) +
  labs(title = "Correlation strength",x = "Top") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_string_prec_overlap_cscore

string_inflation_cscore = inflation_unfil_top
colnames(string_inflation_cscore)[3] = "infla"
string_inflation_cscore = left_join(string_inflation_cscore, inflation_unfil_top2, by=c("Var1", "Var2"))
string_inflation_cscore = string_inflation_cscore[string_inflation_cscore$Var1=="CS-CORE \n(Empirical)",]
string_scaling_factor_infla <- max(string_inflation_cscore$Freq) 

p_string_infla_cscore = ggplot(string_inflation_cscore, aes(x = Var2, group = 1)) +
  geom_col(aes(y = Freq), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = infla * string_scaling_factor_infla), color = "darkred", size = 1.5) +
  geom_point(aes(y = infla * string_scaling_factor_infla), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of true overlaps with STRING",
    sec.axis = sec_axis(trans = ~ . / string_scaling_factor_infla, name = "Prop of misidentified overlaps")) +
  labs(title = "Correlation strength",x = "Top") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_string_infla_cscore


# fixed mis prop-------------------------------------------------------------------
## overlap with biological network with different top threshold -----------------
top_cutoff <- sort(c(seq(1000, 10000, by=1000), 20000,25000,30000,35000,40000,45000,50000))
overlap_string <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_string) <- top_cutoff
rownames(overlap_string) <- colnames(ROSMAP_p)

overlap_true <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_true) <- top_cutoff
rownames(overlap_true) <- colnames(ROSMAP_p)

overlap_prec <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_prec) <- top_cutoff
rownames(overlap_prec) <- colnames(ROSMAP_p)


total_pair <- nrow(estimate_ROSMAP)
for (i in 1:length(top_cutoff)){
  thresh <- top_cutoff[i]
  print(thresh)

  selected_genes_thresh <- list()
  for (j in rownames(overlap_string)){
    cor_ROSMAP <- abs(estimate_ROSMAP[,j])
    deci_ROSMAP <- estimate_ROSMAP[cor_ROSMAP>quantile(cor_ROSMAP, 1-thresh/total_pair),]
    overlap_string[j,i] <- sum(deci_ROSMAP$grp %in% biological_net$grp)
    overlap_true[j,i] <- sum((deci_ROSMAP$grp %in% biological_net$grp)&deci_ROSMAP$true_cor==1)
    overlap_prec[j,i] <- sum(deci_ROSMAP$true_cor==1)/nrow(deci_ROSMAP)
  }
}

overlap_string_long <- reshape2::melt(overlap_true)
colnames(overlap_string_long) <- c("Method", "Top", "Overlap")

inflation_unfil_top <- as.data.frame(as.table((overlap_string-overlap_true)/overlap_string))
colnames(inflation_unfil_top) <- c("Method", "cutoff", "Mis")
inflation_unfil_top$true <- overlap_string_long$Overlap
inflation_unfil_top_sub <- inflation_unfil_top[inflation_unfil_top$Method=="cscore_est",]
inflation_unfil_top_sub$Method <- "cscore_p"
inflation_unfil_top <- rbind(inflation_unfil_top_sub, inflation_unfil_top)

inflation_unfil_top$Method <- recode(inflation_unfil_top$Method,
                                     sct="sctransform", prn="Pearson", spr="Spearman",
                                     propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                                     noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

## overlap with biological network with different p-value threshold --------------------
p_cutoff <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
overlap_string_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_string_p) <- p_cutoff
rownames(overlap_string_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_true_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_true_p) <- p_cutoff
rownames(overlap_true_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_prec_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_prec_p) <- p_cutoff
rownames(overlap_prec_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_cor_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_cor_p) <- p_cutoff
rownames(overlap_cor_p) <- c(colnames(ROSMAP_p), "cscore_p")

estimate_p <- estimate_ROSMAP[,1:9]
estimate_p <- cbind(estimate_p, ROSMAP_p_adj)
estimate_p$cscore_p <- ROSMAP_cscore_p[upper.tri(ROSMAP_cscore_p)]
total_pair <- nrow(estimate_p)

for (i in 1:length(p_cutoff)){
  thresh <- p_cutoff[i]
  print(thresh)

  for (j in rownames(overlap_true_p)){
    ROSMAP_deci <- estimate_p[estimate_p[,j]<thresh,]
    overlap_string_p[j,i] <- sum(ROSMAP_deci$grp %in% biological_net$grp)
    overlap_true_p[j,i] <- sum((ROSMAP_deci$grp %in% biological_net$grp)&ROSMAP_deci$true_cor==1)
    overlap_prec_p[j,i] <- sum(ROSMAP_deci$true_cor==1)/nrow(ROSMAP_deci)
    overlap_cor_p[j,i] <- nrow(ROSMAP_deci)
  }
}


overlap_string_p_long <- melt(overlap_true_p)
colnames(overlap_string_p_long) <- c("Method", "Top", "Overlap")

inflation_p <- as.data.frame(as.table((overlap_string_p-overlap_true_p)/overlap_string_p))
colnames(inflation_p) <- c("Method", "cutoff", "Mis")
inflation_p$true <- overlap_string_p_long$Overlap
inflation_p$Method <- recode(inflation_p$Method,
                             sct="sctransform", prn="Pearson", spr="Spearman",
                             propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                             noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")


i = "CS-CORE"
plot_dat1 <- inflation_p[inflation_p$Method==i,]
plot_dat1$group <- "P-value"
plot_dat2 <- inflation_unfil_top[inflation_unfil_top$Method==i,]
plot_dat2$group <- "Cor-strength"
plot_dat_string <- rbind(plot_dat1, plot_dat2)
fix_mis_cscore_string <- ggplot(plot_dat_string, aes(x=Mis, y=true, color=group))+
    geom_point(size=2)+geom_line(size=1)+labs( x="Prop of \nmisidentified overlaps", y="# of true overlaps with STRING", color="")+
    theme_bw()+guides(color = guide_legend(nrow = 2)) +
    # scale_colour_manual(values = c("P-value"="darkblue", "Cor-strength"="darkred"))+
    theme(text = element_text(size = 14),
         plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1),
        #  legend.position = c(0.98,1.1),
        # legend.justification = c("right", "top"),
        # legend.background = element_blank()
         legend.position = "top", legend.box.margin = margin(b = -10))+xlim(0,0.17)
fix_mis_cscore_string

plot_string_ori <- ggarrange(p_string_prec_overlap_cscore, p_string_prec_overlap_p_cscore,
                             p_string_infla_cscore, 
                             p_string_infla_p_cscore, fix_mis_cscore_string, ncol=5,nrow=1,
                          widths = c(1, 0.9, 1, 0.9,0.8),
          labels = c("F", "G", "H", "I","J"))
p_string = annotate_figure(ggarrange(plot_string_ori, ncol=1, nrow=1),
  top = text_grob("-------------------------------------------------------------------------------------------------------\nOverlap with STRING", 
                  color = "black", face = "bold", size = 16))

pdf('/gpfs/gibbs/pi/zhao/xs282/validation/revision/modify_plot/plot_cscore.pdf', width = 14, height = 8, onefile = T)
ggarrange(p_reproduce, p_string, ncol=1, nrow=2, heights = c(5,5.5))
dev.off()


##############################################################################################################
# overlap with reactome ----------------------------------------------------------------------------------------
##############################################################################################################
# Reactome
hs_filter <- readRDS("Reactome/Processed_Reactome_10_07_2023.rds")
hs_filter <- hs_filter %>%
  mutate(grp = paste(pmax(V1, V2), pmin(V1, V2), sep = "_"))
# number of total overlap ------------------------------------------------------
biological_net <- hs_filter

## overlap with biological network with different top threshold -----------------
top_cutoff <- c(1000,5000, 10000, 50000,100000, sum(estimate_ROSMAP$true_cor))
overlap_string <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_string) <- top_cutoff
rownames(overlap_string) <- colnames(ROSMAP_p)

overlap_true <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_true) <- top_cutoff
rownames(overlap_true) <- colnames(ROSMAP_p)

overlap_prec <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_prec) <- top_cutoff
rownames(overlap_prec) <- colnames(ROSMAP_p)


total_pair <- nrow(estimate_ROSMAP)
for (i in 1:length(top_cutoff)){
  thresh <- top_cutoff[i]
  print(thresh)

  selected_genes_thresh <- list()
  for (j in rownames(overlap_string)){
    cor_ROSMAP <- abs(estimate_ROSMAP[,j])
    deci_ROSMAP <- estimate_ROSMAP[cor_ROSMAP>quantile(cor_ROSMAP, 1-thresh/total_pair),]
    overlap_string[j,i] <- sum(deci_ROSMAP$grp %in% biological_net$grp)
    overlap_true[j,i] <- sum((deci_ROSMAP$grp %in% biological_net$grp)&deci_ROSMAP$true_cor==1)
    overlap_prec[j,i] <- sum(deci_ROSMAP$true_cor==1)/nrow(deci_ROSMAP)
  }
}

overlap_prec_long <- reshape2::melt(overlap_prec)
colnames(overlap_prec_long) <- c("Method", "Top", "Precision")
overlap_string_long <- reshape2::melt(overlap_string)
colnames(overlap_string_long) <- c("Method", "Top", "Overlap")


overlap_prec_overlap <- left_join(overlap_prec_long, overlap_string_long, by=c("Method", "Top"))
overlap_prec_overlap$Top <- as.factor(overlap_prec_overlap$Top)
overlap_prec_overlap$Group <- "Original"
overlap_prec_overlap$Method <- recode(overlap_prec_overlap$Method,
                                           sct="sctransform", prn="Pearson", spr="Spearman",
                                           propr="propr",ana_prn="Analytic PR",
                                           noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

inflation_unfil_top <- as.data.frame(as.table((overlap_string-overlap_true)/overlap_string))
inflation_unfil_top$Var1 <- recode(inflation_unfil_top$Var1,
                                      sct="sctransform", prn="Pearson", spr="Spearman",
                                      propr="propr",ana_prn="Analytic PR",
                                      noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

inflation_unfil_top2 <- as.data.frame(as.table(overlap_true))
inflation_unfil_top2$group <- "True overlaps"
inflation_unfil_top2$Var1 <- recode(inflation_unfil_top2$Var1,
                                    sct="sctransform", prn="Pearson", spr="Spearman",
                                    propr="propr",ana_prn="Analytic PR",
                                    noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")



## overlap with biological network with different p-value threshold --------------------
p_cutoff <- c(0.001, 0.005, 0.01, 0.05, 0.1)
overlap_string_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_string_p) <- p_cutoff
rownames(overlap_string_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_true_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_true_p) <- p_cutoff
rownames(overlap_true_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_prec_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_prec_p) <- p_cutoff
rownames(overlap_prec_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_cor_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_cor_p) <- p_cutoff
rownames(overlap_cor_p) <- c(colnames(ROSMAP_p), "cscore_p")


estimate_p <- estimate_ROSMAP[,1:9]
estimate_p <- cbind(estimate_p, ROSMAP_p_adj)
estimate_p$cscore_p <- ROSMAP_cscore_p[upper.tri(ROSMAP_cscore_p)]
total_pair <- nrow(estimate_p)

for (i in 1:length(p_cutoff)){
  thresh <- p_cutoff[i]
  print(thresh)

  for (j in rownames(overlap_true_p)){
    ROSMAP_deci <- estimate_p[estimate_p[,j]<thresh,]
    overlap_string_p[j,i] <- sum(ROSMAP_deci$grp %in% biological_net$grp)
    overlap_true_p[j,i] <- sum((ROSMAP_deci$grp %in% biological_net$grp)&ROSMAP_deci$true_cor==1)
    overlap_prec_p[j,i] <- sum(ROSMAP_deci$true_cor==1)/nrow(ROSMAP_deci)
    overlap_cor_p[j,i] <- nrow(ROSMAP_deci)
  }
}


overlap_prec_p_long <- melt(overlap_prec_p)
colnames(overlap_prec_p_long) <- c("Method", "Top", "Precision")
overlap_string_p_long <- melt(overlap_string_p)
colnames(overlap_string_p_long) <- c("Method", "Top", "Overlap")

overlap_prec_overlap_p <- left_join(overlap_prec_p_long, overlap_string_p_long, by=c("Method", "Top"))
overlap_prec_overlap_p$Top <- as.factor(overlap_prec_overlap_p$Top)
overlap_prec_overlap_p$Method <- recode(overlap_prec_overlap_p$Method,
                                       sct="sctransform", prn="Pearson", spr="Spearman",
                                       propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                                       noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

inflation_p <- as.data.frame(as.table((overlap_string_p-overlap_true_p)/overlap_string_p))
inflation_p$Var1 <- recode(inflation_p$Var1,
                                 sct="sctransform", prn="Pearson", spr="Spearman",
                                 propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                                 noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

inflation_p2 <- as.data.frame(as.table(overlap_true_p))
inflation_p2$group <- "True overlaps"
inflation_p2$Var1 <- recode(inflation_p2$Var1,
                                  sct="sctransform", prn="Pearson", spr="Spearman",
                                  propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                                  noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

format_supp <- theme(text = element_text(size = 14),
                     legend.position="none")
reactome_prec_overlap_p_sub = overlap_prec_overlap_p[overlap_prec_overlap_p$Method=="CS-CORE",]
reactome_scaling_factor_overlap_p <- max(reactome_prec_overlap_p_sub$Overlap) 
p_reactome_prec_overlap_p_cscore = ggplot(reactome_prec_overlap_p_sub, aes(x = as.factor(Top), group = 1)) +
  geom_col(aes(y = Overlap), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = Precision * reactome_scaling_factor_overlap_p), color = "darkred", size = 1.5) +
  geom_point(aes(y = Precision * reactome_scaling_factor_overlap_p), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of overlaps with Reactome",
    sec.axis = sec_axis(trans = ~ . / reactome_scaling_factor_overlap_p, name = "Precision")) +
  labs(title = "P-value",x = "P-value cutoffs") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_reactome_prec_overlap_p_cscore

reactome_inflation_p_cscore = inflation_p
colnames(reactome_inflation_p_cscore)[3] = "infla"
reactome_inflation_p_cscore = left_join(reactome_inflation_p_cscore, inflation_p2, by=c("Var1", "Var2"))
reactome_inflation_p_cscore = reactome_inflation_p_cscore[reactome_inflation_p_cscore$Var1=="CS-CORE",]
reactome_scaling_factor_infla_p <- max(reactome_inflation_p_cscore$Freq) 

p_reactome_infla_p_cscore = ggplot(reactome_inflation_p_cscore, aes(x = Var2, group = 1)) +
  geom_col(aes(y = Freq), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = infla * reactome_scaling_factor_infla_p), color = "darkred", size = 1.5) +
  geom_point(aes(y = infla * reactome_scaling_factor_infla_p), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of true overlaps with Reactome",
    sec.axis = sec_axis(trans = ~ . / reactome_scaling_factor_infla_p, name = "Prop of misidentified overlaps")) +
  labs(title = "P-value",x = "P-value cutoffs") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_reactome_infla_p_cscore


reactome_prec_overlap_sub = overlap_prec_overlap[overlap_prec_overlap$Method=="CS-CORE \n(Empirical)",]
reactome_scaling_factor_overlap <- max(reactome_prec_overlap_sub$Overlap) 
p_reactome_prec_overlap_cscore = ggplot(reactome_prec_overlap_sub, aes(x = as.factor(Top), group = 1)) +
  geom_col(aes(y = Overlap), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = Precision * reactome_scaling_factor_overlap), color = "darkred", size = 1.5) +
  geom_point(aes(y = Precision * reactome_scaling_factor_overlap), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of overlaps with Reactome",
    sec.axis = sec_axis(trans = ~ . / reactome_scaling_factor_overlap, name = "Precision")) +
  labs(title = "Correlation strength",x = "Top") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_reactome_prec_overlap_cscore

reactome_inflation_cscore = inflation_unfil_top
colnames(reactome_inflation_cscore)[3] = "infla"
reactome_inflation_cscore = left_join(reactome_inflation_cscore, inflation_unfil_top2, by=c("Var1", "Var2"))
reactome_inflation_cscore = reactome_inflation_cscore[reactome_inflation_cscore$Var1=="CS-CORE \n(Empirical)",]
reactome_scaling_factor_infla <- max(reactome_inflation_cscore$Freq) 

p_reactome_infla_cscore = ggplot(reactome_inflation_cscore, aes(x = Var2, group = 1)) +
  geom_col(aes(y = Freq), fill = "darkblue", alpha = 0.7) +
  geom_line(aes(y = infla * reactome_scaling_factor_infla), color = "darkred", size = 1.5) +
  geom_point(aes(y = infla * reactome_scaling_factor_infla), color = "darkred", size = 3) +
  scale_y_continuous(name = "# of true overlaps with Reactome",
    sec.axis = sec_axis(trans = ~ . / reactome_scaling_factor_infla, name = "Prop of misidentified overlaps")) +
  labs(title = "Correlation strength",x = "Top") + theme_bw() +
  theme(
    axis.title.y.left = element_text(color = "darkblue"),
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkred"),
    axis.text.y.right = element_text(color = "darkred"),
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+format_supp
p_reactome_infla_cscore

# fixed mis prop-------------------------------------------------------------------
## overlap with biological network with different top threshold -----------------
top_cutoff <- sort(c(seq(1000, 10000, by=1000), 20000,25000,30000,35000,40000,45000,50000))
overlap_string <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_string) <- top_cutoff
rownames(overlap_string) <- colnames(ROSMAP_p)

overlap_true <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_true) <- top_cutoff
rownames(overlap_true) <- colnames(ROSMAP_p)

overlap_prec <- matrix(NA, ncol = length(top_cutoff), nrow=7)
colnames(overlap_prec) <- top_cutoff
rownames(overlap_prec) <- colnames(ROSMAP_p)


total_pair <- nrow(estimate_ROSMAP)
for (i in 1:length(top_cutoff)){
  thresh <- top_cutoff[i]
  print(thresh)

  selected_genes_thresh <- list()
  for (j in rownames(overlap_string)){
    cor_ROSMAP <- abs(estimate_ROSMAP[,j])
    deci_ROSMAP <- estimate_ROSMAP[cor_ROSMAP>quantile(cor_ROSMAP, 1-thresh/total_pair),]
    overlap_string[j,i] <- sum(deci_ROSMAP$grp %in% biological_net$grp)
    overlap_true[j,i] <- sum((deci_ROSMAP$grp %in% biological_net$grp)&deci_ROSMAP$true_cor==1)
    overlap_prec[j,i] <- sum(deci_ROSMAP$true_cor==1)/nrow(deci_ROSMAP)
  }
}

overlap_string_long <- reshape2::melt(overlap_true)
colnames(overlap_string_long) <- c("Method", "Top", "Overlap")

inflation_unfil_top <- as.data.frame(as.table((overlap_string-overlap_true)/overlap_string))
colnames(inflation_unfil_top) <- c("Method", "cutoff", "Mis")
inflation_unfil_top$true <- overlap_string_long$Overlap
inflation_unfil_top_sub <- inflation_unfil_top[inflation_unfil_top$Method=="cscore_est",]
inflation_unfil_top_sub$Method <- "cscore_p"
inflation_unfil_top <- rbind(inflation_unfil_top_sub, inflation_unfil_top)

inflation_unfil_top$Method <- recode(inflation_unfil_top$Method,
                                     sct="sctransform", prn="Pearson", spr="Spearman",
                                     propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                                     noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")

## overlap with biological network with different p-value threshold --------------------
p_cutoff <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
overlap_string_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_string_p) <- p_cutoff
rownames(overlap_string_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_true_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_true_p) <- p_cutoff
rownames(overlap_true_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_prec_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_prec_p) <- p_cutoff
rownames(overlap_prec_p) <- c(colnames(ROSMAP_p), "cscore_p")

overlap_cor_p <- matrix(NA, ncol = length(p_cutoff), nrow=8)
colnames(overlap_cor_p) <- p_cutoff
rownames(overlap_cor_p) <- c(colnames(ROSMAP_p), "cscore_p")

estimate_p <- estimate_ROSMAP[,1:9]
estimate_p <- cbind(estimate_p, ROSMAP_p_adj)
estimate_p$cscore_p <- ROSMAP_cscore_p[upper.tri(ROSMAP_cscore_p)]
total_pair <- nrow(estimate_p)

for (i in 1:length(p_cutoff)){
  thresh <- p_cutoff[i]
  print(thresh)

  for (j in rownames(overlap_true_p)){
    ROSMAP_deci <- estimate_p[estimate_p[,j]<thresh,]
    overlap_string_p[j,i] <- sum(ROSMAP_deci$grp %in% biological_net$grp)
    overlap_true_p[j,i] <- sum((ROSMAP_deci$grp %in% biological_net$grp)&ROSMAP_deci$true_cor==1)
    overlap_prec_p[j,i] <- sum(ROSMAP_deci$true_cor==1)/nrow(ROSMAP_deci)
    overlap_cor_p[j,i] <- nrow(ROSMAP_deci)
  }
}


overlap_string_p_long <- melt(overlap_true_p)
colnames(overlap_string_p_long) <- c("Method", "Top", "Overlap")
inflation_p <- as.data.frame(as.table((overlap_string_p-overlap_true_p)/overlap_string_p))
colnames(inflation_p) <- c("Method", "cutoff", "Mis")
inflation_p$true <- overlap_string_p_long$Overlap
inflation_p$Method <- recode(inflation_p$Method,
                             sct="sctransform", prn="Pearson", spr="Spearman",
                             propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                             noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")


i = "CS-CORE"
plot_dat1 <- inflation_p[inflation_p$Method==i,]
plot_dat1$group <- "P-value"
plot_dat2 <- inflation_unfil_top[inflation_unfil_top$Method==i,]
plot_dat2$group <- "Cor-strength"
plot_dat_reactome <- rbind(plot_dat1, plot_dat2)
fix_mis_cscore_reactome <- ggplot(plot_dat_reactome, aes(x=Mis, y=true, color=group))+
    geom_point(size=2)+geom_line(size=1)+labs( x="Prop of \nmisidentified overlaps", y="# of true overlaps with Reactome", color="")+
    theme_bw()+guides(color = guide_legend(nrow = 2)) +
    # scale_colour_manual(values = c("P-value"="darkblue", "Cor-strength"="darkred"))+
    theme(text = element_text(size = 14),
         plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1),
        #  legend.position = c(0.98,1.1),
        # legend.justification = c("right", "top"),
        # legend.background = element_blank()
         legend.position = "top", legend.box.margin = margin(b = -10))
fix_mis_cscore_reactome

plot_reactome_ori <- ggarrange(p_reactome_prec_overlap_cscore, p_reactome_prec_overlap_p_cscore,
                             p_reactome_infla_cscore, 
                             p_reactome_infla_p_cscore, fix_mis_cscore_reactome+xlim(0,0.1), ncol=5,nrow=1,
                          widths = c(1, 0.9, 1, 0.9,0.8),
          labels = c("A", "B", "C", "D","E"))
p_reactome = annotate_figure(ggarrange(plot_reactome_ori, ncol=1, nrow=1),
  top = text_grob("Overlap with reactome", 
                  color = "black", face = "bold", size = 16))

pdf('/gpfs/gibbs/pi/zhao/xs282/validation/revision/modify_plot/plot_cscore_reactome.pdf', width = 14, height = 4, onefile = T)
p_reactome
dev.off()
