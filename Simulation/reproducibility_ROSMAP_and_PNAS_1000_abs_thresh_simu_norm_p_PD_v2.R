# reproducibility based if the gene pair is correlated
# /gpfs/gibbs/pi/zhao/xs282/validation/mean_cor/p_value_PNAS_and_ROSMAP_NC_Oli_12_5_2023.R

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggpattern)
library(cowplot)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("/gpfs/gibbs/pi/zhao/xs282/validation/cscore_real_data_function.R")

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
ROSMAP_cscore_est_filter <- ROSMAP_cscore_est
ROSMAP_cscore_est_filter[ROSMAP_cscore_p >= 0.05] <- 0

PNAS_cscore <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_simu_cscore1000_abs_thresh.rds")
PNAS_cscore_p <- MatrixBH(PNAS_cscore$p_value)
PNAS_cscore_p <- PNAS_cscore_p[gene_name, gene_name]
PNAS_cscore_est <- PNAS_cscore$est[gene_name, gene_name]
PNAS_cscore_est_filter <- PNAS_cscore_est
PNAS_cscore_est_filter[PNAS_cscore_p >= 0.05] <- 0

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
count_mat <- matrix(NA, nrow=2, ncol = 8)
colnames(count_mat) <- colnames(ROSMAP_p_adj)
rownames(count_mat) <- c("ROSMAP", "PNAS")
count_mat[1,] <- apply(ROSMAP_p_adj, 2, function(x){sum(x<0.3)})
count_mat[2,] <- apply(PNAS_p_adj, 2, function(x){sum(x<0.3)})
knitr::kable(count_mat)

count_mat <- matrix(NA, nrow=2, ncol = 8)
colnames(count_mat) <- colnames(ROSMAP_p_adj)
rownames(count_mat) <- c("ROSMAP", "PNAS")
count_mat[1,] <- apply(ROSMAP_p_adj, 2, function(x){sum(x<0.01)})
count_mat[2,] <- apply(PNAS_p_adj, 2, function(x){sum(x<0.01)})
knitr::kable(count_mat)

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
# color_setting <- c("CS-CORE \n(Empirical)"="darkgreen", "Noise \nRegularization"="#FF1F5B",
#                    "CS-CORE"="brown", "sctransform"="#AF58BA",
#                    "Pearson"="#FFC61E", "Spearman"="#F28522","Analytic PR"="#a3d0d4",
#                    "propr"="#2166AC")

color_setting <- c("CS-CORE \n(Empirical)"="brown", "Noise \nRegularization"="#AF58BA",
                   "CS-CORE"="#339933", "sctransform"="#ff6699",
                   "Pearson"="#F28522", "Spearman"="#ffff66","Analytic PR"="#99ccff",
                   "propr"="#3366cc")

p_prec_repduc_p <- ggplot(repduc_prec_overlap_p, aes(x=Reproduce, y=ROSMAP, color=Method))+
  geom_point()+ geom_line(size=1)+ylim(c(0.2,1))+theme_bw()+
  # geom_errorbar(aes(ymin = min, ymax = max),width=500)+
  scale_colour_manual(values = color_setting)+
  labs(title="P-value", x="# of reproducible pairs", y="Precision", shape="P-value cutoff")+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_prec_repduc_p

p_prec_repduc_p_pnas <- ggplot(repduc_prec_overlap_p, aes(x=Reproduce, y=PNAS, color=Method))+
  geom_point()+ geom_line(size=1)+ylim(c(0.2,1))+theme_bw()+
  # geom_errorbar(aes(ymin = min, ymax = max),width=500)+
  scale_colour_manual(values = color_setting)+
  labs(title="P-value", x="# of reproducible pairs", y="Precision", shape="P-value cutoff")+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_prec_repduc_p_pnas


inflation_p <- as.data.frame(as.table((repduc_p-repduc_truth_p)/repduc_p))
inflation_p$Var1 <- recode(inflation_p$Var1,
                           sct="sctransform", prn="Pearson", spr="Spearman",
                           propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                           noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")
reproduc_infla_p <- ggplot(inflation_p, aes(x=Var2, y=Freq, color=Var1, group=Var1))+
  geom_point()+geom_line(size=1)+
  labs(y="Misidentified reproducible pairs" , #y="Inflation",
       x="P-value cutoffs", color="Method", title="P-value")+
  theme_bw()+scale_color_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
reproduc_infla_p


inflation_p1 <- as.data.frame(as.table(repduc_p))
inflation_p2 <- as.data.frame(as.table(repduc_truth_p))
inflation_p1$group <- "Identified reproducible pairs"
inflation_p2$group <- "True reproducible pairs"
inflation_p1$Var1 <- recode(inflation_p1$Var1,
                           sct="sctransform", prn="Pearson", spr="Spearman",
                           propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                           noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")
inflation_p2$Var1 <- recode(inflation_p2$Var1,
                            sct="sctransform", prn="Pearson", spr="Spearman",
                            propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                            noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")
inflation_p1_sub <- inflation_p1[inflation_p1$Var2==max(p_cutoff),] %>% arrange(Var2, desc(Freq))
inflation_p1$Var1 <- factor(inflation_p1$Var1, levels = unique(inflation_p1_sub$Var1))
inflation_p2$Var1 <- factor(inflation_p2$Var1, levels = unique(inflation_p1_sub$Var1))


reproduc_count_p <- ggplot()+
  geom_bar(data=inflation_p1, mapping=aes(x=Var2, y=Freq, fill=Var1),
           stat = "identity",position=position_dodge())+
  labs(y="# of reproducible pairs", x="P-value cutoffs", fill="Method", title="P-value")+
  theme_bw()+
  scale_fill_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
reproduc_count_p

reproduc_count_true_p <- ggplot()+
  geom_bar(data=inflation_p2, mapping=aes(x=Var2, y=Freq, fill=Var1),
           stat = "identity",position=position_dodge())+
  labs(y="# of true reproducible pairs", x="P-value cutoffs", fill="Method", title="P-value")+
  theme_bw()+
  scale_fill_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
reproduc_count_true_p

reproduc_count_p+reproduc_count_true_p



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
                            noise="Noise \nRegularization", cscore_est="CS-CORE")
# shape_setting2 <- 1:length(top_cutoff)
# names(shape_setting2) <- top_cutoff

p_prec_repduc_unfil <- ggplot(repduc_prec_overlap, aes(x=Reproduce, y=ROSMAP, color=Method))+
  geom_point()+ geom_line(size=1)+ylim(c(0.2,1))+theme_bw()+
  scale_colour_manual(values = color_setting)+
  labs(title="Correlation strength", x="# of reproducible pairs", y="Precision")+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_prec_repduc_unfil

p_prec_repduc_unfil_pnas <- ggplot(repduc_prec_overlap, aes(x=Reproduce, y=PNAS, color=Method))+
  geom_point()+ geom_line(size=1)+ylim(c(0,1))+theme_bw()+
  scale_colour_manual(values = color_setting)+
  labs(title="Correlation strength", x="# of reproducible pairs", y="Precision", shape="Top")+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_prec_repduc_unfil_pnas

inflation_unfil_top <- as.data.frame(as.table((repduc_top-repduc_truth_top)/repduc_top))
inflation_unfil_top$Var1 <- recode(inflation_unfil_top$Var1,
                           sct="sctransform", prn="Pearson", spr="Spearman",
                           propr="propr",ana_prn="Analytic PR",
                           noise="Noise \nRegularization", cscore_est="CS-CORE")
reproduc_infla_unfil <- ggplot(inflation_unfil_top, aes(x=Var2, y=Freq, color=Var1, group=Var1))+
  geom_point()+geom_line(size=1)+
  labs(y="Misidentified reproducible pairs", #y="Inflation",
       x="Top", color="Method", title="Correlation strength")+
  theme_bw()+scale_color_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
reproduc_infla_unfil


inflation_unfil_top1 <- as.data.frame(as.table(repduc_top))
inflation_unfil_top2 <- as.data.frame(as.table(repduc_truth_top))
inflation_unfil_top1$group <- "Identified reproducible pairs"
inflation_unfil_top2$group <- "True reproducible pairs"
inflation_unfil_top1$Var1 <- recode(inflation_unfil_top1$Var1,
                            sct="sctransform", prn="Pearson", spr="Spearman",
                            propr="propr",ana_prn="Analytic PR",
                            noise="Noise \nRegularization", cscore_est="CS-CORE")
inflation_unfil_top2$Var1 <- recode(inflation_unfil_top2$Var1,
                            sct="sctransform", prn="Pearson", spr="Spearman",
                            propr="propr",ana_prn="Analytic PR",
                            noise="Noise \nRegularization", cscore_est="CS-CORE")
inflation_unfil_top1_sub <- inflation_unfil_top1[inflation_unfil_top1$Var2==max(top_cutoff),] %>% arrange(desc(Freq))
inflation_unfil_top1$Var1 <- factor(inflation_unfil_top1$Var1, levels = unique(inflation_unfil_top1_sub$Var1))
inflation_unfil_top2$Var1 <- factor(inflation_unfil_top2$Var1, levels = unique(inflation_unfil_top1_sub$Var1))

reproduc_count_unfil <- ggplot()+
  geom_bar(data=inflation_unfil_top1, mapping=aes(x=Var2, y=Freq, fill=Var1),
           stat = "identity",position=position_dodge())+
  labs(y="# of reproducible pairs", x="Top", fill="Method", title="Correlation strength")+
  theme_bw()+
  scale_fill_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
reproduc_count_unfil

reproduc_count_true_unfil <- ggplot()+
  geom_bar(data=inflation_unfil_top2, mapping=aes(x=Var2, y=Freq, fill=Var1),
           stat = "identity",position=position_dodge())+
  labs(y="# of true reproducible pairs", x="Top", fill="Method", title="Correlation strength")+
  theme_bw()+
  scale_fill_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
reproduc_count_true_unfil
reproduc_count_unfil+reproduc_count_true_unfil

format_supp <- theme(text = element_text(size = 14),
                     legend.position="none")
plot_rep_ori <- ggarrange(reproduc_count_unfil+format_supp, p_prec_repduc_unfil+format_supp,
                          reproduc_count_true_unfil+format_supp,
                          reproduc_infla_unfil+format_supp+ylim(0,1)+labs(y="Prop of misidentified pairs"),
                          reproduc_count_p+format_supp,p_prec_repduc_p+format_supp,
                          reproduc_count_true_p+format_supp,
                          reproduc_infla_p+format_supp+ylim(0,1)+labs(y="Prop of misidentified pairs"), ncol=4,nrow=2,
                          widths = c(1, 0.9, 1, 0.9),
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
leg <- get_legend(reproduc_count_p+
                    scale_fill_manual(values = color_setting,
                                       breaks = sort(names(color_setting)),
                                       labels = c("Analytic PR", "CS-CORE", "CS-CORE (Empirical)", "Noise Regularization",
                                                  "Pearson", "propr", "sctransform", "Spearman"))+
                    theme(legend.position = "bottom",
                           text = element_text(size = 17),
                           legend.spacing.x = unit(0.5, 'cm'))+
                    guides(fill=guide_legend(nrow=2,byrow=TRUE)))


pdf('mean_cor/semi_PD_sparse/figures/reproduc_v2.pdf', width = 12, height = 7.5, onefile = T)
ggarrange(plot_rep_ori, leg, ncol=1, nrow=2, heights = c(10,1))
dev.off()

plot_string_ori <- ggarrange(string_count_unfil+format_supp, p_prec_string_unfil+format_supp,
                             string_count_true_unfil+format_supp,
                             string_infla_unfil+format_supp+ylim(0,0.5)+labs(y="Prop of misidentified overlaps"),
                             string_count_p+format_supp,p_prec_string_p+format_supp,
                             string_count_true_p+format_supp,
                             string_infla_p+format_supp+ylim(0,0.5)+labs(y="Prop of misidentified overlaps"), ncol=4,nrow=2,
                             widths = c(1, 0.95, 1, 0.95),labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
pdf('mean_cor/semi_PD_sparse/figures/string_v2.pdf', width = 12, height = 8, onefile = T)
ggarrange(plot_string_ori, leg, ncol=1, nrow=2, heights = c(10,1))
dev.off()



leg <- get_legend(reproduc_infla_p+theme(legend.position = "right",
                                         text = element_text(size = 17),
                                         legend.spacing.x = unit(0.5, 'cm')))
pdf('mean_cor/semi_PD_sparse/figures/prec_pnas_v2.pdf', width = 12, height = 3, onefile = T)
ggarrange(NULL, ggarrange(p_prec_repduc_unfil_pnas+format_supp+labs(y="Precision (PNAS)"),
                    p_prec_repduc_p_pnas+format_supp+labs(y="Precision (PNAS)"), ncol=2, nrow = 1, labels = c("A", "B")),
          leg, NULL, ncol=4, nrow=1, widths = c(1.5,5,2,1))
dev.off()



# fixed mis prop --------------------------------------------------------------
# based on p-value--------------------------------------------------------------
count_mat <- matrix(NA, nrow=2, ncol = 8)
colnames(count_mat) <- colnames(ROSMAP_p_adj)
rownames(count_mat) <- c("ROSMAP", "PNAS")
count_mat[1,] <- apply(ROSMAP_p_adj, 2, function(x){sum(x<0.3)})
count_mat[2,] <- apply(PNAS_p_adj, 2, function(x){sum(x<0.3)})
knitr::kable(count_mat)

count_mat <- matrix(NA, nrow=2, ncol = 8)
colnames(count_mat) <- colnames(ROSMAP_p_adj)
rownames(count_mat) <- c("ROSMAP", "PNAS")
count_mat[1,] <- apply(ROSMAP_p_adj, 2, function(x){sum(x<0.01)})
count_mat[2,] <- apply(PNAS_p_adj, 2, function(x){sum(x<0.01)})
knitr::kable(count_mat)

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


repro_p_ls <- list()
for (i in unique(inflation_p$Method)){
  plot_dat1 <- inflation_p[inflation_p$Method==i,]
  plot_dat1$group <- "P-value"
  plot_dat2 <- inflation_unfil_top[inflation_unfil_top$Method==i,]
  plot_dat2$group <- "Cor-strength"
  plot_dat <- rbind(plot_dat1, plot_dat2)
  repro_p_ls[[i]] <- ggplot(plot_dat, aes(x=Mis, y=true, color=group))+
    geom_point(size=2)+geom_line(size=1)+labs(title=i, x="",
                                  y="", color="")+
    theme_bw()+
    scale_colour_manual(values = c("P-value"="darkblue", "Cor-strength"="darkred"))+
    theme(legend.position = "bottom")

}

repro_p_ls[["Pearson"]] <- repro_p_ls[["Pearson"]]+xlim(0,0.1)+ylim(0,1500)
repro_p_ls[["propr"]] <- repro_p_ls[["propr"]]+xlim(0,0.4)+ylim(0,5000)
repro_p_ls[["Spearman"]] <- repro_p_ls[["Spearman"]]+ylim(0,2500)+xlim(0,0.25)
repro_p_ls[["sctransform"]] <- repro_p_ls[["sctransform"]]+xlim(0,0.2)
repro_p_ls[["Analytic PR"]] <- repro_p_ls[["Analytic PR"]]+xlim(0,0.175)


ggarrange(plotlist = repro_p_ls, nrow=2, ncol=4, common.legend = T, legend = "bottom")

plots <- repro_p_ls
legend <- get_legend(plots[[1]]+theme(text = element_text(size = 17)))
adjusted_theme <- theme(legend.position = "none",text = element_text(size = 15),
                        plot.title = element_text(hjust=0.5),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        plot.tag = element_text(size = 14, face = "bold", vjust = 1.3, hjust = -1.5),
                        plot.tag.position = c(0,1),
                        plot.margin = unit(c(0, 0.15, 0, 0), "cm"))

plot_grid <- plot_grid(
  plots[[6]] + adjusted_theme +labs(tag="A"),
  plots[[8]] +adjusted_theme+labs(tag="B"),
  plots[[5]] + adjusted_theme+labs(tag="C"),
  plots[[7]] + adjusted_theme+labs(tag="D"),
  plots[[2]] + adjusted_theme+labs(tag="E"),
  plots[[4]] + adjusted_theme+labs(tag="F"),
  plots[[1]] + adjusted_theme+labs(tag="G"),
  plots[[3]] + adjusted_theme+labs(tag="H"),
  ncol = 4, nrow = 2, align = "hv"
)

final_plot_with_labels <- ggdraw() +
  draw_plot(plot_grid, 0.02, 0.09, 0.98, 0.9, hjust = 0) +
  draw_label("Prop of misidentified reproducible pairs", x = 0.51, y = 0.08, vjust = -0.5, angle = 0, size = 15) +
  draw_label("# of true reproducible pairs", x = 0, y = 0.55, vjust = 1.5, angle = 90, size = 15) +
  draw_plot(legend, 0, 0, 1, 0.1)


# Print the final plot
print(final_plot_with_labels)

pdf('mean_cor/semi_PD_sparse/figures/compare_reproduc_v2.pdf', width = 11, height = 7, onefile = T)
print(final_plot_with_labels)
dev.off()
