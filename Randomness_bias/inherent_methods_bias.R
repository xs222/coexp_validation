# based on the simulated data
library(gridExtra)
library(tidyverse)
library(magrittr)
library(rtracklayer)
library(biomaRt)
library(ggvenn)
library(ggplot2)
library(igraph)
library(venn)
library(data.table)
library(ggpubr)
library(tidyverse)
library(matrixcalc)

set.seed(11272023)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("AFinal/cscore_real_data_function.R")

# prepare biological network----------------------------------------------------
# STRING
hs_filter <- readRDS("STRING/hs_filter_10_4_2023.rds")
hs_filter <- hs_filter %>%
  mutate(grp = paste(pmax(protein1, protein2), pmin(protein1, protein2), sep = "_"))

# prepare coexpression network---------------------------------------------------
ROSMAP_oli_ct <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_sct_cor_NB_simu1000_abs_thresh.rds")

ROSMAP_ori_ests <- readRDS("mean_cor/semi_PD/simu/ROSMAP_NC_Oli_sct1000.rds")
ROSMAP_ori_ests[abs(ROSMAP_ori_ests)<0.015] <- 0
gene_name <- rownames(ROSMAP_ori_ests)
all(rownames(ROSMAP_ori_ests)==colnames(ROSMAP_ori_ests))

# cor estimations
ROSMAP_cscore <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_cscore1000_abs_thresh.rds")
ROSMAP_cscore_p <- MatrixBH(ROSMAP_cscore$p_value)
ROSMAP_cscore_p <- ROSMAP_cscore_p[gene_name, gene_name]
ROSMAP_cscore_est <- ROSMAP_cscore$est[gene_name, gene_name]
ROSMAP_cscore_est_filter <- ROSMAP_cscore_est
ROSMAP_cscore_est_filter[ROSMAP_cscore_p >= 0.05] <- 0
min(abs(ROSMAP_cscore_est_filter)[ROSMAP_cscore_est_filter!=0])


ROSMAP_p <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_simu_norm_p.rds")
ROSMAP_p_adj <- as.data.frame(apply(ROSMAP_p, 2, function(x){p.adjust(x, method = "BH")}))
upper2matrix <- function(est_mat, col_name){
  p_mat <- est_mat-est_mat
  p_mat[upper.tri(p_mat)] <- ROSMAP_p_adj[,col_name]
  p_mat <- p_mat + t(p_mat)
  return(p_mat)
}


ROSMAP_sct <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_sct1000_abs_thresh.rds")
ROSMAP_sct_est <- ROSMAP_sct[gene_name, gene_name]
ROSMAP_sct_p <- upper2matrix(ROSMAP_sct_est, "sct")
ROSMAP_sct_est_filter <- ROSMAP_sct_est
ROSMAP_sct_est_filter[ROSMAP_sct_p >= 0.05] <- 0


ROSMAP_ana_prn <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_ana_prn1000_abs_thresh.rds")
ROSMAP_ana_prn_est <- ROSMAP_ana_prn[gene_name, gene_name]
ROSMAP_ana_prn_p <- upper2matrix(ROSMAP_ana_prn_est, "ana_prn")
ROSMAP_ana_prn_est_filter <- ROSMAP_ana_prn_est
ROSMAP_ana_prn_est_filter[ROSMAP_ana_prn_p >= 0.05] <- 0


ROSMAP_noise <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_noise1000_abs_thresh.rds")
ROSMAP_noise_est <- ROSMAP_noise[gene_name, gene_name]
ROSMAP_noise_est <- apply(ROSMAP_noise_est, c(1,2), as.numeric)
ROSMAP_noise_p <- upper2matrix(ROSMAP_noise_est, "noise")
ROSMAP_noise_est_filter <- ROSMAP_noise_est
ROSMAP_noise_est_filter[ROSMAP_noise_p >= 0.05] <- 0


ROSMAP_propr <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_propr1000_abs_thresh.rds")
ROSMAP_propr_est <- ROSMAP_propr@matrix[gene_name, gene_name]
is.positive.definite((ROSMAP_propr_est+t(ROSMAP_propr_est))/2)
ROSMAP_propr_p <- upper2matrix(ROSMAP_propr_est, "propr")
ROSMAP_propr_est_filter <- ROSMAP_propr_est
ROSMAP_propr_est_filter[ROSMAP_propr_p >= 0.05] <- 0


ROSMAP_prn <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_prn1000_abs_thresh.rds")
ROSMAP_prn_est <- ROSMAP_prn[gene_name, gene_name]
ROSMAP_prn_p <- upper2matrix(ROSMAP_prn_est, "prn")
ROSMAP_prn_est_filter <- ROSMAP_prn_est
ROSMAP_prn_est_filter[ROSMAP_prn_p >= 0.05] <- 0

ROSMAP_spr <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_simu_spr1000_abs_thresh.rds")
ROSMAP_spr_est <- ROSMAP_spr[gene_name, gene_name]
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

estimate_p <- estimate_ROSMAP[,1:9]
ROSMAP_p_adj$cscore_p <- ROSMAP_cscore_p[upper.tri(ROSMAP_cscore_p)]

estimate_p <- cbind(estimate_p, ROSMAP_p_adj<0.05)

estimate_p_long = long_dt <- melt(estimate_p, id.vars = colnames(estimate_p)[c(1:7,9)])

estimate_p_long_sub = estimate_p_long[estimate_p_long$value==1,]
color_setting <- c("CS-CORE (Empirical)"="brown", "Noise Regularization"="#AF58BA",
                   "CS-CORE"="#339933", "sctransform"="#ff6699",
                   "Pearson"="#F28522", "Spearman"="#ffff66","Analytic PR"="#99ccff",
                   "propr"="#3366cc", "True"="black")
estimate_p_long_sub$variable <- recode(estimate_p_long_sub$variable,
                                  sct="sctransform", prn="Pearson", spr="Spearman",
                                  propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                                  noise="Noise Regularization", cscore_est="CS-CORE (Empirical)", true_cor="True")


pdf('/gpfs/gibbs/pi/zhao/xs282/validation/revision/modify_plot/inherent_method_bias.pdf', width = 6, height = 3, onefile = T)

ggplot(estimate_p_long_sub, aes(x=log10mean_mu_ROSMAP, color=variable))+
    geom_density(size=1)+theme_bw()+
    scale_color_manual(values = color_setting)+
    theme(text = element_text(size = 14))+
    labs(x="Gene pair expression level (log10)", y="Density", color="Method")
dev.off()









