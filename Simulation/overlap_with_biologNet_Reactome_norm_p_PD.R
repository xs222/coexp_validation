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
source("/gpfs/gibbs/pi/zhao/xs282/validation/cscore_real_data_function.R")

# prepare biological network----------------------------------------------------
# Reactome
hs_filter <- readRDS("Reactome/Processed_Reactome_10_07_2023.rds")
hs_filter <- hs_filter %>%
  mutate(grp = paste(pmax(V1, V2), pmin(V1, V2), sep = "_"))

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

# color_setting <- c("CS-CORE \n(Empirical)"="darkgreen", "Noise \nRegularization"="#FF1F5B",
#                    "CS-CORE"="brown", "sctransform"="#AF58BA",
#                    "Pearson"="#FFC61E", "Spearman"="#F28522","Analytic PR"="#a3d0d4",
#                    "propr"="#2166AC")
color_setting <- c("CS-CORE \n(Empirical)"="brown", "Noise \nRegularization"="#AF58BA",
                   "CS-CORE"="#339933", "sctransform"="#ff6699",
                   "Pearson"="#F28522", "Spearman"="#ffff66","Analytic PR"="#99ccff",
                   "propr"="#3366cc")


overlap_prec_long <- reshape2::melt(overlap_prec)
colnames(overlap_prec_long) <- c("Method", "Top", "Precision")
overlap_string_long <- reshape2::melt(overlap_string)
colnames(overlap_string_long) <- c("Method", "Top", "Overlap")


# shape_setting <- 1:length(top_cutoff)
# names(shape_setting) <- top_cutoff
overlap_prec_overlap <- left_join(overlap_prec_long, overlap_string_long, by=c("Method", "Top"))
overlap_prec_overlap$Top <- as.factor(overlap_prec_overlap$Top)
overlap_prec_overlap$Group <- "Original"
overlap_prec_overlap$Method <- recode(overlap_prec_overlap$Method,
                                           sct="sctransform", prn="Pearson", spr="Spearman",
                                           propr="propr",ana_prn="Analytic PR",
                                           noise="Noise \nRegularization", cscore_est="CS-CORE")

p_prec_string_unfil <- ggplot(overlap_prec_overlap, aes(x=Overlap, y=Precision, color=Method))+
  geom_point()+ geom_line(size=1)+ylim(c(0.2,1))+theme_bw()+
  # scale_shape_manual(values=shape_setting)+
  scale_colour_manual(values = color_setting)+
  labs(title="Correlation strength", x="Overlaps with Reactome", y="Precision", shape="Top")+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_prec_string_unfil


inflation_unfil_top <- as.data.frame(as.table((overlap_string-overlap_true)/overlap_string))
inflation_unfil_top$Var1 <- recode(inflation_unfil_top$Var1,
                                      sct="sctransform", prn="Pearson", spr="Spearman",
                                      propr="propr",ana_prn="Analytic PR",
                                      noise="Noise \nRegularization", cscore_est="CS-CORE")
string_infla_unfil <- ggplot(inflation_unfil_top, aes(x=Var2, y=Freq, color=Var1, group=Var1))+
  geom_point()+geom_line(size=1)+
  labs(y="Prop of misidentified overlaps", x="Top", color="Method",
       title="Correlation strength")+
  theme_bw()+scale_color_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
string_infla_unfil

inflation_unfil_top1 <- as.data.frame(as.table(overlap_string))
inflation_unfil_top2 <- as.data.frame(as.table(overlap_true))
inflation_unfil_top1$group <- "Identified overlaps"
inflation_unfil_top2$group <- "True overlaps"
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

string_count_unfil <- ggplot()+
  geom_bar(data=inflation_unfil_top1, mapping=aes(x=Var2, y=Freq, fill=Var1),
           stat = "identity",position=position_dodge())+
  labs(y="# of overlaps with Reactome", x="Top", fill="Method", title="Correlation strength")+
  theme_bw()+
  scale_fill_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
string_count_unfil

string_count_true_unfil <- ggplot()+
  geom_bar(data=inflation_unfil_top2, mapping=aes(x=Var2, y=Freq, fill=Var1),
           stat = "identity",position=position_dodge())+
  labs(y="# of true overlaps with Reactome", x="Top", fill="Method", title="Correlation strength")+
  theme_bw()+
  scale_fill_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
string_count_true_unfil


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
p_prec_string_p <- ggplot(overlap_prec_overlap_p, aes(x=Overlap, y=Precision, color=Method))+
  geom_point()+ geom_line(size=1)+ylim(c(0.2,1))+theme_bw()+
  # geom_errorbar(aes(ymin = min, ymax = max),width=500)+
  scale_colour_manual(values = color_setting)+
  labs(title="P-value", x="Overlaps with Reactome", y="Precision", shape="P-value cutoff")+
  theme(plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_prec_string_p


inflation_p <- as.data.frame(as.table((overlap_string_p-overlap_true_p)/overlap_string_p))
inflation_p$Var1 <- recode(inflation_p$Var1,
                                 sct="sctransform", prn="Pearson", spr="Spearman",
                                 propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                                 noise="Noise \nRegularization", cscore_est="CS-CORE \n(Empirical)")
string_infla_p <- ggplot(inflation_p, aes(x=Var2, y=Freq, color=Var1, group=Var1))+
  geom_point()+geom_line(size=1)+
  labs(y="Prop of misidentified overlaps", x="P-value cutoffs", color="Method",  title="P-value")+
  theme_bw()+scale_color_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
string_infla_p


inflation_p1 <- as.data.frame(as.table(overlap_string_p))
inflation_p2 <- as.data.frame(as.table(overlap_true_p))
inflation_p1$group <- "Identified overlaps"
inflation_p2$group <- "True overlaps"
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


string_count_p <- ggplot()+
  geom_bar(data=inflation_p1, mapping=aes(x=Var2, y=Freq, fill=Var1),
           stat = "identity",position=position_dodge())+
  labs(y="# of overlaps with Reactome", x="P-value cutoffs", fill="Method", title="P-value")+
  theme_bw()+
  scale_fill_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
string_count_p

string_count_true_p <- ggplot()+
  geom_bar(data=inflation_p2, mapping=aes(x=Var2, y=Freq, fill=Var1),
           stat = "identity",position=position_dodge())+
  labs(y="# of true overlaps with Reactome", x="P-value cutoffs", fill="Method", title="P-value")+
  theme_bw()+
  scale_fill_manual(values = color_setting)+
  theme(plot.title = element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
string_count_true_p

format_supp <- theme(text = element_text(size = 14),
                     legend.position="none")
leg <- get_legend(string_count_p+
                    scale_fill_manual(values = color_setting,
                                      breaks = sort(names(color_setting)),
                                      labels = c("Analytic PR", "CS-CORE", "CS-CORE (Empirical)", "Noise Regularization",
                                                 "Pearson", "propr", "sctransform", "Spearman"))+
                    theme(legend.position = "bottom",
                          text = element_text(size = 17),
                          legend.spacing.x = unit(0.5, 'cm'))+
                    guides(fill=guide_legend(nrow=2,byrow=TRUE)))

plot_string_ori <- ggarrange(string_count_unfil+format_supp, p_prec_string_unfil+format_supp,
                             string_count_true_unfil+format_supp,
                             string_infla_unfil+format_supp+ylim(0,0.5)+labs(y="Prop of misidentified overlaps"),
                             string_count_p+format_supp,p_prec_string_p+format_supp,
                             string_count_true_p+format_supp,
                             string_infla_p+format_supp+ylim(0,0.5)+labs(y="Prop of misidentified overlaps"), ncol=4,nrow=2,
                             widths = c(1, 0.95, 1, 0.95),labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
pdf('mean_cor/semi_PD_sparse/figures/reactome_v2.pdf', width = 12.4, height = 8.3, onefile = T)
ggarrange(plot_string_ori, leg, ncol=1, nrow=2, heights = c(10,1))
dev.off()


# fixed mis prop

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


string_p_ls <- list()
for (i in unique(inflation_p$Method)){
  plot_dat1 <- inflation_p[inflation_p$Method==i,]
  plot_dat1$group <- "P-value"
  plot_dat2 <- inflation_unfil_top[inflation_unfil_top$Method==i,]
  plot_dat2$group <- "Cor-strength"
  plot_dat <- rbind(plot_dat1, plot_dat2)
  string_p_ls[[i]] <- ggplot(plot_dat, aes(x=Mis, y=true, color=group))+
    geom_point(size=2)+geom_line(size=1)+labs(title=i, x="",
                                  y="", color="")+
    theme_bw()+
    scale_colour_manual(values = c("P-value"="darkblue", "Cor-strength"="darkred"))+
    theme(legend.position = "bottom")

}

string_p_ls[["sctransform"]] <- string_p_ls[["sctransform"]]+xlim(0,0.075)
string_p_ls[["Analytic PR"]] <- string_p_ls[["Analytic PR"]]+xlim(0,0.08)
string_p_ls[["CS-CORE"]] <- string_p_ls[["CS-CORE"]]+xlim(0,0.1)
string_p_ls[["CS-CORE \n(Empirical)"]] <- string_p_ls[["CS-CORE \n(Empirical)"]]+xlim(0,0.13)


ggarrange(plotlist = string_p_ls, nrow=2, ncol=4, common.legend = T, legend = "bottom")

plots <- string_p_ls
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
  draw_label("Prop of misidentified overlaps", x = 0.51, y = 0.08, vjust = -0.5, angle = 0, size = 15) +
  draw_label("# of true overlaps with Reactome", x = 0, y = 0.55, vjust = 1.5, angle = 90, size = 15) +
  draw_plot(legend, 0, 0, 1, 0.1)

# Print the final plot
print(final_plot_with_labels)

pdf('mean_cor/semi_PD_sparse/figures/compare_reactome_v2.pdf', width = 11, height = 7, onefile = T)
print(final_plot_with_labels)
dev.off()
