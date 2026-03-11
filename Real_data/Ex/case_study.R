library(pheatmap)
library(biomaRt)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(pheatmap)
library(WGCNA)
library(org.Hs.eg.db)
library(ggpubr)
# library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(stringr)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source('/gpfs/gibbs/pi/zhao/xs282/AD/AD_functions.r')


gene_name <- readRDS("marginal_fit/ROSMAP_NC_Ex_cor_gene.rds")
ori_ests <- matrix(NA, nrow = length(gene_name), ncol=length(gene_name),
                   dimnames = list(gene_name, gene_name))

marginal_fit_ROSMAP = readRDS('marginal_fit/ROSMAP_NC_Ex_marginal_fit.rds')
mu_ROSMAP <- marginal_fit_ROSMAP[gene_name,]$mu

ncor_gene <- length(gene_name)
mu_col_ROSMAP <- matrix(mu_ROSMAP,ncor_gene,ncor_gene,byrow = T)
mu_row_ROSMAP <- matrix(mu_ROSMAP,ncor_gene,ncor_gene)

tri = upper.tri(ori_ests, diag = FALSE)
idxs = which(tri, arr.ind = T)
estimate_ROSMAP <- data.frame(id1=rownames(ori_ests)[idxs[,1]],
                              id2=colnames(ori_ests)[idxs[,2]],
                              mu_col_ROSMAP=mu_col_ROSMAP[tri],
                              mu_row_ROSMAP=mu_row_ROSMAP[tri])
estimate_ROSMAP$log10mean_mu_ROSMAP <- log10(sqrt(10^estimate_ROSMAP$mu_col_ROSMAP*10^estimate_ROSMAP$mu_row_ROSMAP))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=rownames(ori_ests),mart= mart)
G_list <- G_list %>% group_by(hgnc_symbol) %>%
  dplyr::slice(1) %>% ungroup()
unmapped <- rownames(ori_ests)[!rownames(ori_ests) %in% G_list$hgnc_symbol]
ensembl <- c(unmapped, G_list$ensembl_gene_id)
names(ensembl) <- c(unmapped, G_list$hgnc_symbol)
estimate_ROSMAP$id1 <- ensembl[estimate_ROSMAP$id1]
estimate_ROSMAP$id2 <- ensembl[estimate_ROSMAP$id2]
estimate_ROSMAP <- estimate_ROSMAP %>%
  mutate(grp = paste(pmax(id1, id2), pmin(id1, id2), sep = "_"))


# p-value based approach
p_val <- readRDS("real/ex/ROSMAP_Ex_norm_p.rds")
ROSMAP_est <- readRDS("real/ex/est_cor.rds")
p_val$cscore_p <- ROSMAP_est$cscore_p
ROSMAP_p_adj <- as.data.frame(apply(p_val, 2, function(x){p.adjust(x, method = "BH")}))

p_cor <- ROSMAP_p_adj<0.05
colSums(p_cor)
p_cor_exp <- cbind(p_cor, estimate_ROSMAP)


# cor strength based approach - variable threshold
ROSMAP_est <- readRDS("real/ex/est_cor.rds")
est_sub <- ROSMAP_est[c(10:13, 15:18)]
est_sub$cscore_p <- est_sub$cscore_est

rank_thresh <- colSums(p_cor)
names(rank_thresh) <- colnames(p_cor)
cor_strength_cor <- p_cor
for (i in colnames(p_cor)){
  cor_strength_cor[,i] <- abs(est_sub[,i])>quantile(abs(est_sub[,i]),
                                                    1-rank_thresh[i]/length(est_sub[,i]))
}

colSums(cor_strength_cor)
cor_strength_cor_exp <- cbind(cor_strength_cor, estimate_ROSMAP)

# convert vector to matrix ------------------------------------------------
upp2matrix <- function(data, value_col) {
  wide_df <- data %>%
    dplyr::select(id1, id2, all_of(value_col)) %>%
    pivot_wider(names_from = id2, values_from = all_of(value_col))
  gene_names <- wide_df$id1
  mat <- as.matrix(wide_df[, -1])
  rownames(mat) <- gene_names

  final_mat <- mat[rownames(mat), rownames(mat)]
  final_mat[lower.tri(final_mat)] <- t(final_mat)[lower.tri(final_mat)]
  return(final_mat)
}

unique_gene = unique(c(estimate_ROSMAP$id1,estimate_ROSMAP$id2))
cand_df = data.frame(ana_prn = rep(1, length(unique_gene)),
                     id1 = unique_gene,
                     id2 = unique_gene)
p_cor_exp_ana_prn <- cbind(p_cor[,"ana_prn"], estimate_ROSMAP[,c("id1", "id2")])
colnames(p_cor_exp_ana_prn)[1] = "ana_prn"
p_cor_exp_ana_prn = rbind(p_cor_exp_ana_prn, cand_df)
ana_prn_p = upp2matrix(p_cor_exp_ana_prn, "ana_prn")


cor_strength_cor_exp_ana_prn <- cbind(cor_strength_cor[,"ana_prn"],
                                      estimate_ROSMAP[,c("id1", "id2")])
colnames(cor_strength_cor_exp_ana_prn)[1] = "ana_prn"
cor_strength_cor_exp_ana_prn = rbind(cor_strength_cor_exp_ana_prn, cand_df)
ana_prn_cor_strength = upp2matrix(cor_strength_cor_exp_ana_prn, "ana_prn")


est_ana_prn <- cbind(ROSMAP_est[,"ana_prn"], estimate_ROSMAP[,c("id1", "id2")])
colnames(est_ana_prn)[1] = "ana_prn"
est_ana_prn = rbind(est_ana_prn, cand_df)
est_ana_prn = upp2matrix(est_ana_prn, "ana_prn")

# filtered cor matrix
est_ana_prn_filtered_p = est_ana_prn
est_ana_prn_filtered_p[ana_prn_p==0] = 0

est_ana_prn_filtered_cor_strength = est_ana_prn
est_ana_prn_filtered_cor_strength[ana_prn_cor_strength==0] = 0

### cluster and enrichment analysis ------------------------

p_clu = Clustering(abs(est_ana_prn_filtered_p), colnames(est_ana_prn_filtered_p), power = 1, TOM_clustering = T)
cor_strength_clu = Clustering(abs(est_ana_prn_filtered_cor_strength), colnames(est_ana_prn_filtered_cor_strength), power = 1, TOM_clustering = T)

all_gene = marginal_fit_ROSMAP$gene
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=all_gene,mart= mart)
G_list <- G_list %>% group_by(hgnc_symbol) %>%
  dplyr::slice(1) %>% ungroup()
unmapped <- all_gene[!all_gene %in% G_list$hgnc_symbol]
ensembl <- c(unmapped, G_list$ensembl_gene_id)
names(ensembl) <- c(unmapped, G_list$hgnc_symbol)
saveRDS(ensembl, "BIB_R1/case_study/ex_case_study_gene_map.rds")

all_gene_id = ensembl[all_gene]
p_ego = EnrichGO_ensembl_rm0(p_clu, universe = all_gene_id, pvalueCutoff=0.05)
saveRDS(p_ego, "BIB_R1/case_study/ex_case_study_p_geo.rds")
cor_strength_ego = EnrichGO_ensembl_rm0(cor_strength_clu, universe = all_gene_id, pvalueCutoff=0.05)
saveRDS(cor_strength_ego, "BIB_R1/case_study/ex_case_study_cor_strength_geo.rds")

all_clusters <- names(p_ego[[1]])
plot_data <- lapply(all_clusters, function(clu_name) {
  res <- p_ego[[1]][[clu_name]]
  df <- as.data.frame(res)
  if(nrow(df) > 0) {
    df <- head(df, 5) # Extract top 5
    df$Cluster <- clu_name
    df$Ratio <- sapply(df$GeneRatio, function(x) {
      parts <- as.numeric(unlist(strsplit(x, "/")))
      parts[1] / parts[2]
    })
    return(df)
  } else {
    return(NULL)
  }
}) %>% bind_rows()

plot_data$Cluster <- factor(plot_data$Cluster, levels = all_clusters)
plot_data <- plot_data %>%
  mutate(Description = str_to_sentence(Description)) %>%
  arrange(Cluster, p.adjust) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

common_p_limits <- c(0, 0.05)
common_ratio_limits <- c(0, 0.5)
options(repr.plot.width = 8, repr.plot.height = 10)
p1 = ggplot(plot_data, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Ratio, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue",
                       limits = common_p_limits) +
  scale_size(limits = common_ratio_limits)+
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
  theme_bw(base_size = 13) + # Global text size increase
  theme(
    # Rotate Cluster names (X-axis) to 60 degrees
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),

    # Space between labels and title
    axis.title.y = element_text(margin = margin(r = 20)),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20))+
  labs(title="P-value-based")

all_clusters <- names(cor_strength_ego[[1]])
plot_data <- lapply(all_clusters, function(clu_name) {
  res <- cor_strength_ego[[1]][[clu_name]]
  df <- as.data.frame(res)
  if(nrow(df) > 0) {
    df <- head(df, 5) # Extract top 5
    df$Cluster <- clu_name
    df$Ratio <- sapply(df$GeneRatio, function(x) {
      parts <- as.numeric(unlist(strsplit(x, "/")))
      parts[1] / parts[2]
    })
    return(df)
  } else {
    return(NULL)
  }
}) %>% bind_rows()

plot_data$Cluster <- factor(plot_data$Cluster, levels = all_clusters)
plot_data <- plot_data %>%
  mutate(Description = str_to_sentence(Description)) %>%
  arrange(Cluster, p.adjust) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

options(repr.plot.width = 8, repr.plot.height = 10)
p2 = ggplot(plot_data, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Ratio, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue",
                       limits = common_p_limits) +
  scale_size(limits = common_ratio_limits)+
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
  theme_bw(base_size = 13) + # Global text size increase
  theme(
    # Rotate Cluster names (X-axis) to 60 degrees
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),

    # Space between labels and title
    axis.title.y = element_text(margin = margin(r = 20)),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20))+
  labs(title="Correlation-strength-based")


pdf('BIB_R1/modify_plot/ex_case_study.pdf', width = 16, height = 7.5, onefile = T)
ggarrange(p1,p2, ncol=2, nrow=1, common.legend=T, legend="right")
dev.off()
