library(pheatmap)
library(biomaRt)
library(dplyr)
library(tidyverse)
library(ggpubr)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")

ori_ests <- readRDS("marginal_fit/ROSMAP_NC_Oli_cscore_cor1000.rds")
gene_name <- colnames(ori_ests$est)

marginal_fit_ROSMAP = readRDS('marginal_fit/ROSMAP_NC_Oli_marginal_fit.rds')
mu_ROSMAP <- marginal_fit_ROSMAP[gene_name,]$mu

ncor_gene <- length(gene_name)
mu_col_ROSMAP <- matrix(mu_ROSMAP,ncor_gene,ncor_gene,byrow = T)
mu_row_ROSMAP <- matrix(mu_ROSMAP,ncor_gene,ncor_gene)

tri = upper.tri(ori_ests$est, diag = FALSE)
idxs = which(tri, arr.ind = T)
estimate_ROSMAP <- data.frame(id1=rownames(ori_ests$est)[idxs[,1]],
                              id2=colnames(ori_ests$est)[idxs[,2]],
                              mu_col_ROSMAP=mu_col_ROSMAP[tri],
                              mu_row_ROSMAP=mu_row_ROSMAP[tri])
estimate_ROSMAP$log10mean_mu_ROSMAP <- log10(sqrt(10^estimate_ROSMAP$mu_col_ROSMAP*10^estimate_ROSMAP$mu_row_ROSMAP))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=rownames(ori_ests$est),mart= mart)
G_list <- G_list %>% group_by(hgnc_symbol) %>%
  dplyr::slice(1) %>% ungroup()
unmapped <- rownames(ori_ests$est)[!rownames(ori_ests$est) %in% G_list$hgnc_symbol]
ensembl <- c(unmapped, G_list$ensembl_gene_id)
names(ensembl) <- c(unmapped, G_list$hgnc_symbol)
estimate_ROSMAP$id1 <- ensembl[estimate_ROSMAP$id1]
estimate_ROSMAP$id2 <- ensembl[estimate_ROSMAP$id2]
estimate_ROSMAP <- estimate_ROSMAP %>%
  mutate(grp = paste(pmax(id1, id2), pmin(id1, id2), sep = "_"))


# p-value based approach
p_val <- readRDS("real/oli/ROSMAP_oli_norm_p.rds")
ROSMAP_est <- readRDS("real/oli/est_cor.rds")
p_val$cscore_p <- ROSMAP_est$cscore_p
ROSMAP_p_adj <- as.data.frame(apply(p_val, 2, function(x){p.adjust(x, method = "BH")}))


## count of shared go terms ----------------------------------------------------
go_data <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol", "go_id"),values=rownames(ori_ests$est),mart= mart)
go_list <- split(go_data$go_id, go_data$hgnc_symbol)

tri = upper.tri(ori_ests$est, diag = FALSE)
idxs = which(tri, arr.ind = T)
gene_pair <- data.frame(id1=rownames(ori_ests$est)[idxs[,1]],
                        id2=colnames(ori_ests$est)[idxs[,2]])

thresh <- 0.05
p_cor <- ROSMAP_p_adj<thresh
p_cor_exp <- cbind(p_cor, estimate_ROSMAP)

# cor strength based approach - variable threshold
ROSMAP_est <- readRDS("real/oli/est_cor.rds")
est_sub <- ROSMAP_est[c(10:13, 15:18)]
est_sub$cscore_p <- est_sub$cscore_est

rank_thresh <- colSums(p_cor)
names(rank_thresh) <- colnames(p_cor)
cor_strength_cor <- p_cor
for (i in colnames(p_cor)){
  cor_strength_cor[,i] <- abs(est_sub[,i])>quantile(abs(est_sub[,i]),
                                                    1-rank_thresh[i]/length(est_sub[,i]))
}
cor_strength_cor_exp <- cbind(cor_strength_cor, estimate_ROSMAP)

method <- c("sct", "ana_prn", "noise", "prn", "spr", "propr", "cscore_est", "cscore_p")
method2 <- c("sctransform", "Analytic PR", "Noise \nRegularization", "Pearson", "Spearman", "propr",
             "CS-CORE \n(Empirical)", "CS-CORE")
names(method2) <- c("sct", "ana_prn", "noise", "prn", "spr", "propr", "cscore_est", "cscore_p")
box_ls <- list()
for (i in method){
  common <- intersect(p_cor_exp$grp[p_cor_exp[,i]], cor_strength_cor_exp$grp[cor_strength_cor_exp[,i]])
  go_dat_p <- ROSMAP_est[,i]
  go_dat_p[!p_cor[,i] | p_cor_exp$grp %in% common] <- 0
  p_cor_list <- gene_pair[go_dat_p!=0,]
  p_cor_list$count <- apply(p_cor_list, 1, function(x){
    length(intersect(go_list[[x[1]]], go_list[[x[2]]]))
  })

  go_dat_top <- ROSMAP_est[,i]
  go_dat_top[!cor_strength_cor[,i] | cor_strength_cor_exp$grp %in% common] <- 0
  top_cor_list <- gene_pair[go_dat_top!=0,]
  top_cor_list$count <- apply(top_cor_list, 1, function(x){
    length(intersect(go_list[[x[1]]], go_list[[x[2]]]))
  })

  top_cor_list$group <- "Cor-strength"
  p_cor_list$group <- "P-value"
  comb_dat <-rbind(top_cor_list, p_cor_list)

  box_ls[[i]] <- ggplot(comb_dat, aes(x=group, y=count,fill=group))+
    geom_boxplot(alpha=0.5)+labs(title=method2[i], x="", fill="")+theme_bw()+
    scale_fill_manual(values = c("P-value"="darkblue", "Cor-strength"="darkred"))+
    theme(plot.title = element_text(hjust=0.5),legend.position = "bottom")
}
ggarrange(plotlist = box_ls, ncol=4, nrow = 2, common.legend = T, legend = "bottom")

legend <- get_legend(box_ls[[1]])
adjusted_theme <- theme(text = element_text(size = 17), legend.position = "none",
                        plot.margin = unit(c(0, 0.15, 0, 0), "cm"),
                        title =element_text(size=12),
                        plot.tag = element_text(size = 14, face = "bold", vjust = 1.3, hjust = -1.5),
                        plot.tag.position = c(0,1),)

plot_grid <- plot_grid(box_ls[[2]]+labs(tag="A")+adjusted_theme, box_ls[[8]]+labs(tag="B")+adjusted_theme,
          box_ls[[7]]+labs(tag="C")+adjusted_theme, box_ls[[3]]+labs(tag="D")+adjusted_theme,
          box_ls[[4]]+labs(tag="E")+adjusted_theme, box_ls[[6]]+labs(tag="F")+adjusted_theme,
          box_ls[[1]]+labs(tag="G")+adjusted_theme, box_ls[[5]]+labs(tag="H")+adjusted_theme,
          ncol = 4, nrow = 2, align = "hv")

final_plot_with_labels <- ggdraw() +
  draw_plot(plot_grid)

# Print the final plot
print(final_plot_with_labels)

pdf('real_data/oli_shared_GO_boxplot_05_v2.pdf', width = 12, height = 6.5, onefile = T)
print(final_plot_with_labels)
dev.off()



