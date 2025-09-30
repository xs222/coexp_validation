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
library(grid)
library(ks)
set.seed(11272023)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("AFinal/cscore_real_data_function.R")

# read expression data file
marginal_fit_ROSMAP = readRDS('revision/estimate_cor_lung/lung_Monocyte_marginal_fit.rds')

# reactome ------------------------------------------------------------------------------
hs_filter <- readRDS("Reactome/Processed_Reactome_10_07_2023.rds")
hs_filter <- hs_filter %>%
  mutate(grp = paste(pmax(V1, V2), pmin(V1, V2), sep = "_"))
exp_map <- with(marginal_fit_ROSMAP, setNames(mu, gene))
exp_map[is.infinite(exp_map)] <- min(exp_map[!is.infinite(exp_map)])-1
hs_filter$log10exp1 <- exp_map[hs_filter$V1]
hs_filter$log10exp2 <- exp_map[hs_filter$V2]
hs_filter$log10mean_mu <- log10(sqrt(10^(hs_filter$log10exp1+hs_filter$log10exp2)))
hs_exp_pair <- hs_filter[!(is.na(hs_filter$log10exp1)|is.na(hs_filter$log10exp2)),]


# overall performance ----------------------------------------------------------
### compared with background gene, the gene in biological network are more likely to be highly expressed-------
biological_net <- hs_exp_pair
intersect_gene <- unique(c(biological_net$V1, biological_net$V2))

# The actual bias is based on the gene set
marginal_fit_ROSMAP <- marginal_fit_ROSMAP[order(marginal_fit_ROSMAP$mu, decreasing=T),]

plist <- list()
# top <- c(100, 500, 1000, 1500, 2000, 5000)
top <- c(100, 1000, 5000)
for (i in 1:length(top)){
  gene_set <- marginal_fit_ROSMAP$gene[1:top[i]]
  hs_filter <- biological_net[biological_net$V1 %in% gene_set & biological_net$V2 %in% gene_set,]

  intersect_gene <- unique(c(hs_filter$V1, hs_filter$V2))
  exp_product_matrix <- sqrt(outer(10^(exp_map[intersect_gene]),10^(exp_map[intersect_gene])))
  background_pair_exp <- data.frame(log10mean_mu=log10(exp_product_matrix[upper.tri(exp_product_matrix)]))
  background_pair_exp$group <- "Background"
  string_pair_exp <- data.frame(log10mean_mu=hs_filter$log10mean_mu)
  string_pair_exp$group <- "Reactome"
  pair_exp <- rbind(background_pair_exp, string_pair_exp)
  plist[[i]] <- ggplot(pair_exp, aes(x=log10mean_mu, fill=group, color=group)) +
    geom_histogram(aes(y=..density..), alpha=0.5,
                   position="identity")+
    geom_density(alpha=.2) +
    theme_classic()+
    theme(legend.position = "bottom")+
    labs(title=paste0("Top ", top[i], " \nexpressed genes"), fill="", color="")
}


# The actual bias is based on the gene set highly variable --------------------------
log10mu <- marginal_fit_ROSMAP$mu
mu <- 10^log10mu

log10alpha <- marginal_fit_ROSMAP$alpha
alpha <- 10^log10alpha
beta <- mu/alpha

marginal_fit_ROSMAP$var <- alpha*beta^2
marginal_fit_ROSMAP <- marginal_fit_ROSMAP[order(marginal_fit_ROSMAP$var, decreasing=T),]
marginal_fit_ROSMAP$order_var <- order(marginal_fit_ROSMAP$var, decreasing=T)
marginal_fit_ROSMAP$order_exp <- order(marginal_fit_ROSMAP$mu, decreasing=T)

top <- c(100, 1000, 5000)
for (i in 1:length(top)){
  gene_set <- marginal_fit_ROSMAP$gene[1:top[i]]
  hs_filter <- biological_net[biological_net$V1 %in% gene_set & biological_net$V2 %in% gene_set,]

  intersect_gene <- unique(c(hs_filter$V1, hs_filter$V2))
  exp_product_matrix <- sqrt(outer(10^(exp_map[intersect_gene]),10^(exp_map[intersect_gene])))
  background_pair_exp <- data.frame(log10mean_mu=log10(exp_product_matrix[upper.tri(exp_product_matrix)]))
  background_pair_exp$group <- "Background"
  string_pair_exp <- data.frame(log10mean_mu=hs_filter$log10mean_mu)
  string_pair_exp$group <- "Reactome"
  pair_exp <- rbind(background_pair_exp, string_pair_exp)
  plist[[i+length(top)]] <- ggplot(pair_exp, aes(x=log10mean_mu, fill=group, color=group)) +
    geom_histogram(aes(y=..density..), alpha=0.5,
                   position="identity")+
    geom_density(alpha=.2) +
    theme_classic()+
    theme(legend.position = "bottom")+
    labs(title=paste0("Top ", top[i], " \nvariable genes"), fill="", color="")
}
p1 <- ggarrange(plotlist = plist, ncol = 6, nrow=1, common.legend = T, legend="bottom",
                labels = c("G","H","I", "J", "K", "L"))

# string ----------------------------------------------------------------------------------------------------
hs_filter <- readRDS("STRING/hs_filter_10_4_2023.rds")
hs_filter <- hs_filter %>%
  mutate(grp = paste(pmax(protein1, protein2), pmin(protein1, protein2), sep = "_"))
exp_map <- with(marginal_fit_ROSMAP, setNames(mu, gene))
exp_map[is.infinite(exp_map)] <- min(exp_map[!is.infinite(exp_map)])-1
hs_filter$log10exp1 <- exp_map[hs_filter$protein1]
hs_filter$log10exp2 <- exp_map[hs_filter$protein2]
hs_filter$log10mean_mu <- log10(sqrt(10^(hs_filter$log10exp1+hs_filter$log10exp2)))
hs_exp_pair <- hs_filter[!(is.na(hs_filter$log10exp1)|is.na(hs_filter$log10exp2)),]


# overall performance ----------------------------------------------------------
### compared with background gene, the gene in biological network are more likely to be highly expressed-------
biological_net <- hs_exp_pair[hs_exp_pair$combined_score>500,]
intersect_gene <- unique(c(biological_net$protein1, biological_net$protein2))

# The actual bias is based on the gene set
marginal_fit_ROSMAP <- marginal_fit_ROSMAP[order(marginal_fit_ROSMAP$mu, decreasing=T),]

plist_string <- list()
# top <- c(100, 500, 1000, 1500, 2000, 5000)
top <- c(100, 1000, 5000)
for (i in 1:length(top)){
  gene_set <- marginal_fit_ROSMAP$gene[1:top[i]]
  hs_filter <- biological_net[biological_net$protein1 %in% gene_set & biological_net$protein2 %in% gene_set,]

  intersect_gene <- unique(c(hs_filter$protein1, hs_filter$protein2))
  exp_product_matrix <- sqrt(outer(10^(exp_map[intersect_gene]),10^(exp_map[intersect_gene])))
  background_pair_exp <- data.frame(log10mean_mu=log10(exp_product_matrix[upper.tri(exp_product_matrix)]))
  background_pair_exp$group <- "Background"
  string_pair_exp <- data.frame(log10mean_mu=hs_filter$log10mean_mu)
  string_pair_exp$group <- "STRING"
  pair_exp <- rbind(background_pair_exp, string_pair_exp)
  plist_string[[i]] <- ggplot(pair_exp, aes(x=log10mean_mu, fill=group, color=group)) +
    geom_histogram(aes(y=..density..), alpha=0.5,
                   position="identity")+
    geom_density(alpha=.2) +
    theme_classic()+
    theme(legend.position = "bottom")+
    labs(title=paste0("Top ", top[i], " \nexpressed genes"), fill="", color="")
}

# The actual bias is based on the gene set highly variable --------------------------
log10mu <- marginal_fit_ROSMAP$mu
mu <- 10^log10mu

log10alpha <- marginal_fit_ROSMAP$alpha
alpha <- 10^log10alpha
beta <- mu/alpha

marginal_fit_ROSMAP$var <- alpha*beta^2
marginal_fit_ROSMAP <- marginal_fit_ROSMAP[order(marginal_fit_ROSMAP$var, decreasing=T),]
marginal_fit_ROSMAP$order_var <- order(marginal_fit_ROSMAP$var, decreasing=T)
marginal_fit_ROSMAP$order_exp <- order(marginal_fit_ROSMAP$mu, decreasing=T)

top <- c(100, 1000, 5000)
for (i in 1:length(top)){
  gene_set <- marginal_fit_ROSMAP$gene[1:top[i]]
  hs_filter <- biological_net[biological_net$protein1 %in% gene_set & biological_net$protein2 %in% gene_set,]

  intersect_gene <- unique(c(hs_filter$protein1, hs_filter$protein2))
  exp_product_matrix <- sqrt(outer(10^(exp_map[intersect_gene]),10^(exp_map[intersect_gene])))
  background_pair_exp <- data.frame(log10mean_mu=log10(exp_product_matrix[upper.tri(exp_product_matrix)]))
  background_pair_exp$group <- "Background"
  string_pair_exp <- data.frame(log10mean_mu=hs_filter$log10mean_mu)
  string_pair_exp$group <- "STRING"
  pair_exp <- rbind(background_pair_exp, string_pair_exp)
  plist_string[[i+length(top)]] <- ggplot(pair_exp, aes(x=log10mean_mu, fill=group, color=group)) +
    geom_histogram(aes(y=..density..), alpha=0.5,
                   position="identity")+
    geom_density(alpha=.2) +
    theme_classic()+
    theme(legend.position = "bottom")+
    labs(title=paste0("Top ", top[i], " \nvariable genes"), fill="", color="")
}

p2 <- ggarrange(plotlist = plist_string, ncol = 6, nrow=1, common.legend = T, legend="bottom",
                labels = c("A","B","C", "D", "E", "F"))

fig1 = annotate_figure(ggarrange(p2, ncol=1, nrow=1),
  top = text_grob("STRING", 
                  color = "black", face = "bold", size = 16))
fig2 = annotate_figure(ggarrange(p1, ncol=1, nrow=1),
  top = text_grob("---------------------------------------------------------------------------------------------------\nReactome", 
                  color = "black", face = "bold", size = 16))

pdf('/gpfs/gibbs/pi/zhao/xs282/validation/revision/modify_plot/biological_bias_lung_mono_supp.pdf', width = 13, height = 6, onefile = T)
ggarrange(fig1, fig2, ncol=1, nrow=2, heights = c(5,5.05))
dev.off()






