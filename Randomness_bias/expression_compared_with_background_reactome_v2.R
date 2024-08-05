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
set.seed(11272023)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("/gpfs/gibbs/pi/zhao/xs282/validation/cscore_real_data_function.R")

# read expression data file
marginal_fit_ROSMAP = readRDS('marginal_fit/ROSMAP_NC_Oli_marginal_fit.rds')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=marginal_fit_ROSMAP$gene,mart= mart)
marginal_fit_ROSMAP <- left_join(marginal_fit_ROSMAP, G_list, by=c("gene"="hgnc_symbol"))



hs_filter <- readRDS("Reactome/Processed_Reactome_10_07_2023.rds")
hs_filter <- hs_filter %>%
  mutate(grp = paste(pmax(V1, V2), pmin(V1, V2), sep = "_"))
exp_map <- with(marginal_fit_ROSMAP, setNames(mu, ensembl_gene_id))
exp_map[is.infinite(exp_map)] <- min(exp_map[!is.infinite(exp_map)])-1
hs_filter$log10exp1 <- exp_map[hs_filter$V1]
hs_filter$log10exp2 <- exp_map[hs_filter$V2]
hs_filter$log10mean_mu <- log10(sqrt(10^(hs_filter$log10exp1+hs_filter$log10exp2)))
hs_exp_pair <- hs_filter[!(is.na(hs_filter$log10exp1)|is.na(hs_filter$log10exp2)),]


# overall performance ----------------------------------------------------------
### compared with background gene, the gene in biological network are more likely to be highly expressed-------
biological_net <- hs_exp_pair
intersect_gene <- unique(c(biological_net$V1, biological_net$V2))

gene_exp <- data.frame(log10mu=exp_map, gene=names(exp_map))
gene_exp$BiologicalNet <- ifelse(gene_exp$gene %in% intersect_gene, "Within", "Not within")

set.seed(2252024)
ks_t <- ks.test(gene_exp$log10mu[gene_exp$BiologicalNet=="Within"],
        gene_exp$log10mu[gene_exp$BiologicalNet=="Not within"])
idx1 <- sample(1:sum(gene_exp$BiologicalNet=="Within"), 10000)
idx2 <- sample(1:sum(gene_exp$BiologicalNet=="Not within"), 10000)
kde_t <- kde.test(gene_exp$log10mu[gene_exp$BiologicalNet=="Within"][idx1],
         gene_exp$log10mu[gene_exp$BiologicalNet=="Not within"][idx2])

# Create a text
grob_gene <- grobTree(textGrob("KS test: p-value < 2.2e-16 \nKDE test: p-value < 2.2e-16", x=0.6,  y=0.8,
                          gp=gpar(col="black", fontsize=11, fontface="italic")))
# Plot
p1 <- ggplot(gene_exp, aes(x=log10mu, fill=BiologicalNet, color=BiologicalNet)) +
  geom_histogram(aes(y=..density..), alpha=0.5,
                 position="identity")+
  geom_density(alpha=.2) +
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(title="Gene")+ annotation_custom(grob_gene)
p1

### compared with background gene pairs, the gene pairs in biological network are more likely to be highly expressed-------
exp_product_matrix <- sqrt(outer(10^(exp_map[intersect_gene]),10^(exp_map[intersect_gene])))

background_pair_exp <- data.frame(log10mean_mu=log10(exp_product_matrix[upper.tri(exp_product_matrix)]))
background_pair_exp$group <- "Background"
string_pair_exp <- data.frame(log10mean_mu=biological_net$log10mean_mu)
string_pair_exp$group <- "Reactome"
pair_exp <- rbind(background_pair_exp, string_pair_exp)
# wilcox.test(background_pair_exp$log10mean_mu,
#             string_pair_exp$log10mean_mu,alternative="greater")

ks.test(pair_exp$log10mean_mu[pair_exp$group=="Background"],
        pair_exp$log10mean_mu[pair_exp$group=="Reactome"])
idx1 <- sample(1:sum(pair_exp$group=="Background"), 10000)
idx2 <- sample(1:sum(pair_exp$group=="Reactome"), 10000)
kde.test(pair_exp$log10mean_mu[pair_exp$group=="Background"][idx1],
         pair_exp$log10mean_mu[pair_exp$group=="Reactome"][idx2])

grob_gene_pair <- grobTree(textGrob("KS test: p-value < 2.2e-16 \nKDE test: p-value < 2.2e-16", x=0.6,  y=0.8,
                               gp=gpar(col="black", fontsize=11, fontface="italic")))
p2 <- ggplot(pair_exp, aes(x=log10mean_mu, fill=group, color=group)) +
  geom_histogram(aes(y=..density..), alpha=0.5,
                 position="identity")+
  geom_density(alpha=.2) +
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(title="Gene pair", fill="", color="")+ annotation_custom(grob_gene_pair)
p2

p1_p2 <- ggarrange(p1,p2, ncol = 1, nrow=2, labels = c("A", "B"))

# The actual bias is based on the gene set
marginal_fit_ROSMAP <- marginal_fit_ROSMAP[order(marginal_fit_ROSMAP$mu, decreasing=T),]

plist <- list()
# top <- c(100, 500, 1000, 1500, 2000, 5000)
top <- c(100, 1000, 5000)
for (i in 1:length(top)){
  gene_set <- marginal_fit_ROSMAP$ensembl_gene_id[1:top[i]]
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
    labs(title=paste0("Top ", top[i], " expressed genes"), fill="", color="")
}

p3 <- ggarrange(plotlist = plist, ncol = 3, nrow=1, common.legend = T, legend="none",
                labels = c("C", "D", "E"))



# The actual bias is based on the gene set highly variable --------------------------
log10mu <- marginal_fit_ROSMAP$mu
mu <- 10^log10mu

# use the sampled mean and the trend between mean and alpha to get corresponding alpha
# km_Ex5 <- readRDS('marginal_fit/ROSMAP_NC_Oli_ks_fit_5.rds')
# fitted_trend <- data.frame(mu=km_Ex5$x, alpha=km_Ex5$y)
# log10alpha <- rep(NA,length(mu))
# for (i in 1:length(mu)){
#   idx <- which.min(abs(log10mu[i]-fitted_trend$mu))
#   log10alpha[i] <- fitted_trend$alpha[idx]
# }
log10alpha <- marginal_fit_ROSMAP$alpha
alpha <- 10^log10alpha
beta <- mu/alpha

marginal_fit_ROSMAP$var <- alpha*beta^2
marginal_fit_ROSMAP <- marginal_fit_ROSMAP[order(marginal_fit_ROSMAP$var, decreasing=T),]
marginal_fit_ROSMAP$order_var <- order(marginal_fit_ROSMAP$var, decreasing=T)
marginal_fit_ROSMAP$order_exp <- order(marginal_fit_ROSMAP$mu, decreasing=T)

plist <- list()
top <- c(100, 1000, 5000)
for (i in 1:length(top)){
  gene_set <- marginal_fit_ROSMAP$ensembl_gene_id[1:top[i]]
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
    labs(title=paste0("Top ", top[i], " variable genes"), fill="", color="")
}

p4 <- ggarrange(plotlist = plist, ncol = 3, nrow=1, common.legend = T, legend="bottom",
                labels = c("F","G","H"))


pdf('figures/v2/biological_bias_reactome.pdf', width = 13, height = 5, onefile = T)
ggarrange(p1_p2,ggarrange(p3,p4, ncol = 1, nrow=2,common.legend = T, legend="bottom", heights = c(0.85,1))
, ncol = 2, nrow=1,widths = c(0.5, 1.3), common.legend = T, legend="bottom")
dev.off()
