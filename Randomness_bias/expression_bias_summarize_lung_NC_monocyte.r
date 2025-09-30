library(gridExtra)
library(tidyverse)
library(magrittr)
library(rtracklayer)
library(biomaRt)
library(ggplot2)
library(data.table)
library(ggpubr)
library(tidyverse)
library(grid)
library(ks)
set.seed(11272023)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
source("AFinal/cscore_real_data_function.R")

rand_overlap = function(p, biological_net, md_names){
    p$biological_net = ifelse(p$grp %in% biological_net$grp, 1, 0)
    res = matrix(nrow=length(md_names), ncol=5)
    rownames(res) = md_names
    colnames(res) = c("obs", "random_mean_adjust_exp", "random_mean","random_var_adjust_exp", "random_var")
    for (i in md_names){
        res[i,1] = sum(p[,i] & p$biological_net==1)
        cand = p %>% group_by(exp_gp) %>%
            summarize(s = sum(biological_net), n= n(), c=sum(!!sym(i), na.rm = TRUE))
        res[i,2] = sum(cand$s*cand$c/cand$n)
        sds = cand$s*cand$c/cand$n*(1-cand$s/cand$n)*(cand$n-cand$c)/(cand$n-1)
        sds[is.na(sds)] = 0
        res[i,4] = sum(sds)
    }
    p_sub = p[, md_names]
    res[,3] = colSums(p_sub)*sum(p$biological_net)/nrow(p)
    res[,5] = colSums(p_sub)*sum(p$biological_net)/nrow(p)*(1-sum(p$biological_net)/nrow(p))*(nrow(p)-colSums(p_sub))/(nrow(p)-1)
    res = as.data.frame(res)
    res$adjust_exp = res$obs - res$random_mean_adjust_exp
    res$adjust = res$obs - res$random_mean
    res$rank_exp = rank(-res$adjust_exp)
    res$rank = rank(-res$adjust)
    res$rank_ori = rank(-res$obs)
    res$z_exp = res$adjust_exp/sqrt(res$random_var_adjust_exp)
    res$z = res$adjust/sqrt(res$random_var)
    res$rank_z_exp = rank(-res$z_exp)
    res$rank_z = rank(-res$z)
    res$md = rownames(res)
    res$Method <- recode(res$md,sct="sctransform", prn="Pearson", spr="Spearman",
                       propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                       noise="Noise \nRegularization", cscore_est="CS-CORE Empirical")
    return(res)
}

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

gene_exp <- data.frame(log10mu=exp_map, gene=names(exp_map))
gene_exp$BiologicalNet <- ifelse(gene_exp$gene %in% intersect_gene, "Within", "Not within")

# set.seed(2252024)
# ks.test(gene_exp$log10mu[gene_exp$BiologicalNet=="Within"],
#         gene_exp$log10mu[gene_exp$BiologicalNet=="Not within"])
# idx1 <- sample(1:sum(gene_exp$BiologicalNet=="Within"), 10000)
# idx2 <- sample(1:sum(gene_exp$BiologicalNet=="Not within"), 10000)
# kde.test(gene_exp$log10mu[gene_exp$BiologicalNet=="Within"][idx1],
#          gene_exp$log10mu[gene_exp$BiologicalNet=="Not within"][idx2])

# Create a text
grob_gene <- grobTree(textGrob("KS test: p-value < 2.2e-16 \nKDE test: p-value < 2.2e-16", x=0.5,  y=0.8,
                          gp=gpar(col="black", fontsize=12, fontface="italic")))
# Plot
p1 <- ggplot(gene_exp, aes(x=log10mu, fill=BiologicalNet, color=BiologicalNet)) +
  geom_histogram(aes(y=..density..), alpha=0.5,
                 position="identity")+
  geom_density(alpha=.2) +
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(title="Reactome: gene")+ annotation_custom(grob_gene)
# p1

### compared with background gene pairs, the gene pairs in biological network are more likely to be highly expressed-------
exp_product_matrix <- sqrt(outer(10^(exp_map[intersect_gene]),10^(exp_map[intersect_gene])))

background_pair_exp <- data.frame(log10mean_mu=log10(exp_product_matrix[upper.tri(exp_product_matrix)]))
background_pair_exp$group <- "Background"
string_pair_exp <- data.frame(log10mean_mu=biological_net$log10mean_mu)
string_pair_exp$group <- "Reactome"
pair_exp <- rbind(background_pair_exp, string_pair_exp)
# wilcox.test(background_pair_exp$log10mean_mu,
#             string_pair_exp$log10mean_mu,alternative="greater")

# ks.test(pair_exp$log10mean_mu[pair_exp$group=="Background"],
#         pair_exp$log10mean_mu[pair_exp$group=="Reactome"])
# idx1 <- sample(1:sum(pair_exp$group=="Background"), 10000)
# idx2 <- sample(1:sum(pair_exp$group=="Reactome"), 10000)
# kde.test(pair_exp$log10mean_mu[pair_exp$group=="Background"][idx1],
#          pair_exp$log10mean_mu[pair_exp$group=="Reactome"][idx2])

grob_gene_pair <- grobTree(textGrob("KS test: p-value < 2.2e-16 \nKDE test: p-value < 2.2e-16", x=0.5,  y=0.8,
                               gp=gpar(col="black", fontsize=12, fontface="italic")))
p2 <- ggplot(pair_exp, aes(x=log10mean_mu, fill=group, color=group)) +
  geom_histogram(aes(y=..density..), alpha=0.5,
                 position="identity")+
  geom_density(alpha=.2) +
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(title="Reactome: gene pair", fill="", color="")+ annotation_custom(grob_gene_pair)

# string ---------------------------------------------------------------------
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

gene_exp <- data.frame(log10mu=exp_map, gene=names(exp_map))
gene_exp$BiologicalNet <- ifelse(gene_exp$gene %in% intersect_gene, "Within", "Not within")

set.seed(2252024)
# ks.test(gene_exp$log10mu[gene_exp$BiologicalNet=="Within"],
#         gene_exp$log10mu[gene_exp$BiologicalNet=="Not within"])
# idx1 <- sample(1:sum(gene_exp$BiologicalNet=="Within"), 10000)
# idx2 <- sample(1:sum(gene_exp$BiologicalNet=="Not within"), 10000)
# kde.test(gene_exp$log10mu[gene_exp$BiologicalNet=="Within"][idx1],
#          gene_exp$log10mu[gene_exp$BiologicalNet=="Not within"][idx2])

# Create a text
grob_gene <- grobTree(textGrob("KS test: p-value < 2.2e-16 \nKDE test: p-value < 2.2e-16", x=0.5,  y=0.8,
                          gp=gpar(col="black", fontsize=12, fontface="italic")))
# Plot
p1_string <- ggplot(gene_exp, aes(x=log10mu, fill=BiologicalNet, color=BiologicalNet)) +
  geom_histogram(aes(y=..density..), alpha=0.5,
                 position="identity")+
  geom_density(alpha=.2) +
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(title="STRING: gene")+ annotation_custom(grob_gene)
# p1

### compared with background gene pairs, the gene pairs in biological network are more likely to be highly expressed-------
exp_product_matrix <- sqrt(outer(10^(exp_map[intersect_gene]),10^(exp_map[intersect_gene])))

background_pair_exp <- data.frame(log10mean_mu=log10(exp_product_matrix[upper.tri(exp_product_matrix)]))
background_pair_exp$group <- "Background"
string_pair_exp <- data.frame(log10mean_mu=biological_net$log10mean_mu)
string_pair_exp$group <- "STRING"
pair_exp <- rbind(background_pair_exp, string_pair_exp)
# wilcox.test(background_pair_exp$log10mean_mu,
#             string_pair_exp$log10mean_mu,alternative="greater")

# ks.test(pair_exp$log10mean_mu[pair_exp$group=="Background"],
#         pair_exp$log10mean_mu[pair_exp$group=="STRING"])
# idx1 <- sample(1:sum(pair_exp$group=="Background"), 10000)
# idx2 <- sample(1:sum(pair_exp$group=="STRING"), 10000)
# kde.test(pair_exp$log10mean_mu[pair_exp$group=="Background"][idx1],
#          pair_exp$log10mean_mu[pair_exp$group=="STRING"][idx2])

grob_gene_pair <- grobTree(textGrob("KS test: p-value < 2.2e-16 \nKDE test: p-value < 2.2e-16", x=0.5,  y=0.8,
                               gp=gpar(col="black", fontsize=12, fontface="italic")))
p2_string <- ggplot(pair_exp, aes(x=log10mean_mu, fill=group, color=group)) +
  geom_histogram(aes(y=..density..), alpha=0.5,
                 position="identity")+
  geom_density(alpha=.2) +
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(title="STRING: gene pair", fill="", color="")+ annotation_custom(grob_gene_pair)


# zscore -----------------------------------------------------------------
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
# STRING
hs_filter <- readRDS("STRING/hs_filter_10_4_2023.rds")
hs_filter <- hs_filter %>%
  mutate(grp = paste(pmax(protein1, protein2), pmin(protein1, protein2), sep = "_"))
string <- hs_filter[hs_filter$combined_score>500,]

# reactome
hs_filter <- readRDS("Reactome/Processed_Reactome_10_07_2023.rds")
hs_filter <- hs_filter %>%
  mutate(grp = paste(pmax(V1, V2), pmin(V1, V2), sep = "_"))
reactome <- hs_filter

gene_name = readRDS("revision/estimate_cor_lung/lung_NC_mono_cor_gene.rds")
ori_ests = matrix(ncol=length(gene_name), nrow=length(gene_name))
colnames(ori_ests) = gene_name
rownames(ori_ests) = gene_name

marginal_fit_lung = readRDS('revision/estimate_cor_lung/lung_Monocyte_marginal_fit.rds')
mu_lung <- marginal_fit_lung[gene_name,]$mu

ncor_gene <- length(gene_name)
mu_col_lung <- matrix(mu_lung,ncor_gene,ncor_gene,byrow = T)
mu_row_lung <- matrix(mu_lung,ncor_gene,ncor_gene)

tri = upper.tri(ori_ests, diag = FALSE)
idxs = which(tri, arr.ind = T)
estimate_lung <- data.frame(id1=rownames(ori_ests)[idxs[,1]],
                              id2=colnames(ori_ests)[idxs[,2]],
                              mu_col_lung=mu_col_lung[tri],
                              mu_row_lung=mu_row_lung[tri])
estimate_lung$log10mean_mu_lung <- log10(sqrt(10^estimate_lung$mu_col_lung*10^estimate_lung$mu_row_lung))
estimate_lung <- estimate_lung %>%
  mutate(grp = paste(pmax(id1, id2), pmin(id1, id2), sep = "_"))
# p-value based approach
p_val <- readRDS("revision/estimate_cor_lung/mono/norm_p.rds")
lung_est <- readRDS("revision/estimate_cor_lung/mono/est_cor.rds")
p_val$cscore_p <- lung_est$cscore_p
lung_p_adj <- as.data.frame(apply(p_val, 2, function(x){p.adjust(x, method = "BH")}))

thresh <- 0.05
p_cor <- lung_p_adj<thresh
colSums(p_cor)
p_cor_exp <- cbind(p_cor, estimate_lung)
breaks <- seq(from = min(p_cor_exp$log10mean_mu_lung), to = max(p_cor_exp$log10mean_mu_lung), by = 0.05)
p_cor_exp$exp_gp <- cut(x = p_cor_exp$log10mean_mu_lung, breaks = breaks)

res_string_lung = rand_overlap(p_cor_exp, string, colnames(p_cor))
res_reactome_lung = rand_overlap(p_cor_exp, reactome, colnames(p_cor))

dumbbell_data <- res_string_lung %>%
  mutate(method_label = paste0(Method, " \n(Obs = ", obs, ")")) %>%
  mutate(method_label = fct_reorder(method_label, z_exp))

# Generate the dumbbell plot
dumbbell_string <- ggplot(dumbbell_data, aes(y = method_label)) +
  # Add points for z_exp and z
  geom_point(aes(x = z_exp, color = "Adjust exp"), size = 4) +
  geom_point(aes(x = z, color = "No adjustment"), size = 4) +
  # Add a zero line for reference
  geom_vline(xintercept = qnorm(0.95), linetype = "dashed", color = "grey50") +
  scale_color_manual(name = "Score Type", values = c("Adjust exp" = "#0072B2", "No adjustment" = "#D55E00")) +
  labs(
    title = "STRING",
    x = "Z-Score",
    y = "Method"
  ) +
  theme_minimal()

dumbbell_data <- res_reactome_lung %>%
  mutate(method_label = paste0(Method, " \n(Obs = ", obs, ")")) %>%
  mutate(method_label = fct_reorder(method_label, z_exp))

# Generate the dumbbell plot
dumbbell_reac <- ggplot(dumbbell_data, aes(y = method_label)) +
  # Add points for z_exp and z
  geom_point(aes(x = z_exp, color = "Adjust exp"), size = 4) +
  geom_point(aes(x = z, color = "No adjustment"), size = 4) +
  # Add a zero line for reference
  geom_vline(xintercept = qnorm(0.95), linetype = "dashed", color = "grey50") +
  scale_color_manual(name = "Score Type", values = c("Adjust exp" = "#0072B2", "No adjustment" = "#D55E00")) +
  labs(
    title = "Reactome",
    x = "Z-Score",
    y = "Method"
  ) +
  theme_minimal()

plot_rep_ori1 <- ggarrange(p1_string, p2_string, dumbbell_string, ncol=3,nrow=1,
                           widths = c(1, 1,1.5),
          labels = c("A", "B", "C"))
# fig1 = annotate_figure(ggarrange(plot_rep_ori1, ncol=1, nrow=1),
#   top = text_grob("STRING", 
#                   color = "black", face = "bold", size = 16))

plot_rep_ori2 <- ggarrange(p1,p2,dumbbell_reac, ncol=3,nrow=1,
                           widths = c(1, 1,1.5),
          labels = c("D","E", "F"))
# fig2 = annotate_figure(ggarrange(plot_rep_ori2, ncol=1, nrow=1),
#   top = text_grob("---------------------------------------------------------------------------------------------------\nReactome", 
#                   color = "black", face = "bold", size = 16))

pdf('/gpfs/gibbs/pi/zhao/xs282/validation/revision/modify_plot/biological_bias_lung_mono.pdf', width = 12, height = 7, onefile = T)
ggarrange(plot_rep_ori1, plot_rep_ori2, ncol=1, nrow=2, heights = c(5,5))
dev.off()

