library(pheatmap)
library(biomaRt)
library(dplyr)
library(tidyverse)
library(ggpubr)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")

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

reproduce_overlap = function(p, second_dat, md_names){
    res = matrix(nrow=length(md_names), ncol=5)
    rownames(res) = md_names
    colnames(res) = c("obs", "random_mean_adjust_exp", "random_mean","random_var_adjust_exp", "random_var")
    for (i in md_names){
        biological_net = second_dat[second_dat[,i],]
        p$biological_net = ifelse(p$grp %in% biological_net$grp, 1, 0)
        res[i,1] = sum(p[,i] & p$biological_net==1)
        cand = p %>% group_by(exp_gp) %>%
            summarize(s = sum(biological_net), n= n(), c=sum(!!sym(i), na.rm = TRUE))
        res[i,2] = sum(cand$s*cand$c/cand$n)
        sds = cand$s*cand$c/cand$n*(1-cand$s/cand$n)*(cand$n-cand$c)/(cand$n-1)
        sds[is.na(sds)] = 0
        res[i,4] = sum(sds)
        res[i,3] = sum(p[,i])*sum(p$biological_net)/nrow(p)
        res[i,5] = sum(p[,i])*sum(p$biological_net)/nrow(p)*(1-sum(p$biological_net)/nrow(p))*(nrow(p)-sum(p[,i]))/(nrow(p)-1)
    }
    
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

thresh <- 0.05
p_cor <- ROSMAP_p_adj<thresh
colSums(p_cor)
p_cor_exp <- cbind(p_cor, estimate_ROSMAP)

breaks <- seq(from = min(p_cor_exp$log10mean_mu_ROSMAP), to = max(p_cor_exp$log10mean_mu_ROSMAP), by = 0.05)
p_cor_exp$exp_gp <- cut(x = p_cor_exp$log10mean_mu_ROSMAP, breaks = breaks)

res_string = rand_overlap(p_cor_exp, string, colnames(p_cor))
res_reactome = rand_overlap(p_cor_exp, reactome, colnames(p_cor))

dumbbell_data <- res_string %>%
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

dumbbell_data <- res_reactome %>%
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

# reproducibility -------------------------------------------------------
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
estimate_ROSMAP <- estimate_ROSMAP %>%
  mutate(grp = paste(pmax(id1, id2), pmin(id1, id2), sep = "_"))


# p-value based approach
p_val <- readRDS("real/oli/ROSMAP_oli_norm_p.rds")
ROSMAP_est <- readRDS("real/oli/est_cor.rds")
p_val$cscore_p <- ROSMAP_est$cscore_p
ROSMAP_p_adj <- as.data.frame(apply(p_val, 2, function(x){p.adjust(x, method = "BH")}))

thresh <- 0.05
p_cor <- ROSMAP_p_adj<thresh
# p_cor <- p_val<thresh
colSums(p_cor)
p_cor_exp <- cbind(p_cor, estimate_ROSMAP)
breaks <- seq(from = min(p_cor_exp$log10mean_mu_ROSMAP), to = max(p_cor_exp$log10mean_mu_ROSMAP), by = 0.05)
p_cor_exp$exp_gp <- cut(x = p_cor_exp$log10mean_mu_ROSMAP, breaks = breaks)

setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
ori_ests = matrix(ncol=length(gene_name), nrow=length(gene_name))
colnames(ori_ests) = gene_name
rownames(ori_ests) = gene_name

tri = upper.tri(ori_ests, diag = FALSE)
idxs = which(tri, arr.ind = T)
estimate_pnas <- data.frame(id1=rownames(ori_ests)[idxs[,1]],
                              id2=colnames(ori_ests)[idxs[,2]])
estimate_pnas <- estimate_pnas %>%
  mutate(grp = paste(pmax(id1, id2), pmin(id1, id2), sep = "_"))

# p-value based approach
p_val_pnas <- readRDS("revision/estimate_cor_PNAS_same_gene/oli/norm_p.rds")
pnas_est <- readRDS("revision/estimate_cor_PNAS_same_gene/oli/est_cor.rds")
p_val_pnas$cscore_p <- pnas_est$cscore_p
pnas_p_adj <- as.data.frame(apply(p_val_pnas, 2, function(x){p.adjust(x, method = "BH")}))

thresh <- 0.05
p_cor_pnas <- pnas_p_adj<thresh
# p_cor_pnas <- p_val_pnas<thresh
colSums(p_cor_pnas)
p_cor_exp_pnas <- cbind(p_cor_pnas, estimate_pnas)

res_reproduce=reproduce_overlap(p_cor_exp, p_cor_exp_pnas, colnames(p_cor))

dumbbell_data <- res_reproduce %>%
  mutate(method_label = paste0(Method, " \n(Obs = ", obs, ")")) %>%
  mutate(method_label = fct_reorder(method_label, z_exp))

# Generate the dumbbell plot
dumbbell_reproduce <- ggplot(dumbbell_data, aes(y = method_label)) +
  # Add points for z_exp and z
  geom_point(aes(x = z_exp, color = "Adjust exp"), size = 4) +
  geom_point(aes(x = z, color = "No adjustment"), size = 4) +
  # Add a zero line for reference
  geom_vline(xintercept = qnorm(0.95), linetype = "dashed", color = "grey50") +
  scale_color_manual(name = "Score Type", values = c("Adjust exp" = "#0072B2", "No adjustment" = "#D55E00")) +
  labs(
    title = "Reproducibility",
    x = "Z-Score",
    y = "Method"
  ) +
  theme_minimal()

format = theme(text = element_text(size = 12))
pdf('/gpfs/gibbs/pi/zhao/xs282/validation/revision/modify_plot/adj_biological_bias_ROMSAP_PNAS.pdf', width = 13, height = 4, onefile = T)
ggarrange(dumbbell_string+format, dumbbell_reac+format, dumbbell_reproduce+format, ncol=3, nrow=1, common.legend = T, legend="right")
dev.off()
