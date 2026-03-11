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
    res$z_exp = res$adjust_exp/sqrt(res$random_var_adjust_exp)
    res$z = res$adjust/sqrt(res$random_var)
    res$md = rownames(res)
    res$Method <- recode(res$md,sct="sctransform", prn="Pearson", spr="Spearman",
                       propr="propr",ana_prn="Analytic PR",cscore_p="CS-CORE",
                       noise="Noise \nRegularization", cscore_est="CS-CORE Empirical")
    return(res)
}



# zscore -----------------------------------------------------------------
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
# STRING
hs_filter <- readRDS("STRING/hs_filter_10_4_2023.rds")
hs_filter <- hs_filter %>%
  dplyr::mutate(grp = paste(pmax(protein1, protein2), pmin(protein1, protein2), sep = "_"))
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

intervals = c(0.01,0.05,0.1,0.5)
string_res = list()
reactome_res = list()
for (j in 1:length(intervals)){
  interval = intervals[j]
  breaks <- seq(from = min(p_cor_exp$log10mean_mu_lung), to = max(p_cor_exp$log10mean_mu_lung), by = interval)
  p_cor_exp$exp_gp <- cut(x = p_cor_exp$log10mean_mu_lung, breaks = breaks)

  res_string_lung = rand_overlap(p_cor_exp, string, colnames(p_cor))
  res_reactome_lung = rand_overlap(p_cor_exp, reactome, colnames(p_cor))

  dumbbell_string <- res_string_lung %>%
    mutate(method_label = paste0(Method, " \n(Obs = ", obs, ")")) %>%
    mutate(method_label = fct_reorder(method_label, z_exp))
  dumbbell_string$thresh = interval
  string_res[[j]] = dumbbell_string

  dumbbell_reactome <- res_reactome_lung %>%
    mutate(method_label = paste0(Method, " \n(Obs = ", obs, ")")) %>%
    mutate(method_label = fct_reorder(method_label, z_exp))
  dumbbell_reactome$thresh = interval
  reactome_res[[j]] = dumbbell_reactome
}

string_data = do.call(rbind, string_res)
reactome_data = do.call(rbind, reactome_res)

dumbbell_string <- ggplot(string_data, aes(y = method_label)) +
  # Add points for z_exp and z
  geom_point(aes(x = z_exp, color = "Adjust exp", shape=as.factor(thresh)), size = 4) +
  geom_point(aes(x = z, color = "No adjustment"), size = 4) +
  # Add a zero line for reference
  geom_vline(xintercept = qnorm(0.95), linetype = "dashed", color = "grey50") +
  scale_color_manual(name = "Score Type", values = c("Adjust exp" = "#0072B2", "No adjustment" = "#D55E00")) +
  labs(
    title = "STRING",
    x = "Z-Score",
    y = "Method",
    shape="Interval Size"
  ) +
  theme_minimal()



dumbbell_reac <- ggplot(reactome_data, aes(y = method_label)) +
  # Add points for z_exp and z
  geom_point(aes(x = z_exp, color = "Adjust exp", shape=as.factor(thresh)), size = 4) +
  geom_point(aes(x = z, color = "No adjustment"), size = 4) +
  # Add a zero line for reference
  geom_vline(xintercept = qnorm(0.95), linetype = "dashed", color = "grey50") +
  scale_color_manual(name = "Score Type", values = c("Adjust exp" = "#0072B2", "No adjustment" = "#D55E00")) +
  labs(
    title = "Reactome",
    x = "Z-Score",
    y = "Method",
    shape="Interval Size"
  ) +
  theme_minimal()

plot_rep_ori1 <- ggarrange(dumbbell_string,dumbbell_reac, ncol=2,nrow=1,
                           widths = c(1, 1),common.legend = T, legend="right",
                           labels = c("A", "B"))
plot_rep_ori1

#
# hist(p_cor_exp$log10mean_mu_lung[p_cor_exp$biological_net==1])
# hist(p_cor_exp$log10mean_mu_lung[p_cor_exp$prn==1])

pdf('BIB_R1/modify_plot/robustness_stratification.pdf', width = 8, height = 3.5, onefile = T)
plot_rep_ori1
dev.off()

