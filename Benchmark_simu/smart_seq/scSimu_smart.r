library(Matrix)
library(dplyr)
library(matrixStats)
library(biomaRt)

# INPUT
output_path = "/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/simu/"

count_data = readMM("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/counts.read.txt")
cells = read.table("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/cells.read.new.txt")
genes = read.table("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/genes.read.txt")
meta = read.table("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/meta.txt", sep = "\t", header=T)
colnames(count_data) = cells$V1
rownames(count_data) = genes$V1
meta_sub = meta[meta$Method=="Smart-seq2" & meta$Experiment=="pbmc1",]
meta_sub

ct = "Cytotoxic T cell"
cell_name = meta_sub$NAME[meta_sub$CellType==ct]
count_sub = count_data[,cell_name]
count_sub = as.matrix(count_sub)

# filter genes
gene_idx = which((rowSums(count_sub)>=10) & (rowSums(count_sub>=1)>=10))
nrow(count_sub)-length(gene_idx)
ct_counts = count_sub[gene_idx,]

# determine if there is zero inflation
check_marginal_fit = function(gene){
    library(pscl)
    library(MASS)
    m <- mean(gene)
    v <- var(gene)
    if(m >= v){
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        res = tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue = as.numeric(1 - pchisq(chisq_val, 1))
          c(pvalue,0, ifelse(pvalue<0.05, "ZIP", "Poisson"))
        }, error = function(cond){c(1,0, "Poisson")})
    }else{
        res = tryCatch({
                mle_NB <- glm.nb(gene ~ 1)
                tryCatch({
                    mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
                    chisq_val <- 2 * (logLik(mle_ZINB) - logLik(mle_NB))
                    pvalue = as.numeric(1 - pchisq(chisq_val, 1))
                    c(0, pvalue, ifelse(pvalue<0.05, "ZINB", "NB"))
                }, error = function(cond){c(0, 1, "NB")})
             }, warning = function(w) {
                 tryCatch({
                    mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
                    chisq_val <- 2 * (logLik(mle_ZINB) - logLik(mle_NB))
                    pvalue = as.numeric(1 - pchisq(chisq_val, 1))
                    c(0, pvalue, ifelse(pvalue<0.05, "ZINB", "NB"))
                }, error = function(cond){c(0, 1, "NB")})
             })
    }
    return(res)
}

p_gene = list()
for (i in rownames(ct_counts)){
    gene = ct_counts[i,]
    p_gene[[i]] = check_marginal_fit(gene)
}
p = do.call(rbind, p_gene)
colnames(p) = c("Poisson", "NB","group")
p = as.data.frame(p)
p$gene = rownames(p)

ZINB_genes = p$gene[p$group %in% c("ZINB", "ZIP")]
other_genes = setdiff(p$gene, ZINB_genes)

# estimate the marginal parameter
## estimate non-zero inflated genes
gp_ex <- glm_gp(ct_counts[other_genes,], size_factors = colSums(ct_counts), verbose = T, overdispersion_shrinkage = T, do_cox_reid_adjustment = T)
vanilla <- data.frame(mu = exp(gp_ex$Beta[,1]), alpha = 1/gp_ex$overdispersions, gene=other_genes)

## estimate zero inflated genes
zinb_res = as.data.frame(matrix(nrow=length(ZINB_genes), ncol=4))
colnames(zinb_res) = c("gene", "zero_prop", "alpha", "mu")
rownames(zinb_res) = ZINB_genes
for (g in ZINB_genes){
    df <- data.frame(counts = ct_counts[g,], library_size = log(colSums(ct_counts)))
    fit_ZINB <- zeroinfl(counts ~ 1+offset(library_size)|1, dist = 'negbin', data=df)
    zinb_res[g,] = c(g, plogis(fit_ZINB$coefficients$zero), fit_ZINB$theta, exp(fit_ZINB$coefficients$count))
}
zinb_res[,2:4] = apply(zinb_res[,2:4], 2, as.numeric)
zinb_res$outlier = ifelse(log10(zinb_res$alpha)>(quantile(log10(zinb_res$alpha), 0.75)+2*IQR(log10(zinb_res$alpha))), 1, 0)

marginal_res = list(zinb=zinb_res, nb=vanilla)
saveRDS(marginal_res, paste0(output_path, "/marginal_res.rds"))

# fit the kernel regression
vanilla$up <- ifelse(log10(vanilla$alpha)>2.5, "upper","lower")
marginal_fit_PNAS_sel <- vanilla[vanilla$up=="lower",]
km5 <- ksmooth(log10(marginal_fit_PNAS_sel$mu), log10(marginal_fit_PNAS_sel$alpha),
               kernel="normal", bandwidth = bw.SJ(log10(marginal_fit_PNAS_sel$mu))*5)

mu <- vanilla$mu
gene_name <- vanilla$gene
names(mu) <- gene_name
log10mu = log10(mu)

fitted_trend <- data.frame(mu=km5$x, alpha=km5$y)
log10alpha <- rep(NA,nrow(vanilla))
names(log10alpha) <- vanilla$gene
for (i in 1:nrow(vanilla)){
  idx <- which.min(abs(log10mu[i]-fitted_trend$mu))
  log10alpha[i] <- fitted_trend$alpha[idx]
}
alpha <- 10^log10alpha


marginal_fit_zinb <- zinb_res[zinb_res$outlier==0,]
km5_zinb <- ksmooth(log10(marginal_fit_zinb$mu), log10(marginal_fit_zinb$alpha),
               kernel="normal", bandwidth = bw.SJ(log10(marginal_fit_zinb$mu))*5)

mu_zinb <- zinb_res$mu
gene_name_zinb <- zinb_res$gene
names(mu_zinb) <- gene_name_zinb
log10mu_zinb = log10(mu_zinb)

fitted_trend_zinb <- data.frame(mu=km5_zinb$x, alpha=km5_zinb$y)
log10alpha_zinb <- rep(NA,nrow(zinb_res))
names(log10alpha_zinb) <- zinb_res$gene
for (i in 1:nrow(zinb_res)){
  idx <- which.min(abs(log10mu_zinb[i]-fitted_trend_zinb$mu))
  log10alpha_zinb[i] <- fitted_trend_zinb$alpha[idx]
}
alpha_zinb <- 10^log10alpha_zinb

alpha_comb = c(alpha, alpha_zinb)
mu_comb = c(mu, mu_zinb)

# generate gene cor mat
genes = rownames(ct_counts)
ensg_ids <- sapply(strsplit(genes, "_"), function(x) x[1])

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_lengths <- getBM(attributes = c('ensembl_gene_id', 'transcript_length'),
                      filters = 'ensembl_gene_id',
                      values = ensg_ids,
                      mart = ensembl)

gene_info <- gene_lengths %>%
  group_by(ensembl_gene_id) %>%
  summarise(length = mean(transcript_length))

calculate_tpm <- function(counts, lengths) {
  # Ensure lengths are in kilobases
  lengths_kb <- lengths / 1000
  
  # 1. Normalize for gene length
  rpk <- counts / lengths_kb
  # 2. Normalize for sequencing depth
  per_million_scaling_factor <- colSums(rpk) / 1e6
  tpm <- sweep(rpk, 2, per_million_scaling_factor, "/")
  return(tpm)
}

# mean expression level under each cell type
sel_genes = sapply(strsplit(rownames(ct_counts), "_"), function(x) x[1])
matched_lengths <- gene_info$length[match(sel_genes, gene_info$ensembl_gene_id)]
# Calculate TPM, removing genes with no length information
valid_genes <- !is.na(matched_lengths)
tpm_matrix <- calculate_tpm(ct_counts[valid_genes, ], matched_lengths[valid_genes])
log_tpm_matrix <- log1p(tpm_matrix)
log_tpm_matrix = as.matrix(log_tpm_matrix)

mean_exp = rowMeans(log_tpm_matrix)
mean_exp = mean_exp[order(mean_exp, decreasing=T)]
smart_high_exp = names(mean_exp)[1:1000]
cor_mat = cor(t(log_tpm_matrix[smart_high_exp,]))

saveRDS(cor_mat, paste0(output_path,"/cor_mat.rds"))
                   
# simulate NB data
gene_name = names(mu_comb)
cell_name <- colnames(ct_counts)
seq_depth <- colSums(ct_counts)

source("/gpfs/gibbs/pi/zhao/xs282/validation/AFinal/NB_copula_function.R")
simu_nb <- NB_copula(mu_comb, gene_name, seq_depth, cell_name, alpha_comb,
                     cor_mat, ind=F, seed=11132023)
saveRDS(simu_nb, paste0(output_path,"/simu_NB.rds"))

simu_nb_ind <- NB_copula(mu_comb, gene_name, seq_depth, cell_name, alpha_comb, ind=T, seed=11132023)
saveRDS(simu_nb_ind, paste0(output_path,"/simu_NB_IND.rds"))

# generate zero inflation
set.seed(962025)
zero_inflated_matrix <- simu_nb
for (g in rownames(simu_nb)) {
    num_cells <- ncol(simu_nb)
    if (g %in% zinb_res$gene){
        rate = zinb_res$zero_prop[zinb_res$gene==g]
        zero_mask = rbinom(n = num_cells, size = 1, prob = rate)
    } else{
        zero_mask = rep(0, num_cells)
    }
    
    zero_inflated_matrix[g, ] <- zero_inflated_matrix[g, ] * (1 - zero_mask)
}
saveRDS(zero_inflated_matrix, paste0(output_path,"/simu_ZINB.rds"))

set.seed(962025)
zero_inflated_matrix <- simu_nb_ind
for (g in rownames(simu_nb_ind)) {
    num_cells <- ncol(simu_nb_ind)
    if (g %in% zinb_res$gene){
        rate = zinb_res$zero_prop[zinb_res$gene==g]
        zero_mask = rbinom(n = num_cells, size = 1, prob = rate)
    } else{
        zero_mask = rep(0, num_cells)
    }
    
    zero_inflated_matrix[g, ] <- zero_inflated_matrix[g, ] * (1 - zero_mask)
}
saveRDS(zero_inflated_matrix, paste0(output_path,"/simu_ZINB_IND.rds"))
