library(Matrix)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(biomaRt)
library(RANN)
library(scater)
library(scran)
library(SingleCellExperiment)
library(Seurat)

eval_fun <- function(ct, gene_info){
    eval_metrics <- list()
    sel_genes = sapply(strsplit(rownames(ct), "_"), function(x) x[1])
    matched_lengths <- gene_info$length[match(sel_genes, gene_info$ensembl_gene_id)]
    # Calculate TPM, removing genes with no length information
    valid_genes <- !is.na(matched_lengths)
    tpm_matrix <- calculate_tpm(ct[valid_genes, ], matched_lengths[valid_genes])
    log_tpm_matrix <- log1p(tpm_matrix)
    normalized_data = as.matrix(log_tpm_matrix)
    
    # average of logCPM
    eval_metrics[["gene_mean"]] <- rowMeans(normalized_data)
    
    # variance of logCPM
    eval_metrics[["gene_var"]] <- rowVars(normalized_data)
    
    # coefficient of variation
    eval_metrics[["gene_cv"]] <- sqrt(eval_metrics[["gene_var"]])/eval_metrics[["gene_mean"]]
    
    # fraction of genes with zero counts
    eval_metrics[["gene_frq_zero"]] <- rowMeans(ct==0)
    
    # fraction of cells with zero counts
    eval_metrics[["cell_frq_zero"]] <- colMeans(ct==0)
    
    # library size (total counts)
    eval_metrics[["lib_size"]] <- colSums(ct)
    
    # cell-to-cell correlation
    print("cell_cor")
    if (max(sample_idx)>ncol(normalized_data)){
        new_sample_idx <- sample_idx[sample_idx<=ncol(normalized_data)]
        cell_cor <- cor(as.matrix(normalized_data)[,new_sample_idx], method = "spearman", use = "pairwise.complete.obs")
    } else{
        cell_cor <- cor(as.matrix(normalized_data)[,sample_idx], method = "spearman", use = "pairwise.complete.obs")
    }
    eval_metrics[["cell_cor"]] <- cell_cor[upper.tri(cell_cor)]
    
    # gene-gene correlation
    print("gene_cor")
    genes_selected <- names(sort(eval_metrics[["gene_mean"]], decreasing = T))[1:1000]
    gene_cor = cor(t(normalized_data[genes_selected, ]))
    eval_metrics[["gene_cor"]] <- gene_cor[upper.tri(gene_cor)]

    return(eval_metrics)
}

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

count_data = readMM("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/counts.read.txt")
cells = read.table("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/cells.read.new.txt")
genes = read.table("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/genes.read.txt")
meta = read.table("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/meta.txt", sep = "\t", header=T)
colnames(count_data) = cells$V1
rownames(count_data) = genes$V1
meta_sub = meta[meta$Method=="Smart-seq2" & meta$Experiment=="pbmc1",]

ct = "Cytotoxic T cell"
cell_name = meta_sub$NAME[meta_sub$CellType==ct]
count_sub = count_data[,cell_name]
count_sub = as.matrix(count_sub)

# filter genes
gene_idx = which((rowSums(count_sub)>=10) & (rowSums(count_sub>=1)>=10))
nrow(count_sub)-length(gene_idx)
count_sub = count_sub[gene_idx,]
path <- "/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/simu/cyto_t/"

# gene info
genes = rownames(count_sub)
ensg_ids <- sapply(strsplit(genes, "_"), function(x) x[1])

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_lengths <- getBM(attributes = c('ensembl_gene_id', 'transcript_length'),
                      filters = 'ensembl_gene_id', values = ensg_ids, mart = ensembl)

gene_info <- gene_lengths %>%
  group_by(ensembl_gene_id) %>%
  summarise(length = mean(transcript_length))
                   
                       
# ori ------------------------------------------------------------
set.seed(1272024)
sample_idx <- sample(1:ncol(count_sub), min(500, ncol(count_sub)))
eval_result <- list()
eval_result[["ori"]] <- eval_fun(count_sub, gene_info)


# scDesign3 --------------------------------------------------------------------
print("scDesign3 --------------------------------------------------------------------")               
scdesign3_simu = readRDS(paste0(path,"simu_scDesign3.rds"))   
scdesign3_res = scdesign3_simu$new_count
scdesign3_res = scdesign3_res[rownames(count_sub), colnames(count_sub)]
eval_result[["scDesign3"]] <- eval_fun(scdesign3_res, gene_info)

# scSimu --------------------------------------------------------------------
print("scSimu --------------------------------------------------------------------")                  
scsimu_simu = readRDS(paste0(path,"simu_scSimu.rds"))
scsimu_simu = scsimu_simu[rownames(count_sub), colnames(count_sub)]
eval_result[["scSimu"]] <- eval_fun(scsimu_simu, gene_info)

saveRDS(eval_result, "/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/simu/cyto_t/eval_res.rds")

# KS and KDE-------------------------------------------------------------------------------------------------------------------
library(ks)
simu_method <- c("scDesign3", "scSimu")
eval_metrics <- c("gene_mean", "gene_var", "gene_cv", "gene_frq_zero", "cell_frq_zero",
                  "lib_size", "cell_cor", "gene_cor")

kde_stat <- matrix(nrow=length(simu_method), ncol=length(eval_metrics))
rownames(kde_stat) <- simu_method
colnames(kde_stat) <- c(eval_metrics)

ks_stat <- matrix(nrow=length(simu_method), ncol=length(eval_metrics))
rownames(ks_stat) <- simu_method
colnames(ks_stat) <- c(eval_metrics)


ori_dist <- eval_result[["ori"]]
for (i in simu_method){
  print(paste(i, "----------------------------------------------"))
  simu_dist <- eval_result[[i]]
  for (j in eval_metrics){
    print(j)
    x <- ori_dist[[j]]
    x <- ifelse(is.na(x), 0, x)
    y <- simu_dist[[j]]
    y <- ifelse(is.na(y), 0, y)

    # ks test
    z <- ks.test(x, y)
    ks_stat[i,j] <- z$statistic

    # kde test
    set.seed(1272024)
    if (length(x)>10000 | length(y)>10000){
      if (length(x)==length(y)){
        idx <- sample(1:length(x), 10000)
        x_new <- x[idx]
        y_new <- y[idx]
      } else{
        if (length(x)>10000){
          x_new <- x[sample(1:length(x), 10000)]
        } else{
          x_new <- x
        }

        if (length(y)>10000){
          y_new <- y[sample(1:length(y), 10000)]
        } else{
          y_new <- y
        }
      }
      if (length(unique(x_new))==1){
        x_new[1] <- x_new[1]+0.0001
      }

      if (length(unique(y_new))==1){
        y_new[1] <- y_new[1]+0.0001
      }
    } else{
      x_new <- x
      y_new <- y

      if (length(unique(x_new))==1){
        x_new[1] <- x_new[1]+0.0001
      }

      if (length(unique(y_new))==1){
        y_new[1] <- y_new[1]+0.0001
      }
    }
    kde_stat[i,j] <- kde.test(x_new, y_new)$zstat
  }

}

ks_stat <- as.data.frame(ks_stat)
kde_stat <- as.data.frame(kde_stat)
saveRDS(list(ks=ks_stat, kde=kde_stat), "/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/simu/cyto_t/ks_kde.rds")
