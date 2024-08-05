# library(WGCNA)
# library(org.Hs.eg.db)
# library(clusterProfiler)
# library(enrichplot)
# library(gprofiler2)


prn_p <- function(x, ncell, gene_name, BH=T){
  # calculate the p value for pearson correlation
  sct_t <- x*sqrt((ncell-2)/(1-x^2))
  sct_p <- 2 * (1-pt(q = abs(sct_t), df = ncell - 2))
  if (BH){
    sct_p <- MatrixBH(sct_p)
  }
  sct_p <- sct_p[gene_name,gene_name]
  return(sct_p)
}

post_process_est <- function (est) {
  p <- nrow(est)
  neg_gene_inds <- which(sapply(diag(est), function(x) is.infinite(x) |
                                  is.na(x)))
  if (length(neg_gene_inds) > 0) {
    print(sprintf("%i among %i genes have negative variance estimates. Their co-expressions with other genes were set to 0.",
                  length(neg_gene_inds), p))
  }
  est[neg_gene_inds, ] <- 0
  est[, neg_gene_inds] <- 0
  diag(est) <- 1
  print(sprintf("%.4f%% co-expression estimates were greater than 1 and were set to 1.",
                mean(est[upper.tri(est)] > 1, na.rm = T) * 100))
  print(sprintf("%.4f%% co-expression estimates were smaller than -1 and were set to -1.",
                mean(est[upper.tri(est)] < -1, na.rm = T) * 100))
  est[est > 1] <- 1
  est[est < -1] <- -1
  return(est)
}


MatrixBH = function(p_matrix){
  p_Matrix_BH = p_matrix - p_matrix
  p_Matrix_BH[upper.tri(p_Matrix_BH)] = p.adjust(p_matrix[upper.tri(p_matrix)], method = "BH")
  return(p_Matrix_BH + t(p_Matrix_BH))
}

extract_upp <- function(my_matrix){
  tri <- upper.tri(my_matrix, diag=F)
  idxs <- which(tri, arr.ind = T)
  cor_df <- data.frame(id1 = rownames(my_matrix)[idxs[,1]],
                       id2 = rownames(my_matrix)[idxs[,2]],
                       correlation = my_matrix[tri])
  return(cor_df)
}
