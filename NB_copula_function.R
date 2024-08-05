

copula_fun <- function(cor_mat, ncor_gene, ncell){
  tryCatch(
    {
      cor_mat_up <- chol(cor_mat)
      copula <- matrix(rnorm(ncor_gene*ncell),nrow = ncor_gene)
      copula <- pnorm(t(cor_mat_up)%*%copula)
      return(copula)

    }, error = function(e) {
      print(e)
      print("use eigen decomposition")
      mvnrv <- mvtnorm::rmvnorm(ncell, mean = rep(0, dim(cor_mat)[1]), sigma = cor_mat,
                                checkSymmetry = FALSE, method = "eigen")
      print("prnom")
      copula <- pnorm(t(mvnrv))
      return(copula)
    }
  )
}


NB_copula <- function(mu, gene_name, seq_depth, cell_name,
                      alpha, cor_mat=NULL, ind=T, seed=10312023){
  set.seed(seed)
  ngene <- length(mu)
  ncell <- length(seq_depth)

  # get beta based on mu=alpha*beta
  beta <- mu/alpha

  ###### generate independent expression matrix based on gamma(alpha, beta) ######
  exp_mat <- matrix(rgamma(ngene*ncell, shape = rep(alpha,ncell),
                           scale = rep(beta,ncell)), nrow = ngene)
  colnames(exp_mat) <- cell_name
  rownames(exp_mat) <- gene_name

  ##################### generate correlated gene expression ######################
  if (!ind){
    ncor_gene <- nrow(cor_mat)
    cor_gene_name <- colnames(cor_mat)
    cor_alpha <- alpha[cor_gene_name]
    cor_beta <- beta[cor_gene_name]

    # copula
    copula <- copula_fun(cor_mat, ncor_gene, ncell)
    cor_exp_matrix <- matrix(NA, nrow=ncor_gene, ncol=ncell)
    rownames(cor_exp_matrix) <- cor_gene_name
    for (i in 1:ncor_gene){
      cor_exp_matrix[i,] <- qgamma(copula[i,], shape=cor_alpha[i], scale=cor_beta[i])
    }

    # replace the value for correlated gene in the independent matrix
    exp_mat[cor_gene_name,] <- cor_exp_matrix

  }

  ################ generate count matrix by draw from poisson ####################
  seq_depth_matrix <- matrix(seq_depth, nrow = ngene, ncol = ncell, byrow = T)
  pois_para <- exp_mat*seq_depth_matrix
  count_mat <- matrix(rpois(ngene*ncell, lambda = c(pois_para)), nrow=ngene)
  colnames(count_mat) <- cell_name
  rownames(count_mat) <- gene_name
  return(count_mat)
}
