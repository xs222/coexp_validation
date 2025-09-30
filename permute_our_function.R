
# count_mat: using all cells and all genes(cells in row, genes in column)

permute_our <- function(counts, seq_depth, seed){
  n <- nrow(counts)
  p <- ncol(counts)
  counts_ratio <- counts/seq_depth
  
  counts_permute = MatrixPermutation(counts_ratio, n, p, seed = seed)*seq_depth
  counts_permute = rpois(lambda = c(counts_permute), n = n*p) %>% matrix(nrow = n, ncol = p)
  rownames(counts_permute) <- rownames(counts)
  colnames(counts_permute) <- colnames(counts)
  return(counts_permute)
}

# Permute a matrix by rows
MatrixPermutation = function(mat, n, p, seed){
  ind = matrix(1:(n*p), nrow = n, ncol = p)
  set.seed(seed)
  ind_permute = apply(ind, 2, function(x) sample(x, size = n))
  return(matrix(mat[ind_permute], nrow = n, ncol = p))
}
