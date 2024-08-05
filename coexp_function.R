
sct_cor <- function(sc_obj, sc.sel, sel.gene){
  dimension <- dim(sc.sel)
  sc.sct <- SCTransform(sc_obj,method = "glmGamPoi",ncells=dimension[2],
                        residual.features = sel.gene,return.only.var.genes = F)

  SCT_scaled_res <- sc.sct[["SCT"]]@scale.data
  SCT_scaled_res <- SCT_scaled_res[match(sel.gene, rownames(SCT_scaled_res)),]
  all(rownames(SCT_scaled_res) == sel.gene)

  sct.data <- t(as.matrix(SCT_scaled_res))
  sct_prn <- cor(sct.data, method = "pearson")
  return(sct_prn)
}

ana_prn <- function(count_mat, sel.gene, seq_depth, theta=100){
  sel_count <- count_mat[sel.gene,] %>% as.matrix %>% t
  mu_fit <- outer(seq_depth, colSums(sel_count)/sum(seq_depth))
  resi <- (sel_count-mu_fit)/sqrt(mu_fit+mu_fit^2/theta)

  ncell <- ncol(count_mat)
  resi[resi>sqrt(ncell)] <- sqrt(ncell)
  resi[resi<-sqrt(ncell)] <- -sqrt(ncell)
  cor_est <- cor(resi, method = "pearson")
  return(cor_est)
}

noise_fun <- function(sc_obj, sc.sel, sel.gene, seed, path2){
  dimension <- dim(sc.sel)

  sc.sct <- SCTransform(sc_obj,method = "glmGamPoi",ncells=dimension[2],
                        residual.features = sel.gene,return.only.var.genes = F)

  SCT_scaled_res <- sc.sct[["SCT"]]@scale.data
  SCT_scaled_res <- SCT_scaled_res[match(sel.gene, rownames(SCT_scaled_res)),]

  sc.sct[["my_cluster"]] <- 1
  meta.data <- sc.sct@meta.data
  sct.data <- t(as.matrix(SCT_scaled_res))

  ########################### NoiseRegularization ################################
  # https://github.com/RuoyuZhang/NoiseRegularization/blob/master/notebook/github.liver.example.ipynb
  source("/gpfs/gibbs/pi/zhao/xs282/coexp-sc/simulation/NoiseRegularizationCode.R") # from https://github.com/RuoyuZhang/NoiseRegularization/blob/master/code/code.R

  # Calculate quantile
  # if the min expression of gene<0, then add abs(min(expression)) to calculate the quantile
  sct.quantile.mat <- calculate_quantile_matrix.bygene(m=sct.data,gene.names = NULL,expr.cell.list = NULL,
                                                       quantile.list = c(seq(0.01,0.15,0.01),seq(0.2,0.7,0.1),seq(0.71,1,0.01)))

  # Begin to perform noise regularization
  all.genes <- colnames(sct.data)
  all.cells <- rownames(sct.data)
  # column name in meta data frame to indicate cell clusters, here we use the numbered cluster names
  cluster_column<- 'my_cluster'

  # spearman
  # specify the outdir to save gene gene correlation results
  out.dir.sct.spr <- paste0(path2,"noise/")
  dir.create(out.dir.sct.spr, recursive = T, showWarnings = F)

  # for each cluster, adding different levels of noise (indicated by q.list)
  q.list <- c('1%')

  # begin to run gene gene correlation in cluster ~4min
  gc()
  flush.console()
  # if the 1% quantile is less than max.low=0.1, then unif(0,0.1)
  frac.0.15 <- noise.regularization.cmp(m=sct.data, gene.list=all.genes, max.low=0.1,
                                        quantile.mat = sct.quantile.mat,
                                        quantile.list = q.list,title='sct.bygene',
                                        out.dir=paste0(out.dir.sct.spr,'/cluster',1,'/'),
                                        ncore = 1,nblocks=10,seed = seed)

  cor.dir <- paste0(path2,"noise")
  cor.df.1p <- readRDS(paste0(cor.dir,'/cluster',1, '/', 'sct','.bygene.','1p','.cor.df.rds'))

  # convert result to a correlation matrix
  # create a data contains all values of the correlation matrix
  cor.df.1p.den1 <- separate(cor.df.1p,pair_name, c("gene1", "gene2"), sep = "_")
  cor.df.1p.den1 <- cor.df.1p.den1[c('gene1', 'gene2', 'spearman')]
  cor.df.1p.den2 <- cor.df.1p.den1
  names(cor.df.1p.den2) <- c('gene2', 'gene1', 'spearman')
  cor.df.1p.den2 <- cor.df.1p.den2[c('gene1', 'gene2', 'spearman')]
  cor.df.1p.den3 <- data.frame(gene1=unique(cor.df.1p.den1$gene1),
                               gene2=unique(cor.df.1p.den1$gene1),
                               spearman=rep(1,length(unique(cor.df.1p.den1$gene1))))
  cor.df.1p.den <- rbind(cor.df.1p.den1,cor.df.1p.den2,cor.df.1p.den3)
  noise_cor <- spread(cor.df.1p.den, key = gene2, value = spearman)
  noise_cor <- as.data.frame(noise_cor)
  rownames(noise_cor) <- noise_cor[,1]
  noise_cor <- noise_cor[,-1]
  return(noise_cor)
}

