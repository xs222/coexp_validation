library(dplyr)
library(ggplot2)
library(matrixcalc)
source("/gpfs/gibbs/pi/zhao/xs282/validation/AFinal/NB_copula_function.R")

pre_process_cor_est <- function(cor_est, type = 'cor ests'){
  apply(cor_est, 1, function(ce) sum(is.na(ce))) %>% table %>% print
  
  na_gene_inds <- which(is.na(diag(cor_est)))
  if(anyNA(diag(cor_est))){
    print(sprintf('%i genes have negative variance estimates', length(na_gene_inds)))
    cor_est[na_gene_inds,] <- 0
    cor_est[,na_gene_inds] <- 0    
  }
  
  print(sprintf('#cor est greater than 1: %i', sum(cor_est > 1)))
  print(sprintf('#cor est smaller than -1: %i', sum(cor_est < -1)))
  cor_est[cor_est > 1] <- 1
  cor_est[cor_est < -1] <- -1
  diag(cor_est) <- 1
  return(list(est = cor_est, gene = colnames(cor_est)))
}

est_mat_ROSMAP = readRDS("/gpfs/gibbs/pi/zhao/xs282/validation/marginal_fit/ROSMAP_NC_Oli_cscore_cor1000.rds")
est_cor = est_mat_ROSMAP$est
est_cor = pre_process_cor_est(est_cor)$est

# BH
est_pvalue <- est_mat_ROSMAP$p_value
est_adj_pvalue <- est_pvalue
est_adj_pvalue[est_adj_pvalue!=0] <- 0
est_adj_pvalue[upper.tri(est_pvalue, diag=F)] <- p.adjust(est_pvalue[upper.tri(est_pvalue, diag=F)], method="BH")
est_adj_pvalue <- est_adj_pvalue+t(est_adj_pvalue)

est_cor[est_adj_pvalue>0.05 | abs(est_mat_ROSMAP$est)<0.1] <- 0
est_cor <- (est_cor+t(est_cor))/2
is.positive.semi.definite(est_cor)

# CS-CORE approach
cor_mat_cscore <- est_cor
k <- 0
while(!is.positive.definite(cor_mat_cscore, tol = 1e-05)){
  k <- k+1
  print(sprintf("Matrix is not p.d., rescale iteration %i",k))
  cor_mat_cscore <- est_cor + diag(x = 0.1*k, nrow = nrow(est_cor), ncol = nrow(est_cor))
  cor_mat_cscore <- cor_mat_cscore/cor_mat_cscore[1,1]
}
set.seed(9222025)
mvn_cscore = copula_fun(cor_mat_cscore, ncor_gene=nrow(est_cor), ncell=1000)
mvn_cor_cscore = cor(t(mvn_cscore))

# scSimu approach
set.seed(9222025)
mvn_scsimu = copula_fun(est_cor, ncor_gene=nrow(est_cor), ncell=1000)
mvn_cor_scsimu = cor(t(mvn_scsimu))

# Compare two approaches
comp_res1 = data.frame(true = est_cor[upper.tri(est_cor)],
                     est = mvn_cor_cscore[upper.tri(mvn_cor_cscore)])
comp_res1$group = "CS-CORE"
comp_res2 = data.frame(true = est_cor[upper.tri(est_cor)],
                     est = mvn_cor_scsimu[upper.tri(mvn_cor_scsimu)])
comp_res2$group = "scSimu"

comp_res = rbind(comp_res1, comp_res2)        

png('/gpfs/gibbs/pi/zhao/xs282/validation/revision/modify_plot/scsimu_vs_cscore.png', width = 500, height = 300)
ggplot(comp_res, aes(x=true, y=est, color=group))+
geom_point(alpha=0.2)+
geom_abline(slope=1, intercept=0, color="black", linetype="dashed")+
theme_bw()+labs(x="True correlation", y="Estimated correlation", color="Simumation methods")+
theme(text = element_text(size = 14))
dev.off()