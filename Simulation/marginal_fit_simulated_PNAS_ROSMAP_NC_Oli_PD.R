set.seed(1052023)
library(glmGamPoi)
library(ggplot2)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")

################################# PNAS ####################################
count_sel <- readRDS("mean_cor/semi_PD_sparse/simu/PNAS_NC_Oli_sct_cor_NB_simu1000_abs_thresh.rds")
print(dim(count_sel))
size_factors_sel <- colSums(count_sel)

if(!file.exists('/gpfs/gibbs/pi/zhao/xs282/validation/marginal_fit')){
  dir.create('/gpfs/gibbs/pi/zhao/xs282/validation/marginal_fit')
}

marginal_fit_fn <- 'marginal_fit/PNAS_NC_Oli_PD_sparse_marginal_fit_simulated.rds'

if(!file.exists(marginal_fit_fn)){
  gp_ex <- glm_gp(count_sel, size_factors = size_factors_sel, verbose = T, overdispersion_shrinkage = T, do_cox_reid_adjustment = T)
  vanilla <- data.frame(mu = log10(exp(gp_ex$Beta[,1])),
                        alpha = -log10(gp_ex$overdispersions),
                        deviances = gp_ex$deviances,
                        gene=rownames(count_sel))
  saveRDS(vanilla, marginal_fit_fn)
  print("succeed")
}

vanilla$up <- ifelse(vanilla$alpha>2, "upper","lower")

# plot alpha vs mu
plot1 <- ggplot(vanilla, aes(x=mu, y=alpha))+
  geom_point()+labs(title = "log10(alpha) v.s. log10(mu)")+
  xlab("log10(mu)")+ylab("log10(alpha)")
plot2 <- ggplot(vanilla, aes(x=mu, y=alpha, colour=up))+
  geom_point()+labs(title = "log10(alpha) v.s. log10(mu)")+
  xlab("log10(mu)")+ylab("log10(alpha)")
plot1
plot2

# fit line (only use the lower cluster)
marginal_fit_PNAS_sel <- vanilla[vanilla$up=="lower",]
# kernel smooth
km_ex5_PNAS <- ksmooth(marginal_fit_PNAS_sel$mu, marginal_fit_PNAS_sel$alpha, kernel="normal", bandwidth = bw.SJ(marginal_fit_PNAS_sel$mu)*5)
saveRDS(km_ex5_PNAS, 'marginal_fit/PNAS_NC_Oli_PD_sparse_ks_fit_5_simulated.rds')

################################# ROSMAP ####################################
rm(ls())
count_sel <- readRDS("mean_cor/semi_PD_sparse/simu/ROSMAP_NC_Oli_sct_cor_NB_simu1000_abs_thresh.rds")
print(dim(count_sel))
size_factors_sel <- colSums(count_sel)


marginal_fit_fn <- 'marginal_fit/ROSMAP_NC_Oli_PD_sparse_marginal_fit_simulated.rds'

if(!file.exists(marginal_fit_fn)){
  gp_ex <- glm_gp(count_sel, size_factors = size_factors_sel, verbose = T, overdispersion_shrinkage = T, do_cox_reid_adjustment = T)
  vanilla <- data.frame(mu = log10(exp(gp_ex$Beta[,1])),
                        alpha = -log10(gp_ex$overdispersions),
                        deviances = gp_ex$deviances,
                        gene=rownames(count_sel))
  saveRDS(vanilla, marginal_fit_fn)
  print("succeed")
}

vanilla$up <- ifelse(vanilla$alpha>2, "upper","lower")

# plot alpha vs mu
plot1 <- ggplot(vanilla, aes(x=mu, y=alpha))+
  geom_point()+labs(title = "log10(alpha) v.s. log10(mu)")+
  xlab("log10(mu)")+ylab("log10(alpha)")
plot2 <- ggplot(vanilla, aes(x=mu, y=alpha, colour=up))+
  geom_point()+labs(title = "log10(alpha) v.s. log10(mu)")+
  xlab("log10(mu)")+ylab("log10(alpha)")
plot1
plot2

# fit line (only use the lower cluster)
marginal_fit_ROSMAP_sel <- vanilla[vanilla$up=="lower",]
# kernel smooth
km_ex5_ROSMAP <- ksmooth(marginal_fit_ROSMAP_sel$mu, marginal_fit_ROSMAP_sel$alpha,
                         kernel="normal", bandwidth = bw.SJ(marginal_fit_ROSMAP_sel$mu)*5)
saveRDS(km_ex5_ROSMAP, 'marginal_fit/ROSMAP_NC_Oli_PD_sparse_ks_fit_5_simulated.rds')
