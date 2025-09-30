# marginally fit NB for each gene

set.seed(1052023)
library(glmGamPoi)
library(ggplot2)
library(Seurat)
# data from: https://www.pnas.org/doi/10.1073/pnas.2008762117
sc_obj <- readRDS("/gpfs/gibbs/pi/zhao/xs282/lung_data/lung_nc_myeloid.rds")
table(sc_obj@meta.data$CellType_Category)
sc_obj <- subset(sc_obj, subset = Subclass_Cell_Identity=="Monocyte")


################################# select ####################################
count_sel <- as.matrix(sc_obj[["RNA"]]@counts)
print(dim(count_sel))
size_factors_sel <- colSums(count_sel)

marginal_fit_fn <- '/gpfs/gibbs/pi/zhao/xs282/validation/revision/estimate_cor_lung/lung_Monocyte_marginal_fit.rds'

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
marginal_fit_sel <- vanilla[vanilla$up=="lower",]
# kernel smooth
km_ex5 <- ksmooth(marginal_fit_sel$mu, marginal_fit_sel$alpha, kernel="normal", bandwidth = bw.SJ(marginal_fit_sel$mu)*5)
saveRDS(km_ex5, '/gpfs/gibbs/pi/zhao/xs282/validation/revision/estimate_cor_lung/lung_Monocyte_ks_fit_5.rds')
