
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(glmGamPoi)
library(optparse)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")


# save count matrix
sc_obj <- readRDS("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/MIT_ROSMAP_Multiomics/Gene_Expression/snRNAseq/snRNAseq-10x/processed/Excitatory_neurons_set1.rds")
clic = read.csv("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/ROSMAP/Metadata/ROSMAP_clinical.csv")
sc_obj_meta = sc_obj@meta.data
sc_obj_meta = left_join(sc_obj_meta, clic)
sc_obj$braaksc <- sc_obj_meta$braaksc
count_ROSMAP_ex <- as.matrix(sc_obj[["RNA"]]@counts[, which(sc_obj$braaksc==0)])
saveRDS(count_ROSMAP_ex, "marginal_fit/ROSMAP_NC_Ex_count.rds")

# marginal fit
size_factors_sel <- colSums(count_ROSMAP_ex)
marginal_fit_fn <- 'marginal_fit/ROSMAP_NC_Ex_marginal_fit.rds'

if(!file.exists(marginal_fit_fn)){
  gp_ex <- glm_gp(count_ROSMAP_ex, size_factors = size_factors_sel, verbose = T, overdispersion_shrinkage = T, do_cox_reid_adjustment = T)
  vanilla <- data.frame(mu = log10(exp(gp_ex$Beta[,1])),
                        alpha = -log10(gp_ex$overdispersions),
                        deviances = gp_ex$deviances,
                        gene=rownames(count_ROSMAP_ex))
  saveRDS(vanilla, marginal_fit_fn)
  print("succeed")
}

vanilla$up <- ifelse(vanilla$alpha>2.5, "upper","lower")
# plot alpha vs mu
plot1 <- ggplot(vanilla, aes(x=mu, y=alpha))+
  geom_point()+labs(title = "log10(alpha) v.s. log10(mu)")+
  xlab("log10(mu)")+ylab("log10(alpha)")
plot2 <- ggplot(vanilla, aes(x=mu, y=alpha, colour=up))+
  geom_point()+labs(title = "log10(alpha) v.s. log10(mu)")+
  xlab("log10(mu)")+ylab("log10(alpha)")
plot1+plot2

# fit line (only use the lower cluster)
vanilla_sel <- vanilla[vanilla$up=="lower",]
# kernel smooth
km_ex5_ROSMAP <- ksmooth(vanilla_sel$mu, vanilla_sel$alpha, kernel="normal", bandwidth = bw.SJ(vanilla_sel$mu)*5)
saveRDS(km_ex5_ROSMAP, 'marginal_fit/ROSMAP_NC_Ex_ks_fit_5.rds')

