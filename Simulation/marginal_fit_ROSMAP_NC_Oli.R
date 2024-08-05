set.seed(1052023)
library(glmGamPoi)
library(ggplot2)
library(dplyr)
# data from ROSMAP
setwd("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/MIT_ROSMAP_Multiomics/")
sc_obj <- readRDS("Gene_Expression/snRNAseq/snRNAseq-10x/processed/Oligodendrocytes.rds")

clic = read.csv("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/ROSMAP/Metadata/ROSMAP_clinical.csv")
sc_obj_meta = sc_obj@meta.data
sc_obj_meta = left_join(sc_obj_meta, clic)
length(unique(sc_obj_meta[sc_obj_meta$braaksc==0,]$projid))

sc_obj$braaksc <- sc_obj_meta$braaksc

################################# select Oli ####################################
count_sel <- as.matrix(sc_obj[["RNA"]]@counts[, which(sc_obj$braaksc==0)])
print(dim(count_sel))
size_factors_sel <- colSums(count_sel)

if(!file.exists('/gpfs/gibbs/pi/zhao/xs282/validation/marginal_fit')){
  dir.create('/gpfs/gibbs/pi/zhao/xs282/validation/marginal_fit')
}

marginal_fit_fn <- '/gpfs/gibbs/pi/zhao/xs282/validation/marginal_fit/ROSMAP_NC_Oli_marginal_fit.rds'

if(!file.exists(marginal_fit_fn)){
  gp_ex <- glm_gp(count_sel, size_factors = size_factors_sel, verbose = T, overdispersion_shrinkage = T, do_cox_reid_adjustment = T)
  vanilla <- data.frame(mu = log10(exp(gp_ex$Beta[,1])),
                        alpha = -log10(gp_ex$overdispersions),
                        deviances = gp_ex$deviances,
                        gene=rownames(count_sel))
  saveRDS(vanilla, marginal_fit_fn)
  print("succeed")
}

