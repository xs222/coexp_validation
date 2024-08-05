set.seed(1052023)
library(glmGamPoi)
library(ggplot2)
# data from: https://www.pnas.org/doi/10.1073/pnas.2008762117
sc_obj <- readRDS("/gpfs/gibbs/pi/zhao/cs2629/AD_Nancy/seurat_obj_cell_type_labelled.rds")
table(sc_obj@meta.data$cell_type)

sc_obj$disease = sapply(sc_obj$orig.ident, function(or) {substr(strsplit(or, '_')[[1]][2], 1, 2)})
table(sc_obj$disease)

length(unique(sc_obj$orig.ident[sc_obj$disease=="NC"]))

################################# select Oli ####################################
count_sel <- as.matrix(sc_obj[["RNA"]]@counts[, which(sc_obj$cell_type=="Oli" & sc_obj$disease=="NC")])
print(dim(count_sel))
size_factors_sel <- colSums(count_sel)

if(!file.exists('/gpfs/gibbs/pi/zhao/xs282/validation/marginal_fit')){
  dir.create('/gpfs/gibbs/pi/zhao/xs282/validation/marginal_fit')
}

marginal_fit_fn <- '/gpfs/gibbs/pi/zhao/xs282/validation/marginal_fit/PNAS_NC_Oli_marginal_fit.rds'

if(!file.exists(marginal_fit_fn)){
  gp_ex <- glm_gp(count_sel, size_factors = size_factors_sel, verbose = T, overdispersion_shrinkage = T, do_cox_reid_adjustment = T)
  vanilla <- data.frame(mu = log10(exp(gp_ex$Beta[,1])),
                        alpha = -log10(gp_ex$overdispersions),
                        deviances = gp_ex$deviances,
                        gene=rownames(count_sel))
  saveRDS(vanilla, marginal_fit_fn)
  print("succeed")
}

