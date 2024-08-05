library(SingleCellExperiment)
library(Seurat)
library(dplyr)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")

# PNAS -------------------------------------------------------------------------
PNAS_sc_obj <- readRDS("/gpfs/gibbs/pi/zhao/cs2629/AD_Nancy/seurat_obj_cell_type_labelled.rds")
table(PNAS_sc_obj@meta.data$cell_type)
PNAS_sc_obj$disease = sapply(PNAS_sc_obj$orig.ident, function(or) {substr(strsplit(or, '_')[[1]][2], 1, 2)})
table(PNAS_sc_obj$disease)
length(unique(PNAS_sc_obj$orig.ident[PNAS_sc_obj$disease=="AD"]))
length(unique(PNAS_sc_obj$orig.ident[PNAS_sc_obj$disease=="NC"]))
PNAS_sc_obj_nc <- subset(PNAS_sc_obj, subset = disease=="NC")
length(unique(PNAS_sc_obj_nc$orig.ident))
table(PNAS_sc_obj_nc@meta.data$cell_type)

PNAS_meta <- PNAS_sc_obj_nc@meta.data
write.csv(PNAS_meta, "compare_simulation/compare_simu_multi_data/PNAS_NC/PNAS_NC_meta.csv")

## simulate cell types with larger than 5000 cells: ast, ex, in, mic, oli
## prepare data: ex in control-------------------------------------------------
count_ex <- as.matrix(PNAS_sc_obj_nc[["RNA"]]@counts[, which(PNAS_sc_obj_nc$cell_type=="Ex")])
print(dim(count_ex))
## subset 5000 cells
set.seed(10222023)
cell_idx_ex <- sample(1:ncol(count_ex), size=5000)
ori_ct_ex <- count_ex[, cell_idx_ex]
saveRDS(ori_ct_ex, "compare_simulation/compare_simu_multi_data/PNAS_NC/PNAS_NC_EX_sel5000.rds")

## prepare data: oli in control-------------------------------------------------
count_oli <- as.matrix(PNAS_sc_obj_nc[["RNA"]]@counts[, which(PNAS_sc_obj_nc$cell_type=="Oli")])
print(dim(count_oli))
## subset 5000 cells
set.seed(10222023)
cell_idx_oli <- sample(1:ncol(count_oli), size=5000)
ori_ct_oli <- count_oli[, cell_idx_oli]
saveRDS(ori_ct_oli, "compare_simulation/compare_simu_multi_data/PNAS_NC/PNAS_NC_OLI_sel5000.rds")

## prepare data: ast in control-------------------------------------------------
count_ast <- as.matrix(PNAS_sc_obj_nc[["RNA"]]@counts[, which(PNAS_sc_obj_nc$cell_type=="Ast")])
print(dim(count_ast))
saveRDS(count_ast, "compare_simulation/compare_simu_multi_data/PNAS_NC/PNAS_NC_AST.rds")

## prepare data: In in control-------------------------------------------------
count_in <- as.matrix(PNAS_sc_obj_nc[["RNA"]]@counts[, which(PNAS_sc_obj_nc$cell_type=="In")])
print(dim(count_in))
saveRDS(count_in, "compare_simulation/compare_simu_multi_data/PNAS_NC/PNAS_NC_IN.rds")

## prepare data: Mic in control-------------------------------------------------
count_mic <- as.matrix(PNAS_sc_obj_nc[["RNA"]]@counts[, which(PNAS_sc_obj_nc$cell_type=="Mic")])
print(dim(count_mic))
saveRDS(count_mic, "compare_simulation/compare_simu_multi_data/PNAS_NC/PNAS_NC_MIC.rds")

## prepare data: End in control-------------------------------------------------
count_end <- as.matrix(PNAS_sc_obj_nc[["RNA"]]@counts[, which(PNAS_sc_obj_nc$cell_type=="End")])
print(dim(count_end))
saveRDS(count_end, "compare_simulation/compare_simu_multi_data/PNAS_NC/PNAS_NC_END.rds")

# ROSMAP -------------------------------------------------------------------------
clic <- read.csv("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/ROSMAP/Metadata/ROSMAP_clinical.csv")

## prepare data: oli in control-------------------------------------------------
sc_obj_oli <- readRDS("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/MIT_ROSMAP_Multiomics/Gene_Expression/snRNAseq/snRNAseq-10x/processed/Oligodendrocytes.rds")
sc_obj_meta_oli <- sc_obj_oli@meta.data
sc_obj_meta_oli <- left_join(sc_obj_meta_oli, clic)
sc_obj_oli$braaksc <- sc_obj_meta_oli$braaksc
write.csv(sc_obj_oli@meta.data, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_Oli_sel5000_meta.csv")

count_ROSMAP_oli <- as.matrix(sc_obj_oli[["RNA"]]@counts[, which(sc_obj_oli$braaksc==0)])
dim(count_ROSMAP_oli)
## subset 5000 cells
set.seed(10222023)
cell_idx_oli <- sample(1:ncol(count_ROSMAP_oli), size=5000)
ori_ct_oli <- count_ROSMAP_oli[, cell_idx_oli]
saveRDS(ori_ct_oli, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_Oli_sel5000.rds")

## prepare data: ast in control-------------------------------------------------
sc_obj_ast <- readRDS("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/MIT_ROSMAP_Multiomics/Gene_Expression/snRNAseq/snRNAseq-10x/processed/Astrocytes.rds")
sc_obj_meta_ast <- sc_obj_ast@meta.data
sc_obj_meta_ast <- left_join(sc_obj_meta_ast, clic)
sc_obj_ast$braaksc <- sc_obj_meta_ast$braaksc
write.csv(sc_obj_ast@meta.data, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_Ast_meta.csv")

count_ROSMAP_ast <- as.matrix(sc_obj_ast[["RNA"]]@counts[, which(sc_obj_ast$braaksc==0)])
dim(count_ROSMAP_ast)
saveRDS(count_ROSMAP_ast, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_Ast.rds")

## prepare data: ex in control-------------------------------------------------
sc_obj_ex <- readRDS("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/MIT_ROSMAP_Multiomics/Gene_Expression/snRNAseq/snRNAseq-10x/processed/Excitatory_neurons_set1.rds")
sc_obj_meta_ex <- sc_obj_ex@meta.data
sc_obj_meta_ex <- left_join(sc_obj_meta_ex, clic)
sc_obj_ex$braaksc <- sc_obj_meta_ex$braaksc
write.csv(sc_obj_ex@meta.data, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_ex_meta.csv")

count_ROSMAP_ex <- as.matrix(sc_obj_ex[["RNA"]]@counts[, which(sc_obj_ex$braaksc==0)])
dim(count_ROSMAP_ex)
saveRDS(count_ROSMAP_ex, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_ex.rds")

## prepare data: immune in control-------------------------------------------------
sc_obj_immune <- readRDS("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/MIT_ROSMAP_Multiomics/Gene_Expression/snRNAseq/snRNAseq-10x/processed/Immune_cells.rds")
sc_obj_meta_immune <- sc_obj_immune@meta.data
sc_obj_meta_immune <- left_join(sc_obj_meta_immune, clic)
sc_obj_immune$braaksc <- sc_obj_meta_immune$braaksc
write.csv(sc_obj_immune@meta.data, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_immune_meta.csv")

count_ROSMAP_immune <- as.matrix(sc_obj_immune[["RNA"]]@counts[, which(sc_obj_immune$braaksc==0)])
dim(count_ROSMAP_immune)
saveRDS(count_ROSMAP_immune, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_immune.rds")

## prepare data: in in control-------------------------------------------------
sc_obj_in <- readRDS("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/MIT_ROSMAP_Multiomics/Gene_Expression/snRNAseq/snRNAseq-10x/processed/Inhibitory_neurons.rds")
sc_obj_meta_in <- sc_obj_in@meta.data
sc_obj_meta_in <- left_join(sc_obj_meta_in, clic)
sc_obj_in$braaksc <- sc_obj_meta_in$braaksc
write.csv(sc_obj_in@meta.data, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_in_meta.csv")

count_ROSMAP_in <- as.matrix(sc_obj_in[["RNA"]]@counts[, which(sc_obj_in$braaksc==0)])
dim(count_ROSMAP_in)
saveRDS(count_ROSMAP_in, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_in.rds")

## prepare data: oligodendrocyte precursor cells (OPC) in control-------------------------------------------------
sc_obj_OPC <- readRDS("/gpfs/gibbs/pi/zhao/xs282/ROSMAP/MIT_ROSMAP_Multiomics/Gene_Expression/snRNAseq/snRNAseq-10x/processed/OPCs.rds")
sc_obj_meta_OPC <- sc_obj_OPC@meta.data
sc_obj_meta_OPC <- left_join(sc_obj_meta_OPC, clic)
sc_obj_OPC$braaksc <- sc_obj_meta_OPC$braaksc
write.csv(sc_obj_OPC@meta.data, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_OPC_meta.csv")

count_ROSMAP_OPC <- as.matrix(sc_obj_OPC[["RNA"]]@counts[, which(sc_obj_OPC$braaksc==0)])
dim(count_ROSMAP_OPC)
saveRDS(count_ROSMAP_OPC, "compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_OPC.rds")

