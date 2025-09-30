library(Matrix)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(Seurat)
library(scDesign3)
library(SingleCellExperiment)

count_data = readMM("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/counts.read.txt")
cells = read.table("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/cells.read.new.txt")
genes = read.table("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/genes.read.txt")
meta = read.table("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/meta.txt", sep = "\t", header=T)
colnames(count_data) = cells$V1
rownames(count_data) = genes$V1
meta_sub = meta[meta$Method=="Smart-seq2" & meta$Experiment=="pbmc1",]
meta_sub

ct = "Cytotoxic T cell"
cell_name = meta_sub$NAME[meta_sub$CellType==ct]
count_sub = count_data[,cell_name]
count_sub = as.matrix(count_sub)

# filter genes
gene_idx = which((rowSums(count_sub)>=10) & (rowSums(count_sub>=1)>=10))
nrow(count_sub)-length(gene_idx)
count_sub = count_sub[gene_idx,]

# ---------------------------------------------------------------------------------------------------
sce <- SingleCellExperiment(assays = list(counts = count_sub))
sce$celltype <- ct
colData(sce)$library = colSums(counts(sce))

path <- "/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/simu/cyto_t/"
if(!file.exists(path)){
  dir.create(path,recursive = T)
}

# simulate correlated data
set.seed(10222023)
scDesign3_start_time <- Sys.time()
example_simu <- scdesign3(sce = sce, celltype = "celltype",
  pseudotime = NULL, spatial = NULL,family_use = "zinb",
  other_covariates = "library",
  mu_formula = "offset(log(library))",
  corr_formula = "1",
  parallelization = "pbmcmapply")
scDesign3_end_time <- Sys.time()

saveRDS(example_simu, paste0(path, "simu_scDesign3.rds"))
print(scDesign3_end_time-scDesign3_start_time)

# simulate independent data
set.seed(10222023)
scDesign3_start_time <- Sys.time()
example_simu <- scdesign3(sce = sce, celltype = "celltype",
  pseudotime = NULL, spatial = NULL,family_use = "zinb",
  other_covariates = "library",
  mu_formula = "offset(log(library))",
  corr_formula = "ind",
  parallelization = "pbmcmapply")
scDesign3_end_time <- Sys.time()

saveRDS(example_simu, paste0(path, "simu_scDesign3_ind.rds"))
print(scDesign3_end_time-scDesign3_start_time)
