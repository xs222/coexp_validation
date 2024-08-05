
library(dplyr)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)

pnas_path <- "/gpfs/gibbs/pi/zhao/xs282/validation/compare_simulation/compare_simu_multi_data/PNAS_NC/"
pnas_path <- paste0(list.files(pnas_path, "\\.rds"))
rosmap_path <- "/gpfs/gibbs/pi/zhao/xs282/validation/compare_simulation/compare_simu_multi_data/ROSMAP/"
rosmap_path <- paste0(list.files(rosmap_path, "\\.rds"))
data_path <- c(pnas_path, rosmap_path)
data_path <- sapply(data_path, function(x){gsub(".rds", "", x)})
data_path <- data_path[!(data_path %in% c("ROSMAP_NC_in", "PNAS_NC_MIC"))]

kde_result <- list()
ks_result <- list()

for (data_name in data_path){
  # kde
  kde_stat <- readRDS(paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval_kde.rds"))
  kde_stat_scdesign3 <- readRDS(paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval_kde_scDesign3.rds"))
  kde_stat <- rbind(kde_stat, kde_stat_scdesign3)
  kde_stat$simu_method <- rownames(kde_stat)
  kde_stat_long <- reshape2::melt(kde_stat, id.vars = c("data", "simu_method"))
  kde_stat_long_norm <- kde_stat_long %>% group_by(variable) %>%
    summarise(Norm=(value-min(value))/(max(value)-min(value)),
              simu_method=simu_method,
              data=data, original=value)
  kde_result[[data_name]] <- kde_stat_long_norm

  # ks
  ks_stat <- readRDS(paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval_ks.rds"))
  ks_stat_scdesign3 <- readRDS(paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval_ks_scDesign3.rds"))
  ks_stat <- rbind(ks_stat, ks_stat_scdesign3)
  ks_stat$simu_method <- rownames(ks_stat)
  ks_stat_long <- reshape2::melt(ks_stat, id.vars = c("data", "simu_method"))
  ks_stat_long_norm <- ks_stat_long %>% group_by(variable) %>%
    summarise(Norm=(value-min(value))/(max(value)-min(value)),
              simu_method=simu_method,
              data=data, original=value)
  ks_result[[data_name]] <- ks_stat_long_norm
}

kde_result <- do.call(rbind, kde_result)
ks_result <- do.call(rbind, ks_result)

ind_method <- c("Splat", "powsimR", "SPARSim", "muscat", "SCRIP")
kde_result <- kde_result[kde_result$simu_method %in% ind_method, ]
ks_result <- ks_result[ks_result$simu_method %in% ind_method, ]


# ind data
kde_result_ind <- list()
ks_result_ind <- list()

for (data_name in data_path){
  # kde
  kde_stat <- readRDS(paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval_kde_IND.rds"))
  kde_stat$simu_method <- rownames(kde_stat)
  kde_stat_long <- reshape2::melt(kde_stat, id.vars = c("data", "simu_method"))
  kde_stat_long_norm <- kde_stat_long %>% group_by(variable) %>%
    summarise(Norm=(value-min(value))/(max(value)-min(value)),
              simu_method=simu_method,
              data=data, original=value)
  kde_result_ind[[data_name]] <- kde_stat_long_norm

  # ks
  ks_stat <- readRDS(paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval_ks_IND.rds"))
  ks_stat$simu_method <- rownames(ks_stat)
  ks_stat_long <- reshape2::melt(ks_stat, id.vars = c("data", "simu_method"))
  ks_stat_long_norm <- ks_stat_long %>% group_by(variable) %>%
    summarise(Norm=(value-min(value))/(max(value)-min(value)),
              simu_method=simu_method,
              data=data, original=value)
  ks_result_ind[[data_name]] <- ks_stat_long_norm
}

kde_result_ind <- do.call(rbind, kde_result_ind)
ks_result_ind <- do.call(rbind, ks_result_ind)

kde_result <- rbind(kde_result, kde_result_ind)
ks_result <- rbind(ks_result, ks_result_ind)


kde_mean_across_data <- kde_result %>% group_by(variable, simu_method) %>%
  summarise(mean_norm=mean(Norm))
kde_mean_across_data$simu_method[kde_mean_across_data$simu_method=="scDesign3"] <- "   scDesign3"
kde_mean_across_data <- kde_mean_across_data[!(kde_mean_across_data$variable %in% c("cell_knn", "cell_distance")),]
mean_order_kde <- kde_mean_across_data %>% group_by(simu_method) %>%
  summarise(mean(mean_norm))
mean_order_kde <- mean_order_kde[order(mean_order_kde$`mean(mean_norm)`),]
kde_mean_across_data$simu_method <- factor(kde_mean_across_data$simu_method, ordered=TRUE,
                                            levels = mean_order_kde$simu_method)
kde_mean_across_data$Property <- recode(kde_mean_across_data$variable,
                                        gene_mean='Mean(logCPM)',
                                        gene_var="Var(logCPM)",
                                        gene_cv="Coefficient of variation",
                                        gene_frq_zero="Zero fraction (gene)",
                                        gene_cor = "Gene correlation",
                                        cell_frq_zero = 'Zero fraction (cell)',
                                        cell_cor = 'Cell correlation',
                                        lib_size="Library size")
kde_mean_across_data$simu_method <- recode(kde_mean_across_data$simu_method, NB="scSimu")

ks_mean_across_data <- ks_result %>% group_by(variable, simu_method) %>%
  summarise(mean_norm=mean(Norm))
ks_mean_across_data$simu_method[ks_mean_across_data$simu_method=="scDesign3"] <- "   scDesign3"
ks_mean_across_data <- ks_mean_across_data[!(ks_mean_across_data$variable %in% c("cell_knn", "cell_distance")),]
mean_order_ks <- ks_mean_across_data %>% group_by(simu_method) %>%
  summarise(mean(mean_norm))
mean_order_ks <- mean_order_ks[order(mean_order_ks$`mean(mean_norm)`),]
ks_mean_across_data$simu_method <- factor(ks_mean_across_data$simu_method, ordered=TRUE,
                                           levels = mean_order_ks$simu_method)
ks_mean_across_data$Property <- recode(ks_mean_across_data$variable,
                                       gene_mean='Mean(logCPM)',
                                       gene_var="Var(logCPM)",
                                       gene_cv="Coefficient of variation",
                                       gene_frq_zero="Zero fraction (gene)",
                                       gene_cor = "Gene correlation",
                                       cell_frq_zero = 'Zero fraction (cell)',
                                       cell_cor = 'Cell correlation',
                                       lib_size="Library size")
ks_mean_across_data$simu_method <- recode(ks_mean_across_data$simu_method, NB="scSimu")


kde_ks_mean_across_data <- left_join(ks_result[,1:4], kde_result[,1:4],
                                     by=c("variable", "simu_method", "data"))
kde_ks_mean_across_data$simu_method[kde_ks_mean_across_data$simu_method=="scDesign3"] <- "   scDesign3"
kde_ks_mean_across_data$Norm <- (kde_ks_mean_across_data$Norm.x+kde_ks_mean_across_data$Norm.y)/2
kde_ks_mean_across_data <- kde_ks_mean_across_data %>% group_by(variable, simu_method) %>%
  summarise(mean_norm=mean(Norm))
kde_ks_mean_across_data <- kde_ks_mean_across_data[!(kde_ks_mean_across_data$variable %in% c("cell_knn", "cell_distance")),]
mean_order_ks <- kde_ks_mean_across_data %>% group_by(simu_method) %>%
  summarise(mean(mean_norm))
mean_order_ks <- mean_order_ks[order(mean_order_ks$`mean(mean_norm)`),]
kde_ks_mean_across_data$simu_method <- factor(kde_ks_mean_across_data$simu_method, ordered=TRUE,
                                          levels = mean_order_ks$simu_method)
kde_ks_mean_across_data$Property <- recode(kde_ks_mean_across_data$variable,
                                           gene_mean='Mean(logCPM)',
                                           gene_var="Var(logCPM)",
                                           gene_cv="Coefficient of variation",
                                           gene_frq_zero="Zero fraction (gene)",
                                           gene_cor = "Gene correlation",
                                           cell_frq_zero = 'Zero fraction (cell)',
                                           cell_cor = 'Cell correlation',
                                           lib_size="Library size")
kde_ks_mean_across_data$simu_method <- recode(kde_ks_mean_across_data$simu_method, NB="scSimu")



p_kde <- ggplot(kde_mean_across_data, aes(x=simu_method, y=Property, fill=mean_norm))+
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",
                       na.value = "grey",limits = c(0, 1),breaks = c(0, 1))+
  theme_classic()+
  theme(legend.position="bottom",
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust=0.5))+
  labs(title="KDE Test (IND)", fill="Normalized test statistics")


p_ks <- ggplot(ks_mean_across_data, aes(x=simu_method, y=Property, fill=mean_norm))+
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",
                       na.value = "grey",limits = c(0, 1),breaks = c(0, 1))+
  theme_classic()+
  theme(legend.position="bottom",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust=0.5))+
  labs(title="KS Test (IND)", fill="Normalized test statistics")

p_ks_kde <- ggplot(kde_ks_mean_across_data, aes(x=simu_method, y=Property, fill=mean_norm))+
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",
                       na.value = "grey",limits = c(0, 1),breaks = c(0, 1))+
  theme_classic()+
  theme(legend.position="bottom",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust=0.5))+
  labs(title="Independent", fill="Normalized test statistics")

format_supp1 <- theme(text = element_text(size = 14),
                      legend.spacing.x = unit(2, 'cm'),
                      axis.text.x = element_text(angle = 60, hjust = 1))
format_supp2 <- theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank())

pdf('figures/v2/bench_simu_v2.pdf', width = 12, height = 3.8, onefile = T)
ggarrange(ggarrange(p_ks_kde+format_supp1, p_ks_kde_cor+format_supp1+format_supp2,
                    ncol=2, nrow=1, common.legend = TRUE, legend="bottom",
                    labels = c("A", "B"), widths = c(1, 0.70)),
          ggarrange(p_time_ind+format_supp1,p_time_cor+format_supp1,
                    ncol=2, nrow=1, labels = c("C", "D"), widths = c(0.85, 1)),
          ncol=2, nrow=1, widths = c(1,0.9))
dev.off()

pdf('figures/v2/bench_simu_supp_v2.pdf', width = 12, height = 3.8, onefile = T)
ggarrange(p_ks+format_supp1,
          p_kde+format_supp1+format_supp2,
          p_ks_cor+format_supp1+format_supp2,
          p_kde_cor+format_supp1+format_supp2,
          ncol=4, nrow=1, common.legend = TRUE, legend="bottom",
          labels = c("A", "B", "C", "D"),
          widths = c(0.9,0.6,0.7,0.7))
dev.off()

# correlation structure ---------------------------------------------------
data_name <- "PNAS_NC_EX_sel5000"
simu_ind <- readRDS(paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval_IND.rds"))
simu_cor <- readRDS(paste0("compare_simulation/compare_simu_multi_data/evaluation/", data_name, "_eval.rds"))


vector2matrix <- function(upper_tri_vec, n=1000){
  new_matrix <- matrix(NA, n, n)
  new_matrix[upper.tri(new_matrix)] <- upper_tri_vec
  new_matrix[lower.tri(new_matrix)] <- t(new_matrix)[lower.tri(new_matrix)]
  diag(new_matrix) <- 1
  return(new_matrix)
}

cor_mx_list <- list()
simu_method <- as.character(unique(kde_mean_across_data$simu_method))
simu_method[simu_method=="scSimu"] <- "NB"
for (i in simu_method){
  print(paste(i, "----------------------------------------------"))
  if (i %in% c("Splat", "powsimR", "SPARSim", "muscat", "SCRIP")){
    simu_dist <- simu_cor[[i]]
  } else {
    simu_dist <- simu_ind[[i]]
  }

  cor_mx_list[[i]] <- vector2matrix(simu_dist[["gene_cor"]])
}


for (i in simu_method){
  mat <- cor_mx_list[[i]]
  print(paste0(i, " density: ", mean(mat[upper.tri(mat)]!=0)))
}

my_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(100)
my_breaks <- seq(-1, 1, length.out = length(my_palette) + 1)

heat_list <- list()

for (i in sort(names(cor_mx_list))){
  heat <- pheatmap(cor_mx_list[[i]], show_rownames = F, show_colnames = F,
                   cluster_rows = F, cluster_cols = F, main=i,
                   breaks = my_breaks,color = my_palette, legend=T)
  heat_list[[i]] <- heat[[4]]
}

heat_list[["NB"]] <- pheatmap(cor_mx_list[["NB"]], show_rownames = F, show_colnames = F,
                              cluster_rows = F, cluster_cols = F, main="scSimu",
                              breaks = my_breaks,color = my_palette, legend=T)[[4]]

pdf('compare_simulation/compare_simu_multi_data/plot/simu_ind_cor.pdf', width = 13, height = 5)
grid.arrange(arrangeGrob(grobs= heat_list,ncol=5, nrow=2,
                         common.legend = TRUE, legend="bottom"))
dev.off()
# manually convert to png
