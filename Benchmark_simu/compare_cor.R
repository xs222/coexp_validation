
library(dplyr)
library(ggplot2)

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

kde_mean_across_data <- kde_result %>% group_by(variable, simu_method) %>%
  summarise(mean_norm=mean(Norm))
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



p_kde_cor <- ggplot(kde_mean_across_data, aes(x=simu_method, y=Property, fill=mean_norm))+
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
  labs(title="KDE Test (COR)", fill="Normalized test statistics")


p_ks_cor <- ggplot(ks_mean_across_data, aes(x=simu_method, y=Property, fill=mean_norm))+
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
  labs(title="KS Test (COR)", fill="Normalized test statistics")

p_ks_kde_cor <- ggplot(kde_ks_mean_across_data, aes(x=simu_method, y=Property, fill=mean_norm))+
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
  labs(title="Correlated", fill="Normalized test statistics")

ggarrange(p_kde+theme(legend.spacing.x = unit(2, 'cm')),
          p_ks+theme(legend.spacing.x = unit(2, 'cm')),
          p_ks_kde+theme(legend.spacing.x = unit(2, 'cm')),
          ncol=3, nrow=1, common.legend = TRUE, legend="bottom")





