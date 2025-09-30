
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)

cor_res = readRDS("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/simu/cyto_t/ks_kde.rds")

ks_stat = cor_res$ks
ks_stat$simu_method <- rownames(ks_stat)
ks_stat_long <- reshape2::melt(ks_stat, id.vars = c("simu_method"))

kde_stat = cor_res$kde
kde_stat$simu_method <- rownames(kde_stat)
kde_stat_long <- reshape2::melt(kde_stat, id.vars = c("simu_method"))

mean_order_ks <- rownames(ks_stat)[order(rowMeans(ks_stat[,1:8]))]
ks_stat_long$simu_method <- factor(ks_stat_long$simu_method, ordered=TRUE, levels = mean_order_ks)
ks_stat_long$variable <- recode(ks_stat_long$variable,
                               gene_mean='Mean(logCPM)',
                               gene_var="Var(logCPM)",
                               gene_cv="Coefficient of variation",
                               gene_frq_zero="Zero fraction (gene)",
                               gene_cor = "Gene correlation",
                               cell_frq_zero = 'Zero fraction (cell)',
                               cell_cor = 'Cell correlation',
                               lib_size="Library size")

mean_order_kde <- rownames(kde_stat)[order(rowMeans(kde_stat[,1:8]))]
kde_stat_long$simu_method <- factor(kde_stat_long$simu_method, ordered=TRUE, levels = mean_order_kde)
kde_stat_long$variable <- recode(kde_stat_long$variable,
                               gene_mean='Mean(logCPM)',
                               gene_var="Var(logCPM)",
                               gene_cv="Coefficient of variation",
                               gene_frq_zero="Zero fraction (gene)",
                               gene_cor = "Gene correlation",
                               cell_frq_zero = 'Zero fraction (cell)',
                               cell_cor = 'Cell correlation',
                               lib_size="Library size")


p_kde_cor <- ggplot(kde_stat_long, aes(x=simu_method, y=variable, fill=value))+
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", na.value = "grey",
                       limits = c(min(kde_stat_long$value), max(kde_stat_long$value)),
                       breaks = c(min(kde_stat_long$value), max(kde_stat_long$value)),
                      labels = label_number(accuracy = 1))+
  theme_classic()+
  theme(legend.position="right",
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust=0.5))+
  labs(title="KDE Test (COR)", fill="Test stat")


p_ks_cor <- ggplot(ks_stat_long, aes(x=simu_method, y=variable, fill=value))+
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", na.value = "grey",
                       limits = c(min(ks_stat_long$value), max(ks_stat_long$value)),
                       breaks = c(min(ks_stat_long$value), max(ks_stat_long$value)),
                      labels = label_number(accuracy = 1))+
  theme_classic()+
  theme(legend.position="right",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust=0.5))+
  labs(title="KS Test (COR)", fill="Test stat")

ind_res = readRDS("/gpfs/gibbs/pi/zhao/xs282/validation/revision/smart_seq_data/simu/cyto_t/ks_kde_ind.rds")

ks_stat_ind = ind_res$ks
ks_stat_ind$simu_method <- rownames(ks_stat_ind)
ks_stat_long_ind <- reshape2::melt(ks_stat_ind, id.vars = c("simu_method"))

kde_stat_ind = ind_res$kde
kde_stat_ind$simu_method <- rownames(kde_stat_ind)
kde_stat_long_ind <- reshape2::melt(kde_stat_ind, id.vars = c("simu_method"))

mean_order_ks <- rownames(ks_stat_ind)[order(rowMeans(ks_stat_ind[,1:8]))]
ks_stat_long_ind$simu_method <- factor(ks_stat_long_ind$simu_method, ordered=TRUE, levels = mean_order_ks)
ks_stat_long_ind$variable <- recode(ks_stat_long_ind$variable,
                               gene_mean='Mean(logCPM)',
                               gene_var="Var(logCPM)",
                               gene_cv="Coefficient of variation",
                               gene_frq_zero="Zero fraction (gene)",
                               gene_cor = "Gene correlation",
                               cell_frq_zero = 'Zero fraction (cell)',
                               cell_cor = 'Cell correlation',
                               lib_size="Library size")

mean_order_kde <- rownames(kde_stat_ind)[order(rowMeans(kde_stat_ind[,1:8]))]
kde_stat_long_ind$simu_method <- factor(kde_stat_long_ind$simu_method, ordered=TRUE, levels = mean_order_kde)
kde_stat_long_ind$variable <- recode(kde_stat_long_ind$variable,
                               gene_mean='Mean(logCPM)',
                               gene_var="Var(logCPM)",
                               gene_cv="Coefficient of variation",
                               gene_frq_zero="Zero fraction (gene)",
                               gene_cor = "Gene correlation",
                               cell_frq_zero = 'Zero fraction (cell)',
                               cell_cor = 'Cell correlation',
                               lib_size="Library size")

p_kde_ind <- ggplot(kde_stat_long_ind, aes(x=simu_method, y=variable, fill=value))+
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", na.value = "grey",
                       limits = c(min(kde_stat_long_ind$value), max(kde_stat_long_ind$value)),
                       breaks = c(min(kde_stat_long_ind$value), max(kde_stat_long_ind$value)),
                      labels = label_number(accuracy = 1))+
  theme_classic()+
  theme(legend.position="right",
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust=0.5))+
  labs(title="KDE Test (IND)", fill="Test stat")


p_ks_ind <- ggplot(ks_stat_long_ind, aes(x=simu_method, y=variable, fill=value))+
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", na.value = "grey",
                       limits = c(min(ks_stat_long_ind$value), max(ks_stat_long_ind$value)),
                       breaks = c(min(ks_stat_long_ind$value), max(ks_stat_long_ind$value)),
                      labels = label_number(accuracy = 1))+
  theme_classic()+
  theme(legend.position="right",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust=0.5))+
  labs(title="KS Test (IND)", fill="Test stat")

format_supp1 <- theme(text = element_text(size = 14),
                      axis.text.x = element_text(angle = 60, hjust = 1))
pdf('/gpfs/gibbs/pi/zhao/xs282/validation/revision/modify_plot/simu_smart.pdf', width = 14, height = 3.5, onefile = T)
ggarrange(p_ks_ind+format_supp1, p_kde_ind+format_supp1, p_ks_cor+format_supp1, p_kde_cor+format_supp1,ncol=4, nrow=1,labels = c("A", "B", "C", "D"))
dev.off()



