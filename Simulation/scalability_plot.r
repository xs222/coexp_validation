library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")

ncells = c(8000, 16000, 24000)
ngenes = c(9000, 18000, 27000)
ncor_genes = c(500, 1000, 1500)
res_df = matrix(nrow=27, ncol=23)
colnames(res_df) = c("ncell", "ngene", "ncor_gene", "para_peak_ram_mb", "para_normal_time",
                     "para_user_time", "para_system_time", "sim_peak_ram_mb", "sim_normal_time",
                     "sim_user_time", "sim_system_time", "noise_peak_ram_mb", "noise_normal_time",
                     "noise_user_time", "noise_system_time", "cscore_peak_ram_mb", "cscore_normal_time",
                     "cscore_user_time", "cscore_system_time", "ana_peak_ram_mb", "ana_normal_time",
                     "ana_user_time", "ana_system_time")
i = 0
for (ncell in ncells){
  for (ngene in ngenes){
    res = readRDS(paste0("BIB_R1/scalability/res/ncell", ncell, "_ngene", ngene, ".rds"))
    for (ncor_gene in ncor_genes){
      i = i+1
      bench_noise = res$bench_noise_ls[[paste0("ncor", ncor_gene)]]
      bench_cscore = res$bench_cscore_ls[[paste0("ncor", ncor_gene)]]
      bench_ana = res$bench_ana_ls[[paste0("ncor", ncor_gene)]]
      cpu_noise = res$cpu_noise_ls[[paste0("ncor", ncor_gene)]]
      cpu_cscore = res$cpu_cscore_ls[[paste0("ncor", ncor_gene)]]
      cpu_ana = res$cpu_ana_ls[[paste0("ncor", ncor_gene)]]

      res_df[i, ] = c(ncell, ngene, ncor_gene, res$bench_para[1,"Peak_RAM_Used_MiB"],
                      res$cpu_para["elapsed"], res$cpu_para["user.self"], res$cpu_para["sys.self"],
                      res$bench_sim[1,"Peak_RAM_Used_MiB"],
                      res$cpu_sim["elapsed"], res$cpu_sim["user.self"], res$cpu_sim["sys.self"],
                      bench_noise[1,"Peak_RAM_Used_MiB"], cpu_noise["elapsed"],
                      cpu_noise["user.self"],cpu_noise["sys.self"],bench_cscore[1,"Peak_RAM_Used_MiB"],
                      cpu_cscore["elapsed"], cpu_cscore["user.self"], cpu_cscore["sys.self"],
                      bench_ana[1,"Peak_RAM_Used_MiB"], cpu_ana["elapsed"],
                      cpu_ana["user.self"], cpu_ana["sys.self"])
    }
  }
}



df <- as.data.frame(res_df)
df[] <- lapply(df, function(x) as.numeric(as.character(x)))

ncells_sorted <- sort(unique(df$ncell))
ngenes_sorted <- sort(unique(df$ngene))
ncor_levels <- sort(unique(df$ncor_gene))

plot_data <- df %>%
  pivot_longer(cols = c(noise_peak_ram_mb, cscore_peak_ram_mb, ana_peak_ram_mb),
               names_to = "Tool_Type",
               values_to = "Tool_RAM") %>%
  mutate(Tool_Name = case_when(
    Tool_Type == "noise_peak_ram_mb" ~ "Noise Regularization",
    Tool_Type == "cscore_peak_ram_mb" ~ "CS-CORE",
    Tool_Type == "ana_peak_ram_mb"   ~ "Analytic PR")) %>%
  dplyr::rename(Simulation_RAM = sim_peak_ram_mb) %>%
  pivot_longer(cols = c(Simulation_RAM, Tool_RAM),
               names_to = "Component",
               values_to = "RAM_MiB") %>%
  mutate(ncell_f = factor(paste0("Cells: ", ncell),
                          levels = paste0("Cells: ", ncells_sorted)),
         ngene_f = factor(paste0("Genes: ", ngene),
                          levels = paste0("Genes: ", ngenes_sorted)),
         Component = factor(Component, levels = c("Tool_RAM", "Simulation_RAM"),
                            labels = c("Estimation", "Simulation")),
         Tool_Name = factor(Tool_Name, levels = c("Noise Regularization", "CS-CORE", "Analytic PR")))


plot_data <- plot_data %>%
  mutate(x_idx = as.numeric(factor(ncor_gene, levels = ncor_levels)),
         offset = case_when(
           Tool_Name == "Noise Regularization" ~ -0.25,
           Tool_Name == "CS-CORE"             ~ 0.00,
           Tool_Name == "Analytic PR"          ~ 0.25
         ),
         x_pos = x_idx + offset
  )

p1 = ggplot(plot_data, aes(x = x_pos, y = RAM_MiB / 1024, fill = Tool_Name, alpha = Component)) +
  geom_col(width = 0.2, color = "black", size = 0.1) +
  facet_grid(ncell_f ~ ngene_f, scales = "free_y") +
  scale_x_continuous(breaks = 1:length(ncor_levels), labels = ncor_levels) +
  scale_alpha_manual(values = c("Estimation" = 1, "Simulation" = 0.3)) +
  scale_fill_manual(values = c("Noise Regularization" = "#E41A1C",
                               "CS-CORE" = "#377EB8",
                               "Analytic PR" = "#4DAF4A")) +
  theme_bw() +
  labs(x = "Number of Correlated Genes", y = "Peak RAM (GB)",
       fill = "Estimation Method",
       alpha = "Stack Component") +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray95"),
        panel.grid.minor = element_blank())

# normal time ----------------------------------------------------------
plot_data_time <- df %>%
  pivot_longer(
    cols = c(noise_normal_time, cscore_normal_time, ana_normal_time),
    names_to = "Tool_Type",
    values_to = "Tool_Time"
  ) %>%
  mutate(Tool_Name = case_when(
    Tool_Type == "noise_normal_time" ~ "Noise Regularization",
    Tool_Type == "cscore_normal_time" ~ "CS-CORE",
    Tool_Type == "ana_normal_time"   ~ "Analytic PR"
  )) %>%
  dplyr::rename(Simulation_Time = sim_normal_time) %>%
  pivot_longer(
    cols = c(Simulation_Time, Tool_Time),
    names_to = "Component",
    values_to = "Seconds"
  ) %>%
  mutate(
    # Set sorted factor levels for 3x3 grid
    ncell_f = factor(paste0("Cells: ", ncell), levels = paste0("Cells: ", ncells_sorted)),
    ngene_f = factor(paste0("Genes: ", ngene), levels = paste0("Genes: ", ngenes_sorted)),
    # Stack order: Simulation on bottom, Estimation on top
    Component = factor(Component, levels = c("Tool_Time", "Simulation_Time"),
                       labels = c("Estimation", "Simulation")),
    Tool_Name = factor(Tool_Name, levels = c("Noise Regularization", "CS-CORE", "Analytic PR"))
  )

# 3. Calculate numeric X-positions for dodged stacks
plot_data_time <- plot_data_time %>%
  mutate(
    x_idx = as.numeric(factor(ncor_gene, levels = ncor_levels)),
    offset = case_when(
      Tool_Name == "Noise Regularization" ~ -0.25,
      Tool_Name == "CS-CORE"              ~ 0.00,
      Tool_Name == "Analytic PR"           ~ 0.25
    ),
    x_pos = x_idx + offset
  )

# 4. Create the 3x3 Time Plot (Minutes)
p2 = ggplot(plot_data_time, aes(x = x_pos, y = Seconds / 60, fill = Tool_Name, alpha = Component)) +
  geom_col(width = 0.2, color = "black", size = 0.1) +
  facet_grid(ncell_f ~ ngene_f, scales = "free_y") +
  scale_x_continuous(breaks = 1:length(ncor_levels), labels = ncor_levels) +
  scale_alpha_manual(values = c("Estimation" = 1, "Simulation" = 0.3)) +
  scale_fill_manual(values = c("Noise Regularization" = "#E41A1C",
                               "CS-CORE" = "#377EB8",
                               "Analytic PR" = "#4DAF4A")) +
  theme_bw() +
  labs(x = "Number of Correlated Genes",
       y = "Elapsed Time (Minutes)",
       fill = "Estimation Method",
       alpha = "Stack Component"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95"),
    panel.grid.minor = element_blank()
  )


pdf('BIB_R1/modify_plot/scalability.pdf', width = 10, height = 8, onefile = T)
ggarrange(p2,p1,ncol=1, nrow=2, common.legend = T, legend="bottom",labels = c("A", "B"))
dev.off()
