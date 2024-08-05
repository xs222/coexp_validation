

library(dplyr)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)

setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
pnas_path <- "compare_simulation/compare_simu_multi_data/PNAS_NC/"
pnas_path <- paste0(pnas_path, list.files(pnas_path, "\\.rds"))
rosmap_path <- "compare_simulation/compare_simu_multi_data/ROSMAP/"
rosmap_path <- paste0(rosmap_path, list.files(rosmap_path, "\\.rds"))
data_path <- c(pnas_path, rosmap_path)
data_path <- sapply(data_path, function(x){gsub(".rds", "", x)})
data_path <- data_path[!(data_path %in% c("compare_simulation/compare_simu_multi_data/ROSMAP/ROSMAP_NC_in",
                                          "compare_simulation/compare_simu_multi_data/PNAS_NC/PNAS_NC_MIC"))]


# correlated
time_list <- list()
for(i in data_path){
  time <- readRDS(paste0(i, "/time.rds"))
  time <- time[!is.na(time[,1]),]
  time_SPARSim <- readRDS(paste0(i, "/time_SPARSim.rds"))
  time_hierarchicell <- readRDS(paste0(i, "/time_hierarchicell.rds"))
  time <- rbind(time, time_SPARSim, time_hierarchicell)
  time[,2] <- as.numeric(time[,2])
  time_scdesign3 <- readRDS(paste0(i, "/time_scDesign3.rds"))
  time_scdesign3[,2] <- as.numeric(time_scdesign3[,2])*60*60
  time <- rbind(time, time_scdesign3)
  time <- as.data.frame(time)
  time$data <- i

  time_muscat <- readRDS(paste0(i, "/time_muscat.rds"))
  time[time$V1=="muscat",2] <- as.numeric(time_muscat[1,2])

  time_list[[i]] <- time
}

time <- do.call(rbind, time_list)
time$V1 <- recode(time$V1,NB="scSimu")
time$V2 <- as.numeric(time$V2)
time_order <- time %>% group_by(V1) %>% summarise(mean(V2))
time_order <- time_order[order(time_order$`mean(V2)`),]
time$V1 <- factor(time$V1, ordered=TRUE, levels = time_order$V1)


max(log10(as.numeric(time$V2)))
p_time_cor <- ggplot(time, aes(x=V1, y=log10(as.numeric(V2)), fill=V1))+
  geom_boxplot()+theme_classic()+ylab("log10(Time in seconds)")+
  xlab("Simulation method")+
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust=0.5))+
  labs(title="Correlated", x="")+ylim(0.5,5)

# ind ----------------------------------
ind_method <- c("Splat", "powsimR", "SPARSim", "muscat", "SCRIP")
time_ind1 <- time[time$V1 %in% ind_method, ]

time_list_ind <- list()
for(i in data_path){
  time <- readRDS(paste0(i, "_IND/time.rds"))
  time_permu <- readRDS(paste0(i, "_IND/time_permu.rds"))
  time_scdesign3 <- readRDS(paste0(i, "_IND/time_scDesign3.rds"))

  time <- rbind(time, time_permu, time_scdesign3)
  time <- as.data.frame(time)
  time$data <- i

  time_list_ind[[i]] <- time
}

time <- do.call(rbind, time_list_ind)
time <- rbind(time, time_ind1)
time$V1 <- recode(time$V1,NB="scSimu")
time$V2 <- as.numeric(time$V2)
time_order <- time %>% group_by(V1) %>% summarise(mean(V2))
time_order <- time_order[order(time_order$`mean(V2)`),]
time$V1 <- factor(time$V1, ordered=TRUE, levels = time_order$V1)

p_time_ind <- ggplot(time, aes(x=V1, y=log10(as.numeric(V2)), fill=V1))+
  geom_boxplot()+theme_classic()+ylab("log10(time in seconds)")+
  xlab("Simulation method")+
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust=0.5))+
  labs(title="Independent", x="")+ylim(0.5,5)


p_compu_time <- ggarrange(p_time_ind,p_time_cor, ncol=2, nrow=1, labels = c("C", "D"))
