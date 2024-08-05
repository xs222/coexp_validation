
# fit kernel regression between mu and alpha

library(ggplot2)
setwd("/gpfs/gibbs/pi/zhao/xs282/validation/")
marginal_fit_ROSMAP = readRDS('marginal_fit/ROSMAP_NC_Oli_marginal_fit.rds')
marginal_fit_PNAS = readRDS('marginal_fit/PNAS_NC_Oli_marginal_fit.rds')

################# ROSMAP: fit the kernel regression ##################################
marginal_fit_ROSMAP$up <- ifelse(marginal_fit_ROSMAP$alpha>2, "upper","lower")

# plot alpha vs mu
plot1 <- ggplot(marginal_fit_ROSMAP, aes(x=mu, y=alpha))+
  geom_point()+labs(title = "log10(alpha) v.s. log10(mu)")+
  xlab("log10(mu)")+ylab("log10(alpha)")
plot2 <- ggplot(marginal_fit_ROSMAP, aes(x=mu, y=alpha, colour=up))+
  geom_point()+labs(title = "log10(alpha) v.s. log10(mu)")+
  xlab("log10(mu)")+ylab("log10(alpha)")
plot1
plot2

# fit line (only use the lower cluster)
marginal_fit_ROSMAP_sel <- marginal_fit_ROSMAP[marginal_fit_ROSMAP$up=="lower",]
# kernel smooth
km_ex5_ROSMAP <- ksmooth(marginal_fit_ROSMAP_sel$mu, marginal_fit_ROSMAP_sel$alpha, kernel="normal", bandwidth = bw.SJ(marginal_fit_ROSMAP_sel$mu)*5)
saveRDS(km_ex5_ROSMAP, 'marginal_fit/ROSMAP_NC_Oli_ks_fit_5.rds')

################# PNAS: fit the kernel regression ##################################
marginal_fit_PNAS$up <- ifelse(marginal_fit_PNAS$alpha>2, "upper","lower")

# plot alpha vs mu
plot1 <- ggplot(marginal_fit_PNAS, aes(x=mu, y=alpha))+
  geom_point()+labs(title = "log10(alpha) v.s. log10(mu)")+
  xlab("log10(mu)")+ylab("log10(alpha)")
plot2 <- ggplot(marginal_fit_PNAS, aes(x=mu, y=alpha, colour=up))+
  geom_point()+labs(title = "log10(alpha) v.s. log10(mu)")+
  xlab("log10(mu)")+ylab("log10(alpha)")
plot1
plot2

# fit line (only use the lower cluster)
marginal_fit_PNAS_sel <- marginal_fit_PNAS[marginal_fit_PNAS$up=="lower",]
# kernel smooth
km_ex5_PNAS <- ksmooth(marginal_fit_PNAS_sel$mu, marginal_fit_PNAS_sel$alpha, kernel="normal", bandwidth = bw.SJ(marginal_fit_PNAS_sel$mu)*5)
saveRDS(km_ex5_PNAS, 'marginal_fit/PNAS_NC_Oli_ks_fit_5.rds')
