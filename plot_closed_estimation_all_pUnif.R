#plot_closed_simulation_estimation.R

# Plots results of closed_cmr_simulation_estimation.R

# Tomo Eguchi
# 30 January 2017


rm(list=ls())
source('Cm_SDB_functions.R')
library('ggplot2')
#library('reshape2')
library(tidyverse)

save.plot <- F

RData.files <- c("RData/closed_simulation_pUNIF_k6_2018-01-19.RData",
                 "RData/closed_simulation_pUNIF_k8_2018-01-22.RData",
                 "RData/closed_simulation_pUNIF_k10_2018-01-22.RData",
                 "RData/closed_simulation_pUNIF_k12_2018-01-22.RData",
                 "RData/closed_simulation_pUNIF_k14_2018-01-22.RData",
                 "RData/closed_simulation_pUNIF_k16_2018-01-23.RData")

ks <- seq(from = 6, to = 16, by = 2)
k0 <- 1
# for each capture probability
sample.sizes.list.pUnif <- sim.results.all.list <-
  sim.results.data.list <- p.list <-
  n.list <- mean.n <- vector(mode = 'list',
                             length = length(RData.files))

for (k0 in 1:length(RData.files)){
  load(RData.files[k0])
  sim.results.all.list[[k0]] <- na.omit(do.call(rbind, sim.results.all)) %>%
    filter(., k == ks[k0])
  non.zero.idx <- unlist(lapply(lapply(sample.size, dim), length))
  tmp <- sample.size[non.zero.idx > 0]
  sample.sizes.list.pUnif[[k0]] <- tmp[[length(tmp)]]
  tmp <- sim.results.data[non.zero.idx > 0]

  sim.results.data.list[[k0]] <- tmp[[length(tmp)]]
  p.list[[k0]] <- p
  n.list[[k0]] <- matrix(unlist(lapply(sim.results.data.list[[k0]], colSums)),
                         nrow = 100, ncol = ks[k0], byrow = T)
  mean.n[[k0]] <- apply(n.list[[k0]], MARGIN = 2, FUN = mean)
}

p.list.1 <- list(p.list[[1]][[1]], p.list[[2]][[1]], p.list[[2]][[2]],
                 p.list[[2]][[3]], p.list[[2]][[4]])

all.data <- do.call(rbind, sim.results.all.list) %>%
  mutate(n = unlist(lapply(sample.sizes.list.pUnif,
                           rowSums)))



Nhats <- select(all.data, k, Ndot_hat, Ntime_hat, NHuggins_hat,
                Bayes_Nmed_M0, Bayes_Nmed_Mt, Bayes_Nmed_Mh, n) %>%
  transmute(., k = k, n = n,
            ML_M0 = Ndot_hat,
            ML_Mt = Ntime_hat,
            ML_Huggins = NHuggins_hat,
            Bayes_M0 = Bayes_Nmed_M0,
            Bayes_Mt = Bayes_Nmed_Mt,
            Bayes_Mh = Bayes_Nmed_Mh) %>%
  gather(., Model, Nhat, ML_M0:Bayes_Mh)

Nlow <- select(all.data, k, Ndot_lcl, Ntime_lcl, NHuggins_lcl,
               Bayes_Nlcl_M0, Bayes_Nlcl_Mt, Bayes_Nlcl_Mh, n) %>%
  transmute(., k = k, n = n,
            ML_M0 = Ndot_lcl,
            ML_Mt = Ntime_lcl,
            ML_Huggins = NHuggins_lcl,
            Bayes_M0 = Bayes_Nlcl_M0,
            Bayes_Mt = Bayes_Nlcl_Mt,
            Bayes_Mh = Bayes_Nlcl_Mh) %>%
  gather(., Model, Nlow, ML_M0:Bayes_Mh)

Nhigh <- select(all.data, k, Ndot_ucl, Ntime_ucl, NHuggins_ucl,
                Bayes_Nucl_M0, Bayes_Nucl_Mt, Bayes_Nucl_Mh, n) %>%
  transmute(., k = k, n = n,
            ML_M0 = Ndot_ucl,
            ML_Mt = Ntime_ucl,
            ML_Huggins = NHuggins_ucl,
            Bayes_M0 = Bayes_Nucl_M0,
            Bayes_Mt = Bayes_Nucl_Mt,
            Bayes_Mh = Bayes_Nucl_Mh) %>%
  gather(., Model, Nhigh, ML_M0:Bayes_Mh)

Nestims <- data.frame(k = Nhats$k,
                      n = Nhats$n,
                      model = Nhats$Model,
                      Nhat = Nhats$Nhat,
                      Nlow = Nlow$Nlow,
                      Nhigh = Nhigh$Nhigh) %>%
  filter(., Nhat < 180) %>%
  mutate(Nwidth = Nhigh - Nlow) %>%
  mutate(model_f = factor(model,
                          levels = c("ML_M0", "ML_Mt", "ML_Huggins",
                                     "Bayes_M0", "Bayes_Mt", "Bayes_Mh")))

Nestims_ML <- filter(Nestims, model_f == 'ML_M0' |
                       model_f == "ML_Mt" | model_f == "ML_Huggins")

Nestims_Bayes <- filter(Nestims, model_f == 'Bayes_M0' |
                       model_f == "Bayes_Mt" | model_f == "Bayes_Mh")

p1 <- ggplot(data = Nestims) +
  geom_boxplot(aes(x = k, y = Nhat, group = k)) +
  facet_wrap(~ model_f) +
  scale_x_continuous(name = 'Sample size',
                     breaks = c(6, 8, 10, 12, 14, 16)) +
  scale_y_continuous(name = 'Estimated abundance') +
  ggtitle(paste('Capture probability = UNIF(0.02, 0.1)')) +
  geom_hline(yintercept = 60)  +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

p2 <- ggplot(data = Nestims) +
  geom_boxplot(aes(x = k, y = Nwidth, group = k)) +
  facet_wrap(~ model_f, scales = "free_y") +
  scale_x_continuous(name = 'Sample size',
                     breaks = c(6, 8, 10, 12, 14, 16)) +
  scale_y_continuous(name = 'Width of 95% CI') +
  ggtitle(paste('Capture probability = UNIF(0.02, 0.1)')) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

p3 <- ggplot(data = Nestims) +
  geom_point(aes(x = n, y = Nhat,
                 size = Nwidth),
             shape = 1) +
  facet_wrap(~ model_f, scales = "free_y") +
  scale_x_continuous(name = 'Number of captured turtles') +
  scale_y_continuous(name = 'Estimated abundance') +
  scale_size_continuous(name = 'Width of 95% CI') +
  ggtitle(paste('Capture probability = UNIF(0.02, 0.1)')) +
  geom_hline(yintercept = 60)  +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

p4 <- ggplot(data = Nestims_ML) +
  geom_point(aes(x = n, y = Nwidth, color = model_f, size = k),
             alpha = 0.5) +
  scale_x_continuous(name = 'Number of captured turtles') +
  scale_y_continuous(name = 'Width of 95% CI') +
  scale_color_discrete(name = 'Model') +
  ggtitle(paste('Capture probability = UNIF(0.02, 0.1)')) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.8, 0.5))

p5 <- ggplot(data = Nestims_Bayes) +
  geom_point(aes(x = n, y = Nwidth, color = model_f, size = k),
             alpha = 0.5) +
  scale_x_continuous(name = 'Number of captured turtles') +
  scale_y_continuous(name = 'Width of 95% CI') +
  scale_color_discrete(name = 'Model') +
  ggtitle(paste('Capture probability = UNIF(0.02, 0.1)')) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
#        legend.position = c(0.8, 0.5))

if (save.plot){
  ggsave(plot = p1,
         filename = paste0('figures/Closed_simulation_pUnif_Nhat_',
                           Sys.Date(), '.png'),
         dpi = 600)

  ggsave(plot = p2,
         filename = paste0('figures/Closed_simulation_pUnif_Nwidth_',
                           Sys.Date(), '.png'),
         dpi = 600)

}


