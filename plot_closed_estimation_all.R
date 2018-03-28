#plot_closed_simulation_estimation.R

# Plots results of closed_cmr_simulation_estimation.R

# Tomo Eguchi
# 30 January 2017


rm(list=ls())
source('Cm_SDB_functions.R')
library('ggplot2')
library('reshape2')
library(tidyverse)

save.plot <- T

RData.files <- c("RData/closed_simulation_p2_2017-12-12.RData",
                 "RData/closed_simulation_p4_2017-12-12.RData",
                 "RData/closed_simulation_p6_2017-12-12.RData",
                 "RData/closed_simulation_p8_2017-12-12.RData")

ps <- c(0.02, 0.04, 0.06, 0.08)
ks <- seq(from = 6, to = 16, by = 2)
k0 <- 1
# for each capture probability
for (k0 in 1:length(RData.files)){
  load(RData.files[k0])
  all.data <- do.call(rbind, sim.results.all)
  Nhats <- select(all.data, k, Ndot_hat, Ntime_hat, NHuggins_hat,
                  Bayes_Nmed_M0, Bayes_Nmed_Mt, Bayes_Nmed_Mh) %>%
    transmute(., k = k, ML_M0 = Ndot_hat, ML_Mt = Ntime_hat, ML_Huggins = NHuggins_hat,
              Bayes_M0 = Bayes_Nmed_M0, Bayes_Mt = Bayes_Nmed_Mt,
              Bayes_Mh = Bayes_Nmed_Mh) %>%
    gather(., Model, Nhat, ML_M0:Bayes_Mh)

  Nlow <- select(all.data, k, Ndot_lcl, Ntime_lcl, NHuggins_lcl,
                 Bayes_Nlcl_M0, Bayes_Nlcl_Mt, Bayes_Nlcl_Mh) %>%
    transmute(., k = k, ML_M0 = Ndot_lcl, ML_Mt = Ntime_lcl, ML_Huggins = NHuggins_lcl,
              Bayes_M0 = Bayes_Nlcl_M0, Bayes_Mt = Bayes_Nlcl_Mt,
              Bayes_Mh = Bayes_Nlcl_Mh) %>%
    gather(., Model, Nlow, ML_M0:Bayes_Mh)

  Nhigh <- select(all.data, k, Ndot_ucl, Ntime_ucl, NHuggins_ucl,
                 Bayes_Nucl_M0, Bayes_Nucl_Mt, Bayes_Nucl_Mh) %>%
    transmute(., k = k, ML_M0 = Ndot_ucl, ML_Mt = Ntime_ucl, ML_Huggins = NHuggins_ucl,
              Bayes_M0 = Bayes_Nucl_M0, Bayes_Mt = Bayes_Nucl_Mt,
              Bayes_Mh = Bayes_Nucl_Mh) %>%
    gather(., Model, Nhigh, ML_M0:Bayes_Mh)

  Nestims <- data.frame(k = Nhats$k,
                        model = Nhats$Model,
                        Nhat = Nhats$Nhat,
                        Nlow = Nlow$Nlow,
                        Nhigh = Nhigh$Nhigh) %>%
    filter(., Nhat < 200) %>%
    mutate(Nwidth = Nhigh - Nlow) %>%
    mutate(model_f = factor(model, levels = c("ML_M0", "ML_Mt", "ML_Huggins",
                                              "Bayes_M0", "Bayes_Mt", "Bayes_Mh")))

  # remove outliers:
  p1 <- ggplot(data = Nestims) +
    geom_boxplot(aes(x = k, y = Nhat, group = k)) +
    facet_wrap(~ model_f) +
    scale_x_continuous(name = 'Sample size',
                     breaks = c(6, 8, 10, 12, 14, 16)) +
    scale_y_continuous(name = 'Estimated abundance') +
    ggtitle(paste('Capture probability =', ps[k0])) +
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
    ggtitle(paste('Capture probability =', ps[k0])) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
    if (save.plot){
      ggsave(plot = p1,
             filename = paste0('figures/Closed_simulation_p',
                               ps[k0]*100, '_Nhat_',
                               Sys.Date(), '.png'),
             dpi = 600)

      ggsave(plot = p2,
             filename = paste0('figures/Closed_simulation_p',
                               ps[k0]*100, '_Nwidth_',
                               Sys.Date(), '.png'),
             dpi = 600)

    }
}

