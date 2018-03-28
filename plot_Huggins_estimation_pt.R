#plot_closed_simulation_estimation.R

# Plots results of closed_cmr_simulation_estimation.R

# Tomo Eguchi
# 30 January 2017


rm(list=ls())
source('Cm_SDB_functions.R')
library('ggplot2')

RData.files <- list.files(path = 'RData/RMark/',
                          pattern = '*2017-02-01.RData')
#plot.list <- vector(mode = 'list',
  #                  length = length(RData.files))

k0 <- 1
# for each capture probability
for (k0 in 1:length(RData.files)){
  load(paste0('RData/RMark/', RData.files[k0]))
  # for each sample size
  Nhats.all <- min.n.all <- max.n.all <- mean.n.all <- data.frame(k06 = numeric(0),
                                        k08 = numeric(0),
                                        k10 = numeric(0),
                                        k12 = numeric(0),
                                        k14 = numeric(0),
                                        k16 = numeric(0))

  for (k1 in 1:length(k)){
    for (k2 in 1:length(RMark.out.all[[k1]])){
      # there are 500 estimates
      Nhats.all[k2, k1] <- RMark.out.all[[k1]][[k2]]$NHat
      min.n.all[k2, k1] <- min(sim.results.all[[k1]][[k2]]$sample.n)
      max.n.all[k2, k1] <- max(sim.results.all[[k1]][[k2]]$sample.n)
      mean.n.all[k2, k1] <- mean(sim.results.all[[k1]][[k2]]$sample.n)
    }
  }

  Nhats.all[Nhats.all > 500] <- NA
  Nhats.all$id <- seq(1, 500, by=1)
  tmp <- na.omit(melt(Nhats.all, id = 'id'))

  NAcounts <- colSums(is.na(Nhats.all))
  y.text <- 510

  min.mins <- apply(min.n.all, FUN = min, MARGIN = 2)
  y2.text <- 490

  max.maxs <- apply(max.n.all, FUN = max, MARGIN = 2)
  y3.text <- 470

  mean.means <- apply(max.n.all, FUN = mean, MARGIN = 2)
  y4.text <- 450

  strs <- strsplit(RData.files[[k0]], split = '.RData')
  p1 <- ggplot(data = tmp,
               aes(x = variable, y = value)) +
    geom_boxplot() +
    scale_x_discrete(name = 'Sample size',
                     labels = c('6', '8', '10', '12', '14', '16')) +
    scale_y_continuous(name = 'Estimated abundance') +
    ggtitle('Capture probability = U(0.02, 0.1)') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 60) +
    annotate('text',
             x = c(0.7, 1:6),
             y = rep(y.text, times = 7),
             label = c('>500', NAcounts[1:6])) +
    annotate('text',
             x = c(0.7, 1:6),
             y = rep(y2.text, times = 7),
             label = c('min(n) = ', min.mins)) +
    annotate('text',
             x = c(0.7, 1:6),
             y = rep(y3.text, times = 7),
             label = c('max(n) = ', max.maxs)) +
    annotate('text',
             x = c(0.7, 1:6),
             y = rep(y4.text, times = 7),
             label = c('mean(n) = ', round(mean.means)))

  # ggsave(plot = p1,
  #        filename = paste0('figures/', unlist(strs)[1], '.png'),
  #        width = 8, height = 7, dpi = 600)

  #plot.list[[k0]] <- p1

}
