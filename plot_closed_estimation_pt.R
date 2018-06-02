#plot_closed_simulation_estimation.R

# Plots results of closed_cmr_simulation_estimation.R

# Tomo Eguchi
# 30 January 2017


rm(list=ls())
source('Cm_SDB_functions.R')
library('ggplot2')
library('reshape2')

save.plot <- F

RData.files <- list.files(path = 'RData/RMark/',
                          pattern = '_pt')
#plot.list <- vector(mode = 'list',
  #                  length = length(RData.files))
pts <- c('pt1', 'pt2')
ps <- c('U(0.02, 0.1)', 'U(0.02, 0.06)')
models <- c('_withRMark', '_HugginsWithRMark')
ks.txt <- c('k_6', 'k_8', 'k_10', 'k_12', 'k_14', 'k_16')
k0 <- 1
# for each capture probability
for (k0 in 1:length(pts)){
  files.p <- RData.files[grep(pattern = pts[k0], RData.files)]
  k3 <- 1
  for (k3 in 1:length(models)){
    files.m <- files.p[grep(pattern = models[k3], files.p)]
    # for each sample size
    Nhats.all <- data.frame(k06 = numeric(0),
                            k08 = numeric(0),
                            k10 = numeric(0),
                            k12 = numeric(0),
                            k14 = numeric(0),
                            k16 = numeric(0))

    if (length(files.m) == 1){
      load(paste0('RData/RMark/', files.m))

      for (k1 in 1:length(ks.txt)){
        for (k2 in 1:length(RMark.out.all[[k1]])){
          # there are 500 estimates
          Nhats.all[k2, k1] <- RMark.out.all[[k1]][[k2]]$NHat

        }
      }
    } else {

      k1 <- 1
      for (k1 in 1:length(ks.txt)){
        file.k <- files.m[grep(pattern = ks.txt[k1], files.m)]
        load(paste0('RData/RMark/', file.k))
        for (k2 in 1:length(estims)){
          Nhats.all[k2, k1] <- estims[[k2]]$NHat
        }

      }
    }

    Nhats.all[Nhats.all > 500] <- NA
    Nhats.all$id <- seq(1, 500, by=1)
    tmp <- na.omit(melt(Nhats.all, id = 'id'))

    NAcounts <- colSums(is.na(Nhats.all))
    y.text <- 510

    # min.mins <- apply(min.n.all, FUN = min, MARGIN = 2)
    # y2.text <- 490
    #
    # max.maxs <- apply(max.n.all, FUN = max, MARGIN = 2)
    # y3.text <- 470
    #
    # mean.means <- apply(max.n.all, FUN = mean, MARGIN = 2)
    # y4.text <- 450
    #
    # strs <- strsplit(RData.files[[k0]], split = '.RData')
    p1 <- ggplot(data = tmp,
                 aes(x = variable, y = value)) +
      geom_boxplot() +
      scale_x_discrete(name = 'Sample size',
                       labels = c('6', '8', '10', '12', '14', '16')) +
      scale_y_continuous(name = 'Estimated abundance') +
      ggtitle(paste('Capture probability =', ps[k0])) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_hline(yintercept = 60) +
      annotate('text',
               x = c(0.7, 1:6),
               y = rep(y.text, times = 7),
               label = c('>500', NAcounts[1:6])) #+

    # annotate('text',
    #          x = c(0.7, 1:6),
    #          y = rep(y2.text, times = 7),
    #          label = c('min(n) = ', min.mins)) +
    # annotate('text',
    #          x = c(0.7, 1:6),
    #          y = rep(y3.text, times = 7),
    #          label = c('max(n) = ', max.maxs)) +
    # annotate('text',
    #          x = c(0.7, 1:6),
    #          y = rep(y4.text, times = 7),
    #          label = c('mean(n) = ', round(mean.means)))

    if (save.plot)
      ggsave(plot = p1,
             filename = paste0('figures/Closed_simulation_',
                               pts[k0], models[k3], '_',
                               Sys.Date(), '.png'),
             dpi = 600)
  }


  #plot.list[[k0]] <- p1

}
