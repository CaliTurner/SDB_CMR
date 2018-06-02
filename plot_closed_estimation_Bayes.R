#plot_closed_simulation_estimation.R

# Plots results of closed_cmr_simulation_estimation.R

# Tomo Eguchi
# 30 January 2017


rm(list=ls())
source('Cm_SDB_functions.R')
library('ggplot2')
library('reshape2')

RData.files <- list.files(path = 'RData/Bayes/',
                          pattern = 'withBayes')
#plot.list <- vector(mode = 'list',
  #                  length = length(RData.files))

ps.txt <- c('p2', 'p4', 'p6', 'p8', 'p10')
ks.txt <- c('k_6', 'k_8', 'k_10', 'k_12', 'k_14', 'k_16')
k0 <- 3
# for each capture probability
for (k0 in 1:length(ps.txt)){
  # Get all files for a particular p
  files.p <- RData.files[grep(pattern = ps.txt[k0], RData.files)]

  if (length(files.p) == 1){
    load(paste0('RData/Bayes/', files.p))

    # for each sample size
    Nhats.all <- data.frame(k06 = numeric(0),
                            k08 = numeric(0),
                            k10 = numeric(0),
                            k12 = numeric(0),
                            k14 = numeric(0),
                            k16 = numeric(0))
    for (k1 in 1:length(k)){
      for (k2 in 1:length(bayes.out.all[[k1]])){
        # there are 500 estimates
        Nhats.all[k2, k1] <- bayes.out.all[[k1]][[k2]]$summary$quantiles['N', '50%']

      }
    }

  } else {
    Nhats.all <- data.frame(k06 = numeric(0),
                            k08 = numeric(0),
                            k10 = numeric(0),
                            k12 = numeric(0),
                            k14 = numeric(0),
                            k16 = numeric(0))

    k1 <- 1
    for (k1 in 1:length(ks.txt)){
      file.k <- files.p[grep(pattern = ks.txt[k1], files.p)]
      load(paste0('RData/Bayes/', file.k))
      for (k2 in 1:length(bayes.out)){
        Nhats.all[k2, k1] <- bayes.out[[k2]]$summary$quantiles['N', '50%']
      }

    }

  }

  Nhats.all[Nhats.all > 500] <- NA
  Nhats.all$id <- seq(1, 500, by=1)
  tmp <- na.omit(melt(Nhats.all, id = 'id'))

  NAcounts <- colSums(is.na(Nhats.all))
  y.text <- 510
  p1 <- ggplot(data = tmp,
               aes(x = variable, y = value)) +
    geom_boxplot() +
    scale_x_discrete(name = 'Sample size',
                     labels = c('6', '8', '10', '12', '14', '16')) +
    scale_y_continuous(name = 'Estimated abundance') +
    ggtitle(paste('Capture probability = ', summary.data$p[1])) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 60) +
    annotate('text',
             x = c(0.7, 1:6),
             y = rep(y.text, times = 7),
             label = c('>500', NAcounts[1:6]))
  ggsave(plot = p1,
         filename = paste0('figures/closed_simulation_',
                           ps.txt[k0],
                           '_withRMark_2017-01-30.png'),
         dpi = 600)

  #plot.list[[k0]] <- p1, width = 8, height = 7,

}
