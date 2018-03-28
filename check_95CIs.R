#check_95CIs

# check suspicious lines in sample size vs 95%CIs by extracting
# results - see p4 in plot_closed_estimation_all_pUnif.R

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
                Bayes_Nmed_M0, Bayes_Nmed_Mt, Bayes_Nmed_Mh, n, rep) %>%
  transmute(., k = k, n = n,
            rep = rep,
            ML_M0 = Ndot_hat,
            ML_Mt = Ntime_hat,
            ML_Huggins = NHuggins_hat,
            Bayes_M0 = Bayes_Nmed_M0,
            Bayes_Mt = Bayes_Nmed_Mt,
            Bayes_Mh = Bayes_Nmed_Mh) %>%
  gather(., Model, Nhat, ML_M0:Bayes_Mh)

Nlow <- select(all.data, k, Ndot_lcl, Ntime_lcl, NHuggins_lcl,
               Bayes_Nlcl_M0, Bayes_Nlcl_Mt, Bayes_Nlcl_Mh, n, rep) %>%
  transmute(., k = k, n = n, rep = rep,
            ML_M0 = Ndot_lcl,
            ML_Mt = Ntime_lcl,
            ML_Huggins = NHuggins_lcl,
            Bayes_M0 = Bayes_Nlcl_M0,
            Bayes_Mt = Bayes_Nlcl_Mt,
            Bayes_Mh = Bayes_Nlcl_Mh) %>%
  gather(., Model, Nlow, ML_M0:Bayes_Mh)

Nhigh <- select(all.data, k, Ndot_ucl, Ntime_ucl, NHuggins_ucl,
                Bayes_Nucl_M0, Bayes_Nucl_Mt, Bayes_Nucl_Mh, n, rep) %>%
  transmute(., k = k, n = n, rep,
            ML_M0 = Ndot_ucl,
            ML_Mt = Ntime_ucl,
            ML_Huggins = NHuggins_ucl,
            Bayes_M0 = Bayes_Nucl_M0,
            Bayes_Mt = Bayes_Nucl_Mt,
            Bayes_Mh = Bayes_Nucl_Mh) %>%
  gather(., Model, Nhigh, ML_M0:Bayes_Mh)

Nestims <- data.frame(k = Nhats$k,
                      n = Nhats$n,
                      rep = Nhats$rep,
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

# problems seem to occur in lines... check one at a time:

# most extreme ones:
write.table(filter(Nestims_ML, n == 24 & Nwidth >350),
            file = 'data/series1_CI95.csv',
            sep = ",", quote = F, row.names = F,
            col.names = T, append = F)
write.table(filter(Nestims_ML, n == 23 & Nwidth >250),
            file = 'data/series1_CI95.csv',
            sep = ",", quote = F, row.names = F,
            col.names = F, append = T)
write.table(filter(Nestims_ML, n == 21 & Nwidth >190),
            file = 'data/series1_CI95.csv',
            sep = ",", quote = F, row.names = F,
            col.names = F, append = T)

write.table(filter(Nestims_ML, n == 20 & Nwidth >600),
            file = 'data/series1_CI95.csv',
            sep = ",", quote = F, row.names = F,
            col.names = F, append = T)
write.table(filter(Nestims_ML, n == 19 & Nwidth >600),
            file = 'data/series1_CI95.csv',
            sep = ",", quote = F, row.names = F,
            col.names = F, append = T)
write.table(filter(Nestims_ML, n == 18 & Nwidth >500),
            file = 'data/series1_CI95.csv',
            sep = ",", quote = F, row.names = F,
            col.names = F, append = T)
write.table(filter(Nestims_ML, n == 17 & Nwidth >400),
            file = 'data/series1_CI95.csv',
            sep = ",", quote = F, row.names = F,
            col.names = F, append = T)
write.table(filter(Nestims_ML, n == 16 & Nwidth >300),
            file = 'data/series1_CI95.csv',
            sep = ",", quote = F, row.names = F,
            col.names = F, append = T)
write.table(filter(Nestims_ML, n == 10 | n == 11),
            file = 'data/series1_CI95.csv',
            sep = ",", quote = F, row.names = F,
            col.names = F, append = T)


test.1 <- read.csv(file = 'data/series1_CI95.csv')
ns <- c(10, 11, 16, 17, 18, 19, 20, 21, 23)
ks0 <- seq(6, 16, by=2)
ks.idx <- 1:6
# go through one line at a time of the output and see what we got
# in raw data:
n_turtles.list <- vector(mode = 'list', length = nrow(test.1))
out.df <- data.frame(n=NA, k = NA, rep = NA, Ntime_hat = NA, Ntime_width = NA)
#i2 <- 1
c <- 1
for (i2 in 1:length(ns)){
  tmp1 <- filter(test.1, n == ns[i2])
  ks <- unique(tmp1$k)
  #reps <- unique(tmp1$rep)
 # i3 <- 1
  for (i3 in 1:length(ks)){
    tmp <- filter(all.data, k == ks[i3] & n == ns[i2])
  #  i1 <- 1
    for (i1 in 1:nrow(tmp)){
      n.turtles <- colSums(sim.results.data.list[[ks.idx[ks0 == ks[i3]]]][[tmp$rep[i1]]])
      n_turtles.list[[c]] <- n.turtles
      out.df[c,] <- c(n = ns[i2],
                      k = ks[i3],
                      rep = tmp$rep[i1],
                      Ntime_hat = tmp$Ntime_hat[i1],
                      Ntime_width = tmp$Ntime_ucl[i1] - tmp$Ntime_lcl[i1])
      c <- c + 1
    }

  }
}

#eliminate estimates that are too large:
out.df <- filter(out.df, Ntime_hat < 180)

#out.list <- out.list[unlist(lapply(out.list, length)) > 0]

# It appears that all estimates