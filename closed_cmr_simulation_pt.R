#closed_cmr_simulation_pt.R
# Simulation for closed Pt model

# Tomo Eguchi
# 1 February 2017


rm(list=ls())

source('Cm_SDB_functions.R')


N <- 60
#p <- c(0.02, 0.04, 0.06, 0.08, 0.10)
k <- c(6, 8, 10, 12, 14, 16)
n_sim_reps <- 500

sim.results <- vector(mode = 'list', length = n_sim_reps)

sim.results.all <- vector(mode = 'list', length = length(k))
i <- 1
for (i in 1:length(k)){
  #p <- runif(k[i], 0.02, 0.10)   # pt
  p <- runif(k[i], 0.02, 0.06)    # pt2
  for (k1 in 1:n_sim_reps){
    sim.results[[k1]] <- try(sim.data(N = N,
                                p = p,
                                k = k[i]),
                       silent = TRUE)
  }
  sim.results.all[[i]] <- sim.results

}

save(list = c('N', 'p', 'k', 'sim.results.all'),
     file = paste0('RData/closed_simulation_pt2_',
                   Sys.Date(), '.RData'))



