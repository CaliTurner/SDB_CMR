
# this code was obtained from: https://sites.google.com/site/workshoponcmr/home/sche/13--simulation-using-r-and-rmark


rm(list=ls())

#data.dir<-"C:/Users/mike/Dropbox/teaching/workshop/13"
require(RMark)
library(rjags)

source('Cm_SDB_functions.R')

run.date <- Sys.Date()

N <- 60
p <- c(0.02, 0.04, 0.06, 0.08, 0.10)
k <- c(6, 8, 10, 12, 14, 16)
#p <- runif(length(k), 0.02, 0.10)
n_sim_reps <- 5

n.chains <- 5
n.adapt <- 5000
n.update <- 5000
n.iter <- 50000
nz <- 100   # # augmented individuals

j <- 1
for (j in 1:1){
  summary.data <- data.frame(N = numeric(0),
                             p = numeric(0),
                             k = numeric(0),
                             N.cv = numeric(0),
                             N.hat = numeric(0),
                             N.sd = numeric(0),
                             model.names = character(0))

  sample.size <- vector(mode = 'list', length = length(k))
  sim.results.all <- vector(mode = 'list', length = length(k))
  i <- 1
  for (i in 1:1){

    # sample_sim() function simulates and analyzes. It's in
    # Cm_SDB_functions.R
    sim.results <- try(sample_sim(n_sim_reps = n_sim_reps,
                                  N = N, p = p[j],
                                  k = k[i], nz = nz,
                                  n.chains = n.chains, n.adapt = n.adapt,
                                  n.update = n.update, n.iter = n.iter),
                       silent = TRUE)

    sim.results.all[[i]] <- sim.results$data.all

    new.data <- data.frame(N = N,
                           p = p[j],
                           k = k[i],
                           N.cv = sim.results$N.summary$N.cv,
                           N.hat = sim.results$N.summary$N.avg,
                           N.sd = sim.results$N.summary$N.sd,
                           model.names = names(sim.results$N.summary$N.avg))



    #append to output data
    summary.data <- rbind(summary.data, new.data)
    sample.size[[i]] <- sim.results$sample.size
  }

  save(list = c('summary.data', 'sample.size',
                'N', 'p', 'k',
                'sim.results.all'),
       file = paste0('RData/closed_simulation_p', p[j]*100,
                     '_', run.date, '.RData'))

}




