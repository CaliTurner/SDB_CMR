
# this code was obtained from: https://sites.google.com/site/workshoponcmr/home/sche/13--simulation-using-r-and-rmark


rm(list=ls())

#data.dir<-"C:/Users/mike/Dropbox/teaching/workshop/13"
require(RMark)
library(rjags)
library(R2jags)

source('Cm_SDB_functions.R')
run.date <- Sys.Date()
# n.chains <- 5
# n.adapt <- 1000
# n.update <- 50000
# n.iter <- 5000

# models to consider:
#POPAN (POPAN)
#Closed (Closed)  No explicit N parameter
#Jolly-Seber (Jolly)

####FUNCTION TO SIMULATE CAPTURE HISTORIES UNDER SPECIFIED INPUTS

#simulate capture histories from assumed model and estimated parameter values
#simulate data under homogeneous p
#set up to simulate for specified inputs
#modelName <- "models/Model_M0.txt"

N <- 60
p <- 0.10
#k <- 6
k <- c(6, 8, 10, 12, 14, 16)
#p <- c(0.02, 0.04, 0.06, 0.08, 0.10) #runif(length(k), 0.02, 0.10)
n_sim_reps <- 100

n.chains <- 5
n.adapt <- 5000
n.update <- 5000 # update is not used December 2017
n.iter <- 50000
nz <- 100   # # augmented individuals

j <- i <- 1
for (j in 1:length(p)){

  summary.data <- data.frame(N = numeric(0),
                             p = numeric(0),
                             k = numeric(0),
                             N.cv = numeric(0),
                             N.hat = numeric(0),
                             N.sd = numeric(0),
                             model.names = character(0))

  sample.size <- vector(mode = 'list', length = length(k))
  sim.results.all <- vector(mode = 'list', length = length(k))
  sim.results.data <- vector(mode = 'list', length = length(k))

  for (i in 1:length(k)){
    print(paste('p[', j, '] of', length(p)))
    print(paste('k[', i, '] of', length(k)))

    # sample_sim() function simulates and analyzes. It's in
    # Cm_SDB_functions.R
    sim.results <- sample_sim(n_sim_reps = n_sim_reps,
                              N = N, p = p[j],
                              k = k[i], nz = nz,
                              n.chains = n.chains, n.adapt = n.adapt,
                              n.update = n.update, n.iter = n.iter,
                              parallel = TRUE)

    # estimates and uncertainties for each simulated dataset
    sim.results.all[[i]] <- sim.results$N.output

    # capture histories here:
    sim.results.data[[i]] <- sim.results$data.all

    # Just the averages over all estimates of simularions
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
               'N', 'p', 'k', 'sim.results.data',
               'sim.results.all'),
      file = paste0('RData/closed_simulation_p', p[j]*100,
                    '_', run.date, '.RData'))

}

