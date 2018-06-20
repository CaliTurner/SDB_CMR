
# this code was obtained from: https://sites.google.com/site/workshoponcmr/home/sche/13--simulation-using-r-and-rmark


rm(list=ls())

#data.dir<-"C:/Users/mike/Dropbox/teaching/workshop/13"
require(RMark)
#library(rjags)

source('Cm_SDB_functions.R')
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
#p <- c(0.02, 0.04, 0.06, 0.08, 0.10)
k <- c(6, 8, 10, 12, 14, 16)
p <- runif(length(k), 0.02, 0.10)
n_sim_reps <- 500

summary.data <- data.frame(N = numeric(0),
                           p = numeric(0),
                           k = numeric(0),
                           cv = numeric(0),
                           Nhat = numeric(0),
                           N.sd = numeric(0))

sample.size <- vector(mode = 'list', length = length(k))
sim.results.all <- vector(mode = 'list', length = length(k))
for (i in 1:length(k)){

  # sample_sim() function simulates and analyzes. It's in
  # Cm_SDB_functions.R
  sim.results <- try(sample_sim(n_sim_reps = n_sim_reps,
                                N = N, p = p[i],
                                k = k[i]),
                     silent = TRUE)

  sim.results.all[[i]] <- sim.results$data.all

  new.data <- data.frame(N = N,
                         p = p[i],
                         k = k[i],
                         cv = sim.results$N.summary$N.cv,
                         Nhat = sim.results$N.summary$N.avg,
                         N.sd = sim.results$N.summary$N.sd)

  #append to output data
  summary.data <- rbind(summary.data, new.data)
  sample.size[[i]] <- sim.results$sample.size
}

save(list = c('summary.data', 'sample.size', 'N', 'p', 'k',
              'sim.results.all', 'sim.results.Nhats'),
     file = paste0('RData/closed_simulation_pt_',
                   Sys.Date(), '.RData'))



