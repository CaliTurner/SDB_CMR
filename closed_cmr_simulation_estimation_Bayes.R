
# this code was obtained from: https://sites.google.com/site/workshoponcmr/home/sche/13--simulation-using-r-and-rmark


rm(list=ls())

#data.dir<-"C:/Users/mike/Dropbox/teaching/workshop/13"
#require(RMark)
library(rjags)

source('Cm_SDB_functions.R')
n.chains <- 5
n.adapt <- 1000
n.update <- 50000
n.iter <- 5000
nz <- 50 # the augmented additional rows of individuals

RData.files <- list.files(path = 'RData/',
                          pattern = '*2017-01-30.RData')

k <- 1
for (k in 1:length(RData.files)){
  load(paste0('RData/', RData.files[k]))
  bayes.out.all <- vector(mode = 'list',
                          length = length(sim.results.all))
  k1 <- 1
  for (k1 in 1:length(sim.results.all)){   # different sample occassions
    k2 <- 1
    bayes.out <- vector(mode = 'list',
                        length = length(sim.results.all[[k1]]))
    for (k2 in 1:length(sim.results.all[[k1]])){   # # simulations
      bayes.out[[k2]] <- estim_Bayes(sim.results.all[[k1]][[k2]],
                                     'models/Model_M0.txt',
                                     params = c("N", "p", "Omega", "deviance"),
                                     nz = nz,
                                     n.chains = n.chains,
                                     n.adapt = n.adapt,
                                     n.update = n.update,
                                     n.iter = n.iter)

    }
    bayes.out.all[[k1]] <- bayes.out

  }
  save(list = c('summary.data', 'sample.size', 'N', 'p', 'k',
                'sim.results.all', 'bayes.out.all'),
       file = paste0('RData/closed_simulation_p', p[j]*100,
                     '_withBayes_', Sys.Date(), '.RData'))

}





