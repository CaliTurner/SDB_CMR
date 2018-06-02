# closed_cmr_estimation_Bayes.R

# runs the augmented CMR analysis to estimate abundance
# from simualted data - see closed_cmr_simulation_estimation.R


# Tomo Eguchi
# 30 January 2017


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

k0 <- 1
for (k0 in 2:length(RData.files)){
  load(paste0('RData/', RData.files[k0]))
  bayes.out.all <- vector(mode = 'list',
                          length = length(sim.results.all))
  k1 <- 1
  for (k1 in 1:length(sim.results.all)){   # different sample occassions
    k2 <- 1
    bayes.out <- vector(mode = 'list',
                        length = length(sim.results.all[[k1]]))
    for (k2 in 1:length(sim.results.all[[k1]])){   # # simulations
    	print(paste0('k1 = ', k1, '; k2 = ', k2))
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
       file = paste0('RData/Bayes/', unlist(strsplit(RData.files[k0],
                                               split = '2017-01-30'))[1],
                     'withBayes_', Sys.Date(), '.RData'))

}





