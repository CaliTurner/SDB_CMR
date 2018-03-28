#closed_cmr_estimation_RMark.R
#Using already created capture-recapture data, this runs RMark to estimate
# phats and Nhats.


rm(list=ls())

#data.dir<-"C:/Users/mike/Dropbox/teaching/workshop/13"
library(RMark)
#library(rjags)

source('Cm_SDB_functions.R')

# models to consider:
#POPAN (POPAN)
#Closed (Closed)  No explicit N parameter
#Jolly-Seber (Jolly)


p.dot <- list(formula = ~1)
# RData.files <- list.files(path = 'RData/',
#                           pattern = '*_2017-01-30.RData')
RData.files <- list.files(path = 'RData/',
                          pattern = '*_2017-01-30.RData')

k0 <- 1
for (k0 in 1:length(RData.files)){
  load(paste0('RData/', RData.files[k0]))
  #RMark.out.all <- vector(mode = 'list',
  #                        length = length(sim.results.all))
  k1 <- 1
  for (k1 in 1:length(sim.results.all)){   # different sample occassions
    k2 <- 1
    estims <- vector(mode = 'list', length = length(sim.results.all[[k1]]))
    for (k2 in 1:length(sim.results.all[[k1]])){
      print(paste('k1 = ', k1, '; k2 = ', k2))
      y <- sim.results.all[[k1]][[k2]]
      dp <- process.data(data = y$capt.hist,
                         model = 'Closed')
      ddl <- make.design.data(dp)

      #y2 <- y[rowSums(y$y.full) > 0, ]
      #capt.hist <- data.frame(ch = pasty(y2[,1:ncol(y2)]),
      #                        ind = 1)
      m0 <- mark(data = dp,
                 ddl = ddl,
                 model.parameters = list(p = p.dot),
                 silent = T,
                 output = F,
                 delete = T)
      p.hat <- get.real(m0, "p")
      N.hat <- get.real(m0, 'f0') + nrow(capt.hist)

      estims[[k2]] <- list(pHat = p.hat, NHat = N.hat)

    }
    #RMark.out.all[[k1]] <- estims
    save(list = c('summary.data', 'sample.size', 'N', 'p',
                  'sim.results.all', 'estims'),
         file = paste0('RData/RMark/', unlist(strsplit(RData.files[k0],
                                                       split = '2017-01-30'))[1],
                       'k_', k[k1],  '_withRMark_2017-01-30.RData'))

  }


}


