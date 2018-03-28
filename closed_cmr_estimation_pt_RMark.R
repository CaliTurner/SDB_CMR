#closed_cmr_estimation_pt_RMark.R
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

#Share = T sets p = c
p.time <- list(formula = ~time, share=TRUE)
RData.files <- list.files(path = 'RData/',
                          pattern = '*pt2_2017-02-02.RData')

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

      m0 <- mark(data = dp,
                 ddl = ddl,
                 model.parameters = list(p = p.time),
                 silent = T,
                 output = F,
                 delete = T)
      p.hat <- get.real(m0, "p")
      N.hat <- get.real(m0, 'f0') + nrow(dp$data)

      estims[[k2]] <- list(pHat = p.hat, NHat = N.hat)

    }
    #RMark.out.all[[k1]] <- estims
    save(list = c('N', 'p', 'k', 'sim.results.all', 'estims'),
         file = paste0('RData/RMark/',
                       unlist(strsplit(RData.files[k0],
                                       split = '2017-02-02'))[1],
                       'k_', k[k1],
                       '_withRMark_', Sys.Date(), '.RData'))

  }


}


