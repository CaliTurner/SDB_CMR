#Analysis-1


# Tomo Eguchi
# 1 September 2017

rm(list=ls())
#library(RMark)
#library(rjags)
#library(R2jags)
library(tidyverse)
library(jagsUI)

source('Cm_SDB_functions.R')

n.chains <- 5
n.adapt <- 10000
n.update <- 10000
n.iter <- 50000
nz <- 50 # the augmented additional rows of individuals

#dat0 <- read.csv(file = 'data/SDB_CMR_2017.csv')
dat0 <- read.csv(file = "data/CMR2019_size_data_2019-09-23.csv")
y.full <- select(dat0, -c("ID", "CCL_cm"))

#y.full <- dat0[1:nrow(dat0), 2:ncol(dat0)]

bayes.out.M0 <- estim_Bayes(y.full,
                            'M0',
                            params = c("N", "p", "Omega", "deviance"),
                            nz = nz,
                            n.chains = n.chains,
                            n.adapt = n.adapt,
                            n.update = n.update,
                            n.iter = n.iter)

bayes.out.Mt <- estim_Bayes(y.full,
                            'Mt',
                            params = c("N", "p", "Omega", "deviance"),
                            nz = nz,
                            n.chains = n.chains,
                            n.adapt = n.adapt,
                            n.update = n.update,
                            n.iter = n.iter)

bayes.out.Mb <- estim_Bayes(y.full,
                            'Mb',
                            params = c("N", "p", "c", "Omega", "deviance"),
                            nz = nz,
                            n.chains = n.chains,
                            n.adapt = n.adapt,
                            n.update = n.update,
                            n.iter = n.iter)

# for individual heterogeneity models, more nz is needed to capture the
# tail of the distribution:
nz <- 100
bayes.out.Mh <- estim_Bayes(y.full,
                            'Mh',
                            params = c("N", "mean.p", "sd", "Omega", "deviance"),
                            nz = nz,
                            n.chains = n.chains,
                            n.adapt = n.adapt * 3,
                            n.update = n.update,
                            n.iter = n.iter * 2)

bayes.out.Mth <- estim_Bayes(y.full,
                             'Mth',
                             params = c("N", "mean.p", "mean.lp",
                                        "sd", "Omega", "deviance"),
                             nz = nz,
                             n.chains = n.chains,
                             n.adapt = n.adapt * 3,
                             n.update = n.update,
                             n.iter = n.iter * 2)

bayes.out.Mtbh <- estim_Bayes(y.full,
                              'Mtbh',
                              params = c("N", "mean.p", "gamma",
                                         "sd", "Omega", "deviance"),
                              nz = nz,
                              n.chains = n.chains,
                              n.adapt = n.adapt * 3,
                              n.update = n.update,
                              n.iter = n.iter * 2)

data.frame(models = c('M0', 'Mt', 'Mb',
                      'Mh', 'Mth', 'Mtbh'),
           DIC = c(bayes.out.M0$DIC, bayes.out.Mt$DIC,
                   bayes.out.Mb$DIC, bayes.out.Mh$DIC,
                   bayes.out.Mth$DIC, bayes.out.Mtbh$DIC),
           N = c(bayes.out.M0$summary$statistics['N', 'Mean'],
                 bayes.out.Mt$summary$statistics['N', 'Mean'],
                 bayes.out.Mb$summary$statistics['N', 'Mean'],
                 bayes.out.Mh$summary$statistics['N', 'Mean'],
                 bayes.out.Mth$summary$statistics['N', 'Mean'],
                 bayes.out.Mtbh$summary$statistics['N', 'Mean']),
           SE = c(bayes.out.M0$summary$statistics['N', 'SD'],
                  bayes.out.Mt$summary$statistics['N', 'SD'],
                  bayes.out.Mb$summary$statistics['N', 'SD'],
                  bayes.out.Mh$summary$statistics['N', 'SD'],
                  bayes.out.Mth$summary$statistics['N', 'SD'],
                  bayes.out.Mtbh$summary$statistics['N', 'SD']),
           N_lcl = c(bayes.out.M0$summary$quantiles['N', '2.5%'],
                     bayes.out.Mt$summary$quantiles['N', '2.5%'],
                     bayes.out.Mb$summary$quantiles['N', '2.5%'],
                     bayes.out.Mh$summary$quantiles['N', '2.5%'],
                     bayes.out.Mth$summary$quantiles['N', '2.5%'],
                     bayes.out.Mtbh$summary$quantiles['N', '2.5%']),
           N_ucl = c(bayes.out.M0$summary$quantiles['N', '97.5%'],
                     bayes.out.Mt$summary$quantiles['N', '97.5%'],
                     bayes.out.Mb$summary$quantiles['N', '97.5%'],
                     bayes.out.Mh$summary$quantiles['N', '97.5%'],
                     bayes.out.Mth$summary$quantiles['N', '97.5%'],
                     bayes.out.Mtbh$summary$quantiles['N', '97.5%'])) %>%
  arrange(., DIC) -> Bayes.results.table

save(list = ls(),
     file = paste0('RData/Bayes_Analysis_output_',
                   Sys.Date(), '.RData'))

