#Cm_SCB_functions


# 27 January 2017
# Tomo Eguchi

sysInfo <- Sys.info()
ifelse(sysInfo[1] == 'Linux',
       source('~/Documents/R/TomosFunctions.R'),
       source('~/R/TomosFunctions.R'))

##FIRST DEFINE SOME FUNCTIONS NEEDED
## FUNCTION TO CREATE CAPTURE HISTORY CHARACTER STRINGS
#https://sites.google.com/site/cmrsoftware/lecture-lab-schedule/3-intro-to-mark-rmark/formatting-data-for-mark-and-rmark
pasty<-function(x)
{
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {
    out[i]<-paste(x[i,],collapse="")
  }
  return(out)
}

#Simulating data over n_sim_reps (ideally 1000 times)
sim.data<-function(N, p, k){
  #simulate capture histories
  y <- array(dim=c(N,k))
  ind <- array(dim=N)

  for(i in 1:N) {
    y[i, 1:k] <- rbinom(k, 1, p)
    ind[i] <- sum(y[i,]) > 0
  }
  #capture history data frame
  capt.hist <- data.frame(ch = pasty(y[,1:k]),
                          ind = ind)
  capt.hist <- subset(capt.hist, ind==T, select=c(ch))

  sample.n <- colSums(y)
  #end of function

  out <- list(capt.hist = capt.hist,
              sample.n = sample.n,
              y.full = y)
  return(out)
}

sample_sim <- function(n_sim_reps,N,p,k,nz = 50,
                       n.chains = 5, n.adapt = 1000,
                       n.update = 5000, n.iter = 10000) {
  # WE CAN DEFINE SOME THINGS OUTSIDE THE LOOP FOR EXAMPLE
  # WE WILL USE THE SAME TIME CONSTANT MODEL EACH SIMU;LLATION
  # Define parameters
  p.dot <- list(formula = ~1, share = TRUE)
  p.time <- list(formula = ~time, share = TRUE)

  #set up an empty data frame to store all the simulation results
  output.data <- data.frame(k=numeric(0),
                            rep=numeric(0),
                            Ndot_hat=numeric(0),
                            Ntime_hat = numeric(0),
                            Bayes_Nmed_M0 = numeric(0),
                            Bayes_Nmed_Mt = numeric(0),
                            Bayes_Nmed_Mh = numeric(0))

  sample.n.out <- matrix(nrow = n_sim_reps,
                         ncol = k)

  data.out <- bayes.diag <- bayes.DIC <- bayes.model <- vector(mode = 'list',
                                                               length = n_sim_reps)
  r <- 1   # this helps when running line by line
  for(r in 1:n_sim_reps){
    print(paste('Simulation', r, 'of', n_sim_reps))
    #function to create capture histories from simulated data
    sim.out <- sim.data(N=N, p=p, k=k)
    data.out[[r]] <- sim.out$y.full

    #additional parameters to supress output each simulation
    # note that "Closed" with get.real(m0, "N") as shown in this example
    # does not work because Closed models have no N parameter.
    M.dot <- mark(sim.out$capt.hist,
                  model = "Closed",
                  model.parameters = list(p = p.dot),
                  silent = T,
                  output = F,
                  delete = T)

    M.t <- mark(sim.out$capt.hist,
                model = "Closed",
                model.parameters = list(p = p.time),
                silent = T,
                output = F,
                delete = T)

    dp.Huggins <- process.data(data = sim.out$capt.hist,
                               model = 'Huggins')
    ddl.Huggins <- make.design.data(dp.Huggins)

    Mt.Huggins <- mark(data = dp.Huggins,
                       ddl = ddl.Huggins,
                       model.parameters = list(p = p.time),
                       silent = T,
                       output = F,
                       delete = T)

    #pull off just the estimate of N (you can select other parameters if you want)
    p.dot.hat <- get.real(M.dot, "p")
    p.dot.se <- M.dot[["results"]][["real"]][["se"]][1]
    p.dot.lcl<- M.dot[["results"]][["real"]][["lcl"]][1]
    p.dot.ucl<- M.dot[["results"]][["real"]][["ucl"]][1]

    N.dot.hat <- M.dot[["results"]][["derived"]][["N Population Size"]][["estimate"]]
    N.dot.se <- M.dot[["results"]][["derived"]][["N Population Size"]][["se"]]
    N.dot.lcl <- M.dot[["results"]][["derived"]][["N Population Size"]][["lcl"]]
    N.dot.ucl <- M.dot[["results"]][["derived"]][["N Population Size"]][["ucl"]]

    p.time.hat <- get.real(M.t, "p")
    p.time.se <- M.t[["results"]][["real"]][["se"]][1:k]
    p.time.lcl <- M.t[["results"]][["real"]][["lcl"]][1:k]
    p.time.ucl <- M.t[["results"]][["real"]][["ucl"]][1:k]

    N.time.hat <- M.t[["results"]][["derived"]][["N Population Size"]][["estimate"]]
    N.time.se <-  M.t[["results"]][["derived"]][["N Population Size"]][["se"]]
    N.time.lcl <-  M.t[["results"]][["derived"]][["N Population Size"]][["lcl"]]
    N.time.ucl <-  M.t[["results"]][["derived"]][["N Population Size"]][["ucl"]]

    p.Huggins.hat <- get.real(Mt.Huggins, "p")
    p.Huggins.se <- Mt.Huggins[["results"]][["real"]][["se"]]
    p.Huggins.lcl <- Mt.Huggins[["results"]][["real"]][["lcl"]]
    p.Huggins.ucl <- Mt.Huggins[["results"]][["real"]][["ucl"]]

    N.Huggins.hat <- Mt.Huggins[["results"]][["derived"]][["N Population Size"]][["estimate"]]
    N.Huggins.se <- Mt.Huggins[["results"]][["derived"]][["N Population Size"]][["se"]]
    N.Huggins.lcl <- Mt.Huggins[["results"]][["derived"]][["N Population Size"]][["lcl"]]
    N.Huggins.ucl <- Mt.Huggins[["results"]][["derived"]][["N Population Size"]][["ucl"]]

    # Run Bayesian analysis on the same data:
    bayes.out.Mt <- estim_Bayes(sim.out$y.full,
                                'Mt',
                                params = c("N", "p", "Omega", "deviance"),
                                nz = nz,
                                n.chains = n.chains,
                                n.adapt = n.adapt,
                                n.update = n.update,
                                n.iter = n.iter)

    bayes.out.M0 <- estim_Bayes(sim.out$y.full,
                                'M0',
                                params = c("N", "p", "Omega", "deviance"),
                                nz = nz,
                                n.chains = n.chains,
                                n.adapt = n.adapt,
                                n.update = n.update,
                                n.iter = n.iter)

    bayes.out.Mh <- estim_Bayes(sim.out$y.full,
                                'Mh',
                                params = c("N", "mean.p", "sd", "Omega", "deviance"),
                                nz = 100,
                                n.chains = n.chains,
                                n.adapt = n.adapt * 3,
                                n.update = n.update,
                                n.iter = n.iter * 2)

    bayes.diag[[r]] <- c(bayes.out.M0$diag, bayes.out.Mt$diag, bayes.out.Mh$diag)
    bayes.DIC[[r]] <- c(bayes.out.M0$DIC, bayes.out.Mt$DIC, bayes.out.Mh$DIC)
    bayes.model[[r]] <- c(bayes.out.M0$model, bayes.out.Mt$model, bayes.out.Mh$model)

    #put in data frame
    new.data <- data.frame(k = k,
                           rep = r,
                           Ndot_hat = N.dot.hat,
                           Ndot_se = N.dot.se,
                           Ndot_lcl = N.dot.lcl,
                           Ndot_ucl = N.dot.ucl,
                           Ntime_hat = N.time.hat,
                           Ntime_se = N.time.se,
                           Ntime_lcl = N.time.lcl,
                           Ntime_ucl = N.time.ucl,
                           NHuggins_hat = N.Huggins.hat,
                           NHuggins_se = N.Huggins.se,
                           NHuggins_lcl = N.Huggins.lcl,
                           NHuggins_ucl = N.Huggins.ucl,
                           Bayes_Nmed_M0 = bayes.out.M0$summary$quantiles['N', '50%'],
                           Bayes_Nse_M0 = bayes.out.M0$summary$statistics['N', 'SD'],
                           Bayes_Nlcl_M0 = bayes.out.M0$summary$quantiles['N', '2.5%'],
                           Bayes_Nucl_M0 = bayes.out.M0$summary$quantiles['N', '97.5%'],
                           Bayes_Nmed_Mt = bayes.out.Mt$summary$quantiles['N', '50%'],
                           Bayes_Nse_Mt = bayes.out.Mt$summary$statistics['N', 'SD'],
                           Bayes_Nlcl_Mt = bayes.out.Mt$summary$quantiles['N', '2.5%'],
                           Bayes_Nucl_Mt = bayes.out.Mt$summary$quantiles['N', '97.5%'],
                           Bayes_Nmed_Mh = bayes.out.Mh$summary$quantiles['N', '50%'],
                           Bayes_Nse_Mh = bayes.out.Mh$summary$statistics['N', 'SD'],
                           Bayes_Nlcl_Mh = bayes.out.Mh$summary$quantiles['N', '2.5%'],
                           Bayes_Nucl_Mh = bayes.out.Mh$summary$quantiles['N', '97.5%'])

    #append to output data
    output.data<-rbind(output.data, new.data)

    sample.n.out[r,] <- sim.out$sample.n
  }
  #create  summay statistics
  #empirical mean, sd, and cv
  N.avg <- apply(output.data[, c('Ndot_hat',
                                 'Ntime_hat',
                                 'NHuggins_hat',
                                 'Bayes_Nmed_M0',
                                 'Bayes_Nmed_Mt',
                                 'Bayes_Nmed_Mh')],
                 MARGIN = 2, FUN = mean)
  N.sd <- apply(output.data[, c('Ndot_hat',
                                'Ntime_hat',
                                'NHuggins_hat',
                                'Bayes_Nmed_M0',
                                'Bayes_Nmed_Mt',
                                'Bayes_Nmed_Mh')],
                MARGIN = 2, FUN = sd)

  N.cv<-N.sd/N.avg
  #return simulated data and summary stats in list objects
  N.summary<-list(N.avg = N.avg,
                  N.sd = N.sd,
                  N.cv = N.cv)

  output <- list(N.summary=N.summary,
                 N.output=output.data,
                 sample.size = sample.n.out,
                 data.all = data.out,
                 bayes.diag = bayes.diag,
                 bayes.DIC = bayes.DIC,
                 bayes.model = bayes.model)
  return(output)
}

estim_RMark_Closed <- function(cap.hist, model.parameters){
  m0 <- mark(capt.hist,
               model = 'Closed',
               model.parameters = model.parameters,
               silent = T,
               output = F,
               delete = T)
  #pull off just the estimate of N (you can select other parameters if you want)
  p.hat <- get.real(m0, "p")
  N.hat <- get.real(m0, 'f0') + nrow(capt.hist)

  output <- list(pHat = p.hat, NHat = N.hat)
  return(output)
}

estim_Bayes <- function(y.full, modelName,
                        params = c("N", "p", "Omega", "deviance"),
                        nz = 50,
                        n.chains = 5, n.adapt = 1000,
                        n.update = 5000, n.iter = 10000){
  # also use Bayesian model:
  y.data <- y.full[rowSums(y.full) > 0,]
  yobs <- as.matrix(y.data)
  yaug <- rbind(yobs,
                array(0, dim = c(nz, dim(yobs)[2])))

  bugs.data <- list(yaug = yaug,
                    M = nrow(yaug),
                    T = ncol(yaug))

  if (modelName == 'M0'){
    inits <- function() list(z = rep(1, nrow(yaug)),
                             p = runif(1, 0, 1))

  } else if (modelName == 'Mt'){
    inits <- function() list(z = rep(1, nrow(yaug)),
                             p = runif(ncol(y.full), 0, 1))
  } else if (modelName == 'Mb'){
    inits <- function() list(z = rep(1, nrow(yaug)),
                             p = runif(1, 0, 1))
  } else if (modelName == 'Mh'){
    # for model Mh, capture histories are converted into capture frequencies:
    y <- sort(apply(y.data, 1, sum), decreasing = TRUE)
    yaug <- c(y, rep(0, nz))
    bugs.data <- list(yaug = yaug,
                      M = length(yaug),
                      T = ncol(y.data))
    inits <- function() list(z = rep(1, length(yaug)),
                             mean.p = runif(1, 0, 1))

  } else if (modelName == 'Mth'){
    inits <- function() list(z = rep(1, nrow(yaug)),
                             mean.p = runif(ncol(yaug), 0, 1),
                             sd = runif(1, 0.1, 0.9))
  } else if (modelName == 'Mtbh'){
    inits <- function() list(z = rep(1, nrow(yaug)),
                             mean.p = runif(ncol(yaug), 0, 1),
                             sd = runif(1, 0.1, 0.9))

  }

  model.file <- paste0('models/Model_', modelName, '.txt')

  jm <- jags.model(model.file,
                   data = bugs.data,
                   inits,
                   n.chains = n.chains,
                   n.adapt = n.adapt)

  #update(jm, n.iter = n.update)

  load.module("dic")
  zm <- coda.samples(jm,
                     variable.names = params,
                     n.iter = n.iter)

  g.diag <- gelman.diag(zm)
  h.diag <- heidel.diag(zm)
  r.diag <- raftery.diag(zm)

  # dicOut <- dic.samples(jm,
  #                       n.iter = n.iter,
  #                       type = "pD")
  # dicOut2 <- dic.samples(jm,
  #                        n.iter = n.iter,
  #                        type = "popt")

  sum_zm <- summary(zm)
  meanDev <- sum_zm$statistics[rownames(sum_zm$statistics) == "deviance",
                               colnames(sum_zm$statistics) == "Mean"]
  sdDev <- sum_zm$statistics[rownames(sum_zm$statistics) == "deviance",
                             colnames(sum_zm$statistics) == "SD"]
  DIC <- 0.5*(sdDev^2) + meanDev
  out <- list(diag = list(g.diag = g.diag,
                          h.diag = h.diag,
                          r.diag = h.diag),
              summary = sum_zm,
              DIC = DIC,
              model = model.file,
              sample = runjags::combine.mcmc(zm))
  return(out)

}

estim_Bayes_parallel <- function(y.full,
                                 modelName,
                                 params = c("N", "p", "Omega", "deviance"),
                                 nz = 50,
                                 n.chains = 5,
                                 n.adapt = 1000,
                                 n.update = 5000,
                                 n.iter = 10000){
  # also use Bayesian model:
  y.data <- y.full[rowSums(y.full) > 0,]
  yobs <- as.matrix(y.data)
  yaug <- rbind(yobs,
                array(0, dim = c(nz, dim(yobs)[2])))

  bugs.data <- list(yaug = yaug,
                    M = nrow(yaug),
                    T = ncol(yaug))

  if (modelName == 'M0'){
    inits <- function(yaug){
      A <- list(z = rep(1, nrow(yaug)),
                p = runif(1, 0, 1))
      return(A)}

  } else if (modelName == 'Mt'){
    inits <- function(yaug){
      A <- list(z = rep(1, nrow(yaug)),
                p = runif(ncol(yaug), 0, 1),
                .RNG.name = "lecuyer::RngStream")
      return(A)}
  } else if (modelName == 'Mb'){
    inits <- function(yaug) list(z = rep(1, nrow(yaug)),
                             p = runif(1, 0, 1))
  } else if (modelName == 'Mh'){
    # for model Mh, capture histories are converted into capture frequencies:
    y <- sort(apply(y.data, 1, sum), decreasing = TRUE)
    yaug <- c(y, rep(0, nz))
    bugs.data <- list(yaug = yaug,
                      M = length(yaug),
                      T = ncol(y.data))
    inits <- function(yaug) list(z = rep(1, length(yaug)),
                             mean.p = runif(1, 0, 1))

  } else if (modelName == 'Mth'){
    inits <- function(yaug) list(z = rep(1, nrow(yaug)),
                             mean.p = runif(ncol(yaug), 0, 1),
                             sd = runif(1, 0.1, 0.9))
  } else if (modelName == 'Mtbh'){
    inits <- function(yaug) list(z = rep(1, nrow(yaug)),
                                 mean.p = runif(ncol(yaug), 0, 1),
                                 sd = runif(1, 0.1, 0.9))

  }

  model.file <- paste0('models/Model_', modelName, '.txt')

  # jm <- jags.model(model.file,
  #                  data = bugs.data,
  #                  inits,
  #                  n.chains = n.chains,
  #                  n.adapt = n.adapt)
  #
  # #update(jm, n.iter = n.update)
  #
  # load.module("dic")
  # zm <- coda.samples(jm,
  #                    variable.names = params,
  #                    n.iter = n.iter)
  jags.fit <- jags.parallel(data = bugs.data,
                            parameters.to.save = params,
                            model.file = model.file,
                            n.burnin = 10000,
                            n.chains = n.chains,
                            n.iter = n.iter,
                            n.thin = 1,
                            jags.module = c('dic', 'glm', 'lecuyer'))

  zm <- as.mcmc(jags.fit)

  g.diag <- gelman.diag(zm)
  h.diag <- heidel.diag(zm)
  r.diag <- raftery.diag(zm)

  # dicOut <- dic.samples(jm,
  #                       n.iter = n.iter,
  #                       type = "pD")
  # dicOut2 <- dic.samples(jm,
  #                        n.iter = n.iter,
  #                        type = "popt")

  sum_zm <- summary(zm)
  meanDev <- sum_zm$statistics[rownames(sum_zm$statistics) == "deviance",
                               colnames(sum_zm$statistics) == "Mean"]
  sdDev <- sum_zm$statistics[rownames(sum_zm$statistics) == "deviance",
                             colnames(sum_zm$statistics) == "SD"]
  DIC <- 0.5*(sdDev^2) + meanDev
  out <- list(diag = list(g.diag = g.diag,
                          h.diag = h.diag,
                          r.diag = h.diag),
              summary = sum_zm,
              DIC = DIC,
              model = model.file,
              sample = runjags::combine.mcmc(zm))
  return(out)

}
