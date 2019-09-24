#Analysis-1


# Tomo Eguchi
# 1 September 2017

rm(list=ls())
library(RMark)
library(tidyverse)
library(lubridate)
#library(geosphere)

source('Cm_SDB_functions.R')
#dat0 <- read.csv(file = 'data/SDB_CMR_2019.csv')

# The input file for this was created by extract_data_TurtleDB_2019.R
# which access the turtle database, extracts data for the CMR study,
# saves into a .csv file.
dat.size <- read_csv(file = 'data/CMR2019_size_data_2019-09-23.csv')

CH <- select(dat.size, -c("ID", "CCL_cm"))
# extract the capture history portion:
capt.hist <- data.frame(ch = pasty(CH))

# effort data:  not entered into databaes as of 2019-09-23 so removing it
# bring in the net data - created by Katie O'Laughlin for the tide project
# dat.net <- read_csv(file = 'data/Standardized Net Data.csv') %>%
#   mutate(., Date = with_tz(as.Date(Field_Date, format = '%m/%d/%Y'),
#                            'america/los_angeles')) %>%
#   filter(., Date >= as.Date('2017-05-02') &
#            Date <= as.Date('2017-07-06')) %>%
#   mutate(., NetDist = distGeo(cbind(Net_Lon_End1, Net_Lat_End1),
#                               cbind(Net_Lon_End2, Net_Lat_End2)))
# 
# # missing values are filled in with the average of others
# dat.net[is.na(dat.net$NetDist), 'NetDist'] <- mean(dat.net$NetDist, na.rm = T)
# 
# # Compute the number of net deployment hours, multiply by lengths
# # to get the total effort:
# dat.net %>% mutate(., NetHrs = difftime(Net_Retrieval_Time,
#                                         Net_Deployment_Time,
#                                         units = 'hours')) %>%
#   mutate(., Effort = NetDist * as.numeric(NetHrs)) -> dat.net
# 
# # total effort per day then remove June 13 - no capture:
# dat.net.summary <- dat.net %>%
#   group_by(., Date) %>%
#   summarise(., TotalEffort = sum(Effort)) %>%
#   filter(., Date != '2017-06-13')
# 
# effort.1000 <- dat.net.summary$TotalEffort/1000

# CMR analyses start here. Closed models can't have individual covaraites

# set up models:
# share = TRUE makes p = c. p is the dominant parameter so that needs to be
# defined first.
time <- list(formula = ~ time, share=TRUE)
dot <- list(formula = ~ 1, share = TRUE)
#effort <- list(formula = ~ effort, share = TRUE)
#Time <- list(formula = ~ Time, share = TRUE)

# closed model
dp.closed <- process.data(data = capt.hist,
                          model = 'Closed')

ddl.closed <- make.design.data(dp.closed)
# ddl.closed$p$effort <- 0
# ddl.closed$p$effort <- effort.1000
# 
# ddl.closed$c$effort <- 0
# ddl.closed$c$effort <- effort.1000[2:length(effort.1000)]

# capture probability
p.dot <- dot
p.time <- time
#p.effort <- effort

f0.dot <- dot  # # never caught

models.closed <- create.model.list("Closed")

model.list.closed <- mark.wrapper(models.closed,
                                  data = dp.closed,
                                  ddl = ddl.closed,
                                  silent = F)

# do the median c-hat stuff here:
export.chdata(dp.closed, filename = 'data/closed_data', replace = T)

# Fletcher's c-hat looks okay.

# Fletcher's chat seems to be just about 1.0 so probably don't have to
# do anything here
med.chat <- 1.0  # adjust according to what comes out of the test
model.list.closed <- adjust.chat(med.chat, model.list.closed)

best.model.closed <- as.numeric(row.names(model.list.closed$model.table)[1])
best.closed.summary <- summary(model.list.closed[[best.model.closed]])

model.list.closed$model.table

# extract the Nhat and its SEs
Nhat.closed <- eval(parse(text = paste0('model.list.closed$',
                                 names(model.list.closed[best.model.closed]),
                                 '$results$derived')))

# Add heterogeneity: (not sure if this is useful - see Huggins' below)
# pi.dot <- list(formula = ~ 1)
# rm(list = c('p.effort', 'p.time'))
# dp.closed.het <- process.data(data = capt.hist,
#                           model = 'HetClosed')
#
# ddl.closed.het <- make.design.data(dp.closed.het)
# models.closed.het <- create.model.list("HetClosed")
# model.list.closed.het <- mark.wrapper(models.closed.het,
#                                       data = dp.closed.het,
#                                       ddl = ddl.closed.het,
#                                       silent = F)


# Huggins - should be able to incorporate individual covariates; CCL in this case
# may change them to size classes?
dp.Huggins <- process.data(data = capt.hist,
                           model = 'Huggins')

dp.Huggins$data$size <- dat.size$CCL_cm
dp.Huggins$data$size_cat1 <- ifelse(dat.size$CCL_cm > 90, 1, 0)
dp.Huggins$data$size_cat2 <- ifelse(dat.size$CCL_cm > 80, 1, 0)
dp.Huggins$data$size_cat3 <- ifelse(dat.size$CCL_cm > 70, 1, 0)

size <- list(formula = ~size, share = TRUE)
size_cat1 <- list(formula = ~ size_cat1, share = TRUE)
size_cat2 <- list(formula = ~ size_cat2, share = TRUE)
size_cat3 <- list(formula = ~ size_cat3, share = TRUE)

ddl.Huggins <- make.design.data(dp.Huggins)
#ddl.Huggins$p$effort <- 0
#ddl.Huggins$p$effort <- effort.1000

#ddl.Huggins$c$effort <- 0
#ddl.Huggins$c$effort <- effort.1000[2:length(effort.1000)]

# capture probability
p.dot <- dot
p.time <- time
#p.effort <- effort
p.size <- size
p.size_cat1 <- size_cat1
p.size_cat2 <- size_cat2
p.size_cat3 <- size_cat3

models.Huggins <- create.model.list("Huggins")

model.list.Huggins <- mark.wrapper(models.Huggins,
                                   data = dp.Huggins,
                                   ddl = ddl.Huggins,
                                   silent = F)

model.list.Huggins$model.table
best.model.Huggins <- as.numeric(row.names(model.list.Huggins$model.table)[1])
best.Huggins.summary <- summary(model.list.Huggins[[best.model.Huggins]])


Nhat.Huggins <- eval(parse(text = paste0('model.list.Huggins$',
                                         names(model.list.Huggins[best.model.Huggins]),
                                         '$results$derived')))


# rm(list = c('p.effort', 'p.size', 'p.size_cat1',
#             'p.size_cat2', 'p.size_cat3'))
# pi.dot <- list(formula = ~ 1)
# p.dot <- list(formula = ~ 1)
#
# models.Huggins.het <- create.model.list('HugHet')
# dp.Huggins.het <- process.data(data = capt.hist,
#                            model = 'HugHet')

# dp.Huggins.het$data$size <- dat0.size$CCL
# dp.Huggins.het$data$size_cat1 <- ifelse(dat0.size$CCL > 90, 1, 0)
# dp.Huggins.het$data$size_cat2 <- ifelse(dat0.size$CCL > 80, 1, 0)
# dp.Huggins.het$data$size_cat3 <- ifelse(dat0.size$CCL > 70, 1, 0)

# ddl.Huggins.het <- make.design.data(dp.Huggins.het)
# model.list.Huggins <- mark.wrapper(models.Huggins.het,
#                                    data = dp.Huggins.het,
#                                    ddl = ddl.Huggins.het,
#                                    silent = F)
#
#

save(list = ls(),
     file = paste0('RData/Mark_Analysis_output_',
                   Sys.Date(), '.RData'))
