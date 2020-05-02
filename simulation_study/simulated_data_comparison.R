setwd("~/MS_Project/simulation_study/")
simdata <- readRDS("../simulated_Y.rds")
dat <- read.csv("../ForestDataFuzzed.csv")

source("../ZIGFunctions.R")
source("preprocess.R")
source("copula_func.R")
source("kriging_func.R")
source("rfsp_func.R")
source("metrics.R")

library(gstat)
library(tidyverse)
library(rgdal)
library(sp)
library(automap)

alb.xy <- project(with(dat,cbind(lon_fuzzed,lat_fuzzed)),
                  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
colnames(alb.xy) <- c("x","y")
dat <- cbind(alb.xy,dat)

train.n <- 300
sim.sets.n <- 1000

obs_responses <- vector(mode = "list",
                        length = sim.sets.n)
gauscop_predictions_alg <- vector(mode = "list",
                              length =sim.sets.n)
gauscop_predictions_krige <- vector(mode = "list",
                              length =sim.sets.n)
rfsp_predictions <- vector(mode = "list",
                          length = sim.sets.n)
ok_predictions <- vector(mode = "list",
                         length = sim.sets.n)

gauscop_times <- vector(mode = "numeric", length = sim.sets.n)
rfsp_times <- vector(mode = "numeric", length = sim.sets.n)
ok_times <- vector(mode = "numeric", length = sim.sets.n)

epsilons <- vector(mode = "numeric", length = sim.sets.n)

start <- Sys.time()

for (i in 1:sim.sets.n) {
    
    set.seed(182)
    full_simdata <- data.frame(
        x = dat$x,
        y = dat$y,
        resp = simdata[i,]
    )
    
    theta <- NULL
    rfsp_model <- NULL
    
    tmp <- split_data(full_simdata, n = train.n)
    
    # Create spatial random forest model
    rfsp_start <- Sys.time()
    rfsp_model <- calculate_rfsp_predictions_rasterpd(training_data = tmp$train, test_points = tmp$test)
    rfsp_time <- Sys.time() - rfsp_start
    
    # while(is.null(rfsp_model)){
    #     tmp <- split_data(full_simdata, n = train.n)
    #     try(
    #         calculate_rfsp_predictions(training_data = tmp$train,
    #                                    test_points = tmp$test) -> rfsp_model
    #     )
    # }
    
    # Create spatial Gaussian copula model
    gauscop_start <- Sys.time()
    transformed.training.set <- cube_root_fix_zeros(tmp$train)
    epsilons[i] <- transformed.training.set$epsilon
    
    theta <- calculate_spatial_params(simdata = transformed.training.set$transformed_data,
                                      testdata = tmp$test,
                                      zeros = transformed.training.set$zeros)
    krige_zhat_gauscop(resp = transformed.training.set$transformed_data$resp,
                       mu = theta$mu,
                       beta = theta$beta,
                       epsilon = transformed.training.set$epsilon,
                       Pi = transformed.training.set$Pi,
                       training_data = tmp$train,
                       test_data = tmp$test) -> zhat.zigs.df.krige
    calculate_zhat_gauscop(resp = transformed.training.set$transformed_data$resp,
                           mu = theta$mu,
                           beta = theta$beta,
                           epsilon = transformed.training.set$epsilon,
                           Pi = transformed.training.set$Pi,
                           Sigma_obs = theta$Sigma_obs,
                           Sigma_obs_noobs = theta$Sigma_obs_noobs) -> zhat.zigs.df.alg
    
    yhat_gauscop_krige <- backtransform_gauscop(zhat.zigs.df.krige$zigs, transformed.training.set$epsilon)
    yhat_gauscop_alg <- backtransform_gauscop(zhat.zigs.df.alg$zigs, transformed.training.set$epsilon)
    gauscop_time <- Sys.time() - gauscop_start
    
    # Create ordinary kriging model
    ok_start <- Sys.time()
    training_data.cpy <- tmp$train
    test_data.cpy <- tmp$test
    coordinates(training_data.cpy) <- ~x+y
    coordinates(test_data.cpy) <- ~x+y
    
    ak_model <- autoKrige(resp ~ 1, training_data.cpy, test_data.cpy)
    ok_time <- Sys.time() - ok_start
    
    obs_responses[[i]] <- tmp$test$resp
    gauscop_predictions_alg[[i]] <- yhat_gauscop_alg
    gauscop_predictions_krige[[i]] <- yhat_gauscop_krige
    rfsp_predictions[[i]] <- rfsp_model$predictions
    ok_predictions[[i]] <- ak_model$krige_output$var1.pred
    
    gauscop_times[i] <- gauscop_time
    rfsp_times[i] <- rfsp_time
    ok_times[i] <- ok_time
    
    print(paste0("Completed set: ", toString(i)))
}

Sys.time() - start

sim.final.df <- data.frame(
    obs = unlist(obs_responses),
    gauscop_krige = unlist(gauscop_predictions_krige),
    gauscop_alg = unlist(gauscop_predictions_alg),
    rfsp = unlist(rfsp_predictions),
    ok = unlist(ok_predictions)
)
saveRDS(sim.final.df,
        file = "sim_final_df.rds")

# list(
#     obs = obs_responses,
#     gauscop = gauscop_predictions,
#     rfsp = rfsp_predictions,
#     ok = ok_predictions
# ) -> all.output.data
# 
# saveRDS(all.output.data,
#         file = "sim_final_lists.rds")

summary(gauscop_times)
summary(ok_times)
summary(rfsp_times)

# RMSPE of the three methods
rmspe(sim.final.df$gauscop_krige, sim.final.df$obs)
rmspe(sim.final.df$gauscop_alg, sim.final.df$obs)
rmspe(sim.final.df$rfsp, sim.final.df$obs)
rmspe(sim.final.df$ok, sim.final.df$obs)

# Histograms of predictions vs observed values
par(mfrow=c(2,2))
hist(sim.final.df$gauscop_krige,
     main = "Gaussian Copula (kriging) predictions",
     breaks = seq(0, 3000, 100),
     xaxt = 'n')
hist(sim.final.df$gauscop_alg[sim.final.df$gauscop_alg <= 3000],
     main = "Gaussian Copula (matrix algebra) predictions",
     breaks = seq(0, 3000, 100),
     xaxt = 'n')
hist(sim.final.df$ok,
     main = "Ordinary Kriging predictions",
     breaks = seq(min(sim.final.df$ok), 3000, 100),
     xaxt = 'n')
hist(sim.final.df$rfsp,
     main = "Random forest predictions",
     breaks = seq(0, 3000, 100),
     xaxt = 'n')
# hist(sim.final.df$obs,
#      main = "Observed values",
#      breaks = seq(0, max(sim.final.df$obs)+100, 100))
par(mfrow=c(1,1))

# How are zeros showing up in the predictions
sum(sim.final.df$gauscop == 0) / nrow(sim.final.df)
sum(sim.final.df$obs == 0) / nrow(sim.final.df)

# Let's do some residual plots
# par(mfrow=c(1,3))
# plot(x = sim.final.df$obs,
#      y = (sim.final.df$gauscop - sim.final.df$obs),
#      main = "Residual plot of Gaussian Copula",
#      xlab = "Observed Values",
#      ylab = "Residuals")
# plot(x = sim.final.df$obs,
#      y = (sim.final.df$ok - sim.final.df$obs),
#      main = "Residual plot of Ordinary Kriging",
#      xlab = "Observed Values",
#      ylab = "Residuals")
# plot(x = sim.final.df$obs,
#      y = (sim.final.df$rfsp - sim.final.df$obs),
#      main = "Residual plot of Random Forest",
#      xlab = "Observed Values",
#      ylab = "Residuals")
# par(mfrow=c(1,1))