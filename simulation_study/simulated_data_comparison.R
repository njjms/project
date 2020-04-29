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
gauscop_predictions <- vector(mode = "list",
                              length =sim.sets.n)
rfsp_predictions <- vector(mode = "list",
                          length = sim.sets.n)
ok_predictions <- vector(mode = "list",
                         length = sim.sets.n)

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
    
    # Create spatial random forest model
    while(is.null(rfsp_model)){
        tmp <- split_data(full_simdata, n = train.n)
        try(
            calculate_rfsp_predictions(training_data = tmp$train,
                                       test_points = tmp$test) -> rfsp_model
        )
    }
    
    # Create spatial Gaussian copula model
    transformed.training.set <- cube_root_fix_zeros(tmp$train)
    while(is.null(theta)) {
        try(
            theta <- calculate_spatial_params(simdata = transformed.training.set$transformed_data,
                                              testdata = tmp$test,
                                              zeros = transformed.training.set$zeros)
        )
    }
    krige_zhat_gauscop(resp = transformed.training.set$transformed_data$resp,
                       mu = theta$mu,
                       beta = theta$beta,
                       epsilon = transformed.training.set$epsilon,
                       Pi = transformed.training.set$Pi,
                       training_data = tmp$train,
                       test_data = tmp$test) -> zhat.zigs.df
    
    # calculate_zhat_gauscop(resp = transformed.training.set$transformed_data$resp,
    #                        mu = theta$mu,
    #                        beta = theta$beta,
    #                        epsilon = transformed.training.set$epsilon,
    #                        Pi = transformed.training.set$Pi,
    #                        Sigma_obs = theta$Sigma_obs,
    #                        Sigma_obs_noobs = theta$Sigma_obs_noobs) -> zhat.zigs.df
    
    yhat_gauscop <- backtransform_gauscop(zhat.zigs.df$zigs, transformed.training.set$epsilon)
    
    # Create ordinary kriging model
    training_data.cpy <- tmp$train
    test_data.cpy <- tmp$test
    coordinates(training_data.cpy) <- ~x+y
    coordinates(test_data.cpy) <- ~x+y
    
    ak_model <- autoKrige(resp ~ 1, training_data.cpy, test_data.cpy)
    
    # Finally, get RMSPE for all three models
    
    # rmspe(rfsp_model$predictions, tmp$test$resp)
    # rmspe(ak_model$krige_output$var1.pred, tmp$test$resp)
    # rmspe(yhat_gauscop, tmp$test$resp)
    
    obs_responses[[i]] <- tmp$test$resp
    gauscop_predictions[[i]] <- yhat_gauscop
    rfsp_predictions[[i]] <- rfsp_model$predictions
    ok_predictions[[i]] <- ak_model$krige_output$var1.pred
    
    print(paste0("Completed set: ", toString(i)))
    
}

Sys.time() - start

sim.final.df <- data.frame(
    obs = unlist(obs_responses),
    gauscop = unlist(gauscop_predictions),
    rfsp = unlist(rfsp_predictions),
    ok = unlist(ok_predictions)
)
saveRDS(sim.final.df,
        file = "sim_final_df.rds")

list(
    obs = obs_responses,
    gauscop = gauscop_predictions,
    rfsp = rfsp_predictions,
    ok = ok_predictions
) -> all.output.data
str(all.output.data)

saveRDS(all.output.data,
        file = "sim_final_lists.rds")

rmspe(sim.final.df$gauscop, sim.final.df$obs)
rmspe(sim.final.df$rfsp, sim.final.df$obs)
rmspe(sim.final.df$ok, sim.final.df$obs)

par(mfrow=c(2,2))
hist(sim.final.df$gauscop)
hist(sim.final.df$ok)
hist(sim.final.df$rfsp)
hist(sim.final.df$obs)
par(mfrow=c(1,1))
# 
sum(sim.final.df$gauscop == 0) / nrow(sim.final.df)
sum(sim.final.df$obs == 0) / nrow(sim.final.df)
