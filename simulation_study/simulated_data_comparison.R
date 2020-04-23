setwd("~/MS_Project/simulation_study/")
simdata <- readRDS("../simulated_Y.rds")
dat <- read.csv("../ForestDataFuzzed.csv")

source("../ZIGFunctions.R")
source("preprocess.R")
source("copula_func.R")
source("kriging_func.R")
source("rfsp_func.R")

require(gstat)
require(tidyverse)
require(rgdal)
library(sp)

alb.xy <- project(with(dat,cbind(lon_fuzzed,lat_fuzzed)),
                  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
colnames(alb.xy) <- c("x","y")
dat <- cbind(alb.xy,dat)
full_simdata <- data.frame(
    x = dat$x,
    y = dat$y,
    resp = simdata[2,]
)

set.seed(182)
tmp <- split_data(full_simdata, n = 300)
transformed.training.set <- cube_root_fix_zeros(tmp$train)
# str(transformed.training.set)

theta <- calculate_spatial_params(simdata = transformed.training.set$transformed_data,
                                  testdata = tmp$test,
                                  zeros = transformed.training.set$zeros)
str(theta)

calculate_zhat_gauscop(resp = transformed.training.set$transformed_data$resp,
                       mu = theta$mu,
                       beta = theta$beta,
                       epsilon = transformed.training.set$epsilon,
                       Pi = transformed.training.set$Pi,
                       Sigma_obs = theta$Sigma_obs,
                       Sigma_obs_noobs = theta$Sigma_obs_noobs) -> zhat.zigs.df

hist(zhat.zigs.df$zigs)

yhat <- backtransform_gauscop(zhat.zigs.df$zigs, transformed.training.set$epsilon)
plot(tmp$test$resp, yhat)

# Compare this to ordinary kriging

ord_kriging(training_data = tmp$train,
            test_points = tmp$test) -> ok_model


plot(tmp$test$resp, ok_model$var1.pred,
     main = "Ordinary Kriging",
     xlab = "Actual Responses",
     ylab = "Predicted Responses")

par(mfrow=c(1,3))
hist(tmp$test$resp)
hist(yhat)
hist(ok_model$var1.pred)
par(mfrow=c(1,1))

# Compare this to RFSP

set.seed(1)
tmp1 <- split_data(full_simdata, n = 300)
training_data <- tmp1$train
test_points <- tmp1$test

calculate_rfsp_predictions(training_data = tmp1$train,
                           test_points = tmp1$test) -> rfsp_model
hist(rfsp_model$predictions)
plot(tmp1$test$resp, rfsp_model$predictions)


