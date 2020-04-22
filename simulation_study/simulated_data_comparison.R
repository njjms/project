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

tmp <- split_data(full_simdata, n = 300, cube.root.transform = TRUE)
transformed_tmp <- cube_root_fix_zeros(tmp$train)
tmp$train <- transformed_tmp$transformed_data
tmp$zeros <- transformed_tmp$zeros
tmp$epsilon <- transformed_tmp$epsilon
tmp$Pi <- transformed_tmp$Pi


hist(spt_tmp$Y)
hist(pzig(y = spt_tmp$Y, mu = spt_tmp$mu, beta = spt_tmp$beta, epsilon = spt_tmp$epsilon, Pi = spt_tmp$Pi))
hist(z)
qqnorm(z)
qqline(z)




training_data <- spt_tmp$train
test_points <- spt_tmp$test
hist(training_data$resp)
ok <- ord_kriging(training_data, test_points)
yhat <- backtransform(ok$var1.pred, spt_tmp$epsilon)
hist(yhat)
plot(spt_tmp$test$resp, yhat)

training_data$z <- z
ok_model <- krige(z ~ 1, training_data, test_points, model = z.fit)
yhat <- ok_model$var1.pred
hist(yhat)
hist(pz)
pz <- pnorm(yhat)
pz[which(pz==1)] <- 1-.Machine$double.eps
pz[which(pz==0)] <- .Machine$double.eps

zigs <- qzig(u = pz, mu = spt_tmp$mu, beta = spt_tmp$beta, epsilon = spt_tmp$epsilon, Pi = spt_tmp$Pi[1])  
hist(zigs)
hist(zigs)

plot(spt_tmp$test$resp, zigs)

class(full_simdata)
head(full_simdata)
coordinates(full_simdata) <- ~x+y
lzn.vgm <- variogram(resp ~ 1, full_simdata)
lzn.fit <- fit.variogram(lzn.vgm, vgm("Gau"))
plot(lzn.vgm, lzn.fit)

