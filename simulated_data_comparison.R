setwd("~/MS_Project/code")
simdata <- readRDS("simulated_Y.rds")
dat <- read.csv("ForestDataFuzzed.csv")

source("ZIGFunctions.R")
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

split_data <- function(full_data, n = 300, cube.root.transform = TRUE, correct.zeros = TRUE) {
  
  colnames(full_data) <- c("x", "y", "resp")
  training_idx <- sample(1:nrow(full_data), n, replace=FALSE)
  training_data <- full_data[training_idx,]
  test_data <- full_data[-training_idx,]
  output <- list(train = training_data, 
                 test = test_data)
  return(output)
}

cube_root_fix_zeros <- function(data, cube.root.transform = TRUE, correct.zeros = TRUE) {
  if (!("resp" %in% colnames(data))) {
    stop("Needs response column in dataset")
  }
  
  if (cube.root.transform) {
    data$resp <- data$resp^(1/3)
  }
  zeros <- which(data$resp == 0)
  epsilon <- min(data$resp[-zeros])
  Pi <- rep(length(zeros)/n, n)
  
  if (correct.zeros) {
    data$resp[zeros] <- runif(length(zeros), 0, epsilon)
  }
  
  output <- list(transformed_data = data,
                 zeros = zeros,
                 epsilon = epsilon,
                 Pi = Pi)
  return(output)
}

calculate_spatial_params <- function(full_simdata, n = 300, cube.root.transform = TRUE, plot.var = FALSE) {
    #' Calculates Method of Moments estimates for Gamma model and the spatial covariance parameters for ordinary kriging
    #' For use with the simulated datasets - no covariate information used here
    
    tmp <- split_data(full_simdata, n = 300, cube.root.transform = TRUE)
    transformed_tmp <- cube_root_fix_zeros(tmp$train)
    tmp$train <- transformed_tmp$transformed_data
    tmp$zeros <- transformed_tmp$zeros
    tmp$epsilon <- transformed_tmp$epsilon
    tmp$Pi <- transformed_tmp$Pi
    
    H <- as.matrix(dist(tmp$train[,c("x", "y")]))
    
    # Use method of moment estimator for gamma distribution 
    mu <- mean(tmp$train$resp[-tmp$zeros])
    beta <- var(tmp$train$resp[-tmp$zeros])/mu
    mu <- rep(mu, times=nrow(tmp$train))
    
    #Calculate covariance parameters using variogram
    emp_var <- variogram(resp~1,
                         loc=~x+y,
                         data = tmp$train)
    v_fit <- tryCatch({
        fit.variogram(emp_var,
                         vgm("Gau"))
    }, warning = function(war) {
        warning("Variogram did not converge.")
        stop("Variogram did not converge.")
    }, message = function(mes) {
        message("Hmm something weird.")
    })
    if (plot.var) {
        plot(emp_var,v_fit)
    }
    
    alphaN<-v_fit$psill[2]/sum(v_fit$psill) # nugget parameter
    alphaR<-v_fit$range[2] # decay parameter
   
    #Calculate spatial correlation matrix w. exponential form for observed points!
    Sigma_obs<-alphaN*exp(-(H/alphaR)^2)
    diag(Sigma)<-1
    
    #Calculate spatial correlation matrix for observed/unobserved points! 
    H_noobs <- as.matrix(rdist(tmp$train[,c("x", "y")], tmp$test[,c("x", "y")]))
    Sigma_noobs <- alphaN*exp(-(H_noobs/alphaR)^2)
    
    output <- list(
      Y = tmp$train$resp,
      Sigma = Sigma_obs,
      Sigma_noobs = Sigma_noobs,
      alphaN = alphaN,
      alphaR = alphaR,
      epsilon = tmp$epsilon,
      mu = mu,
      beta = beta,
      Pi = tmp$Pi,
      test = tmp$test,
      train = tmp$train
    )
    
    return(output)
}

calculate_yhat <- function(full_simdata, n = 300, cube.root.transform = TRUE, plot.var = FALSE) {
  #' Calculate yhat using spatial gaussian copula
  
  spt_tmp <- calculate_spatial_params(full_simdata, n = n, cube.root.transform = cube.root.transform, plot.var = FALSE)
  z <- qnorm(pzig(y = spt_tmp$Y, mu = spt_tmp$mu, beta = spt_tmp$beta, epsilon = spt_tmp$epsilon, Pi = spt_tmp$Pi))
  S_obs_inv <- solve(spt_tmp$Sigma, tol = 1e-22)
  S_noobs_tr <- t(spt_tmp$Sigma_noobs)
  zhat <- S_noobs_tr %*% S_obs_inv %*% z
  hist(pzig(y = spt_tmp$Y, mu = spt_tmp$mu, beta = spt_tmp$beta, epsilon = spt_tmp$epsilon, Pi = spt_tmp$Pi)) 
  pz <- pnorm(zhat)
  pz[which(pz==1)] <- 1-.Machine$double.eps
  pz[which(pz==0)] <- .Machine$double.eps
  
  zigs <- qzig(u = pz, mu = spt_tmp$mu, beta = spt_tmp$beta, epsilon = spt_tmp$epsilon, Pi = spt_tmp$Pi[1])  
}

hist(spt_tmp$Y)
hist(pzig(y = spt_tmp$Y, mu = spt_tmp$mu, beta = spt_tmp$beta, epsilon = spt_tmp$epsilon, Pi = spt_tmp$Pi))
hist(z)
qqnorm(z)
qqline(z)


ord_kriging_z <- function(training_data, test_points, test_z) {
  
  if (!("x" %in% colnames(training_data) && "y" %in% colnames(training_data))){
    stop("Need to specify coordinates using x and y")
  }
  if (!("x" %in% colnames(test_points) && "y" %in% colnames(test_points))){
    stop("Need to specify coordinates using x and y")
  }
  
  training_data$z <- test_z
  test_points <- test_points[,c("x", "y")]
  coordinates(training_data) <- ~x+y
  coordinates(test_points) <- ~x+y
  
  z.vgm <- variogram(z ~ 1, training_data)
  tryCatch({
    z.fit <- fit.variogram(z.vgm, vgm("Gau"))
  }, warning = function(war) {
    stop("Variogram did not converge.")
  })
  
  tryCatch({
    invisible(ok_model_z <- krige(z ~ 1, training_data, test_points, model = z.fit))
  }, warning = function(war) {
    stop("Ordinary Kriging not successful.")
  })
  
  return(ok_model_z)
}

ord_kriging <- function(training_data, test_points) {
  
  if (!("x" %in% colnames(training_data) && "y" %in% colnames(training_data))){
    stop("Need to specify coordinates using x and y in training set")
  }
  if (!("x" %in% colnames(test_points) && "y" %in% colnames(test_points))){
    stop("Need to specify coordinates using x and y in test set")
  }
  
  test_points <- test_points[,c("x", "y")]
  coordinates(training_data) <- ~x+y
  coordinates(test_points) <- ~x+y
  
  resp.vgm <- variogram(resp ~ 1, training_data)
  tryCatch({
    resp.fit <- fit.variogram(z.vgm, vgm("Gau"))
  }, warning = function(war) {
    stop("Variogram did not converge.")
  })
  
  tryCatch({
    invisible(ok_model <- krige(resp ~ 1, training_data, test_points, model = resp.fit))
  }, warning = function(war) {
    stop("Ordinary Kriging not successful.")
  })
  
  return(ok_model)
}

backtransform <- function(preds, epsilon) {
  preds[preds < epsilon] <- 0
  return(preds^3)
}

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

calculate_rfsp <- function(training_set, test_points) {
  
  if (!("x" %in% colnames(training_data) && "y" %in% colnames(training_data))){
    stop("Need to specify coordinates using x and y in training set")
  }
  if (!("x" %in% colnames(test_points) && "y" %in% colnames(test_points))){
    stop("Need to specify coordinates using x and y in test set")
  }
  
  test_points <- test_points[,c("x", "y")]
  coordinates(training_data) <- ~x+y
  coordinates(test_points) <- ~x+y
}

library(raster)
library(GSIF)
training_data <- spt_tmp$train
test_data <- spt_tmp$test

points <- SpatialPointsDataFrame(training_data[,c("x", "y")], as.data.frame(training_data[,c("resp")]))
training_points_spdf <- SpatialPixelsDataFrame(points = points, data = points@data, tolerance = .61)

test_points <- SpatialPointsDataFrame(test_data[,c("x", "y")], data.frame(resp = test_data[,c("resp")]))
test_points_spdf <- SpatialPixelsDataFrame(points = test_points, data = test_points@data, tolerance = .61)

coordinates(training_data) <- ~x+y
grid.dist0 <- GSIF::buffer.dist(training_data, training_points_spdf, as.factor(1:nrow(training_data)))
dn0 <- paste(names(grid.dist0), collapse="+") 
fm0 <- as.formula(paste("resp ~ ", dn0))   
ov.resp <- over(training_data, grid.dist0)

rm.resp <- cbind(data.frame(resp = training_data$resp), ov.resp) 
resp_ranger <- ranger(fm0, rm.resp, quantreg=TRUE, num.trees=150) 
resp_ranger

resp.predict <- predict(resp_ranger, grid.dist0@data)
hist(resp.predict$predictions)

# need to predict values at the test points?

coordinates(training_data) <- ~x+y
grid.dist1 <- GSIF::buffer.dist(training_data, test_points_spdf, as.factor(1:nrow(training_data)))
dn1 <- paste(names(grid.dist1), collapse="+") 
fm1 <- as.formula(paste("resp ~ ", dn0))   
ov.resp <- over(training_data, grid.dist1)

rm.resp <- cbind(data.frame(resp = training_data$resp), ov.resp) 
summary(rm.resp)
resp_ranger <- ranger(fm0, rm.resp, quantreg=TRUE, num.trees=150) 
resp_ranger

resp.predict <- predict(resp_ranger, grid.dist1@data)
plot(test_data$resp, resp.predict$predictions)
