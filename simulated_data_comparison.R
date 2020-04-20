setwd("~/MS_Project/code")
simdata <- readRDS("simulated_Y.rds")
dat <- read.csv("ForestDataFuzzed.csv")

source("ZIGFunctions.R")
require(gstat)
require(tidyverse)
require(rgdal)

alb.xy <- project(with(dat,cbind(lon_fuzzed,lat_fuzzed)),
                  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
colnames(alb.xy) <- c("x","y")
dat <- cbind(alb.xy,dat)
full_simdata <- data.frame(
    x = dat$x,
    y = dat$y,
    resp = simdata[2,]
)

split_data <- function(full_data, n = 300, cube.root.transform = TRUE) {
  
  colnames(full_data) <- c("x", "y", "resp")
  training_idx <- sample(1:nrow(full_data), n, replace=FALSE)
  training_data <- full_data[training_idx,]
  if (cube.root.transform) {
    training_data$resp <- training_data$resp^(1/3)
  }
  zeros <- which(training_data$resp == 0)
  epsilon <- min(training_data$resp[-zeros])
  Pi <- rep(length(zeros)/n, n)
  training_data$resp[zeros] <- runif(length(zeros), 0, epsilon)
  
  test_data <- full_data[-training_idx,]
  output <- list(train = training_data, test = test_data,
                 zeros = zeros,
                 epsilon = epsilon,
                 Pi = Pi) 
  return(output)
}

calculate_spatial_params <- function(full_simdata, n = 300, cube.root.transform = TRUE, plot.var = FALSE) {
    #' Calculates Method of Moments estimates for Gamma model and the spatial covariance parameters for ordinary kriging
    #' For use with the simulated datasets - no covariate information used here
    
    tmp <- split_data(full_simdata, n = 300, cube.root.transform = TRUE)
    
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
      Pi = tmp$Pi
    )
    
    return(output)
}

calculate_yhat <- function(full_simdata, n = 300, cube.root.transform = TRUE, plot.var = FALSE) {
  spt_tmp <- calculate_spatial_params(full_simdata, n = n, cube.root.transform = cube.root.transform, plot.var = FALSE)
  z <- qnorm(pzig(y = spt_tmp$Y, mu = spt_tmp$mu, beta = spt_tmp$beta, epsilon = spt_tmp$epsilon, Pi = spt_tmp$Pi))
  S_obs_inv <- solve(spt_tmp$Sigma, tol = 1e-22)
  S_noobs_tr <- t(spt_tmp$Sigma_noobs)
  zhat <- S_noobs_tr %*% S_obs_inv %*% z
  
  pz <- pnorm(zhat)
  pz[which(pz==1)] <- 1-.Machine$double.eps
  pz[which(pz==0)] <- .Machine$double.eps
  
  zigs <- qzig(u = pz, mu = spt_tmp$mu, beta = spt_tmp$beta, epsilon = spt_tmp$epsilon, Pi = spt_tmp$Pi[1])  
}
