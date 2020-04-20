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
    cov1 = dat$annpre,
    resp = simdata[1,]
)

calculate_theta <- function(full_data, n = 300, plot.var = FALSE, echo.time=FALSE, optim.method = "BFGS") {
    #' Calculates the MLEs by optimizing the negative log likelihood of the spatial Gaussian copula model
    
    training_n <- n
    training_idx <- sample(1:nrow(full_data), training_n, replace=FALSE)
    training_data <- full_data[training_idx,]
    colnames(training_data) <- c("x", "y", "cov1", "resp")
    
    training_data$resp <- (training_data$resp^(1/3))
    
    if(ncol(training_data) < 4) {
        stop("Input should have x, y, covariate, and response values")
    }
    
    
    H <- as.matrix(dist(training_data[,c("x", "y")]))
    
    zeros <- which(training_data$resp == 0) 
    epsilon <- min(training_data$resp[-zeros])
    Pi <- rep(length(zeros)/nrow(training_data), nrow(training_data))
    
    #Replace 0s with "small" uniform r.v.
    training_data$resp[zeros] <- runif(length(zeros), 0, epsilon)
    
    # Use method of moment estimator for gamma distribution 
    mu <- mean(training_data$resp[-zeros])
    beta <- var(training_data$resp[-zeros])/mu
    mu <- rep(mu, times=nrow(training_data))
    
    #Calculate covariance parameters using variogram
    emp_var <- variogram(resp~1,
                       loc=~x+y,
                       training_data)
    v_fit <- tryCatch({
        fit.variogram(emp_var,
                         vgm("Gau"))
    }, warning = function(war) {
        warning("Variogram did not converge.")
        stop("Variogram did not converge.")
    }, message = function(mes) {
        warning("Variogram did not converge.")
        stop("Variogram did not converge.")
    })
    
    if (plot.var) {
        plot(emp_var,v_fit)
    }
    
    alphaN<-v_fit$psill[2]/sum(v_fit$psill) # nugget parameter
    alphaR<-v_fit$range[2] # decay parameter
   
    #Calculate spatial correlation matrix w. exponential form 
    Sigma<-alphaN*exp(-(H/alphaR)^2)
    diag(Sigma)<-1
    
    #Calculate logistic regression coefficients
    X <- matrix(cbind(rep(1, nrow(training_data)),
                      training_data$cov1),
                nrow = nrow(training_data))
    zero_resp_train <- ifelse(training_data$resp < epsilon, 1, 0)
    zero_glm_train <- glm(zero_resp_train ~ training_data$cov1,
                          family=binomial)
    xi_Z <- zero_glm_train$coefficients
   
    #Calculate gamma regression coefficients 
    nonzero_resp <- which(training_data$resp >= epsilon)
    nonzero_glm_train <- glm(resp ~ cov1,
                             data = training_data[nonzero_resp,],
                             family=Gamma("log"))
    xi_NZ <- nonzero_glm_train$coefficients
    
    #Calculate MLEs using optim()
    theta <- c(logit(alphaN),
               log(alphaR),
               log(1/beta),
               xi_Z,
               xi_NZ)
    names(theta) <- NULL

    start <- Sys.time()
    result <- optim(theta, 
                    NL, 
                    method = optim.method,
                    Y = training_data$resp,
                    epsilon = epsilon,
                    H = H,
                    X = X)
    end <- Sys.time() - start

    data.frame(
      initial = c(alphaN, alphaR, beta, xi_Z, xi_NZ),
      mle = c(alogit(result$par[1]), exp(result$par[2]), 1/exp(result$par[3]), result$par[4:length(result$par)])
    ) -> output
    
    if (echo.time) {
        print(paste0("Calculation of MLEs took ", toString(end), " seconds."))
    }
    
    return(output)
}

calculate_theta(full_simdata, n = 300, plot.var = TRUE, echo.time=TRUE, optim.method = "BFGS")
