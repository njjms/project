# MSc Project 2020: Comparison of Spatial Gaussian Copula and Random Forests in Spatial Prediction
# Sun, N
#
# Code for using spatial random forests to predict values

library(ranger)
library(raster)
library(GSIF)

calculate_rfsp_predictions <- function(training_data, test_points, num.trees = 150, tolerance = .7) {
    
    #' Function for using RFSP to calculate predictions for unobserved data on untransformed data
    #' 
    #' @param training_set dataframe of observed locations, including observed responses
    #' @param test_points dataframe of unobserved locations
  
    if (!("x" %in% colnames(training_data) && "y" %in% colnames(training_data))){
        stop("Need to specify coordinates using x and y in training set")
    }
    if (!("x" %in% colnames(test_points) && "y" %in% colnames(test_points))){
        stop("Need to specify coordinates using x and y in test set")
    }
    
    # coordinates(training_data) <- ~x+y
    # coordinates(test_points) <- ~x+y
  
    training_points <- SpatialPointsDataFrame(training_data[,c("x", "y")], as.data.frame(training_data[,c("resp")]))
    training_points_spdf <- SpatialPixelsDataFrame(points = training_points, data = training_points@data, tolerance = tolerance)
    
    test_points <- SpatialPointsDataFrame(test_points[,c("x", "y")], data.frame(resp = test_points[,c("resp")]))
    test_points_spdf <- SpatialPixelsDataFrame(points = test_points, data = test_points@data, tolerance = tolerance)
    
    grid.dist0 <- GSIF::buffer.dist(training_points, training_points_spdf, as.factor(1:nrow(training_data)))
    dn0 <- paste(names(grid.dist0), collapse="+") 
    fm0 <- as.formula(paste("resp ~ ", dn0))   
    ov.resp <- over(training_points, grid.dist0)
    
    grid.dist1 <- GSIF::buffer.dist(training_points, test_points_spdf, as.factor(1:nrow(training_data)))
    
    rm.resp <- cbind(data.frame(resp = training_data$resp), ov.resp) 
    resp_ranger <- ranger(fm0,
                          rm.resp,
                          quantreg=TRUE,
                          num.trees=num.trees) 
    
    resp.predict <- predict(resp_ranger, grid.dist1@data)
    return(resp.predict)
}

