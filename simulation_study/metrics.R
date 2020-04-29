# MSc Project 2020: Comparison of Spatial Gaussian Copula and Random Forests in Spatial Prediction
# Sun, N
#
# Code for calculating prediction metrics

rmspe <- function(predicted_values, observed_values) {
    
    #' Computes root mean squared prediction error 
    #' 
    #' @param predicted_values vector of predicted values
    #' @param observed_values vector of observed values
    #' 
    #' @return RMSPE single value
    
    if(length(predicted_values) != length(observed_values)) {
        stop("Vectors are of unequal length.")
    }
    
    sq.diffs <- (predicted_values - observed_values)^2
    return(sqrt(mean(sq.diffs)))
     
}

mspe <- function(predicted_values, observed_values) {
    #' Computes sum of squared residuals
    #' 
    #' @param predicted_values vector of predicted values
    #' @param observed_values vector of observed values
    #' 
    #' @return RMSPE single value
    
    if(length(predicted_values) != length(observed_values)) {
        stop("Vectors are of unequal length.")
    }
    
    sum.sq.diffs <- sum((predicted_values - observed_values)^2)
    return((sum.sq.diffs))
}

make_plot <- function(data) {
    
    ggplot(data, aes(x=x/1000,y=y/1000,color=totvol)) +
      geom_point(size=2) +
      scale_color_gradient(low = "white", high = "red") +
      scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
      theme(axis.text.x=element_text(angle=90,hjust=1))
    
}