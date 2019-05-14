#' Plot the trace plots of the accepted samples.
#'
#' This function takes the results from an inference function, and plots the trace plot of all the accepted samples,
#'  with a red line to represent the true value if it falls within the limits.
#'  If a parameter was fixed it will plot a straight line.
#'
#' @param res The results matrix from the MCMC functions.
#' @param params A named vector of the parameters, with names "beta1", "beta2", "gamma", "d", "p", and "llh".
#' @param burnin The ammount of burnin that wants to be removed (defaults to NA to plot the full trace plots).
#'
#' @keywords description results acceptance probabilities
#' @export
#'
#' @return A description of the average acceptance probabilities.
#'
#' @examples
#' 
#' Plot_trace(res, params, burnin = NA)
#' 
#' 


Plot_trace <- function(res, params, burnin = NA){
  
  if(is.na(burnin)){samp <- 1:length(res[,1]); sub <- "Full"}else{samp <- -(1:burnin); sub <- "Removed burn-in"}
  
  
    par(mfrow = c(3,2))
    
    {plot(res[samp ,2], type = "l", main = "beta_1", sub = sub, ylab = "Value")
    abline(h = params["beta1"], col = "red")}
    
    {plot(res[samp ,3], type = "l", main = "beta_2", sub = sub, ylab = "Value")
    abline(h = params["beta2"], col = "red")}
    
    {plot(res[samp ,4], type = "l", main = "gamma", sub = sub, ylab = "Value")
    abline(h = params["gamma"], col = "red")}
    
    {plot(res[samp ,5], type = "l", main = "d", sub = sub, ylab = "Value")
    abline(h = params["d"], col = "red")}
    
    {plot(res[samp ,6], type = "l", main = "p", sub = sub, ylab = "Value")
      abline(h = params["p"], col = "red")}
    
    {plot(res[samp ,7], type = "l", main = "Log-likelihood", sub = sub, ylab = "Value")
      abline(h = params["llh"], col = "red")}
  
}