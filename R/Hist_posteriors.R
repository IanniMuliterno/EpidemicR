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


Hist_posteriors <- function(res, params, burnin = NA){
  
  if(is.na(burnin)){samp <- 1:length(res[,1]); sub <- "Full"}else{samp <- -(1:burnin); sub <- "Removed burn-in"}
  
  par(mfrow = c(3,2))
  
  {hist(res[samp,2], main = "beta_1", sub=sub, breaks = 100, xlab =  "Value")
    abline(v = params["beta1"], col = "red")}
  
  {hist(res[samp,3], main = "beta_2", sub = sub, breaks = 100, xlab =  "Value")
    abline(v = params["beta2"], col = "red")}
  
  {hist(res[samp,4], main = "gamma", sub = sub, breaks = 100, xlab =  "Value")
    abline(v = params["gamma"], col = "red")}
  
  {hist(res[samp,5], main = "d", sub =sub, breaks = 100, xlab =  "Value")
    abline(v = params["d"], col = "red")}
  
  {hist(res[samp,6], main = "p", sub = sub, breaks = 100, xlab =  "Value")
    abline(v = params["p"], col = "red")}
  
}