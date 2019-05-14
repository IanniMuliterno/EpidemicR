#' Calculate the acceptance probability of a Metropolis-Hastings step (ratio of likelihoods case only.)
#'
#' This function calculates alpha, the acceptance probability of a MH step, for a random walk or similar where the proposal is
#'  symmetric so cancels and one is left with a ratio of likelihoods.
#'
#' @param llh The value of the log-likelihood with the current value of the paramter in question.
#' @param llh.prime The value of the log-likelihood with the newly drawn value of the paramter in question.
#'
#' @keywords alpha acceptance probability MH Metropolis-Hastings Metropolis Hastings
#' @export
#'
#' @return Returns the min(1, ratio of likelihoods).
#'
#' @examples
#' MH_accept(llh, llh.prime)

MH_accept <- function(llh, llh.prime){
  
  alpha <- exp( llh.prime - llh )
  
  return(min(c(1, alpha)))
}