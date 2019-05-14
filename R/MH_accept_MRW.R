#' Calculate the acceptance probability of a Metropolis-Hastings step for a multiplicative random walk.
#'
#' This function calculates alpha, the acceptance probability of a MH step, for a multiplicative random walk, which is where proposals
#'  are made on the log scale. It assumes a Gamma prior on the parameter of interest.
#'
#' @param llh The value of the log-likelihood with the current value of the paramter in question.
#' @param llh.prime The value of the log-likelihood with the newly drawn value of the paramter in question.
#' @param nu.prior Assuming a Gamma prior on the parameter in question, the value of the parameter nu in that prior.
#' @param lambda.prior Assuming a Gamma prior on the parameter in question, the value of the parameter lambda in that prior.
#' @param param.cur The current value of the parameter of interest.
#' @param param.prime The proposed value of the parameter of interest.
#'
#' @keywords alpha acceptance probability MH Metropolis-Hastings Metropolis Hastings multiplicative random walk
#' @export
#'
#' @return Returns the min(1, alpha).
#'
#' @examples
#' MH_accept_MRW(llh, llh.prime, nu.beta, lambda.beta, beta.cur, beta.prime)

MH_accept_MRW <- function(llh, llh.prime, nu.prior, lambda.prior, param.cur, param.prime){
  
  alpha <- exp( llh.prime + (nu.prior)*log(param.prime) - lambda.prior*param.prime -
                  ( llh + (nu.prior)*log(param.cur) - lambda.prior*param.cur)  )
  
  return(min(c(1, alpha)))
}