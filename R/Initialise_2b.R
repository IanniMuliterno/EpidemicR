#' Initialise the parameters and infection times.
#'
#' This function initialises an epidemic (checking it is valid) which has two infection rates.
#'
#' @param N The size of the total population.
#' @param rem.times A vector of the removal times.
#' @param beta1.true The true value of beta1, will fix beta1 if given, defaults to NA.
#' @param beta2.true The true value of beta2, will fix beta2 if given, defaults to NA.
#' @param p.true The true value of p, will fix p if given, defaults to NA.
#' @param dist.true The true value of d, will fix d if given, defaults to NA.
#' @param inf.times A vector of the infection times.
#' @param d.lower Lower bound on the proposal of d.
#' @param d.upper Upper bound on the proposal of d.
#' @param dist.mat An NxN distance matrix.
#' @param reparam TRUE/FALSE statement of whether we are using the reparameterisation beta2 = p * beta1 (p in (0,1)).
#'
#' @keywords Initialise epidemic two betas
#' @export
#'
#' @return Returns a list with values (beta1.cur, beta2.cur, p.cur, dist.cur, inf.times).
#'
#' @examples
#' Initilise_2b(beta1.true = 0.004, dist.mat = distance.mat, reparam = T)

# Randomly generate viable initial choices for the infection times
# If the generated infection times are valid then the likelihood will be non-zero.

Initialise_2b <- function(N=N, rem.times = rem.times, 
                         beta1.true = NA, beta2.true = NA, p.true = NA, dist.true = NA, inf.times = NA,
                         d.lower = d.lower, d.upper = d.upper, dist.mat, 
                         reparam = FALSE){

  repeat{
    # Initialise parameters
    if(is.na(beta1.true)){  beta1.cur <- runif(1,0,0.4) }else{ beta1.cur <- beta1.true } 
    if(is.na(beta2.true)){  beta2.cur <- runif(1,0,beta1.cur) }else{ beta2.cur <- beta2.true }
    if(is.na(p.true)){  p.cur <- runif(1,0,1) }else{ p.cur <- p.true }
    if(is.na(dist.true)){  dist.cur <- runif(1,d.lower,d.upper) }else{ dist.cur <- dist.true }
       
    # Initialise the beta matrix
    if(reparam == T){beta.mat <- Beta_mat_form(dist.mat, c(beta1.cur, (p.cur*beta1.cur)), dist.cur)}else{
                     beta.mat <- Beta_mat_form(dist.mat, c(beta1.cur, beta2.cur), dist.cur)
    }
    
    # Initialise the infection times
    if(anyNA(inf.times)){
      inf.periods <- rexp(N, 0.15)
      inf.times.cur <- rem.times - inf.periods
    }else{
      inf.times.cur <- inf.times
    }
    
    #Calculate the log-likelihood
    if(log_likelihood(inf.times.cur, rem.times, beta.mat) != -Inf){break}
  }
  
  print("initialised")
  
  return(list(beta1.cur, beta2.cur, p.cur, dist.cur, inf.times.cur))
}