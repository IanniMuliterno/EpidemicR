#' Make inference on a homogeneous General Stochastic Epidemic.
#'
#' This function uses realistically available information from a homogeneous GSE to perform Bayesian inference and attempt to
#'  recover the parameters.
#'
#' The function has the functionality to allow the fixing of some (or all) parameters. For the function parameters that begin inc.
#'  a list can be fed to them of the form (true parameter value(s), TRUE/FALSE), where the TRUE dictates that inference should be
#'  made for that parameter. Possible options for this list are; list(parameters, F) in which case the function will fix that parameter
#'  value, list(NA, T) in which case the function will randomly generate a valid initial value for the parameter and then make inference,
#'  list(parameter, T) in which case the function will initialise the parameter at its true value and then make inference, or (NA, F) 
#'  which we do not suggest using as it will fix the parameter at a random value.
#'  
#'  Gamma is always assumed to be initialised at the value 0.15.
#'
#' @param N.its The number of desired iterations of the MCMC algorithm.
#' @param N The total size of the population.
#' @param inf.ids A vectors of the IDs of the infected individuals.
#' @param rem.times A vector of the removal times, ordered by individual ID.
#' @param dist.mat An NxN distance matrix.
#' @param lambda.b The rate parameter for beta, assuming a Gamma prior.
#' @param nu.b The shape parameter for beta, assuming a Gamma prior.
#' @param lambda.g The rate parameter for gamma, assuming a Gamma prior.
#' @param nu.g The shape parameter for gamma, assuming a Gamma prior.
#' @param inc.beta A list object of 2 levels, the true value of beta1, and T/F binary value that says 
#'  whether to make inference on beta1 (T = make inference, see details for the different options).
#' @param inc.inf.times A list object of 2 levels; a vector of the true infection times, and T/F binary value that says
#'  whether to make inference on the infection times (T = make inference, see details for the different options).
#' @param inc.gamma See inc.beta, but for the removal rate gamma.
#' @param infupdate Tuning parameter. The number of infection times that should be updated at each iteration.
#'
#' @keywords MCMC Gibbs MH Metropolis Hastings inference reparameterisation
#' @export
#'
#' @return This function returns a list object with elements; a matrix of results (which included all the accepted samples and
#'  the log-likelihood at the end of each iteration), the acceptance rate of the infection times, the acceptance rate of p, and
#'  the acceptance rate of d.
#' 
#' @examples
#' infernce_r <- MCMC_homo(N.its, N, inf.ids, rem.times, dist.mat, 
#' lambda.b = 0.001, nu.b = 1, lambda.g = 0.001 , nu.g = 1,
#' inc.beta = list(NA, T), inc.inf.times = list(NA, T), inc.gamma = list(NA, T),
#' infupdate = 1)


MCMC_homo <- function(N.its, N, inf.ids, rem.times, dist.mat, 
                         lambda.b = 0.001, nu.b = 1, lambda.g = 0.001 , nu.g = 1,
                         inc.beta = list(NA, T), inc.inf.times = list(NA, T), inc.gamma = list(NA, T),
                         infupdate = 1){
  
  ##############################################
  ### Initialise the parameters and epidemic ###
  ##############################################
  
  InitialiseEpi <- Initialise_2b(N=N, rem.times,
                                 beta1.true = inc.beta[[1]], beta2.true = NA, p.true = 1,
                                 dist.true = 1, inf.times = inc.inf.times[[1]],
                                 d.lower = 1, d.upper = 2, dist.mat, 
                                 reparam = TRUE)
  beta.cur <- InitialiseEpi[[1]]
  beta2.cur <- InitialiseEpi[[2]]
  p.cur <- InitialiseEpi[[3]]
  dist.cur <- InitialiseEpi[[4]]
  inf.times <- InitialiseEpi[[5]]
  beta.mat <- Beta_mat_form(dist.mat, c(beta.cur, beta.cur), dist.cur)
  
  ######################    
  ### Results matrix ###
  ######################
  
  res <- matrix(ncol = 4, nrow = N.its)
  colnames(res) <- c("sample", "beta", "gamma", "llh")
  
  ##########################
  ### Functional objects ###
  ##########################
  
  it = 1
  n.I = length(inf.ids)
  acc.sum.I = 0
  llh <- log_likelihood(inf.times, rem.times, beta.mat)
  
  ###########################
  ### ~~ THE ALGORITHM ~~ ###
  ###########################
  
  while(it <= N.its){
    
    ###############################
    ### Gibbs sampler for gamma ###
    ###############################
    
    if(inc.gamma[[2]] == T){
      
      #Calculate the removal integral
      g.ints <- sum(rem.times[inf.ids] - inf.times[inf.ids])
      
      # Draw gamma
      gamma.cur<- rgamma(n=1, shape = (n.I+nu.g), rate = (lambda.g + g.ints))
      
    }
    
    
    ###############################
    ### Gibbs sampler for beta ###
    ###############################
    
    if(inc.beta[[2]] == T){
      
      # Form the p matrix
      p.mat <- P_mat_form(dist.mat, ps = c(1, 1), 1)
      
      # Calculate the infection integral
      b.ints <- integral_part(inf.times, rem.times, p.mat, inf.ids)
      
      # Draw beta
      beta.cur <- rgamma(n=1, shape = (n.I-1+nu.b), rate = (lambda.b + b.ints))
      
      # Calculate functional objects
      beta.mat <- Beta_mat_form(dist.mat, c(beta.cur, beta.cur), 1)
      llh <- log_likelihood(inf.times, rem.times, beta.mat)
      
    }
    
    
    ###################################
    ### MH step for infection times ###
    ###################################
    
    if(inc.inf.times[[2]] == T){
      
      # Which infection time is being replaced
      Ireplace <- sample(inf.ids, infupdate)
      
      # Draw new infection time
      Qdraw <- rexp(infupdate , gamma.cur) 
      inf.times.prime <- inf.times
      inf.times.prime[Ireplace] <- rem.times[Ireplace] - Qdraw
      
      # Calculate functional objects
      llh.prime <- log_likelihood(inf.times.prime, rem.times, beta.mat)
      
      # MH acceptance probability
      alpha.I = MH_accept(llh, llh.prime)
      
      # Do we accept the new time(s) or not?
      accept.test <- runif(1,0,1)
      if(accept.test < alpha.I){ #If yes:
        
        inf.times <- inf.times.prime
        llh <- llh.prime
        
        acc.sum.I = acc.sum.I+1
      }
    }
    
    
    ##########################
    ### Record the results ###
    ##########################
    
    # Record parameters
    res[it,] <- c(it, beta.cur, gamma.cur, llh)
    
    # Update count
    it = it + 1
  }
  
  return(list(res, (acc.sum.I/N.its)))
}
