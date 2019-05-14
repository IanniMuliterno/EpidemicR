#' Make inference on a General Stochastic Epidemic
#'
#' This function uses realistically available information from a GSE to perform Bayesian inference and attempt to recover the parameters.
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
#' @param lambda.b1 The rate parameter for beta1, assuming a Gamma prior.
#' @param nu.b1 The shape parameter for beta1, assuming a Gamma prior.
#' @param lambda.b2 The rate parameter for beta2, assuming a Gamma prior.
#' @param nu.b2 The shape parameter for beta2, assuming a Gamma prior.
#' @param lambda.g The rate parameter for gamma, assuming a Gamma prior.
#' @param nu.g The shape parameter for gamma, assuming a Gamma prior.
#' @param inc.beta1 A list object of 2 levels, the true value of beta1, and T/F binary value that says 
#'  whether to make inference on beta1 (T = make inference, see details for the different options).
#' @param inc.beta2 See inc.beta1, but for beta2.
#' @param inc.dist See inc.beta1, but for the distance d.
#' @param inc.inf.times A list object of 2 levels; a vector of the true infection times, and T/F binary value that says
#'  whether to make inference on the infection times (T = make inference, see details for the different options).
#' @param inc.gamma See inc.beta1, but for the removal rate gamma.
#' @param d.upper The upper bound on the distance d.
#' @param sigmab1 Tuning parameter. The standard deviation of the multiplicative random walk for beta1.
#' @param sigmab2 Tuning parameter. The standard deviation of the multiplicative random walk for beta1.
#' @param sigmad Tuning parameter. The standard deviation of the folded random walk for d.
#' @param infupdate Tuning parameter. The number of infection times that should be updated at each iteration.
#'
#' @keywords MCMC Gibbs MH Metropolis Hastings inference reparameterisation
#' @export
#'
#' @return This function returns a list object with elements; a matrix of results (which included all the accepted samples and
#'  the log-likelihood at the end of each iteration), the acceptance rate of the infection times, the acceptance rate of beta1,
#'  the acceptance rate of beta2, and the acceptance rate of d.
#' 
#' @examples
#' inference <- MCMC(N.its = 100000, N = 10, inf.ids, rem.times, dist.mat, 
#' lambda.b1 = 0.001, nu.b1 = 1, lambda.b2 = 0.001, nu.b2 = 1, lambda.g = 0.001 , nu.g = 1,
#' inc.beta1 = list(0.004, T), inc.beta2 = list(NA, T), inc.dist = list(NA, T),
#' inc.inf.times = list(inf.times, F), inc.gamma = list(NA, T),
#' d.upper = 1.5, sigmab1 = 1, sigmab2 = 1, sigmad = 2, infupdate = 1)



MCMC <- function(N.its, N, inf.ids, rem.times, dist.mat, 
                 lambda.b1 = 0.001, nu.b1 = 1, lambda.b2 = 0.001, nu.b2 = 1, lambda.g = 0.001 , nu.g = 1,
                 inc.beta1 = list(NA, T), inc.beta2 = list(NA, T), inc.dist = list(NA, T),
                 inc.inf.times = list(NA, T), inc.gamma = list(NA, T),
                 d.upper, sigmab1, sigmab2, sigmad, infupdate = 1){
  
  ############################################    
  ### Calculate upper and lower bound on d ###
  ############################################
      
      inf.dist.mat <- dist.mat[inf.ids, inf.ids]
      d.lower <- min(inf.dist.mat[inf.dist.mat > 0])
  
  ##############################################
  ### Initialise the parameters and epidemic ###
  ##############################################

      InitialiseEpi <- Initialise_2b(N=N, rem.times,
                                     beta1.true = inc.beta1[[1]], beta2.true = inc.beta2[[1]], p.true = NA,
                                     dist.true = inc.dist[[1]], inf.times = inc.inf.times[[1]],
                                     d.lower = d.lower, d.upper = d.upper, dist.mat, 
                                     reparam = FALSE)
      beta1.cur <- InitialiseEpi[[1]]
      beta2.cur <- InitialiseEpi[[2]]
      p.cur <- InitialiseEpi[[3]]
      dist.cur <- InitialiseEpi[[4]]
      inf.times <- InitialiseEpi[[5]]
      beta.mat <- Beta_mat_form(dist.mat, c(beta1.cur, beta2.cur), dist.cur)
  
  ######################    
  ### Results matrix ###
  ######################
      
      res <- matrix(ncol = 7, nrow = N.its)
      colnames(res) <- c("sample", "beta1", "beta2", "gamma", "d", "p", "llh")
  
  ##########################
  ### Functional objects ###
  ##########################
      
      it = 1
      n.I = length(inf.ids)
      acc.sum.I = acc.sum.B1 = acc.sum.B2 = acc.sum.d = 0
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
          
        
        #########################
        ### MH step for beta1 ###
        #########################
        
        if(inc.beta1[[2]] == T){
          
          # Draw new beta1 value
          beta1.draw <- exp(rnorm( 1, log(beta1.cur), sigmab1))
          
          # Calculate functional objects
          beta.mat.prime <- Beta_mat_form(dist.mat, c(beta1.draw, beta2.cur), dist.cur)
          llh.prime <- log_likelihood(inf.times, rem.times, beta.mat.prime)
          
          # MH acceptance probability
          alpha.B1 = MH_accept_MRW(llh, llh.prime, nu.b1, lambda.b1, beta1.cur, beta1.draw)
          
          # Do we accept the new p or not?
          accept.test <- runif(1,0,1)
          if(accept.test < alpha.B1){ #If yes:
            
            beta1.cur <- beta1.draw
            beta.mat <- beta.mat.prime
            llh <- llh.prime
            
            acc.sum.B1 = acc.sum.B1+1
          }
        }

       
        #########################
        ### MH step for beta2 ###
        #########################
        
        if(inc.beta2[[2]] == T){
          
          # Draw new beta1 value
          beta2.draw <- exp(rnorm( 1, log(beta2.cur), sigmab2))
          
          # Calculate functional objects
          beta.mat.prime <- Beta_mat_form(dist.mat, c(beta1.cur, beta2.draw), dist.cur)
          llh.prime <- log_likelihood(inf.times, rem.times, beta.mat.prime)
          
          # MH acceptance probability
          alpha.B2 = MH_accept_MRW(llh, llh.prime, nu.b2, lambda.b2, beta2.cur, beta2.draw)
          
          # Do we accept the new p or not?
          accept.test <- runif(1,0,1)
          if(accept.test < alpha.B2){ #If yes:
            
            beta2.cur <- beta2.draw
            beta.mat <- beta.mat.prime
            llh <- llh.prime
            
            acc.sum.B2 = acc.sum.B2+1
          }
        }  
          
          
        #####################
        ### MH step for d ###
        #####################
        
        if(inc.dist[[2]] == T){  
          
            # Draw new d value
            dist.draw <- Folded_draw(param.cur = dist.cur, param.sigma = sigmad, lower = d.lower, upper = d.upper)
            
            # Calculate functional objects
            beta.mat.prime <- Beta_mat_form(dist.mat, c(beta1.cur, (p.cur*beta1.cur)), dist.draw)
            llh.prime <- log_likelihood(inf.times, rem.times, beta.mat.prime)
            
            # MH acceptance probability
            alpha.d = MH_accept(llh, llh.prime)
            
            # Do we accept the new d or not?
            accept.test <- runif(1,0,1)
            if(accept.test < alpha.d){ #If yes:
              
              dist.cur <- dist.draw
              beta.mat <- beta.mat.prime
              llh <- llh.prime
              
              acc.sum.d = acc.sum.d+1
            }     
        }
            
        ##########################
        ### Record the results ###
        ##########################
        
            # Record parameters
            res[it,] <- c(it, beta1.cur, beta2.cur, gamma.cur, dist.cur, (beta2.cur/beta1.cur) , llh)
            
            # Update count
            it = it + 1
      }
  
  return(list(res, (acc.sum.I/N.its), (acc.sum.B1/N.its), (acc.sum.B2/N.its), (acc.sum.d/N.its)))
}
