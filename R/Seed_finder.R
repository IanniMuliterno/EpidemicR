#' Find a seed for the desired outbreak size
#'
#' Runs a series of independent simulations for a series of seeds, and outputs a seed(s) that give a desired final size,
#'  as well as a histogram of all final sizes.
#'
#' @param model Takes values "Unif", "Clust", "Edge", "Latt", "Lines", or "Trials", corresponding to the appropriate distance matrix function.
#' @param N The population size.
#' @param ... The parameters of the specifief distance matrix function (see example).
#' @param betas A vector of infection rates, used for the Beta_mat_form function.
#' @param ds A vector of distances at which the infection rate changes, used for the Beta_mat_form function.
#' @param size The final size of the epidemic that is desired (total number infected).
#' @param error The allowed error on the final size above and below.
#' @param seeds The number of seeds desired to be tested. Always iterates from 1:seeds.
#'
#' @keywords Histogram seeds seed final size
#' @export
#'
#' @return A histogram of the final size of all tested simulations, and a table of seeds and final sizes that fit the requirements.
#'
#' @examples
#' Seed_finder("Clust", 
#'             N=100, xlim = 20, ylim = 20, centres = 5, spread = 1, 
#'             betas = c(0.1, 0.001), ds = 0.25, 
#'             size = 25, error = 2, seeds = 250)

Seed_finder <- function(model, N, ..., betas, ds, size, error, seeds){
  
  if(model == "Unif"){
    distfunc <- Dist_mat_unif
  }
  if(model == "Clust"){
    distfunc <- Dist_mat_clust
  }
  if(model == "Edge"){
    distfunc <- Dist_mat_edge
  }
  if(model == "Latt"){
    distfunc <- Dist_mat_latt
  }
  if(model == "Lines"){
    distfunc <- Dist_mat_lines
  }
  if(model == "Trials"){
    distfunc <- Dist_mat_trials
  }
  
  final.size <- rep(NA, seeds)
  
  for(k in 1:seeds){
    
    set.seed(k)
    
    # Generate a distance matrix
    distance_mat <- distfunc(N, ...)[[2]]
    
    # Generate an associated infection rate rate matrix
    rate_mat <- Beta_mat_form(distance_mat, betas, ds)
    
    # Generate a simulated epidemic
    Hetero_sim <- GSE_sim(N, beta.mat = rate_mat, gamma = 0.15)
    
    # Total number of infected
    final.size[k] <- sum(Hetero_sim[,3] != Inf)
  }
  
  print(summary(final.size))
  
  hist(final.size, breaks = (seeds/5))
  
  whichseeds <- which(final.size > (size-error) & final.size < (size+error))
  rbind("seed" = whichseeds, "size" = final.size[whichseeds])
}












