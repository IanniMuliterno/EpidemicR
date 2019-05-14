#' Generate an infection rate matrix.
#'
#' This function generates a beta matrix given a distance function, values of beta, and the distances at which
#'  the infection rate changes. It works for a heterogeneous population with a near vs. far spatial kernel.
#'
#' @param dist.mat An NxN matrix of the pairwise distances between individuals.
#' @param betas A descending vector of at least length 2 of the different infection rates at each level.
#' @param ds An acsending vector of length length(betas)-1 which contains the distances at which the infection rate changes.
#'
#' @keywords beta infection rate matrix
#' @export
#'
#' @return The function returns an NxN matrix of the pairwise infection rates between individuals.
#'
#' @examples
#' distance_mat <- Dist_mat_unif(N=100, xlim = 20, ylim = 20)[[2]]
#' Beta_mat_form(distance_mat, c(0.004, 0.002), 10)

Beta_mat_form <- function(dist.mat, betas, ds){
  
  levels <- length(betas)
  
  beta.mat <- dist.mat
  
  for(i in (levels-1):1){
    beta.mat[which(dist.mat < ds[i])] <- betas[i]
  }
  
  beta.mat[which(dist.mat >= ds[(levels-1)])] <- betas[levels]
  diag(beta.mat) <- 0
  
  return(beta.mat)
}
