#' Generate an proportion of maximum infection rate matrix.
#'
#' This function generates a p matrix given a distance function, values of beta, and the distances at which
#'  the infection rate changes. It works for a heterogeneous population with a near vs. far spatial kernel.
#'  It is essentially just a wrapper for the Beta_mat_form function, but used to calculate p.
#'
#' @param dist.mat An NxN matrix of the pairwise distances between individuals.
#' @param ps A descending vector of at least length 2 of the different proportions of maximum infection rate at each level.
#'           The first element is assumed to be 1.
#' @param ds An acsending vector of length length(betas)-1 which contains the distances at which the infection rate changes.
#'
#' @keywords p proportion infection rate matrix
#' @export
#'
#' @return The function returns an NxN matrix of the pairwise proportion of maximum infection rates between individuals.
#'
#' @examples
#' distance_mat <- Dist_mat_unif(N=100, xlim = 20, ylim = 20)[[2]]
#' P_mat_form(distance_mat, c(1, 0.1), 10)

P_mat_form <- function(dist.mat, ps = 1, ds){
  
  levels <- length(ps)
  
  p.mat <- dist.mat
  
  for(i in (levels-1):1){
    p.mat[which(dist.mat < ds[i])] <- ps[i]
  }
  
  p.mat[which(dist.mat >= ds[(levels-1)])] <- ps[levels]
  diag(p.mat) <- 0
  
  return(p.mat)
}
