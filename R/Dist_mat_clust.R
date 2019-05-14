
#' Clustered population distribution.
#'
#' This function generates a distance matrix for a clustered population, across a specified rectangular plane.
#' 
#' The function designates a number of individuals specified to be clusters, around each it then generates an appropriate
#'  number of "children" using an isotrpoic distribution with variance as specificed (defaults to identity). 
#'  
#'  Attempts to ensure no points have negative coordinates by not allowing centres to have coordinate values below two.
#'
#' @param N The total number of individuals in the population.
#' @param xlim The width of the plane on which individuals will be generated (defaults to 20 units wide).
#' @param ylim The height of the place on which individuals will be generated (defaults to 20 units high).
#' @param centres The number of centres (clusters). N must be divisable by the number of centres.
#' @param spread The variance of the isotropic Gaussian distribution around each cluster (defaults to identity).
#'
#' @keywords Clustered population distribution distance matrix
#' @export
#'
#' @return The function returns a list. The first object in the list is a Nx2 matrix of the coordinates of the individuals.
#'   The second object is an NxN distance matrix.
#'
#' @examples
#' xy.coords <- Dist_mat_clust(N=100, xlim = 20, ylim = 20, centres = 5, spread = 1)[[1]]
#' distance_mat <- Dist_mat_clust(N=100, xlim = 20, ylim = 20, centres = 5, spread = 1)[[2]]
#' plot(xy.coords[,1], xy.coords[,2], type  = "p")
#' 

Dist_mat_clust <- function(N, xlim = 20, ylim = 20, centres = 5, spread = 1){
  
  x.clust.coords <- runif(centres, 2, xlim)
  y.clust.coords <- runif(centres, 2, ylim)
  clust.coords = cbind(x.clust.coords, y.clust.coords)
  
  std.dev <- matrix(c(spread,0,0,spread),2,2)
  
  clustees <- (N-centres)/centres
  
  child.coords = NA
  xy.coords = matrix(nrow =  N, ncol = 2)
  xy.coords[1:centres, ] <- clust.coords 
  for(j in 1:centres){
    child.coords <- MASS::mvrnorm(n = clustees, clust.coords[j,], std.dev)
    xy.coords[((1:clustees) + centres + (j-1)*clustees), ] <- child.coords
  }
  
  dist.mat <- as.matrix(dist(xy.coords, upper = T))
  
  return(list(xy.coords, dist.mat))
}
