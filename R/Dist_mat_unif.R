
#' Uniform population distribution.
#'
#' This function generates a distance matrix with uniform distribution of individuals across a specified rectangle.
#'
#' @param N The total number of individuals in the population.
#' @param xlim The width of the plane on which individuals will be generated (defaults to 20 units wide).
#' @param ylim The height of the place on which individuals will be generated (defaults to 20 units high).
#'
#' @keywords Uniform population distribution distance matrix
#' @export
#'
#' @return The function returns a list. The first object in the list is a Nx2 matrix of the coordinates of the individuals.
#'   The second object is an NxN distance matrix.
#'
#' @examples
#' xy.coords <- Dist_mat_unif(N=100, xlim = 20, ylim = 20)[[1]]
#' distance_mat <- Dist_mat_unif(N=100, xlim = 20, ylim = 20)[[2]]
#' plot(xy.coords[,1], xy.coords[,2], type  = "p")
#' 

Dist_mat_unif <- function(N, xlim = 20, ylim = 20){
  
  x.coords = runif(N, 0, xlim)
  y.coords = runif(N, 0, ylim)
  
  xy.coords <- cbind(x.coords, y.coords)
  
  dist.mat <- as.matrix(dist(cbind(x.coords, y.coords), upper = T))
  
  return(list(xy.coords, dist.mat))
}
