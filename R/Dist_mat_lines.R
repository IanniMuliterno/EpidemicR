
#' Spreading lines population distribution.
#'
#' This function generates a distance matrix where the population is arranged on a lattice
#'  where the gaps between rows and columns increases by one unit with each successive individual.
#'
#' @param N The total number of individuals in the population.(Must be a square number).
#'
#' @keywords Uniform population distribution distance matrix
#' @export
#'
#' @return The function returns a list. The first object in the list is a Nx2 matrix of the coordinates of the individuals.
#'   The second object is an NxN distance matrix.
#'
#' @examples
#' xy.coords <- Dist_mat_lines(N=100)[[1]]
#' distance_mat <- Dist_mat_lines(N=100)[[2]]
#' plot(xy.coords[,1], xy.coords[,2], type  = "p")
#' 

Dist_mat_lines <- function(N){
  
  nrow = sqrt(N)
  ncol = N/nrow
  
  line.coords <- 1
  for(i in 2:ncol){
    line.coords[i] <- line.coords[i-1] + (i-1)
  }
  
  x.coords <- rep(line.coords/10, each = nrow)
  y.coords <- rep(line.coords/10, times = nrow)
  
  xy.coords <- cbind(x.coords, y.coords)
  
  dist.mat <- as.matrix(dist(cbind(x.coords, y.coords), upper = T))
  
  return(list(xy.coords, dist.mat))
}
