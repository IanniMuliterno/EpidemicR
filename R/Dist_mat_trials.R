
#' Multiple lines population distribution.
#'
#' This function generates a distance matrix where the population is arranged on seperate lines
#'  where the gaps between rows increase by one unit with each successive individual. The multiple lines are far apart and 
#'  can be thought of as multiple trials of a 2D epidemic, if the first individual in every row is an initial infective.
#'
#' @param N The total number of individuals in the population.
#' @param nline The number of lines(trials). Must divide N.
#'
#' @keywords Uniform population distribution distance matrix
#' @export
#'
#' @return The function returns a list. The first object in the list is a Nx2 matrix of the coordinates of the individuals.
#'   The second object is an NxN distance matrix.
#'
#' @examples
#' xy.coords <- Dist_mat_trials(N=100, nline = 10)[[1]]
#' distance_mat <- Dist_mat_trials(N=100, nline = 10)[[2]]
#' plot(xy.coords[,1], xy.coords[,2], type  = "p")
#' 

Dist_mat_trials <- function(N, nline){
  
  perrow <- N/nline
  
  line.coords <- 1
  for(i in 2:perrow){
    line.coords[i] <- line.coords[i-1] + (i-1)
  }
  line.coords <- line.coords/10
  
  linegap <- 2*line.coords[perrow]
  
  x.coords <- rep(line.coords, each = nline)
  y.coords <- rep(1, nline) + seq(0, (nline-1))*linegap
  
  xy.coords <- cbind(x.coords, y.coords)
  
  dist.mat <- as.matrix(dist(cbind(x.coords, y.coords), upper = T))
  
  return(list(xy.coords, dist.mat))
}
