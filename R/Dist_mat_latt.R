
#' Lattice population distribution (with buddy).
#'
#' This function generates a distance matrix for individuals on a lattice, across a specified rectangular plane. 
#'     There is also the option to generate a "buddy" for each individual.
#'
#' @param N The total number of individuals in the population.
#' @param nrow The number of rows of the lattice. N must be divisible by nrow.
#' @param xlim The width of the plane on which individuals will be generated (defaults to 20 units wide).
#' @param ylim The height of the place on which individuals will be generated (defaults to 20 units high).
#' @param buddy Dicatates whether to generate a "buddy" for each point, generated from an isotropic
#'  distribution around each point with variance (xlim/200, ylim/200).
#'
#' @keywords Lattice population distribution distance matrix clustered
#' @export
#'
#' @return The function returns a list. The first object in the list is a Nx2 matrix of the coordinates of the individuals.
#'   The second object is an NxN distance matrix.
#'
#' @examples
#' xy.coords <- Dist_mat_latt(N = 200, nrow = 10, xlim = 20, ylim = 20, buddy = TRUE)[[1]]
#' distance_mat <- Dist_mat_latt(N = 200, nrow = 10, xlim = 20, ylim = 20, buddy = TRUE)[[2]]
#' plot(xy.coords[,1], xy.coords[,2], type  = "p")
#' 

Dist_mat_latt <- function(N, nrow, xlim = 20, ylim = 20, buddy = TRUE){
  
  if(buddy == TRUE){
    N.lat  <- N/2
  }else{N.lat = N}
  
  n_per_col <- N.lat/nrow 
  
  lat.x.coords <- rep(seq(2,xlim,length.out = nrow), each = n_per_col)
  lat.y.coords <- rep(seq(2, xlim, length.out = n_per_col), times = nrow)
  lattice.coords = cbind(lat.x.coords, lat.y.coords)
  
  if(buddy == TRUE){
    
    library(MASS)
    
    std.dev <- matrix(c((xlim/200),0,0,(ylim/200)),2,2)
  
    child.coords = NA
    xy.coords = matrix(nrow =  N, ncol = 2)
    xy.coords[1:N.lat, ] <- lattice.coords 
    for(j in (N.lat+1):N){
      child.coords <- MASS::mvrnorm(n = 1, lattice.coords[(j - N.lat),], std.dev)
      xy.coords[j,] <- child.coords
    }
  }else{
    xy.coords <- lattice.coords 
  }

  dist.mat <- as.matrix(dist(xy.coords, upper = T))
  
  return(list(xy.coords, dist.mat))
}
