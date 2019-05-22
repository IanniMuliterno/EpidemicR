
#' Clustered edge population distribution.
#'
#' This function generates a distance matrix for a population that is clustered on the left and right edges of the space only.
#'
#' @param N The total number of individuals in the population.
#' @param xlim The width of the plane on which individuals will be generated (defaults to 20 units wide).
#' @param ylim The height of the place on which individuals will be generated (defaults to 20 units high).
#' @param clusters The total number of clusters. Half the clusters will be on the right edge, the other half on the left.
#'                  (Must be even, greater than 2, and divide N). Clusters are generated on blocks of size (xlim/10, ylim/10).
#'
#' @keywords Clustered edge population distribution distance matrix
#' @export
#'
#' @return The function returns a list. The first object in the list is a Nx2 matrix of the coordinates of the individuals.
#'   The second object is an NxN distance matrix.
#'
#' @examples
#' xy.coords <- Dist_mat_edge(N=100, xlim = 20, ylim = 20, clusters = 10)[[1]]
#' distance_mat <- Dist_mat_edge(N=100, xlim = 20, ylim = 20, clusters = 10)[[2]]
#' plot(xy.coords[,1], xy.coords[,2], type  = "p")
#' 

Dist_mat_edge <- function(N, xlim = 20, ylim = 20, clusters = 8){
  
  N.clust <- N/clusters # Number of individuals per cluster
  n.side.clust <- clusters/2 # Number of clusters on each side
  n.gaps <- n.side.clust - 1 # Number of gaps on each side
  
  clust.width <- xlim/10 # Cluster width
  clust.height <- ylim/10 # Cluster height
  
  gaps.y <- (ylim - clust.height*n.side.clust)/n.gaps # Height of gaps
  
  clust.start <- 0
  for(i in 1:n.gaps){clust.start[i+1] <- i*(clust.height + gaps.y)}
  
  Left.y <- NA
  Right.y <- NA
  for(i in 1:n.side.clust){
    Left.y[1:N.clust + (i-1)*N.clust]<- runif(N.clust, clust.start[i], (clust.start[i] + clust.height))
    Right.y[1:N.clust + (i-1)*N.clust]<- runif(N.clust, clust.start[i], (clust.start[i] + clust.height))
  }
  
  Left.x <- runif(n.side.clust*N.clust, clust.start[1], clust.start[1]+clust.width)
  Right.x <- runif(n.side.clust*N.clust, clust.start[n.side.clust], clust.start[n.side.clust]+clust.width)
  
  xy.coords <- cbind(c(Left.x, Right.x), c(Left.y, Right.y))
  
  dist.mat <- as.matrix(dist(xy.coords, upper = T))
  
  return(list(xy.coords, dist.mat))
}































