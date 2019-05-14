#' Calculate the interval intersect.
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates the interval
#'  intersect used in the calculation of the intergral for the infection likelihood. 
#'  This is the interval \eqn{[ (R_i ^ I_j) - (I_i) ^ (I_j) ]} for all i and j.
#'
#' @param interval_i Is the matrix cbind(t_inf, t_rem)[infected, ]. This is a 2 column matrix of the infection ]
#'                   and removal times for the infected individuals.
#' @param interval_j Is the matrix cbind(0, t_inf). This is a 2 column matrix of 0s in the first column and
#'                   the infection times in the second.
#'
#' @keywords interval intersection infection GSE likelihood
#' @export
#'
#' @return Returns a matrix where each column is a vector for each j.
#'
#' @examples
#' This function is utilised by the log_likelihood function.

# == Components ==
# int_start <- calculates the maximum of 0 and the infection time of infected individual i
#              since i can only infect j after it has become infected itself
#              this is (0 max I_i)
# int_end <- calculates the minimum of the infection times of the every individual j in the population
#            and the removal times of the each infected individual
#            this is (I_j min R_i)

# == Functions ==
# "pmax()" calculates the "parallel" maximum of two vectors ie. pmax(c(1,2,3), c(3,2,1)) = c(3,2,3)
# similarly for "pmin()"

interval_intersect = function(inf_times, rem_times, infected_inds){
  interval_j <- cbind(inf_times, rem_times)
  interval_i <- interval_j[infected_inds, ]
  
  int_start <- sapply(interval_j[,1], function(x) pmin(x, interval_i[,1]))
  int_end <- sapply(interval_j[,1], function(x) pmin(x, interval_i[,2]))
  int_end - int_start
  #returns a matrix where each column is a vector for each j
}