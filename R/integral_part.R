#' Calculate the integral part of the infection likelihood.
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates the integral part
#'  of the infection likelihood, using the alternative double-sum parameterisation.
#'
#' @param inf_times A vector of the infection times of all individuals (Inf if not infected).
#' @param rem_times A vector of the removal times of all infected individuals.
#' @param B The infection rate matrix.
#' @param infected_inds A vector of the infected individuals.
#'
#' @keywords Integral infection GSE likelihood
#' @export
#'
#' @return Returns the value of the integral detailed above.
#'
#' @examples
#' This function is utilised by the log_likelihood function.

# == Components ==
# i_infected = vector of individuals which are infected
# E = matrix where e_ij is defined as in the interval intersect function.
# sum(integral) calculates the double sum but in matrix notation

integral_part = function(inf_times, rem_times, B, infected_inds){
  E = interval_intersect(inf_times, rem_times, infected_inds)
  integral = E * B[infected_inds,]
  sum(integral)
}