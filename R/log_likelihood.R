#' Calculate the log-likelihood for an epidemic, up to a constant of proportionality.
#'
#' This function is used tp calculate the likelihood of a General Stochastic Epidemic, and brings together all the previously 
#'  specified functions; prod_part, interval_intersect, and integral_part.
#'
#' @param inf_times A vector of the infection times of all individuals (Inf if not infected).
#' @param rem_times A vector of the removal times of all infected individuals.
#' @param B The infection rate matrix.
#'
#' @keywords likelihood log loglikelihood infection GSE 
#' @export
#'
#' @return Returns the value of the log likelihood, up to a constant of proportionality.
#'
#' @examples
#' This function is utilised in the inference functions.

log_likelihood = function(inf_times, rem_times, B){
  infected_inds <- inf_times < Inf
  prod = prod_part(inf_times, rem_times, B, infected_inds)
  integral = integral_part(inf_times, rem_times, B, infected_inds)
  prod - integral
}