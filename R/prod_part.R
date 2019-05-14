#' Calculate the product part of the infection likelihood
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates 
#'  \eqn{\prod_{j != kappa}^{n_I} [ \sum_{ i in I_{n_{j-}} } [ beta_{i,j} ] ]}.
#'
#' @param inf_times A vector of the infection times of all individuals (Inf if not infected).
#' @param rem_times A vector of the removal times of all infected individuals.
#' @param B The infection rate matrix.
#' @param infected_inds A vector of the infected individuals.
#'
#' @keywords Product infection GSE likelihood
#' @export
#'
#' @return Returns the value of the product detailed above.
#'
#' @examples
#' This function is utilised by the log_likelihood function.

# == Components ==
# is_infected = 0/1 vector which says if each individual became infected or not.
# waifw = "who acquired infection from whom" matrix. True false matrix where the value of element (i,j) is 1 if $i$ could 
#          have infected $j$ and 0 otherwise. Susceptibles just have vector of 0s.
# lambda_j = vector of the sums of the infectious pressure exererted on each individual who became infected
#          = in other words it is the sum of the beta_ij for i in 
#            {the set of infected individuals just before individual j becomes infected}.
# I0 = the index of the initial infective

prod_part <- function(inf_times, rem_times, B, infected_inds) {
  waifw <- sapply(inf_times, function(t) inf_times < t & t < rem_times)
  lambdaj <- colSums(B[,infected_inds] * waifw[, infected_inds])
  I0 <- which.min(inf_times)
  sum(log(lambdaj[-I0]))
  # this calculates what we called sumlogsum
}