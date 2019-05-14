#' Draw a parameter from a folded Normal distribution
#'
#' This function draws a single value from a folded Normal distribution with mean the current value of the parameter,
#'  and variance as required.
#'
#' @param param.cur The current value of the parameter.
#' @param param.sigma The standard deviation of the folded Normal distribution.
#' @param lower The lower bound of the distribution.
#' @param upper The upper bound of the distribution.
#'
#' @keywords folded Normal distribution proposal
#' @export
#'
#' @return Returns a draw from the specified folded Normal distribution.
#'
#' @examples
#' Folded_draw(param.cur = 0.4, param.sigma = 0.1, lower = 0, upper = 1)

Folded_draw <- function(param.cur, param.sigma, lower, upper){
  
  p.draw <- rnorm(1, param.cur, param.sigma)
  
  while(p.draw > upper | p.draw < lower){
    if(p.draw > upper){p.draw <- p.draw - 2*(p.draw - upper)}
    if(p.draw < lower){ p.draw <- p.draw + 2*(lower - p.draw)}
  }
  
  return(p.draw)
}