#' Describe the acceptance probabilities.
#'
#' This function takes the results from an inference function, and describes in words the acceptance probabilities of the parameters.
#'  Will show 0% in the case where the parameter is fixed.
#'
#' @param res.list The list of results outputted by the MCMC functions.
#' @param reparam Whether we are using the reparameterisation to make inference or not (defaults to FALSE).
#'
#' @keywords descirption results acceptance probabilities
#' @export
#'
#' @return A description of the average acceptance probabilities.
#'
#' @examples
#' 
#' Avg_acc_probs(res.list, reparam = F)
#' 
#' 


Avg_acc_probs <- function(res.list, reparam = F){
  
  cat(c("The average acceptance probability for the infection times was ", res.list[[2]] * 100, "%.\n"), sep = "")
  
  if(reparam == F){
    cat(c("The average acceptance probability for beta_1 was ", res.list[[3]] * 100, "%.\n"), sep = "")
    cat(c("The average acceptance probability for beta_2 was ", res.list[[4]] * 100, "%.\n"), sep = "")
    cat(c("The average acceptance probability for d was ", res.list[[5]] * 100, "%.\n\n"), sep = "")
  }else{
    cat(c("The average acceptance probability for p was ", res.list[[3]] * 100, "%.\n"), sep = "")
    cat(c("The average acceptance probability for d was ", res.list[[4]] * 100, "%.\n\n"), sep = "")
  }
  
}