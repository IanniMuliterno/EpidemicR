
#' An epidemic simulation
#'
#' This function allows you to simulate a stochastic SIR epidemic, with a heterogeneous contact matrix.
#' 
#' The infectious life history of each individual is comprised of two parts;
#' \itemize{
#'  \item \eqn{Q_i}, the length of individuals \eqn{i}'s infectious period.
#'  \item \eqn{W_{i,j}}, the length of time, after individual \eqn{i} becomes infected, before they make contact with individual \eqn{j}.
#' }
#' The distributions of which are given by;
#' \itemize{
#'  \item \eqn{Q_i \overset{iid}{\sim}}{Q_i ~} Exp(\eqn{\gamma}{gamma}).
#'  \item \eqn{W_{i,j} \overset{iid}{\sim}}{W_{i,j} ~} Exp(\eqn{\beta_{i,j}}{beta_{i,j}}).
#' }
#' The simulation begins with an initial infective infectious at time 0. 
#' 
#' \eqn{I_i} is the time at which individual \eqn{i} becomes 
#' infected, and their removal time is given by \eqn{R_i = I_i + Q_i}. 
#' 
#' If \eqn{W_{i,j} < Q_i}, then individual \eqn{i} makes infectious contact with individual \eqn{j} at time \eqn{I_i + W_{i,j}}. 
#' 
#' If individual \eqn{j} is susceptible when this happens, then they become infected, otherwise nothing happens. 
#' The simulation assumes that \eqn{W_{i,j} != W_{j,i}}.
#' 
#' The simulation continues until there are no infected individuals left in the population.
#'
#' @param N The size of the population.
#' @param beta.mat An NxN matrix of the pair wise infection rates between all individuals.
#' @param gamma The removal rate.
#' @keywords General Stochastic Epidemic GSE
#' @return The function returns a matrix containing 4 columns; the id number of each individual, their infection time (NA if not infected),
#'  their removal time (NA if not infected), and their generated infectious period length.
#' @export
#' @examples
#' 
#' # Generate a distance matrix
#'    distance_mat <- Dist_mat_unif(N=100, xlim = 20, ylim = 20)[[2]]
#' 
#' # Generate an associated infection rate rate matrix
#'    rate_mat <- Beta_mat_form(distance_mat, c(0.004, 0.002), 10)
#'    
#' # Generate a simulated epidemic
#'    Hetero_sim <- GSE_sim(N = 100, beta.mat = rate_mat, gamma = 0.15)

GSE_sim <- function(N, beta.mat, gamma){
  
  ### Initialise the "Time until contact" matrix ###
  
      # Is one or more of the infection rates 0?
        whichzero <- which(beta.mat == 0)
        
        madness = 0
        while(madness == 0){
          
          suppressWarnings(
            W <- rexp(N^2, beta.mat) 
          )
          W[is.nan(W)] <- Inf
          Wmat = matrix(W, N, N, byrow = T) # "Time until contact" matrix
          
          # Matrix has;
                    # Infected individuals indexed on the rows
                    # Susceptible individuals indexed on the columns
          
          uniquecount <- N*N - length(whichzero) + 1 # No. of elements - Number of zeros in rate matrix + 1 to count elements with value 0
          if(length(unique(as.vector(Wmat))) == uniquecount){
            madness = 1
          } # If all the non-diagonal values of the matrix are unique (as would be expected in reality, 
            # and is necessary for the functioning of the code), then accept this matrix as an initial value
          
        } # End of while loop
          

  
  ### Initialise the "Infectious period" matrix ### 
      
      Q = rexp(N, gamma) 
      Qmat <- matrix(rep(Q, N), N, N) # "Infectious period" matrix
      # Matrix has;
              # Top row: infectious period of first indiviual (one number repeated)
              # Second row: infectious period of second individual (")
              # and so on...
   
  ### Initialise the "Infection times" vector ### 
      
      inf.times.full <- c( 0, rep(Inf, (N-1)) ) # Vector of the true infection times
                                # Initially only know infection time of intital infective
      
  ### Initialise functional objects ### 
      
      pop = seq(1,N) # Vector of population ids
      infect = c(1) # Vector of the previously infected individuals
      sus = pop[-infect] # Vector of the currently susceptible individuals
  
      
  ### The simulation ###
      
      STOP =  0
      
      while(STOP == 0){
        
        ### Calculate the current submatrices that are relevant ###
        
            Wcur <- Wmat[infect, sus, drop = F] # "Time, after infection, until contact" matrix
            Qcur <- Qmat[infect, sus, drop=F] # "Infectious period" matrix
            
        ### Calculate which individuals can become infected ###
        
            within <- Wcur<Qcur # Which pairs have contact within the infectious period of the infected individual
            which.inf <- which(within == T, arr.ind = T)
            # A new infection occurs only if Q>W for atleast one pair
            # This also ensure that it is not a recovered individual infecting a new person
            
            track <- rbind(infect[which.inf[,1]], sus[which.inf[,2]]) # A table to track who gets infected
            
        ### Calculate the infection times of all the potential new infectees ###
            
            inf.times <- inf.times.full[track[1, ]] + Wcur[which.inf]
        
        ### Recover the index of the newly infected individual ###
        
            #The individual who becomes infected next will be the one with the smallest infection time
            suppressWarnings(
              newinftime <- min(inf.times)
            )
            up <- which.min(inf.times) #index of smallest infection time
            newinf <- track[,up][2]
        
        ### Update population ###
            
            #If a new person was infected then
            if(is.na(newinf) == F){
              
              inf.times.full[newinf] <- newinftime #update infection times
              
              infect <- c(infect, newinf) #update infectious statuses
              sus <- pop[-infect]
              
            #If no one new was infected then  
            }else{
              STOP = 1
            }
        
      } #end of while loop
      
  ### Record the results ###
      
      rem.times.full <- inf.times.full + Q # Removal times
      
      res <- cbind(1:N,    # Individual ids ordered by infection status
                   inf.times.full,    # Infection times, with NA for those not infected
                   rem.times.full,    # Removal times, with NA for those not infected
                   Q)                 # Infectious periods
      
      colnames(res) <- c("id", "inf.time", "rem.time", "inf.period")
  
  return(res)
}
