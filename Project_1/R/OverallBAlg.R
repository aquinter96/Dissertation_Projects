#' Wrapper function for full X matrix estimation
#' 
#' The function is a wrapper for the EMAlgBAdLassoCV
#' function with the additional layer of estimating
#' the optimal number of latent factors for X using grid
#' search. Specifically, the optimal model for each
#' input number of latent factors is estimated and saved,
#' after which all estimates are compared by BIC and
#' the estimate of s and corresponding model estimates
#' are returned
#' 
#' @param X The X dataset
#' @param tuningpB A vector of numbers containing the candidate tuning
#' parameter values to be compared via grid search for all models fit
#' @param s_seq A vector of positive integers containing the different
#' numbers of latent factors for X that are to be compared
#' @param nfolds Deprecated parameter.
#' 
#' @return The optimal number of latent parameters estimated and the
#' corresponding optimal model estimates, BIC, and the optimal tuning
#' value chosen for said optimal model
#' 
#' @examples
#' OberallBAlg(data$X, seq(0, 5, 0.1), 1:4)
#' 
#' 
#' 
#' @export
OverallBAlg <- function(X, tuningpB = seq(0,15,0.1), s_seq = 1:4, nfolds=10){
  Blist <- replicate(length(s_seq),NA,simplify=F)
  BIClist <- rep(0, length(s_seq))
  for(i in 1:length(s_seq)){
    Binit <- EMAlgB(X, s_seq[i])
    Blist[[i]] <- EMAlgBAdLassoCV(X, s_seq[i], Binit$B, Binit$Phi1, tuningpB, weights = Binit$B)
    BIClist[i] <- Blist[[i]]$BICopt
  }
  sopt <- s_seq[which(BIClist == min(BIClist))]
  optB <- Blist[[which(BIClist == min(BIClist))]]
  return(list("Bopt" = optB, "optimal s" = sopt, "BIC" = optB$BICopt, "lambda" = optB$`optimal lambda`))
}
