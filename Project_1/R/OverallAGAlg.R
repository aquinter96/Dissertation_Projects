#' Wrapper function for full Y matrix estimation
#' 
#' The function is a wrapper for the EMAlgAGammaAdLassoCV
#' function with the additional layer of estimating
#' the optimal number of latent factors for X using grid
#' search. Specifically, the optimal model for each
#' input number of latent factors is estimated and saved,
#' after which all estimates are compared by BIC and
#' the estimate of s and corresponding model estimates
#' are returned
#' 
#' @param X The X dataset
#' @param Y The Y dataset
#' @param Bobj The direct output from OverakllBAlg, containing the
#' various estimates from the X model needed for estimation in the Y model
#' @param tuningpA A vector of numbers containing the candidate tuning
#' parameter values to be compared via grid search for all models fit
#' @param m_seq A vector of positive integers containing the different
#' numbers of latent factors for Y that are to be compared
#' @param nfolds Deprecated parameter.
#' 
#' @return The optimal number of latent parameters estimated and the
#' corresponding optimal model estimates, BIC, the optimal tuning
#' value chosen for said optimal model, and the various estimates
#' from Bobj
#' 
#' @examples
#' OberallAGAlg(data$X,, data$Y, OptBests, seq(0, 5, 0.1), 1:4)
#' 
#' 
#' 
#' @export
OverallAGAlg <- function(X, Y, Bobj,tuningpA = seq(0,15,0.1), m_seq = 1:4, nfolds=10){
  Ainit <- NA
  s <- Bobj$`optimal s`
  B <- Bobj$Bopt$B
  B0 <- Bobj$Bopt$B0
  Phi1 <- Bobj$Bopt$Phi1
  Alist <- replicate(length(m_seq),NA,simplify=F)
  BIClist <- rep(0, length(m_seq))
  for(i in 1:length(m_seq)){
    inits <- initvalcalc(X, Y, B, s, m_seq[i], nrow(Y))
    Aweight <- EMAlgAGammaAdLassoCV(X, Y, s, m_seq[i], B, B0, Phi1, inits$Ainit, inits$Ginit, diag(diag(inits$Phi2init)), inits$Phi3init, 0, weights = matrix(rep(0,ncol(Y)*m_seq[i]),ncol=m_seq[i]))
    Alist[[i]] <- EMAlgAGammaAdLassoCV(X, Y, s, m_seq[i], B, B0, Phi1, inits$Ainit, inits$Ginit, diag(diag(inits$Phi2init)), inits$Phi3init, tuningpA, weights = Aweight$A)
    BIClist[i] <- Alist[[i]]$BICopt
  }
  mopt <- m_seq[which(BIClist == min(BIClist))]
  optA <- Alist[[which(BIClist == min(BIClist))]]
  return(list("model results" = optA, "X" = X, "B" = B, "B0" = B0, "Phi1" = Phi1, "optimal s BIC" = s, "optimal m BIC" = mopt, "ABIC" = optA$BICopt, "lambdaA" = optA$`optimal lambda`, "BBIC" = Bobj$BIC, "lambdaB" = Bobj$lambda))
}
