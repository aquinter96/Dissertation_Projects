## OverallAGAlg.R

OverallAGAlg <- function(X, Y, Bobj,tuningpA = seq(0,15,0.1), m_seq = 1:4, nfolds=10){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of A (and corresponding estimates
  ## of A0/Gamma/Phi2/Phi3) and corresponding BIC for each # of latent factors for Y in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.) Also extract the optimal estimates
  ## calculated in OverallBAlg.R of B/B0/Phi1 for use in the estimation of A/A0/Gamma/Phi2/Phi3
  #################################################################################
  
  Ainit <- NA
  s <- Bobj$`optimal s`
  B <- Bobj$Bopt$B
  B0 <- Bobj$Bopt$B0
  Phi1 <- Bobj$Bopt$Phi1
  Alist <- replicate(length(m_seq),NA,simplify=F)
  BIClist <- rep(0, length(m_seq))
  
  #################################################################################
  ## for each # of latent factors in the grid search, determine the optimal estimate
  ## of A/A0/GammaPhi2/Phi3, then save the corresponding estimates and BIC in the respective lists
  #################################################################################
  
  for(i in 1:length(m_seq)){
    
    #################################################################################
    ## for each # of latent factors in the grid search, determine the optimal estimate
    ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
    #################################################################################
    
    inits <- initvalcalc(X, Y, B, s, m_seq[i], nrow(Y))
    Aweight <- EMAlgAGammaAdLassoCV(X, Y, s, m_seq[i], B, B0, Phi1, inits$Ainit, inits$Ginit, diag(diag(inits$Phi2init)), inits$Phi3init, 0, weights = matrix(rep(0,ncol(Y)*m_seq[i]),ncol=m_seq[i]))
    Alist[[i]] <- EMAlgAGammaAdLassoCV(X, Y, s, m_seq[i], B, B0, Phi1, inits$Ainit, inits$Ginit, diag(diag(inits$Phi2init)), inits$Phi3init, tuningpA, weights = Aweight$A)
    BIClist[i] <- Alist[[i]]$BICopt
  }
  
  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  mopt <- m_seq[which(BIClist == min(BIClist))]
  optA <- Alist[[which(BIClist == min(BIClist))]]
  return(list("model results" = optA, "X" = X, "B" = B, "B0" = B0, "Phi1" = Phi1, "optimal s BIC" = s, "optimal m BIC" = mopt, "ABIC" = optA$BICopt, "lambdaA" = optA$`optimal lambda`, "BBIC" = Bobj$BIC, "lambdaB" = Bobj$lambda))
}

## end of code