## OverallBAlg.R

OverallBAlg <- function(X, tuningpB = seq(0,15,0.1), s_seq = 1:4, nfolds=10){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  Blist <- replicate(length(s_seq),NA,simplify=F)
  BIClist <- rep(0, length(s_seq))
  for(i in 1:length(s_seq)){
    
    #################################################################################
    ## for each # of latent factors in the grid search, determine the optimal estimate
    ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
    #################################################################################
    
    Binit <- EMAlgB(X, s_seq[i])
    Blist[[i]] <- EMAlgBAdLassoCV(X, s_seq[i], Binit$B, Binit$Phi1, tuningpB, weights = Binit$B)
    BIClist[i] <- Blist[[i]]$BICopt
  }
  
  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  sopt <- s_seq[which(BIClist == min(BIClist))]
  optB <- Blist[[which(BIClist == min(BIClist))]]
  return(list("Bopt" = optB, "optimal s" = sopt, "BIC" = optB$BICopt, "lambda" = optB$`optimal lambda`))
}

## end of code