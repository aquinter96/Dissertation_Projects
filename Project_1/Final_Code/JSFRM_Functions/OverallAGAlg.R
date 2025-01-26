## OverallAGAlg.R

OverallAGAlg <- function(Data, Best, tuningpA = seq(0, 15, 0.1), m_seq = 1:4){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of A (and corresponding estimates
  ## of A0/Gamma/Phi2/Phi3) and corresponding BIC for each # of latent factors for Y in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.) Also extract the optimal estimates
  ## calculated in OverallBAlg.R of B/B0/Phi1 for use in the estimation of A/A0/Gamma/Phi2/Phi3
  #################################################################################
  
  Aestlist <- replicate(length(m_seq), NULL, simplify=F)
  Alist <- replicate(length(tuningpA), list(), simplify=F)
  Aest <- NULL
  A_final <- NULL
  All_final <- list()
  BICs <- rep(0, length(tuningpA))
  BIC_opt <- rep(0, length(m_seq))
  
  #################################################################################
  ## for each # of latent factors in the grid search, determine the optimal estimate
  ## of A/A0/GammaPhi2/Phi3, then save the corresponding estimates and BIC in the respective lists
  #################################################################################
  
  for(i in 1:length(m_seq)){
    #################################################################################
    ## for each # of latent factors in the grid search, determine the optimal estimate
    ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
    #################################################################################
    
    Ainit <- A_inits(Data, Best, m_seq[i])
    Alist <- parLapply(cl, tuningpA, EMAlgAGammaAdLassoCV, Data, Best, Ainit, Ainit$A)

    for(j in 1:length(tuningpA)){
      BICs[j] <- Alist[[j]]$BIC
    }
    
    BIC_opt[i] <- min(BICs)
    
    Aestlist[[i]] <-Alist[[min(which(BICs == BIC_opt[i]))]]
    
  }

  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  Aopt <- Aestlist[[which(BIC_opt == min(BIC_opt))]]
  
  A_final$A_final_pars <- Aopt$est_Model_param
  A_final$m_opt <- ncol(Aopt$est_Model_param$A)
  A_final$A_log.Lik <- Aopt$log.Lik
  A_final$diffList <- Aopt$diffList
  A_final$A_BIC <- Aopt$BIC
  A_final$tuningpA <- Aopt$tuningpA
  A_final$A_time_diff <- Aopt$time.diff

  All_final$B_estimates <- Best
  All_final$A_estimates <- A_final
  
  return(All_final)
}

## end of code