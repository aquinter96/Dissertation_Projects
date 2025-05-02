## OverallAGAlg.R

OverallMAlg <- function(Data, Xest, tuningpA = seq(0, 15, 0.1), m_seq = 1:4){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of A (and corresponding estimates
  ## of A0/Gamma/Phi2/Phi3) and corresponding BIC for each # of latent factors for Y in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.) Also extract the optimal estimates
  ## calculated in OverallBAlg.R of B/B0/Phi1 for use in the estimation of A/A0/Gamma/Phi2/Phi3
  #################################################################################
  
  Aestlist <- replicate(length(m_seq), NULL, simplify=F)
  Alist <- replicate(length(tuningpA), list(), simplify=F)
  Aest <- NULL
  M_final <- NULL
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
    
    Ainit <- A_inits(Data, Xest, m_seq[i])
    Alist <- parLapply(cl, tuningpA, Singular_ErrorM, Data, Xest, Ainit, Ainit$A)
    # Alist <- parLapply(cl, tuningpA, EMAlgAGammaAdLassoCV, Data, Best, Ainit, Ainit$A)
    
    Alist[sapply(Alist, is.null)] <- NULL
    length(BICs) <- length(Alist)
    
    if(length(Alist) >= 1){
      for(j in 1:length(Alist)){
        BICs[j] <- Alist[[j]]$BIC
      }
      
      Alist <- Alist[which(!is.nan(BICs))]
      BICs <- BICs[!is.nan(BICs)]
      
      BIC_opt[i] <- min(BICs)
      
      Aestlist[[i]] <- Alist[[min(which(BICs == BIC_opt[i]))]]
    }
    else{
      BIC_opt[i] <- NA
      Aestlist[[i]] <- NULL
    }
    
  }

  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  Aopt <- Aestlist[[which(BIC_opt == min(BIC_opt, na.rm = T))]]
  
  M_final$M_final_pars <- Aopt$est_Model_param
  M_final$m_opt <- ncol(Aopt$est_Model_param$A)
  M_final$M_log.Lik <- Aopt$log.Lik
  M_final$diffList <- Aopt$diffList
  M_final$M_BIC <- Aopt$BIC
  M_final$tuningpA <- Aopt$tuningpA
  M_final$M_time_diff <- Aopt$time.diff

  All_final$X_estimates <- Xest
  All_final$M_estimates <- M_final
  
  return(All_final)
}

## end of code