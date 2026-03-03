## OverallAGAlg.R

OverallMAlg <- function(Data, Xest, m_seq = 1:4){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of A (and corresponding estimates
  ## of A0/Gamma/Phi2/Phi3) and corresponding BIC for each # of latent factors for Y in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.) Also extract the optimal estimates
  ## calculated in OverallBAlg.R of B/B0/Phi1 for use in the estimation of A/A0/Gamma/Phi2/Phi3
  #################################################################################
  
  Mestlist <- replicate(length(m_seq), NULL, simplify=F)
  Mest <- NULL
  M_final <- NULL
  All_final <- list()
  BIC_opt <- rep(0, length(m_seq))
  tuningvals <- 0
  
  #################################################################################
  ## for each # of latent factors in the grid search, determine the optimal estimate
  ## of A/A0/GammaPhi2/Phi3, then save the corresponding estimates and BIC in the respective lists
  #################################################################################
  
    #################################################################################
    ## for each # of latent factors in the grid search, determine the optimal estimate
    ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
    #################################################################################
    
    Minits <- lapply(m_seq, NMWrapperMinit, Data, Xest)
    Minits <- Minits[!sapply(Minits,is.null)]
    Mestlist <- lapply(Minits, NMWrapperM, Data, Xest, tuningvals)
    Mestlist <- Mestlist[!sapply(Mestlist,is.null)]

  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  if(length(Mestlist) == 0){
    return(print("All models failed to converge"))
  }else{
    
  BIC_opt <- unlist(lapply(Mestlist, function(x)x$value))
  Minit_opt <- Minits[[which(BIC_opt == min(BIC_opt, na.rm = T))]]
  tune_opt <- Mestlist[[which(BIC_opt == min(BIC_opt, na.rm = T))]]$par

  Mopt <- EMAlgM(tune_opt, Data, Xest, Minit_opt, NM = F)
  
  M_final$M_final_pars <- Mopt$est_Model_param
  M_final$m_opt <- ncol(Mopt$est_Model_param$A)
  M_final$M_log.Lik <- Mopt$log.Lik
  M_final$diffList <- Mopt$diffList
  M_final$M_BIC <- Mopt$BIC
  M_final$tuningpA <- Mopt$tuningpA
  M_final$M_time_diff <- Mopt$time.diff

  All_final$X_estimates <- Xest
  All_final$M_estimates <- M_final
  
  return(All_final)
  }
}

## end of code