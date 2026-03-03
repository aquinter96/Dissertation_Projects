## OverallAGAlg.R

OverallYAlg <- function(Data, Mest, t_seq = 1:4){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of A (and corresponding estimates
  ## of A0/Gamma/Phi2/Phi3) and corresponding BIC for each # of latent factors for Y in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.) Also extract the optimal estimates
  ## calculated in OverallBAlg.R of B/B0/Phi1 for use in the estimation of A/A0/Gamma/Phi2/Phi3
  #################################################################################
  
  Yestlist <- replicate(length(t_seq), NULL, simplify=F)
  Yest <- NULL
  Y_final <- NULL
  All_final <- list()
  BIC_opt <- rep(0, length(t_seq))
  tuningvals <- rep(0, 2)
  
  #################################################################################
  ## for each # of latent factors in the grid search, determine the optimal estimate
  ## of A/A0/GammaPhi2/Phi3, then save the corresponding estimates and BIC in the respective lists
  #################################################################################
  
    #################################################################################
    ## for each # of latent factors in the grid search, determine the optimal estimate
    ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
    #################################################################################
  
  Yinits <- lapply(t_seq, NMWrapperYinit, Data, Mest)
  Yinits <- Yinits[!sapply(Yinits,is.null)]
  Yestlist <- lapply(Yinits, NMWrapperY, Data, Mest, tuningvals)
  Yestlist <- Yestlist[!sapply(Yestlist,is.null)]
  
  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################

  if(length(Yestlist) == 0){
    return(print("All models failed to converge"))
  }else{
    
  BIC_opt <- unlist(lapply(Yestlist, function(x)x$value))
  Yinit_opt <- Yinits[[which(BIC_opt == min(BIC_opt, na.rm = T))]]
  tune_opt <- Yestlist[[which(BIC_opt == min(BIC_opt, na.rm = T))]]$par
  
  Yopt <- EMAlgY(tune_opt, Data, Mest, Yinit_opt, NM = F)
    
  Y_final$Y_final_pars <- Yopt$est_Model_param
  Y_final$t_opt <- ncol(Yopt$est_Model_param$C)
  Y_final$Y_log.Lik <- Yopt$log.Lik
  Y_final$diffList <- Yopt$diffList
  Y_final$Y_BIC <- Yopt$BIC
  Y_final$tuningpC <- Yopt$tuningpC
  Y_final$tuningpD <- Yopt$tuningpD
  Y_final$Y_time_diff <- Yopt$time.diff
  
  All_final$X_estimates <- Mest$X_estimates
  All_final$M_estimates <- Mest$M_estimates
  All_final$Y_estimates <- Y_final
  
  return(All_final)
  }
}

## end of code