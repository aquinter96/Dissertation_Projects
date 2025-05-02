## OverallAGAlg.R

OverallYAlg <- function(Data, Mest, tuningpC = seq(0, 15, 0.1), tuningpD = seq(0, 15, 0.1), t_seq = 1:4){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of A (and corresponding estimates
  ## of A0/Gamma/Phi2/Phi3) and corresponding BIC for each # of latent factors for Y in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.) Also extract the optimal estimates
  ## calculated in OverallBAlg.R of B/B0/Phi1 for use in the estimation of A/A0/Gamma/Phi2/Phi3
  #################################################################################
  
  tuning_list <- list(tuningpC, tuningpD)
  tuning_grid <- expand.grid(tuning_list)
  tuning_grid <- lapply(seq_len(nrow(tuning_grid)), function(x) tuning_grid[x,])

  Cestlist <- replicate(length(t_seq), NULL, simplify=F)
  Clist <- replicate(length(tuningpC)*length(tuningpD), list(), simplify=F)
  Cest <- NULL
  Y_final <- NULL
  All_final <- list()
  BICs <- rep(0, length(tuningpC)*length(tuningpD))
  BIC_opt <- rep(0, length(t_seq))
  
  #################################################################################
  ## for each # of latent factors in the grid search, determine the optimal estimate
  ## of A/A0/GammaPhi2/Phi3, then save the corresponding estimates and BIC in the respective lists
  #################################################################################
  
  for(i in 1:length(t_seq)){
    #################################################################################
    ## for each # of latent factors in the grid search, determine the optimal estimate
    ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
    #################################################################################
    Cinit <- C_inits(Data, Mest, t_seq[i])
    Clist <- parLapply(cl, tuning_grid, EMAlgY, Data, Mest, Cinit)
    #Clist <- parLapply(cl, tuning_grid, EMAlgY, Data, Model_Params, Model_Params)

    Clist[sapply(Clist, is.null)] <- NULL
    length(BICs) <- length(Clist)

    if(length(Clist) >= 1){
      for(j in 1:length(Clist)){
        BICs[j] <- Clist[[j]]$BIC
      }
      
      Clist <- Clist[which(!is.nan(BICs))]
      BICs <- BICs[!is.nan(BICs)]
      
      BIC_opt[i] <- min(BICs)
      
      Cestlist[[i]] <- Clist[[min(which(BICs == BIC_opt[i]))]]
    }
    else{
      BIC_opt[i] <- NA
      Cestlist[[i]] <- NULL
    }
  }
  
  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################

  #Copt <- Cestlist[[which(BIC_opt == min(BIC_opt, na.rm = T))]]
  Copt <- Cestlist[[1]]
  
  Y_final$Y_final_pars <- Copt$est_Model_param
  Y_final$m_opt <- ncol(Copt$est_Model_param$C)
  Y_final$Y_log.Lik <- Copt$log.Lik
  Y_final$diffList <- Copt$diffList
  Y_final$Y_BIC <- Copt$BIC
  Y_final$tuningpC <- Copt$tuningpC
  Y_final$tuningpD <- Copt$tuningpD
  Y_final$Y_time_diff <- Copt$time.diff
  
  All_final$X_estimates <- Mest$X_estimates
  All_final$M_estimates <- Mest$M_estimates
  All_final$Y_estimates <- Y_final
  
  return(All_final)
}

## end of code