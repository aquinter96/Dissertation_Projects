## OverallBAlg.R

OverallXAlg <- function(X, tuningpB = seq(0, 15, 0.1), s_seq = 1:4){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  Bestlist <- replicate(length(s_seq), NULL, simplify=F)
  Blist <- replicate(length(tuningpB), list(), simplify=F)
  Best <- NULL
  X_final <- NULL
  BICs <- rep(0, length(tuningpB))
  BIC_opt <- rep(0, length(s_seq))
  
  for(i in 1:length(s_seq)){
    #################################################################################
    ## for each # of latent factors in the grid search, determine the optimal estimate
    ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
    #################################################################################
    
    Binit <- B_inits(X, s_seq[i])
    Blist <- parLapply(cl, tuningpB, Singular_ErrorX, X, Binit, Binit$B)
    
    #Blist <- parLapply(cl, tuningpB, EMAlgBAdLassoCV, X, Binit, Binit$B)
    Blist[sapply(Blist, is.null)] <- NULL
    length(BICs) <- length(Blist)
    
    if(length(Blist) >= 1){
      for(j in 1:length(Blist)){
        BICs[j] <- Blist[[j]]$BIC
      }
      
      Blist <- Blist[which(!is.nan(BICs))]
      BICs <- BICs[!is.nan(BICs)]
      
      BIC_opt[i] <- min(BICs)
      
      Bestlist[[i]] <- Blist[[min(which(BICs == BIC_opt[i]))]]
    }
    else{
      BIC_opt[i] <- NA
      Bestlist[[i]] <- NULL
    }
  }
  
  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  Bopt <- Bestlist[[which(BIC_opt == min(BIC_opt, na.rm = T))]]
  
  X_final$X_final_pars <- Bopt$est_Model_param
  X_final$s_opt <- ncol(Bopt$est_Model_param$B)
  X_final$X_log.Lik <- Bopt$log.Lik
  X_final$diffList <- Bopt$diffList
  X_final$X_BIC <- Bopt$BIC
  X_final$tuningpB <- Bopt$tuningpB
  X_final$X_time_diff <- Bopt$time.diff
  
  return(X_final)
}

## end of code