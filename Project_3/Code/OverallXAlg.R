## OverallBAlg.R

OverallXAlg <- function(X, s_seq = 1:4){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  Xestlist <- replicate(length(s_seq), NULL, simplify=F)
  Xest <- NULL
  X_final <- NULL
  BIC_opt <- rep(0, length(s_seq))
  tuningvals <- 0
  
    #################################################################################
    ## for each # of latent factors in the grid search, determine the optimal estimate
    ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
    #################################################################################
    
    Xinits <- lapply(s_seq, NMWrapperXinit, X)
    Xinits <- Xinits[!sapply(Xinits,is.null)]
    Xestlist <- lapply(Xinits, NMWrapperX, X, tuningvals)
    Xestlist <- Xestlist[!sapply(Xestlist,is.null)]

  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  if(length(Xestlist) == 0){
    return(print("All models failed to converge"))
  }else{

  BIC_opt <- unlist(lapply(Xestlist, function(x)x$value))
  Xinit_opt <- Xinits[[which(BIC_opt == min(BIC_opt, na.rm = T))]]
  tune_opt <- Xestlist[[which(BIC_opt == min(BIC_opt, na.rm = T))]]$par
  
  Xopt <- EMAlgX(tune_opt, X, Xinit_opt, NM = F)

  X_final$X_final_pars <- Xopt$est_Model_param
  X_final$s_opt <- ncol(Xopt$est_Model_param$B)
  X_final$X_log.Lik <- Xopt$log.Lik
  X_final$diffList <- Xopt$diffList
  X_final$X_BIC <- Xopt$BIC
  X_final$tuningpB <- Xopt$tuningpB
  X_final$X_time_diff <- Xopt$time.diff
  
  return(X_final)
  }
}

## end of code