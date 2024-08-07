## OverallXAlg.R

OverallXAlg <- function(Params, Data, tuningpA1 = seq(0, 15, 0.1), tuningpB1 = seq(0, 15, 0.1), tuningpA2 = seq(0, 15, 0.1), tuningpB2 = seq(0, 15, 0.1)){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  tuning_list <- list(tuningpA1, tuningpB1, tuningpA2, tuningpB2)
  
  tuning_grid <- expand.grid(tuning_list)
  tuning_grid <- lapply(seq_len(nrow(tuning_grid)), function(x) tuning_grid[x,])
  
  Xlist <- replicate(length(tuningpA1)*length(tuningpB1)*length(tuningpA2)*length(tuningpB2), list(), simplify=F)
  Xest <- NULL
  X_final <- NULL
  BICs <- rep(0, length(tuningpA1)*length(tuningpB1)*length(tuningpA2)*length(tuningpB2))
  
  #################################################################################
  ## for each # of latent factors in the grid search, determine the optimal estimate
  ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
  #################################################################################

  Xlist <- parLapply(cl, tuning_grid, EMAlg_X, Data, Params, Params)

  for(i in 1:length(tuning_grid)){
    BICs[i] <- Xlist[[i]]$BIC
  }

  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  Xopt <- Xlist[[which(BICs == min(BICs))]]
  
  X_final$X_final_pars <- Xopt$est_Model_param
  X_final$m_opt <- ncol(Xopt$est_Model_param$A1)
  X_final$u1_opt <- ncol(Xopt$est_Model_param$B1)
  X_final$u2_opt <- ncol(Xopt$est_Model_param$B2)
  X_final$X_log.Lik <- Xopt$log.Lik
  X_final$diffList <- Xopt$diffList
  X_final$X_BIC <- Xopt$BIC
  X_final$tuningpA1 <- Xopt$tuningpA1
  X_final$tuningpB1 <- Xopt$tuningpB1
  X_final$tuningpA2 <- Xopt$tuningpA2
  X_final$tuningpB2 <- Xopt$tuningpB2
  X_final$X_time_diff <- Xopt$time.diff
  
  return(X_final)
  
}

## end of code