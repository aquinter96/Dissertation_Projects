## OverallBAlg.R

OverallYAlg <- function(Params, Data, Xest, tuningpG1 = seq(0, 15, 0.1), tuningpD1 = seq(0, 15, 0.1), tuningpG2 = seq(0, 15, 0.1), tuningpD2 = seq(0, 15, 0.1)){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  tuning_list <- list(tuningpG1, tuningpD1, tuningpG2, tuningpD2)
  
  tuning_grid <- expand.grid(tuning_list)
  tuning_grid <- lapply(seq_len(nrow(tuning_grid)), function(x) tuning_grid[x,])

  Ylist <- replicate(length(tuningpG1)*length(tuningpD1)*length(tuningpG2)*length(tuningpD2), list(), simplify=F)
  Y_final <- NULL
  All_final <- NULL
  BICs <- rep(0, length(tuningpG1)*length(tuningpD1)*length(tuningpG2)*length(tuningpD2))

  Ylist <- parLapply(cl, tuning_grid, EMAlg_Y, Data, Xest, Params, Params)
      
  for(i in 1:length(tuning_grid)){
    BICs[i] <- Ylist[[i]]$BIC
  }

  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  Yopt <- Ylist[[which(BICs == min(BICs))]]
  
  Y_final$Y_final_pars <- Yopt$est_Model_param
  Y_final$j_opt <- ncol(Yopt$est_Model_param$G1)
  Y_final$z1_opt <- ncol(Yopt$est_Model_param$D1)
  Y_final$z2_opt <- ncol(Yopt$est_Model_param$D2)
  Y_final$Y_log.Lik <- Yopt$log.Lik
  Y_final$diffList <- Yopt$diffList
  Y_final$Y_BIC <- Yopt$BIC
  Y_final$tuningpG1 <- Yopt$tuningpG1
  Y_final$tuningpD1 <- Yopt$tuningpD1
  Y_final$tuningpG2 <- Yopt$tuningpG2
  Y_final$tuningpD2 <- Yopt$tuningpD2
  Y_final$Y_time_diff <- Yopt$time.diff
  
  All_final$X_estimates <- Xest
  All_final$Y_estimates <- Y_final
  
  return(All_final)
}

## end of code