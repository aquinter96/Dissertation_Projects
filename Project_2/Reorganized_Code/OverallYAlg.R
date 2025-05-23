## OverallBAlg.R

OverallYAlg <- function(Data, Xest, tuningpG1 = seq(0, 15, 0.1), tuningpD1 = seq(0, 15, 0.1), tuningpG2 = seq(0, 15, 0.1), tuningpD2 = seq(0, 15, 0.1), j_seq = 1:4, z_seq = list(1:4, 1:4)){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  tuning_list <- list(tuningpG1, tuningpD1, tuningpG2, tuningpD2)
  
  tuning_grid <- expand.grid(tuning_list)
  
  z_combs <- expand.grid(z_seq)
  z_grid <- lapply(seq_len(nrow(z_combs)), function(x) z_combs[x,])
  
  Yestlist <- replicate(length(j_seq)*length(z_grid), NULL, simplify=F)
  Ylist <- replicate(length(tuningpG1)*length(tuningpD1)*length(tuningpG2)*length(tuningpD2), list(), simplify=F)
  Yest <- NULL
  Y_final <- NULL
  BICs <- rep(0, length(tuningpG1)*length(tuningpD1)*length(tuningpG2)*length(tuningpD2))
  BIC_opt <- rep(0, length(j_seq)*length(z_grid))
  
  c <- 0

  for(a in 1:length(j_seq)){
    for(b in 1:length(z_grid)){
      
      c <- c + 1
      
      Yinit <- Y_inits(Data, Xest, j_seq[a], as.numeric(z_grid[[b]]))
      
      Ylist <- parLapply(cl, tuning_grid, EMAlg_Y, Data, Xest, Yinit, Yinit)
      
      for(d in 1:length(z_grid)){
        BICs <- Ylist[[d]]$BIC
      }
      
      BIC_opt[c] <- min(BICs)
      
      Xestlist[[c]] <-Xlist[[which(BICs == BIC_opt[c])]]
      
    }
  }
  
  
  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  Yopt <- Yestlist[[which(BIC_opt == min(BIC_opt))]]
  
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