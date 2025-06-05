## OverallXAlg.R

OverallXAlg <- function(Data, m_seq = 1:4, u_seq = list(1:4, 1:4)){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  lat_combs <- expand.grid(c(list(m_seq), u_seq))
  lat_grid <- lapply(seq_len(nrow(lat_combs)), function(x) lat_combs[x,])
  
  Xestlist <- replicate(length(lat_grid), NULL, simplify=F)
  Xest <- NULL
  X_final <- NULL
  BIC_opt <- rep(0, length(lat_grid))
  k <- 0
  tuningvals <- rep(0, 4)
  
  # for(i in 1:length(m_seq)){
  #   for(j in 1:length(u_grid)){
  #     
  #     k <- k + 1
      #################################################################################
      ## for each # of latent factors in the grid search, determine the optimal estimate
      ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
      #################################################################################
      
      Xinits <- parLapply(cl, lat_grid, X_inits, Data)
      Xestlist <- parLapply(cl, Xinits, NMWrapperX, Data, tuningvals)

      # Xinit <- X_inits(Data, m_seq[i], as.numeric(u_grid[[j]]))
      # Xest <- neldermead(tuningvals, EMAlg_X, lower = rep(0, 4), Data = Data, initial_Model_Params = Xinit)
      
  #   }
  # }
  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  BIC_opt <- unlist(lapply(Xestlist, function(x)x$value))
  Xinit_opt <- Xinits[[which(BIC_opt == min(BIC_opt, na.rm = T))]]
  tune_opt <- Xestlist[[which(BIC_opt == min(BIC_opt, na.rm = T))]]$par
  
  Xopt <- EMAlg_X(tune_opt, Data, Xinit_opt, NM = F)
  
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