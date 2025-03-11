## OverallXAlg.R

OverallXAlg <- function(Data, tuningpA1 = seq(0, 15, 0.1), tuningpB1 = seq(0, 15, 0.1), tuningpA2 = seq(0, 15, 0.1), tuningpB2 = seq(0, 15, 0.1), m_seq = 1:4, u_seq = list(1:4, 1:4)){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  tuning_list <- list(tuningpA1, tuningpB1, tuningpA2, tuningpB2)
  
  tuning_grid <- expand.grid(tuning_list)
  tuning_grid <- lapply(seq_len(nrow(tuning_grid)), function(x) tuning_grid[x,])
  u_combs <- expand.grid(u_seq)
  u_grid <- lapply(seq_len(nrow(u_combs)), function(x) u_combs[x,])
  lat_combs <- expand.grid(c(list(m_seq), u_seq))
  lat_grid <- lapply(seq_len(nrow(lat_combs)), function(x) lat_combs[x,])
  
  Xestlist <- replicate(length(tuningpA1)*length(tuningpB1)*length(tuningpA2)*length(tuningpB2), list(), simplify=F)
  Xlist <- replicate(length(m_seq)*length(u_grid), NULL, simplify=F)
  Xest <- NULL
  X_final <- NULL
  BICs <- rep(0, length(m_seq)*length(u_grid))
  BIC_opt <- rep(0, length(tuningpA1)*length(tuningpB1)*length(tuningpA2)*length(tuningpB2))
  m <- 0
  
  for(i in 1:length(tuningpA1)){
    for(j in 1:length(tuningpB1)){
      for(k in 1:length(tuningpA2)){
        for(l in 1:length(tuningpB2)){
          
      m <- m + 1
      #################################################################################
      ## for each # of latent factors in the grid search, determine the optimal estimate
      ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
      #################################################################################
      Xinits <- parLapply(cl, lat_grid, X_inits, Data)
      #Xinits <- lapply(lat_grid, X_inits, Data)
      
      #Xlist <- lapply(Xinits, EMAlg_X, Data, tuningpA1[i], tuningpB1[j], tuningpA2[k], tuningpB2[l])
      Xlist <- parLapply(cl, Xinits, Singular_ErrorX2, Data, tuningpA1[i], tuningpB1[j], tuningpA2[k], tuningpB2[l])
      Xlist[sapply(Xlist, is.null)] <- NULL
      length(BICs) <- length(Xlist)
      if(length(Xlist) >= 1){
        for(h in 1:length(Xlist)){
          BICs[h] <- Xlist[[h]]$BIC
        }
        
        Xlist <- Xlist[which(!is.nan(BICs))]
        BICs <- BICs[!is.nan(BICs)]
        
        BIC_opt[k] <- min(BICs)
        
        Xestlist[[k]] <- Xlist[[min(which(BICs == BIC_opt[k]))]]
      }
      else{
        BIC_opt[k] <- NA
        Xestlist[[k]] <- NULL
      }
        }
      }
    }
  }
  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  Xopt <- Xestlist[[which(BIC_opt == min(BIC_opt, na.rm = T))]]
  
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