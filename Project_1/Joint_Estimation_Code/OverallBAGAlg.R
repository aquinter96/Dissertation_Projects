## OverallAGAlg.R

OverallBAGAlg <- function(Data, tuningpB = seq(0, 15, 0.1), tuningpA = seq(0, 15, 0.1), s_seq = 1:4, m_seq = 1:4, Ainit){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of A (and corresponding estimates
  ## of A0/Gamma/Phi2/Phi3) and corresponding BIC for each # of latent factors for Y in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.) Also extract the optimal estimates
  ## calculated in OverallBAlg.R of B/B0/Phi1 for use in the estimation of A/A0/Gamma/Phi2/Phi3
  #################################################################################
  tuning_list <- list(tuningpB, tuningpA)
  
  tuning_grid <- expand.grid(tuning_list)
  tuning_grid <- lapply(seq_len(nrow(tuning_grid)), function(x) tuning_grid[x,])
  ABestlist <- replicate(length(s_seq)*length(m_seq), NULL, simplify=F)
  ABest <- NULL
  ABlist <- replicate(length(tuningpB)*length(tuningpA), list(), simplify=F)
  AB_final <- list()
  BICBs <- rep(0, length(tuningpB))
  BICs <- rep(0, length(tuningpB)*length(tuningpA))
  BIC_opt <- rep(0, length(s_seq)*length(m_seq))
  k <- 0
  #################################################################################
  ## for each # of latent factors in the grid search, determine the optimal estimate
  ## of A/A0/GammaPhi2/Phi3, then save the corresponding estimates and BIC in the respective lists
  #################################################################################
  
  # for(i in 1:length(s_seq)){
  #   for(j in 1:length(m_seq)){
      #################################################################################
      ## for each # of latent factors in the grid search, determine the optimal estimate
      ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
      #################################################################################

      k <- k + 1
      
      # Binit <- B_inits(Data$X, s_seq[i])
      # Blist <- parLapply(cl, tuningpB, EMAlgBAdLassoCV, Data$X, Binit, Binit$B)
      # 
      # for(l in 1:length(tuningpB)){
      #   BICBs[l] <- Blist[[l]]$BIC
      # }
      # 
      # BICB_opt <- min(BICBs)
      # 
      # Best <-Blist[[which(BICBs == BICB_opt)]]
      # inits <- AB_inits(Data, Best, s_seq[i], m_seq[j])
      
      inits <- list()
      
      inits$B0 <- Ainit$B_estimates$B_final_pars$B0
      inits$B <- Ainit$B_estimates$B_final_pars$B
      inits$Psi <- Ainit$B_estimates$B_final_pars$Psi
      inits$Phi1 <- Ainit$B_estimates$B_final_pars$Phi1
      inits$A0 <- Ainit$A_estimates$A_final_pars$A0
      inits$A <- Ainit$A_estimates$A_final_pars$A
      inits$Gamma <- Ainit$A_estimates$A_final_pars$Gamma
      inits$Phi2 <- Ainit$A_estimates$A_final_pars$Phi2
      inits$Phi3 <- Ainit$A_estimates$A_final_pars$Phi3
      
      ABlist <- lapply(tuning_grid, EMAlgBAGammaAdLassoCV, Data, inits, inits$B, inits$A)
      
      for(l in 1:length(tuning_grid)){
        BICs[l] <- ABlist[[l]]$BIC
      }

      ABest <- ABlist[[min(which(BICs == min(BICs)))]]

      # for(l in 1:length(tuning_grid)){
      #   BICs[l] <- ABlist[[l]]$BIC
      # }
      # 
      # BIC_opt[k] <- min(BICs)
      # 
      # ABestlist[[k]] <-ABlist[[which(BICs == BIC_opt[k])]]
      # 
  #   }
  # }

  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  ABopt <- ABest
  
  AB_final$final_pars <- ABopt$est_Model_param
  AB_final$s_opt <- ncol(ABopt$est_Model_param$B)
  AB_final$m_opt <- ncol(ABopt$est_Model_param$A)
  AB_final$log.Lik <- ABopt$log.Lik
  AB_final$diffList <- ABopt$diffList
  AB_final$BIC <- ABopt$BIC
  AB_final$tuningpB <- ABopt$tuningpB
  AB_final$tuningpA <- ABopt$tuningpA
  AB_final$time_diff <- ABopt$time.diff

  return(AB_final)
}

## end of code