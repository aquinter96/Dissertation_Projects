## OverallBAlg.R

OverallXAlg <- function(Data, tuningpA1 = seq(0, 15, 0.1), tuningpB1 = seq(0, 15, 0.1), tuningpA2 = seq(0, 15, 0.1), tuningpB2 = seq(0, 15, 0.1), m_seq = 1:4, u_seq = list(1:4, 1:4)){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  Xestlist <- replicate(length(m_seq)*length(u_seq[1])*length(u_seq[2]), NULL, simplify=F)
  Xest <- NULL
  X_final <- NULL
  BIC_previous <- 0
  BIC_opt <- rep(0, length(m_seq)*length(u_seq[1])*length(u_seq[2]))
  ulist <- c()
  h <- 0

  # param_grid <- list(tuningpA, tuningpB, m_seq)
  # for(i in 1:length(u_seq)){
  #   param_grid[[3+i]] <- u_seq[[i]]
  # }
  
  for(a in 1:length(m_seq)){
    for(b in 1:length(u_seq[[1]])){
      for(c in 1:length(u_seq[[2]])){
        
        ulist[1] <- u_seq[[1]][b]
        ulist[2] <- u_seq[[2]][c]
        
        Xinit <- X_inits(Data, m_seq[a], ulist)
        
        h <- h + 1
        
        for(d in 1:length(tuningpA1)){
          for(e in 1:length(tuningpB1)){
            for(f in 1:length(tuningpA2)){
              for(g in 1:length(tuningpB2)){
                

      #################################################################################
      ## for each # of latent factors in the grid search, determine the optimal estimate
      ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
      #################################################################################
      
            # for(f in 1:length(u_seq)){
            #   ulist[f] <- u_seq[f][]
            # }
                
                Xest <- EMAlg_X(Data, Xinit, tuningpA1[d], tuningpB1[e], tuningpA2[f], tuningpB2[g], weights = Xinit)
      
                if(d == 1 | Xest$BIC < BIC_previous){
                  Xestlist[[h]] <- Xest
                  BIC_previous <- Xestlist[[h]]$BIC
                }
              }
            }
          }
        }
        
        BIC_opt[h] <- BIC_previous
        
      }
    }
  }
  
  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  Xopt <- Xestlist[[which(BIC_opt == min(BIC_opt))]]
  
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