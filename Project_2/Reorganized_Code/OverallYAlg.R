## OverallBAlg.R

OverallYAlg <- function(Data, Xest, tuningpG1 = seq(0, 15, 0.1), tuningpD1 = seq(0, 15, 0.1), tuningpG2 = seq(0, 15, 0.1), tuningpD2 = seq(0, 15, 0.1), j_seq = 1:4, z_seq = list(1:4, 1:4)){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  Yestlist <- replicate(length(j_seq)*length(z_seq[1])*length(z_seq[2]), NULL, simplify=F)
  Yest <- NULL
  Y_final <- NULL
  All_final <- list()
  BIC_previous <- 0
  BIC_opt <- rep(0, length(j_seq)*length(z_seq[1])*length(z_seq[2]))
  zlist <- c()
  h <- 0
  
  # param_grid <- list(tuningpA, tuningpB, m_seq)
  # for(i in 1:length(u_seq)){
  #   param_grid[[3+i]] <- u_seq[[i]]
  # }
  
  for(a in 1:length(j_seq)){
    for(b in 1:length(z_seq[1])){
      for(c in 1:length(z_seq[2])){
        
        zlist[1] <- z_seq[[1]][b]
        zlist[2] <- z_seq[[2]][c]
        
        Yinit <- Y_inits(Data, Xest, j_seq[a], zlist)
        
        h <- h + 1
        
        for(d in 1:length(tuningpG1)){
          for(e in 1:length(tuningpD1)){
            for(f in 1:length(tuningpG2)){
              for(g in 1:length(tuningpD2)){
                
                
                #################################################################################
                ## for each # of latent factors in the grid search, determine the optimal estimate
                ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
                #################################################################################
                
                # for(f in 1:length(u_seq)){
                #   ulist[f] <- u_seq[f][]
                # }
                
                Yest <- EMAlg_Y(Data, Xest, Yinit, tuningpG1[d], tuningpD1[e], tuningpG2[f], tuningpD2[g], weights = Yinit)
                
                if(d == 1 | Yest$BIC < BIC_previous){
                  Yestlist[[h]] <- Yest
                  BIC_previous <- Yestlist[[h]]$BIC
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
  
  Yopt <- Yestlist[[which(BIC_opt == min(BIC_opt))]]
  
  Y_final$Y_final_pars <- Yopt$est_Model_param
  Y_final$m_opt <- ncol(Yopt$est_Model_param$A1)
  Y_final$u1_opt <- ncol(Yopt$est_Model_param$B1)
  Y_final$u2_opt <- ncol(Yopt$est_Model_param$B2)
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