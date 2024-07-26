## OverallBAlg.R

OverallBAlg <- function(Data, tuningpA = seq(0, 15, 0.1), tuningpB = seq(0, 15, 0.1), m_seq = 1:4, u_seq = list(1:4, 1:4)){
  
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
  ulist <- list()
  g <- 0
  param_grid <- list(tuningpA, tuningpB, m_seq)
  for(i in 1:length(u_seq)){
    param_grid[[3+i]] <- u_seq[[i]]
  }
  
  for(a in 1:length(m_seq)){
    for(b in 1:length(u_seq[1])){
      for(c in 1:length(u_seq[2])){
        
        g <- g + 1
        
        for(d in 1:length(tuningpA)){
          for(e in 1:length(tuningpB)){
      
      #################################################################################
      ## for each # of latent factors in the grid search, determine the optimal estimate
      ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
      #################################################################################
      
            # for(f in 1:length(u_seq)){
            #   ulist[f] <- u_seq[f][]
            # }
            
            ulist[1] <- u_seq[1][b]
            ulist[2] <- u_seq[2][c]
            
            Xinit <- X_inits(Data, m_seq[a], ulist)
            Xest <- EMAlg_X(Data, Xinit, tuningpA[d], tuningpB[e], weights = Xinit)
      
            if(d == 1 | Xest$BIC < BIC_previous){
              Xestlist[[g]] <- Xest
              BIC_previous <- Xestlist[[i]]$BIC
            }
      
    }
    
    BIC_opt[i] <- BIC_previous
    
  }
  
  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  
  Bopt <- Bestlist[[which(BIC_opt == min(BIC_opt))]]
  
  B_final$B_final_pars <- Bopt$est_Model_param
  B_final$s_opt <- ncol(Bopt$est_Model_param$B)
  B_final$B_log.Lik <- Bopt$log.Lik
  B_final$diffList <- Bopt$diffList
  B_final$B_BIC <- Bopt$BIC
  B_final$tuningpB <- Bopt$tuningpB
  B_final$B_time_diff <- Bopt$time.diff
  
  return(B_final)
}

## end of code