## OverallBAlg.R

OverallBAlg <- function(X, tuningpB = seq(0, 15, 0.1), s_seq = 1:4){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  Bestlist <- replicate(length(s_seq), NULL, simplify=F)
  Best <- NULL
  B_final <- NULL
  BIC_previous <- 0
  BIC_opt <- rep(0, length(s_seq))
  
  for(i in 1:length(s_seq)){
    for(j in 1:length(tuningpB)){
    #################################################################################
    ## for each # of latent factors in the grid search, determine the optimal estimate
    ## of B/B0/Phi1, then save the corresponding estimates and BIC in the respective lists
    #################################################################################
    
      Binit <- B_inits(X, s_seq[i])
      Best <- EMAlgBAdLassoCV(X, Binit, tuningpB[j], weights = Binit$B)

      if(j == 1 | Best$BIC < BIC_previous){
        Bestlist[[i]] <- Best
        BIC_previous <- Bestlist[[i]]$BIC
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