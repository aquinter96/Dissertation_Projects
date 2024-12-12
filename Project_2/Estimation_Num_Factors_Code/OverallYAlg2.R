## OverallBAlg.R

OverallYAlg <- function(Data, Xest, tuningpG1 = seq(0, 15, 0.1), tuningpD1 = seq(0, 15, 0.1), tuningpG2 = seq(0, 15, 0.1), tuningpD2 = seq(0, 15, 0.1), j_seq = 1:4, z_seq = list(1:4, 1:4)){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################
  
  tuning_list <- list(tuningpG1, tuningpD1, tuningpG2, tuningpD2)
  
  tuning_grid <- expand.grid(tuning_list)
  tuning_grid <- lapply(seq_len(nrow(tuning_grid)), function(x) tuning_grid[x,])
  z_combs <- expand.grid(z_seq)
  z_grid <- lapply(seq_len(nrow(z_combs)), function(x) z_combs[x,])
  lat_combs <- expand.grid(c(list(j_seq), z_seq))
  lat_grid <- lapply(seq_len(nrow(lat_combs)), function(x) lat_combs[x,])
  
  Yestlist <- replicate(length(tuningpG1)*length(tuningpD1)*length(tuningpG2)*length(tuningpD2), NULL, simplify=F)
  Ylist <- replicate(length(j_seq)*length(z_grid), list(), simplify=F)
  Yest <- NULL
  Y_final <- NULL
  BICs <- rep(0, length(j_seq)*length(z_grid))
  BIC_opt <- rep(0, length(tuningpG1)*length(tuningpD1)*length(tuningpG2)*length(tuningpD2))
  
  All_final <- list()
  
  c <- 0
  
  for(i in 1:length(tuningpG1)){
    for(j in 1:length(tuningpD1)){
      for(k in 1:length(tuningpG2)){
        for(l in 1:length(tuningpD2)){
          
          c <- c + 1
          
          #Yinits <- lapply(lat_grid, Y_inits, Data, Xest)
          Yinits <- parLapply(cl, lat_grid, Y_inits, Data, Xest)
          print("init done")
          #Ylist <- lapply(Yinits, EMAlg_Y, Data, Xest, tuningpG1[i], tuningpD1[j], tuningpG2[k], tuningpD2[l])
          Ylist <- parLapply(cl, Yinits, EMAlg_Y, Data, Xest, tuningpG1[i], tuningpD1[j], tuningpG2[k], tuningpD2[l])
          print("em done")
          for(d in 1:length(Ylist)){
            BICs[[d]] <- Ylist[[d]]$BIC
          }
          
          BIC_opt[c] <- min(BICs)
          
          Yestlist[[c]] <-Ylist[[which(BICs == BIC_opt[c])]]
        }
      }
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