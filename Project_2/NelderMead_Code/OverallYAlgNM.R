## OverallBAlg.R

OverallYAlg <- function(Data, Xest, j_seq = 1:4, z_seq = list(1:4, 1:4)){
  
  #################################################################################
  ## create empty lists to store the optimal estimate of B (and corresponding estimates
  ## of B0/Phi1) and corresponding BIC for each # of latent factors for X in the grid search 
  ## (i.e. 1 latent factor, 2 latent factors, etc.)
  #################################################################################

  lat_combs <- expand.grid(c(list(j_seq), z_seq))
  lat_grid <- lapply(seq_len(nrow(lat_combs)), function(x) lat_combs[x,])
  
  Yestlist <- replicate(length(lat_grid), NULL, simplify=F)
  Yest <- NULL
  Y_final <- NULL
  All_final <- list()
  BIC_opt <- rep(0, length(lat_grid))
  c <- 0
  tuningvals <- rep(0, 4)
  
  # for(a in 1:length(j_seq)){
  #   for(b in 1:length(z_grid)){
      
      # c <- c + 1
      # 

      Yinits <- parLapply(cl, lat_grid, Y_inits, Data, Xest)
      Yestlist <- parLapply(cl, Yinits, NMWrapperY, Data, Xest, tuningvals)
      
      # Yinits <- Y_inits(lat_grid[[1]], Data, Xest)
      # Yestlist <- NMWrapperY(Yinits, Data, Xest, tuningvals)

  #   }
  # }
  #return(Yest)
  
  #################################################################################
  ## extract the # of latent factors that had the lowest BIC and the corresponding
  ## model estimates
  #################################################################################
  #return(nm)
  
  BIC_opt <- unlist(lapply(Yestlist, function(x)x$value))
  Yinit_opt <- Yinits[[which(BIC_opt == min(BIC_opt, na.rm = T))]]
  tune_opt <- Yestlist[[which(BIC_opt == min(BIC_opt, na.rm = T))]]$par

  Yopt <- EMAlg_Y(tune_opt, Data, Xest, Yinit_opt, NM = F)
      
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