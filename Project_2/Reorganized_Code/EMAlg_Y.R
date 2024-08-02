## EMAlg_X.R

EMAlg_Y <- function(Data, Xest, initial_Model_Params, tuningpG1, tuningpD1, tuningpG2, tuningpD2, weights){
  
  alg.start <- Sys.time()
  
  Old_Par <- initial_Model_Params
  Old_Par$A1 <- Xest$X_final_pars$A1
  Old_Par$B1 <- Xest$X_final_pars$B1
  Old_Par$A2 <- Xest$X_final_pars$A2
  Old_Par$B2 <- Xest$X_final_pars$B2
  Old_Par$PhiR <- Xest$X_final_pars$PhiR
  Old_Par$PhiS1 <- Xest$X_final_pars$PhiS1
  Old_Par$PhiS2 <- Xest$X_final_pars$PhiS2
  Old_Par$Phi11 <- Xest$X_final_pars$Phi11
  Old_Par$Phi12 <- Xest$X_final_pars$Phi12

  continue_status <- 1
  logLikList <- NULL 
  myDiffList <- NULL
  niter <- 0
  
  #define initial values for B0/B/Phi1
  
  ############################################################################ 
  ## The EM iterations: the while loop
  ############################################################################ 
  
  while(continue_status){
    ## Here grid search is used to find optimal tuning parameters for B in the LASSO part. 
    #for the first iteration of the algorithm ONLY, use grid search to calculate the optimal tuning parameter combination
    #for B
    E_estimates <- Estep_Y(Data, Old_Par)
    
    M_estimates <- Mstep_Y(Data, Old_Par, E_estimates, tuningpG1, tuningpD1, tuningpG2, tuningpD2, weights)
    
    log.lik <- logLik_Y(Data, E_estimates) 
    myDiff <- Convergence_check(Old_Par[-c(16:27)], M_estimates)
    diff_criteria <-  max(myDiff)
    
    if (niter > 300 | diff_criteria< 0.001 ){
      continue_status <- 0
    }
    
    Old_Par <- M_estimates
    
    logLikList <- c(logLikList, log.lik) 
    myDiffList <- c(myDiffList,  diff_criteria)
    
    niter <- niter + 1
  }
  
  BICval <- BIC_Y(Data, Old_Par, log.lik)
  
  time.diff <- Sys.time() - alg.start
  
  ## to save results 
  return( list(est_Model_param = Old_Par,  log.Lik = logLikList, diffList = myDiffList, BIC = BICval, tuningpG1 = tuningpG1, tuningpD1 = tuningpD1, tuningpG2 = tuningpG2, tuningpD2 = tuningpD2, time.diff = time.diff) )
}

## end of code