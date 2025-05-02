## EMAlgAGammaAdLassoCV.R

EMAlgY <- function(tuningList, Data, Mest, initial_Model_Params){
  
  alg.start <- Sys.time()
  
  tuningpC <- as.numeric(tuningList[1])
  #tuningpC <- 0
  tuningpD <- as.numeric(tuningList[2])
  #tuningpD <- 0
  
  #define initial values for intercept vector C0 and coefficient matrix C
  
  Old_Par <- initial_Model_Params
  Old_Par$A <- Mest$M_estimates$M_final_pars$A
  Old_Par$A0 <- Mest$M_estimates$M_final_pars$A0
  Old_Par$Gamma <- Mest$M_estimates$M_final_pars$Gamma
  Old_Par$Phi2 <- Mest$M_estimates$M_final_pars$Phi2
  Old_Par$Phi3 <- Mest$M_estimates$M_final_pars$Phi3
  Old_Par$B <- Mest$X_estimates$X_final_pars$B
  Old_Par$B0 <- Mest$X_estimates$X_final_pars$B0
  Old_Par$Phi1 <- Mest$X_estimates$X_final_pars$Phi1
  Old_Par$Psi <- Mest$X_estimates$X_final_pars$Psi
  weights <- Old_Par

  # Old_Par <- initial_Model_Params
  # # Old_Par$A <- Mest$M_estimates$M_final_pars$A
  # # Old_Par$A0 <- Mest$M_estimates$M_final_pars$A0
  # # Old_Par$Gamma <- Mest$M_estimates$M_final_pars$Gamma
  # # Old_Par$Phi2 <- Mest$M_estimates$M_final_pars$Phi2
  # # Old_Par$Phi3 <- Mest$M_estimates$M_final_pars$Phi3
  # # Old_Par$B <- Mest$X_estimates$X_final_pars$B
  # # Old_Par$B0 <- Mest$X_estimates$X_final_pars$B0
  # # Old_Par$Phi1 <- Mest$X_estimates$X_final_pars$Phi1
  # # Old_Par$Psi <- Mest$X_estimates$X_final_pars$Psi
  # weights <- Old_Par
  
  continue_status <- 1
  logLikList <- NULL 
  myDiffList <- NULL
  niter <- 0
  ############################################################################ 
  ## The EM iterations: the while loop
  ############################################################################ 
  
  while(continue_status){
    
    E_estimates <- Estep_Y(Data, Old_Par)
    
    M_estimates <- Mstep_Y(Data, Old_Par, E_estimates, tuningpC, tuningpD, weights)
    
    log.lik <- logLik_Y(Data, E_estimates)
    myDiff <- Convergence_check(Old_Par[-c(7:15)], M_estimates)
    #myDiff <- Convergence_check(Old_Par[-c(1:8, 10:15)], M_estimates[-c(1:8, 10:15)])
    
    diff_criteria <-  max(myDiff)
    
    if (niter > 1000 | diff_criteria< 0.001 ){
      continue_status <- 0
    }
    
    Old_Par <- M_estimates
    
    logLikList <- c(logLikList, log.lik) 
    myDiffList <- c(myDiffList,  diff_criteria) 
    
    niter <- niter + 1
  }
  
  Old_Par <- Old_Par[-c(7:15)]
  #Old_Par <- Old_Par[9]
  
  BICval <- BIC_Y(Data, Old_Par, log.lik)
  
  time.diff <- Sys.time() - alg.start
  
  return( list(est_Model_param = Old_Par,  log.Like = logLikList, diffList = myDiffList, BIC = BICval, tuningpC = tuningpC, tuningpD = tuningpD, time.diff = time.diff) )
}

## end of code
