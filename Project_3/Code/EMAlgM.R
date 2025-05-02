## EMAlgAGammaAdLassoCV.R

EMAlgM <- function(tuningpA, Data, Xest, initial_Model_Params, weights){
  
  alg.start <- Sys.time()
  
  #define initial values for intercept vector B0 and coefficient matrix B
  Old_Par <- initial_Model_Params
  Old_Par$B <- Xest$X_final_pars$B
  Old_Par$B0 <- Xest$X_final_pars$B0
  Old_Par$Phi1 <- Xest$X_final_pars$Phi1
  Old_Par$Psi <- Xest$X_final_pars$Psi
  
  continue_status <- 1
  logLikList <- NULL 
  myDiffList <- NULL
  niter <- 0
  ############################################################################ 
  ## The EM iterations: the while loop
  ############################################################################ 

  while(continue_status){
    
    E_estimates <- Estep_M(Data, Old_Par)
    
    M_estimates <- Mstep_M(Data, Old_Par, E_estimates, tuningpA, weights)
    
    log.lik <- logLik_M(Data, E_estimates)
    myDiff <- Convergence_check(Old_Par[-c(6, 7, 8, 9)], M_estimates)
    diff_criteria <-  max(myDiff)
    
    if (niter > 1000 | diff_criteria< 0.001 ){
      continue_status <- 0
    }
    
    Old_Par <- M_estimates
    
    logLikList <- c(logLikList, log.lik) 
    myDiffList <- c(myDiffList,  diff_criteria) 
    
    niter <- niter + 1
  }
  
  Old_Par <- Old_Par[-c(6, 7, 8, 9)]
  
  BICval <- BIC_M(Data, Old_Par, log.lik)
  
  time.diff <- Sys.time() - alg.start
  
  return( list(est_Model_param = Old_Par,  log.Like = logLikList, diffList = myDiffList, BIC = BICval, tuningpA = tuningpA, time.diff = time.diff) )
}

## end of code
