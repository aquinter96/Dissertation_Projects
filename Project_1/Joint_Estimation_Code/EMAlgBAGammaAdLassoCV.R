## EMAlgBAGammaAdLassoCV.R

EMAlgBAGammaAdLassoCV <- function(tuningplist, Data, initial_Model_Params, weightsB, weightsA){
  alg.start <- Sys.time()
  
  #define initial values for intercept vector B0 and coefficient matrix B

  tuningpB <- as.numeric(tuningplist[1])
  tuningpA <- as.numeric(tuningplist[2])
  
  Old_Par <- initial_Model_Params
  continue_status <- 1
  logLikList <- NULL 
  myDiffList <- NULL
  niter <- 0
  ############################################################################ 
  ## The EM iterations: the while loop
  ############################################################################ 

  while(continue_status){
    E_estimates <- Estep_XY(Data, Old_Par)
    M_estimates <- Mstep_XY(Data, Old_Par, E_estimates, tuningpB, tuningpA, weightsB, weightsA)
    log.lik <- logLik_XY(Data, E_estimates)
    myDiff <- Convergence_check(Old_Par, M_estimates)
    diff_criteria <-  max(myDiff)
    if (niter > 1000 | diff_criteria< 0.00001 ){
      continue_status <- 0
    }
    
    Old_Par <- M_estimates
    
    logLikList <- c(logLikList, log.lik) 
    myDiffList <- c(myDiffList,  diff_criteria) 
    
    niter <- niter + 1
  }

  BICval <- BIC_XY(Data, Old_Par, log.lik)
  time.diff <- Sys.time() - alg.start
  
  return( list(est_Model_param = Old_Par,  log.Like = logLikList, diffList = myDiffList, BIC = BICval, tuningpB = tuningpB, tuningpA = tuningpA, time.diff = time.diff) )
}

## end of code
