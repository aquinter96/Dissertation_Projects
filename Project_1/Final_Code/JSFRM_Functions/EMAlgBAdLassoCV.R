## EMAlgBAdLassoCV.R

EMAlgBAdLassoCV <- function(tuningpB, Data, initial_Model_Params, weights){
  
  alg.start <- Sys.time()

  Old_Par <- initial_Model_Params
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
    E_estimates <- Estep_X(Data, Old_Par)

    M_estimates <- Mstep_X(Data, Old_Par, E_estimates, tuningpB, weights)
    
    log.lik <- logLik_X(Data, E_estimates) 
    myDiff <- Convergence_check(Old_Par, M_estimates)
    diff_criteria <-  max(myDiff)
    
    if (niter > 300 | diff_criteria< 0.001 ){
      continue_status <- 0
    }
    
    Old_Par <- M_estimates
    
    logLikList <- c(logLikList, log.lik) 
    myDiffList <- c(myDiffList,  diff_criteria) 

    niter <- niter + 1
  }

  BICval <- BIC_X(Data, Old_Par, log.lik)

  time.diff <- Sys.time() - alg.start
  
  ## to save results 
  return( list(est_Model_param = Old_Par,  log.Lik = logLikList, diffList = myDiffList, BIC = BICval, tuningpB = tuningpB, time.diff = time.diff) )
}

## end of code