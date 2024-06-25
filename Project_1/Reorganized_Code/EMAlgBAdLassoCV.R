## EMAlgBAdLassoCV.R

EMAlgBAdLassoCV <- function(Data, initial_Model_Params, tuningpB = seq(5.1,5.5,0.1), weights){
  alg.start <- Sys.time()
  
  tuningpBinit <- tuningpB
  
  #define initial values for B0/B/Phi1
  niter <- 0
  n <- nrow(X)
  q <- ncol(X)
  B0 <- matrix(0, nrow = q, ncol = 1, byrow = T)
  B0new <-  matrix(colMeans(X), nrow = q, ncol = 1, byrow = T)
  B <- matrix(0, nrow = q, ncol = xk_sig)
  Bnew <- Binit
  Phi1 <- matrix(0, nrow = q, ncol = q)
  Phi1new <- Phi1init

  ############################################################################ 
  ## The EM iterations: the while loop
  ############################################################################ 
  
  while(( (norm(B0 - B0new, type = "F") > 0.00001) | (norm(B - Bnew, type = "F") > 0.00001) | (norm(Phi1 - Phi1new, type = "F") > 0.00001) ) & (niter < 15000)){
    B0 <- B0new
    B <- Bnew
    Phi1 <- Phi1new
    
    ## Here grid search is used to find optimal tuning parameters for B in the LASSO part. 
    #for the first iteration of the algorithm ONLY, use grid search to calculate the optimal tuning parameter combination
    #for B
    if(niter == 0){
      
      deltaijB <- 0 #initialize coordinate descent starting values
      thetaijB <- 0
      omegaijB <- 0
      lambdaB <- NULL
      bbar <- NULL
      
      B0Cv <- rep(0, q)
      BCV <- matrix(0, nrow = n, ncol = q)
      Phi1CV <- matrix(0, nrow = q, ncol = q)
      B0start <- B0
      Bstart <- B
      Phi1start <- Phi1
      
      BICmat <- rep(0,length(tuningpB))
      
      cv.start <- Sys.time()
      maxit <- 0
      
      for(a in 1:length(tuningpB)){
        
        B0CV <- B0start
        BCV <- Bstart
        Phi1CV <- Phi1start
        
        B0 <- rep(0, q)
        B <- matrix(0, nrow = q, ncol = xk_sig)
        Phi1 <- matrix(0, nrow = q, ncol = q)
        
        maxit <- 0
        
        ##########################################
        ## This while loop is for calculating estimates of B0/B/Phi1
        ## for the current tuning parameter value a
        ##########################################
        while(( (norm(B0 - B0CV, type = "F") > 0.00001) | (norm(B - BCV, type = "F") > 0.00001) | (norm(Phi1 - Phi1CV, type = "F") > 0.00001) ) & (maxit < 15000)){
          
          B0 <- B0CV
          B <- BCV
          Phi1 <- Phi1CV
          
          E_estimates <- Estep_X(X, B0, B, Phi1)

          M_estimates <- Mstep_X(X, E_estimates$EW, E_estimates$condvar, B0, B, Phi1, tuningpB[a], weights)
          
          B0CV <- M_estimates$B0
          BCV <- M_estimates$B
          PhiCV <- M_estimates$Phi1
          
          maxit <- maxit + 1
        }
        #once the coordinate descent algorithm has converged, use converged B estimate to calculate predicted values of test set
        invmodcv <- Matrix::solve(Phi1CV + tcrossprod(BCV))
        log.lik <- -(n/2)*(sum(diag(invmodcv%*%(cov(X)+(colMeans(X)-colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
                             q*log((2*pi))+log(det(Phi1CV +tcrossprod(BCV))))
        
        BICmat[a] <- log(n)*(sum(BCV != 0) - xk_sig) - 2*log.lik
        
        #calculate CV MSE for each tuning parameter
        #using a running sum over each fold 
      }
      if(length(BICmat) > 1){
        sigmaNAs <- is.na(BICmat)
        tuningpB <- tuningpB[!sigmaNAs]
        BICmat <- BICmat[!sigmaNAs]
      }
      if(length(BICmat) == 0){
        if(length(tuningpBinit) == 1){
          tuningpBhat <- tuningpBinit
        }
        if(length(tuningpBinit) > 1){
          tuningpBhat <- min(tuningpBinit)
        }
      }else{
        tuningpBhat <- tuningpB[which(BICmat == min(BICmat))]
      }
      
      B0 <- B0start
      B <- Bstart
      Phi1 <- Phi1start
      
      #once BIC has been calculated for every tuning parameter, extract tuning parameter
      #that corresponds to the smallest CV SE (i.e. the optimal tuning parameter)
      cv.end <- Sys.time()
      cv.time <- cv.end - cv.start
    }
    
    
    E_estimates <- Estep_X(X, B0, B, Phi1)
    
    M_estimates <- Mstep_X(X, E_estimates$EW, E_estimates$condvar, B0, B, Phi1, tuningpBhat, weights)
    
    B0new <- M_estimates$B0
    Bnew <- M_estimates$B
    Phinew <- M_estimates$Phi1
    
    niter <- niter + 1
  }
  B0 <- B0new
  B <- Bnew
  Phi1 <- Phi1new
  
  invmodv <- Matrix::solve(Phi1 + tcrossprod(B))
  condvar <- diag(xk_sig) - crossprod(B,invmodv)%*%B
  EW <- crossprod(B,invmodv)%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))
  
  log.lik <- -(n/2)*(sum(diag(invmodv%*%(cov(X)+(colMeans(X)-colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
                       q*log((2*pi))+log(det(Phi1 +tcrossprod(B))))
  
  BICopt <- log(n)*(sum(B != 0) - xk_sig) - 2*log.lik
  
  alg.end <- Sys.time()
  alg.time <- alg.end - alg.start - cv.time
  output <- list("B" = B, "B0" = B0, "Phi1" = Phi1, "iterations" = niter, "optimal lambda" = tuningpBhat, "X" = X, "EW" = EW, "computation time" = alg.time, "CV time" = cv.time, "BICopt" = BICopt)
  return(output)
}

## end of code