## EMAlgAGammaAdLassoCV.R

EMAlgAGammaAdLassoCV <- function(X, Y, xk_sig, k_sig, B, B0, Phi1, Ainit, Ginit, Phi2init, Phi3init, tuningpA = seq(0.1,10,0.1), nfolds=10, weights){
  
  alg.start <- Sys.time()
  
  #define initial values for intercept vector B0 and coefficient matrix B
  niter <- 0
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  
  if(xk_sig == 1){
    B <- as.matrix(B)
    Phi3init <- as.matrix(Phi3init)
  }
  
  tuningpAinit <- tuningpA
  
  A0 <- matrix(0, nrow = p, ncol = 1, byrow = T) #use PC for initial estimates of A and B
  A0new <-  matrix(colMeans(Y), nrow = p, ncol = 1, byrow = T)
  A <- matrix(0, nrow = p, ncol = k_sig)
  Anew <- Ainit
  
  Gamma <- matrix(0, nrow = k_sig, ncol = xk_sig)
  Gammanew <- Ginit
  
  Phi2 <- matrix(0, nrow = p, ncol = p)
  Phi2new <- Phi2init
  Phi3 <- matrix(0, nrow = k_sig, ncol = k_sig)
  Phi3new <- Phi3init

  ############################################################################ 
  ## The EM iterations: the while loop
  ############################################################################ 
  
  while(((norm(A0 - A0new) > 0.00001) | (norm(A - Anew) > 0.00001) | (norm(Gamma - Gammanew) > 0.00001) | (norm(Phi2 - Phi2new) > 0.00001) | (norm(Phi3 - Phi3new) > 0.00001)) & (niter < 15000)){
    
    A0 <- A0new
    A <- Anew
    Gamma <- Gammanew
    Phi2 <- Phi2new
    Phi3 <- Phi3new
    
    ## Here grid search is used to find optimal tuning parameters for A in the LASSO part. 
    #for the first iteration of the algorithm ONLY, use grid search to calculate the optimal tuning parameter combination
    #for A
    if(niter == 0){
      
      deltaijA <- 0 #initialize coordinate descent starting values
      thetaijA <- 0
      omegaijA <- 0
      lambdaA <- NULL
      abar <- NULL
      A0CV <- rep(0, p)
      ACV <- matrix(0, nrow = p, ncol = k_sig)
      GammaCV <- matrix(0, nrow = k_sig, ncol = xk_sig)
      Phi2CV <- matrix(0, nrow = p, ncol = p)
      Phi3CV <- matrix(0, nrow = k_sig, ncol = k_sig)
      A0start <- A0
      Astart <- A
      Gammastart <- Gamma
      Phi2start <- Phi2
      Phi3start <- Phi3
      
      BICmat <- rep(0, length(tuningpA))
      
      cv.start <- Sys.time()
      for(a in 1:length(tuningpA)){
        
        ACV <- Astart
        A0CV <- A0start
        Phi2CV <- Phi2start
        Phi3CV <- Phi3start
        GammaCV <- Gammastart
        
        A <- matrix(0, nrow = p, ncol = k_sig)
        A0 <- rep(0, p)
        Gamma <- matrix(0, nrow = k_sig, ncol = xk_sig)
        Phi2 <- matrix(0, nrow = p, ncol = p)
        Phi3 <- matrix(0, nrow = k_sig, ncol = k_sig)
        
        maxit <- 0
        
        ##########################################
        ## This while loop is for calculating estimates of A0/A/Gamma/Phi2/Phi3
        ## for the current tuning parameter value a
        ##########################################
        while(((norm(A - ACV) > 0.0001) | (norm(A0 - A0CV) > 0.0001) | (norm(Gamma - GammaCV) > 0.0001) | (norm(Phi2 - Phi2CV) > 0.0001) | (norm(Phi3 - Phi3CV) > 0.0001)) & (maxit < 15000)){
          
          A0 <-A0CV
          A <- ACV
          Gamma <- GammaCV
          Phi2 <- Phi2CV
          Phi3 <- Phi3CV
          
          E_estimates <- Estep_Y(X, Y, A0, B0, A, B, Gamma, Phi1, Phi2, Phi3)
          
          M_estimates <- Mstep_Y(X, Y, E_estimates$EW, E_estimates$EZ, E_estimates$condvarWZ, E_estimates$condvar, A0, A, Gamma, Phi2, Phi3, tuningpA[a], weights)
          
          A0CV <- M_estimates$A0
          ACV <- M_estimates$A
          GammaCV <- M_estimates$Gamma
          Phi2CV <- M_estimates$Phi2
          Phi3CV <- M_estimates$Phi3
          
          maxit <- maxit + 1
        }
        
        sigma21 <-  rbind( matrix( cbind( t(B), t(GammaCV)%*% t(ACV)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(GammaCV%*%t(B),(Phi3CV + GammaCV%*%t(GammaCV))%*%t(ACV)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = GammaCV W + E
        
        sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), ACV%*%GammaCV%*%t(B)), rbind(B%*%t(ACV%*%GammaCV),
                                                                                           (Phi2CV + ACV%*%(Phi3CV + tcrossprod(GammaCV))%*%t(ACV)))), nrow = q+p, ncol = q+p)
        
        sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(GammaCV) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( GammaCV  , Phi3CV + GammaCV%*%t(GammaCV) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
        
        inv11 <- Matrix::solve(sigma11)
        condvarWZ <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
        condvar  <- condvarWZ[xk_sig+(1:k_sig), xk_sig+(1:k_sig) ]  ## m x m    
        
        log.lik <- -(n/2)*(sum(diag(inv11%*%(cov(cbind(X,Y))+(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))%*%t(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))))) + (q+p)*log((2*pi))+log(det(sigma11)))
        
        BICmat[a] <- log(n)*(sum(ACV != 0) - k_sig) - 2*log.lik
        
      }
      if(length(BICmat) > 1){
        sigmaNAs <- is.na(BICmat)
        tuningpA <- tuningpA[!sigmaNAs]
        BICmat <- BICmat[!sigmaNAs]
      }
      if(length(BICmat) == 0){
        if(length(tuningpAinit) == 1){
          tuningpAhat <- tuningpAinit
        }
        if(length(tuningpAinit) > 1){
          tuningpAhat <- min(tuningpAinit)
        }
      }else{
        tuningpAhat <- tuningpA[which(BICmat == min(BICmat))]
      }
      A <- Astart
      A0 <- A0start
      Phi2 <- Phi2start
      Phi3 <- Phi3start
      Gamma <- Gammastart
      
      #once CV MSE has been calculated for every tuning parameter, extract tuning parameter
      #that corresponds to the smallest CV SE (i.e. the optimal tuning parameter)
      if(length(tuningpAhat) > 1){
        tuningpAhat <- min(tuningpAhat)
      }
      cv.end <- Sys.time()
      cv.time <- cv.end - cv.start
    }
    
    E_estimates <- Estep_Y(X, Y, A0, B0, A, B, Gamma, Phi1, Phi2, Phi3)
    
    M_estimates <- Mstep_Y(X, Y, E_estimates$EW, E_estimates$EZ, E_estimates$condvarWZ, E_estimates$condvar, A0, A, Gamma, Phi2, Phi3, tuningpAhat, weights)
    
    A0new <- M_estimates$A0
    Anew <- M_estimates$A
    Gammanew <- M_estimates$Gamma
    Phi2new <- M_estimates$Phi2
    Phi3new <- M_estimates$Phi3
    
    niter <- niter + 1
  }
  
  A <- Anew
  A0 <- A0new
  Gamma <- Gammanew
  Phi2 <- Phi2new
  Phi3 <- Phi3new
  
  alg.end <- Sys.time()
  alg.time <- alg.end - alg.start - cv.time
  
  sigma21 <-  rbind( matrix( cbind( t(B), t(Gamma)%*% t(A)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(Gamma%*%t(B),(Phi3 + Gamma%*%t(Gamma))%*%t(A)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = Gamma W + E
  
  sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), A%*%Gamma%*%t(B)), rbind(B%*%t(A%*%Gamma),
                                                                                 (Phi2 + A%*%(Phi3 + tcrossprod(Gamma))%*%t(A)))), nrow = q+p, ncol = q+p)
  
  sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(Gamma) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( Gamma  , Phi3 + Gamma%*%t(Gamma) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
  
  inv11 <- Matrix::solve(sigma11)
  
  EWZ <- sigma21%*%inv11%*%(rbind(t(X),t(Y))-matrix(rbind(B0,A0), nrow = q+p, ncol = n, byrow = F)) ## (m + s) x n
  EW <-  matrix(EWZ[1:xk_sig,], nrow = xk_sig)   ## q x n
  EZ <-  matrix(EWZ[xk_sig+(1:k_sig), ], nrow = k_sig)  ## p x n
  
  log.lik <- -(n/2)*(sum(diag(inv11%*%(cov(cbind(X,Y))+(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))%*%t(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))))) + (q+p)*log((2*pi))+log(det(sigma11)))
  
  BICopt <- log(n)*(sum(A != 0) - k_sig) - 2*log.lik
  output <- list("A" = A, "A0" = A0, "Gamma" = Gamma, "Phi2" = Phi2, "Phi3" = Phi3, "iterations" = niter, "optimal lambda" = tuningpAhat, "Y" = Y, "EZ" = EZ, "EW" = EW, "computation time" = alg.time, "CV time" = cv.time, "BICopt" = BICopt)
  return(output)
}

## end of code
