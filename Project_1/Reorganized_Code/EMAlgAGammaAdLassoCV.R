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
  ## assume sparsity for A, and thus, consider LASSO penalty  
  ## Here, we define a lasso function: a soft-thresholding function
  ############################################################################ 
  lasso <- function(x,y){ #initialize soft-thresholding function
    result <- NULL
    if(abs(x) <= y){
      result <- 0
    } else{
      result <- x - y*sign(x)
    }
    return(result)
  }
  
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
          
          sigma21 <-  rbind( matrix( cbind( t(B), t(Gamma)%*% t(A)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(Gamma%*%t(B),(Phi3 + Gamma%*%t(Gamma))%*%t(A)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = Gamma W + E
          
          sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), A%*%Gamma%*%t(B)), rbind(B%*%t(A%*%Gamma),
                                                                                         (Phi2 + A%*%(Phi3 + tcrossprod(Gamma))%*%t(A)))), nrow = q+p, ncol = q+p)
          
          sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(Gamma) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( Gamma  , Phi3 + Gamma%*%t(Gamma) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
          
          inv11 <- Matrix::solve(sigma11)
          condvarWZ <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
          condvar  <- condvarWZ[xk_sig+(1:k_sig), xk_sig+(1:k_sig) ]  ## m x m    
          
          EWZ <- sigma21%*%inv11%*%(rbind(t(X),t(Y))-matrix(rbind(B0,A0), nrow = q+p, ncol = n, byrow = F)) ## (m + s) x n
          EW <-  matrix(EWZ[1:xk_sig,], nrow = xk_sig)   ## q x n
          EZ <-  matrix(EWZ[xk_sig+(1:k_sig), ], nrow = k_sig)  ## p x n
          
          ##################################              
          ## for each element in the first part of A: mxm  
          ## S is p*m
          ##################################              
          for(i in 1:k_sig){
            for(j in 1:k_sig){
              deltaijA <- 0
              thetaijA <- 0
              omegaijA <- 0
              #calculate coordinate descent updates for the first sxs submatrix of B, while ensuring that said sxs matrix is lower triangular
              if(j < i){
                deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
                thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0CV[i,], n)))
                omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%ACV[i,-j]
                abar <- ((thetaijA-omegaijA))/(deltaijA)
                if(weights[i,j] == 0){
                  lambdaA <- (tuningpA[a]*Phi2CV[i,i])/deltaijA*(1/1e-5)
                }
                else{
                  lambdaA <- (tuningpA[a]*Phi2CV[i,i])/deltaijA*(1/abs(weights[i,j]))
                }
                ACV[i,j] <- lasso(abar, lambdaA)
              }
            }
          }
          for(i in (k_sig+1):p){
            for(j in 1:k_sig){
              deltaijA <- 0
              thetaijA <- 0
              omegaijA <- 0
              deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
              thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0CV[i,], n)))
              omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%ACV[i,-j]
              abar <- ((thetaijA-omegaijA))/(deltaijA)
              if(weights[i,j] == 0){
                lambdaA <- (tuningpA[a]*Phi2CV[i,i])/deltaijA*(1/1e-5)
              }
              else{
                lambdaA <- (tuningpA[a]*Phi2CV[i,i])/deltaijA*(1/abs(weights[i,j]))
              }
              ACV[i,j] <- lasso(abar, lambdaA)
            }
          }
          
          #################################################################################
          ## back to the while loop to update A0/Gamma/Phi2/Phi3
          #################################################################################
          
          A0CV <- (1/n)*as.matrix(colSums(Y - t(ACV%*%EZ)))
          
          Phi2CV <- (1/n)*diag(diag(t(Y)%*%Y - 2*t(Y)%*%matrix(A0, nrow = n, ncol = p, byrow = T) - 2*t(Y)%*%t(EZ)%*%t(ACV) + matrix(A0, nrow = p, ncol = n, byrow = F)%*%t(matrix(A0, nrow = p, ncol = n, byrow = F)) + 
                                      2*matrix(A0, nrow = p, ncol = n, byrow = F)%*%t(EZ)%*%t(ACV) + n*ACV%*%condvar%*%t(ACV) + ACV%*%EZ%*%t(EZ)%*%t(ACV)))
          
          GammaCV <- as.matrix((n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)] + EZ%*%t(EW))%*%solve(n*condvarWZ[(1:xk_sig),(1:xk_sig)] + EW%*%t(EW)))
          
          #Phi3 Update
          if(k_sig == 1){
            Phi3CV <- (1/n)*as.matrix(n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),((xk_sig+1):(xk_sig+k_sig))] + EZ%*%t(EZ) - 2*n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)]%*%t(GammaCV) -
                                        2*EZ%*%t(EW)%*%t(GammaCV) + n*GammaCV%*%condvarWZ[(1:xk_sig),(1:xk_sig)]%*%t(GammaCV) + GammaCV%*%EW%*%t(EW)%*%t(GammaCV))
          }else{
            Phi3CV <- (1/n)*as.matrix(diag(diag(n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),((xk_sig+1):(xk_sig+k_sig))] + EZ%*%t(EZ) - 2*n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)]%*%t(GammaCV) -
                                                  2*EZ%*%t(EW)%*%t(GammaCV) + n*GammaCV%*%condvarWZ[(1:xk_sig),(1:xk_sig)]%*%t(GammaCV) + GammaCV%*%EW%*%t(EW)%*%t(GammaCV))))
          }
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
    
    sigma21 <-  rbind( matrix( cbind( t(B), t(Gamma)%*% t(A)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(Gamma%*%t(B),(Phi3 + Gamma%*%t(Gamma))%*%t(A)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = Gamma W + E
    
    sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), A%*%Gamma%*%t(B)), rbind(B%*%t(A%*%Gamma),
                                                                                   (Phi2 + A%*%(Phi3 + tcrossprod(Gamma))%*%t(A)))), nrow = q+p, ncol = q+p)
    
    sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(Gamma) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( Gamma  , Phi3 + Gamma%*%t(Gamma) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
    
    inv11 <- Matrix::solve(sigma11)
    condvarWZ <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
    condvar  <- condvarWZ[xk_sig+(1:k_sig), xk_sig+(1:k_sig) ]  ## m x m    
    EZ <-  matrix(EWZ[xk_sig+(1:k_sig), ], nrow = k_sig)  ## p x n
    EW <-  matrix(EWZ[1:xk_sig,], nrow = xk_sig)   ## q x n
    
    #now that the optimal tuning parameter has been chosen, run the EM algorithm as normal, using the optimal tuning parameter for all subsequent iterations
    
    for(i in 1:k_sig){
      for(j in 1:k_sig){
        deltaijA <- 0
        thetaijA <- 0
        omegaijA <- 0
        if(j < i){
          deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
          thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0[i,], n)))
          omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%Anew[i,-j]
          abar <- ((thetaijA-omegaijA))/(deltaijA)
          if(weights[i,j] == 0){
            lambdaA <- (tuningpAhat*Phi2[i,i])/deltaijA*(1/1e-5)
          }
          else{
            lambdaA <- (tuningpAhat*Phi2[i,i])/deltaijA*(1/abs(weights[i,j]))
          }
          Anew[i,j] <- lasso(abar, lambdaA)
        }
      }
    }
    for(i in (k_sig+1):p){
      for(j in 1:k_sig){
        deltaijA <- 0
        thetaijA <- 0
        omegaijA <- 0
        deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
        thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0[i,], n)))
        omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%Anew[i,-j]
        abar <- ((thetaijA-omegaijA))/(deltaijA)
        if(weights[i,j] == 0){
          lambdaA <- (tuningpAhat*Phi2[i,i])/deltaijA*(1/1e-5)
        }
        else{
          lambdaA <- (tuningpAhat*Phi2[i,i])/deltaijA*(1/abs(weights[i,j]))
        }
        Anew[i,j] <- lasso(abar, lambdaA)
      }
    }
    
    #A0 Update
    A0new <- (1/n)*as.matrix(colSums(Y - t(Anew%*%EZ)))
    
    #Phi2 Update
    Phi2new <- (1/n)*diag(diag(t(Y)%*%Y - 2*t(Y)%*%matrix(A0new, nrow = n, ncol = p, byrow = T) - 2*t(Y)%*%t(EZ)%*%t(Anew) + matrix(A0new, nrow = p, ncol = n, byrow = F)%*%t(matrix(A0new, nrow = p, ncol = n, byrow = F)) + 
                                 2*matrix(A0new, nrow = p, ncol = n, byrow = F)%*%t(EZ)%*%t(Anew) + n*Anew%*%condvar%*%t(Anew) + Anew%*%EZ%*%t(EZ)%*%t(Anew)))
    
    #conditional covariance matrix of W and Z
    condvarWZ <- matrix(cbind(rbind(diag(xk_sig), Gamma), rbind(t(Gamma),(Phi3 + Gamma%*%t(Gamma)))), nrow = xk_sig+k_sig, ncol = xk_sig+k_sig) -
      matrix(cbind(rbind(t(B), Gamma%*%t(B)), rbind(t(Gamma)%*%t(Anew),(Phi3 + Gamma%*%t(Gamma))%*%t(Anew))), nrow = xk_sig+k_sig, ncol = q+p)%*%solve(matrix(cbind(rbind((Phi1 + tcrossprod(B)), Anew%*%Gamma%*%t(B)), rbind(B%*%t(Anew%*%Gamma),
                                                                                                                                                                                                                              (Phi2new + Anew%*%(Phi3 + tcrossprod(Gamma))%*%t(Anew)))), nrow = q+p, ncol = q+p))%*%t(matrix(cbind(rbind(t(B), Gamma%*%t(B)), rbind(t(Gamma)%*%t(Anew),(Phi3 + Gamma%*%t(Gamma))%*%t(Anew))), nrow = xk_sig+k_sig, ncol = q+p))
    
    #Gamma Update
    Gammanew <- as.matrix((n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)] + EZ%*%t(EW))%*%solve(n*condvarWZ[(1:xk_sig),(1:xk_sig)] + EW%*%t(EW)))
    
    #Phi3 Update
    if(k_sig == 1){
      Phi3new <- as.matrix((1/n)*(n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),((xk_sig+1):(xk_sig+k_sig))] + EZ%*%t(EZ) - 2*n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)]%*%t(Gammanew) -
                                    2*EZ%*%t(EW)%*%t(Gammanew) + n*Gammanew%*%condvarWZ[(1:xk_sig),(1:xk_sig)]%*%t(Gammanew) + Gammanew%*%EW%*%t(EW)%*%t(Gammanew)))
    }else{
      Phi3new <- as.matrix((1/n)*diag(diag((n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),((xk_sig+1):(xk_sig+k_sig))] + EZ%*%t(EZ) - 2*n*condvarWZ[((xk_sig+1):(xk_sig+k_sig)),(1:xk_sig)]%*%t(Gammanew) -
                                              2*EZ%*%t(EW)%*%t(Gammanew) + n*Gammanew%*%condvarWZ[(1:xk_sig),(1:xk_sig)]%*%t(Gammanew) + Gammanew%*%EW%*%t(EW)%*%t(Gammanew)))))
    }
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
  log.lik <- -(n/2)*(sum(diag(inv11%*%(cov(cbind(X,Y))+(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))%*%t(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))))) + (q+p)*log((2*pi))+log(det(sigma11)))
  
  BICopt <- log(n)*(sum(A != 0) - k_sig) - 2*log.lik
  output <- list("A" = A, "A0" = A0, "Gamma" = Gamma, "Phi2" = Phi2, "Phi3" = Phi3, "iterations" = niter, "optimal lambda" = tuningpAhat, "Y" = Y, "EZ" = EZ, "EW" = EW, "computation time" = alg.time, "CV time" = cv.time, "BICopt" = BICopt)
  return(output)
}

## end of code
