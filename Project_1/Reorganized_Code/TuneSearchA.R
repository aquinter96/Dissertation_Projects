TuneSearchA.R

tunesearchA <- function(tuningpA, X, B0, B, Y, A0, A, Gamma, Phi1, Phi2, Phi3, weights){
  
  ACV <- A
  A0CV <- A0
  Phi2CV <- Phi2
  Phi3CV <- Phi3
  GammaCV <- Gamma
  
  q <- nrow(B)
  p <- nrow(A)
  xk_sig <- ncol(Gamma)
  k_sig <- nrow(Gamma)
  n <- nrow(Y)
  
  A <- matrix(0, nrow = p, ncol = k_sig)
  A0 <- rep(0, p)
  Gamma <- matrix(0, nrow = k_sig, ncol = xk_sig)
  Phi2 <- matrix(0, nrow = p, ncol = p)
  Phi3 <- matrix(0, nrow = k_sig, ncol = k_sig)
  
  deltaijA <- 0
  thetaijA <- 0
  omegaijA <- 0
  
  maxit <- 0
  
  ##########################################
  ## This while loop is for calculating estimates of 
  ## A0/A/Gamma/Phi2/Phi3 for the given tuning
  ## parameter value tuningpA
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
    EW <- matrix(EWZ[1:xk_sig,], nrow = xk_sig)   ## q x n
    EZ <- matrix(EWZ[xk_sig+(1:k_sig), ], nrow = k_sig)  ## p x n
    
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
            lambdaA <- (tuningpA*Phi2CV[i,i])/deltaijA*(1/1e-5)
          }
          else{
            lambdaA <- (tuningpA*Phi2CV[i,i])/deltaijA*(1/abs(weights[i,j]))
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
          lambdaA <- (tuningpA*Phi2CV[i,i])/deltaijA*(1/1e-5)
        }
        else{
          lambdaA <- (tuningpA*Phi2CV[i,i])/deltaijA*(1/abs(weights[i,j]))
        }
        ACV[i,j] <- lasso(abar, lambdaA)
      }
    }
    
    print(ACV)
    #once new updates for every entry of B are calculated, calculate Frobenius norm of new B matrix for comparison to norm of previous B matrix
    #to check for convergence
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
  
  sigma21 <-  rbind( matrix( cbind( t(B), t(GammaCV)%*% t(ACV)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(GammaCV%*%t(B),(Phi3 + GammaCV%*%t(GammaCV))%*%t(ACV)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = GammaCV W + E
  
  sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), ACV%*%GammaCV%*%t(B)), rbind(B%*%t(ACV%*%GammaCV),
                                                                                     (Phi2CV + ACV%*%(Phi3 + tcrossprod(GammaCV))%*%t(ACV)))), nrow = q+p, ncol = q+p)
  
  sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(GammaCV) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( GammaCV  , Phi3 + GammaCV%*%t(GammaCV) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
  
  inv11 <- Matrix::solve(sigma11)
  condvarWZ <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
  condvar  <- condvarWZ[xk_sig+(1:k_sig), xk_sig+(1:k_sig) ]  ## m x m    
  
  log.lik <- -(n/2)*(sum(diag(inv11%*%(cov(cbind(X,Y))+(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))%*%t(colMeans(cbind(X,Y)) - colMeans(cbind(X,Y)))))) + (q+p)*log((2*pi))+log(det(sigma11)))
  
  BICs <- log(n)*(sum(ACV != 0) - k_sig) - 2*log.lik
  
  ##########################################
  ## Return BIC corresponding to given tuning parameter
  ##########################################
  
  BICs
  
}

## end of code