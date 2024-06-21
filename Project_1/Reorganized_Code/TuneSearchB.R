# TuneSearchB.R

tunesearchB <- function(tuningpB, X, B0, B, Phi1, weights){
  
  B0CV <- B0
  BCV <- B
  Phi1CV <- Phi1
  
  q <- nrow(B)
  xk_sig <- ncol(B)
  n <- nrow(X)
  
  B0 <- rep(0, q)
  B <- matrix(0, nrow = q, ncol = xk_sig)
  Phi1 <- matrix(0, nrow = q, ncol = q)
  
  deltaijB <- 0
  thetaijB <- 0
  omegaijB <- 0
  
  maxit <- 0
  
  ##########################################
  ## This while loop is for calculating estimates of B0/B/Phi1
  ## for the given tuning parameter value tuningpB
  ##########################################
  while(( (norm(B0 - B0CV, type = "F") > 0.00001) | (norm(B - BCV, type = "F") > 0.00001) | (norm(Phi1 - Phi1CV, type = "F") > 0.00001) ) & (maxit < 15000)){
    
    B0 <- B0CV
    B <- BCV
    Phi1 <- Phi1CV
    
    invmodcv <- Matrix::solve(Phi1 + tcrossprod(B))
    condvar <- diag(xk_sig) - crossprod(B,invmodcv)%*%B
    EW <- crossprod(B,invmodcv)%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))
    
    for(i in 1:xk_sig){
      for(j in 1:xk_sig){
        deltaijB <- 0
        thetaijB <- 0
        omegaijB <- 0
        #calculate coordinate descent updates for the first sxs submatrix of B, while ensuring that said sxs matrix is lower triangular
        if(j < i){
          deltaijB <- (tcrossprod(EW) + n*condvar)[j,j]
          thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0CV[i,], n)))
          omegaijB <- (tcrossprod(EW) + n*condvar)[-j,j]%*%BCV[i,-j]
          bbar <- ((thetaijB-omegaijB))/(deltaijB)
          if(weights[i,j] == 0){
            lambdaB <- (tuningpB*Phi1CV[i,i])/deltaijB*(1/1e-5)
          }
          else{
            lambdaB <- (tuningpB*Phi1CV[i,i])/deltaijB*(1/abs(weights[i,j]))
          }
          BCV[i,j] <- lasso(bbar, lambdaB)
        }
      }
    }
    #calculate coordinate descent updates for remaining (q-(s+1))xs submatrix of B
    for(i in (xk_sig+1):q){
      for(j in 1:xk_sig){
        deltaijB <- 0
        thetaijB <- 0
        omegaijB <- 0
        deltaijB <- (tcrossprod(EW) + n*condvar)[j,j]
        thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0CV[i,], n)))
        omegaijB <- (tcrossprod(EW) + n*condvar)[-j,j]%*%BCV[i,-j]
        bbar <- ((thetaijB-omegaijB))/(deltaijB)
        if(weights[i,j] == 0){
          lambdaB <- (tuningpB*Phi1CV[i,i])/deltaijB*(1/1e-5)
        }
        else{
          lambdaB <- (tuningpB*Phi1CV[i,i])/deltaijB*(1/abs(weights[i,j]))
        }
        BCV[i,j] <- lasso(bbar, lambdaB)
      }
    }
    #once new updates for every entry of B are calculated, calculate Frobenius norm of new B matrix for comparison to norm of previous B matrix
    #to check for convergence
    B0CV <- (1/n)*as.matrix(colSums(X - t(BCV%*%EW)))
    Phi1CV <- (1/n)*diag(diag(t(X)%*%X - 2*t(X)%*%matrix(B0CV, nrow = n, ncol = q, byrow = T) - 2*t(X)%*%t(EW)%*%t(BCV) + matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(matrix(B0CV, nrow = q, ncol = n, byrow = F)) + 
                                2*matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(EW)%*%t(BCV) + n*BCV%*%condvar%*%t(BCV) + BCV%*%EW%*%t(EW)%*%t(BCV)))
    maxit <- maxit + 1
  }
  
  invmodcv <- Matrix::solve(Phi1CV + tcrossprod(BCV))
  log.lik <- -(n/2)*(sum(diag(invmodcv%*%(cov(X)+(colMeans(X)-colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
                       q*log((2*pi))+log(det(Phi1CV +tcrossprod(BCV))))
  
  BICs <- log(n)*(sum(BCV != 0) - xk_sig) - 2*log.lik
  
  ##########################################
  ## Return BIC corresponding to given tuning parameter
  ##########################################
  
  BICs
  
}

## end of code