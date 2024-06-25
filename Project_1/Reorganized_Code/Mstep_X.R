## Mstep_X.R

Mstep_X <- function(X, EW, condvar, B0, B, Phi1, tuningpB, weights){
  
  n <- nrow(X)
  q <- nrow(B)
  xk_sig <- ncol(B)
  
  BCV <- B

  ##################################              
  ## for each element in the first part of B: sxs  
  ## B is q*s
  ##################################              
  for(i in 1:xk_sig){
    for(j in 1:xk_sig){
      deltaijB <- 0
      thetaijB <- 0
      omegaijB <- 0
      #calculate coordinate descent updates for the first sxs submatrix of B, while ensuring that said sxs matrix is lower triangular
      if(j < i){
        deltaijB <- (tcrossprod(EW) + n*condvar)[j,j]
        thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0[i,], n)))
        omegaijB <- (tcrossprod(EW) + n*condvar)[-j,j]%*%B[i,-j]
        bbar <- ((thetaijB-omegaijB))/(deltaijB)
        if(weights[i,j] == 0){
          lambdaB <- (tuningpB*Phi1[i,i])/deltaijB*(1/1e-5)
        }
        else{
          lambdaB <- (tuningpB*Phi1[i,i])/deltaijB*(1/abs(weights[i,j]))
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
      thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0[i,], n)))
      omegaijB <- (tcrossprod(EW) + n*condvar)[-j,j]%*%B[i,-j]
      bbar <- ((thetaijB-omegaijB))/(deltaijB)
      if(weights[i,j] == 0){
        lambdaB <- (tuningpB*Phi1[i,i])/deltaijB*(1/1e-5)
      }
      else{
        lambdaB <- (tuningpB*Phi1[i,i])/deltaijB*(1/abs(weights[i,j]))
      }
      BCV[i,j] <- lasso(bbar, lambdaB)
    }
  }
  
  #################################################################################
  ## back to the while loop to update B0/Phi1
  #################################################################################
  
  B0CV <- (1/n)*as.matrix(colSums(X - t(BCV%*%EW)))
  Phi1CV <- (1/n)*diag(diag(t(X)%*%X - 2*t(X)%*%matrix(B0CV, nrow = n, ncol = q, byrow = T) - 2*t(X)%*%t(EW)%*%t(BCV) + matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(matrix(B0CV, nrow = q, ncol = n, byrow = F)) + 
                              2*matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(EW)%*%t(BCV) + n*BCV%*%condvar%*%t(BCV) + BCV%*%EW%*%t(EW)%*%t(BCV)))
  
  ests <- list("B0" = B0CV, "B" = BCV, "Phi1" = Phi1CV)

  return(ests)
}

## end of code