## Mstep_X.R

Mstep_X <- function(Data, Old_Par, E_estimates, tuningpB, weights){
  
  B <- Old_Par$B
  B0 <- Old_Par$B0
  Phi1 <- Old_Par$Phi1
  EW <- E_estimates$EW
  condvar <- E_estimates$condvar
  X <- Data
  
  n <- nrow(X)
  q <- nrow(B)
  s <- ncol(B)
  
  BCV <- B
  
  ests <- list()

  ##################################              
  ## for each element in the first part of B: sxs  
  ## B is q*s
  ##################################              
  for(i in 1:s){
    for(j in 1:s){
      deltaijB <- 0
      thetaijB <- 0
      omegaijB <- 0
      #calculate coordinate descent updates for the first sxs submatrix of B, while ensuring that said sxs matrix is lower triangular
      if(j < i){
        deltaijB <- (tcrossprod(EW) + n*condvar)[j,j]
        thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0[i], n)))
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
  for(i in (s+1):q){
    for(j in 1:s){
      deltaijB <- 0
      thetaijB <- 0
      omegaijB <- 0
      deltaijB <- (tcrossprod(EW) + n*condvar)[j,j]
      thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0[i], n)))
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
  
  ests$B0 <- B0CV
  ests$B <- BCV
  ests$Phi1 <- Phi1CV

  return(ests)
}

## end of code