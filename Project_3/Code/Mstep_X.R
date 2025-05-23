## Mstep_X.R

Mstep_X <- function(Data, Old_Par, E_estimates, tuningpB, weights){
  
  ests <- list()
  B <- as.matrix(Old_Par$B)
  B0 <- as.matrix(Old_Par$B0)
  Psi <- as.matrix(Old_Par$Psi)
  Phi1 <- as.matrix(Old_Par$Phi1)
  X <- Data
  
  n <- nrow(X)
  q <- nrow(B)
  s <- ncol(B)
  
  EV <- E_estimates$EV
  condvar <- E_estimates$condvar
  
  BCV <- B
  
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
        deltaijB <- (tcrossprod(EV) + n*condvar)[j,j]
        thetaijB <- crossprod(EV[j,],(X[,i]-rep(B0[i], n)))
        omegaijB <- (tcrossprod(EV) + n*condvar)[-j,j]%*%B[i,-j]
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
      deltaijB <- (tcrossprod(EV) + n*condvar)[j,j]
      thetaijB <- crossprod(EV[j,],(X[,i]-rep(B0[i], n)))
      omegaijB <- (tcrossprod(EV) + n*condvar)[-j,j]%*%B[i,-j]
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
  
  B0CV <- (1/n)*as.matrix(colSums(X - t(BCV%*%EV)))
  if(s == 1){
    PsiCV <- (1/n)*(n*condvar + EV%*%t(EV))
  }else{
    PsiCV <- (1/n)*diag(diag(n*condvar + EV%*%t(EV)))
  }
  Phi1CV <- (1/n)*diag(diag(t(X)%*%X - 2*t(X)%*%matrix(B0CV, nrow = n, ncol = q, byrow = T) - 2*t(X)%*%t(EV)%*%t(BCV) + matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(matrix(B0CV, nrow = q, ncol = n, byrow = F)) + 
                              2*matrix(B0CV, nrow = q, ncol = n, byrow = F)%*%t(EV)%*%t(BCV) + n*BCV%*%condvar%*%t(BCV) + BCV%*%EV%*%t(EV)%*%t(BCV)))
  
  ests$B0 <- B0CV
  ests$B <- BCV
  ests$Psi <- PsiCV
  ests$Phi1 <- Phi1CV
  
  return(ests)
}

## end of code