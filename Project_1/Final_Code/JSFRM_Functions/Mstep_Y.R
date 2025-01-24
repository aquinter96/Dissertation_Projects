## Mstep_Y.R

Mstep_Y <- function(Data, Old_Par, E_estimates, tuningpA, weights){
  
  X <- Data$X
  Y <- Data$Y
  A0 <- Old_Par$A0
  A <- Old_Par$A
  Gamma <- Old_Par$Gamma
  Phi2 <- Old_Par$Phi2
  Phi3 <- Old_Par$Phi3
  EW <- E_estimates$EW
  EZ <- E_estimates$EZ
  condvar <- E_estimates$condvar
  condvarWZ <- E_estimates$condvarWZ
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  s <- ncol(Gamma)
  m <- ncol(A)
  
  ACV <- A
  
  ests <- list()
  
  ##################################              
  ## for each element in the first part of A: mxm  
  ## S is p*m
  ##################################              
  for(i in 1:m){
    for(j in 1:m){
      deltaijA <- 0
      thetaijA <- 0
      omegaijA <- 0
      #calculate coordinate descent updates for the first sxs submatrix of B, while ensuring that said sxs matrix is lower triangular
      if(j < i){
        deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
        thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0[i,], n)))
        omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%A[i,-j]
        abar <- ((thetaijA-omegaijA))/(deltaijA)
        if(weights[i,j] == 0){
          lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/1e-5)
        }
        else{
          lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/abs(weights[i,j]))
        }
        ACV[i,j] <- lasso(abar, lambdaA)
      }
    }
  }
  for(i in (m+1):p){
    for(j in 1:m){
      deltaijA <- 0
      thetaijA <- 0
      omegaijA <- 0
      deltaijA <- (tcrossprod(EZ) + n*condvar)[j,j]
      thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0[i,], n)))
      omegaijA <- (tcrossprod(EZ) + n*condvar)[-j,j]%*%A[i,-j]
      abar <- ((thetaijA-omegaijA))/(deltaijA)
      if(weights[i,j] == 0){
        lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/1e-5)
      }
      else{
        lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/abs(weights[i,j]))
      }
      ACV[i,j] <- lasso(abar, lambdaA)
    }
  }
  
  #################################################################################
  ## back to the while loop to update A0/Gamma/Phi2/Phi3
  #################################################################################
  
  A0CV <- (1/n)*as.matrix(colSums(Y - t(ACV%*%EZ)))
  
  Phi2CV <- (1/n)*diag(diag(t(Y)%*%Y - 2*t(Y)%*%matrix(A0CV, nrow = n, ncol = p, byrow = T) - 2*t(Y)%*%t(EZ)%*%t(ACV) + matrix(A0CV, nrow = p, ncol = n, byrow = F)%*%t(matrix(A0CV, nrow = p, ncol = n, byrow = F)) + 
                              2*matrix(A0CV, nrow = p, ncol = n, byrow = F)%*%t(EZ)%*%t(ACV) + n*ACV%*%condvar%*%t(ACV) + ACV%*%EZ%*%t(EZ)%*%t(ACV)))
  
  GammaCV <- as.matrix((n*condvarWZ[((s+1):(s+m)),(1:s)] + EZ%*%t(EW))%*%solve(n*condvarWZ[(1:s),(1:s)] + EW%*%t(EW)))
  
  #Phi3 Update
  if(m == 1){
    Phi3CV <- (1/n)*as.matrix(n*condvarWZ[((s+1):(s+m)),((s+1):(s+m))] + EZ%*%t(EZ) - 2*n*condvarWZ[((s+1):(s+m)),(1:s)]%*%t(GammaCV) -
                                2*EZ%*%t(EW)%*%t(GammaCV) + n*GammaCV%*%condvarWZ[(1:s),(1:s)]%*%t(GammaCV) + GammaCV%*%EW%*%t(EW)%*%t(GammaCV))
  }else{
    Phi3CV <- (1/n)*as.matrix(diag(diag(n*condvarWZ[((s+1):(s+m)),((s+1):(s+m))] + EZ%*%t(EZ) - 2*n*condvarWZ[((s+1):(s+m)),(1:s)]%*%t(GammaCV) -
                                          2*EZ%*%t(EW)%*%t(GammaCV) + n*GammaCV%*%condvarWZ[(1:s),(1:s)]%*%t(GammaCV) + GammaCV%*%EW%*%t(EW)%*%t(GammaCV))))
  }
  
  ests$A0 <- A0CV
  ests$A <- ACV
  ests$Gamma <- GammaCV
  ests$Phi2 <- Phi2CV
  ests$Phi3 <- Phi3CV
  ests$B <- Old_Par$B
  ests$B0 <- Old_Par$B0
  ests$Phi1 <- Old_Par$Phi1
  ests$Psi <- Old_Par$Psi
  
  return(ests)
  
}
  
## end of code