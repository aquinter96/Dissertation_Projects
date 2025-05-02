## Mstep_M.R

Mstep_M <- function(Data, Old_Par, E_estimates, tuningpA, weights){
  
  X <- Data$X
  M <- Data$M
  A0 <- Old_Par$A0
  A <- Old_Par$A
  Gamma <- Old_Par$Gamma
  Phi2 <- Old_Par$Phi2
  Phi3 <- Old_Par$Phi3
  EV <- E_estimates$EV
  EW <- E_estimates$EW
  condvar <- E_estimates$condvar
  condvarVW <- E_estimates$condvarVW
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(M)
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
        deltaijA <- (tcrossprod(EW) + n*condvar)[j,j]
        thetaijA <- crossprod(EW[j,],(M[,i]-rep(A0[i,], n)))
        omegaijA <- (tcrossprod(EW) + n*condvar)[-j,j]%*%A[i,-j]
        abar <- ((thetaijA-omegaijA))/(deltaijA)
        if(weights[i,j] == 0){
          lambdaA <- (tuningpA*Phi3[i,i])/deltaijA*(1/1e-5)
        }
        else{
          lambdaA <- (tuningpA*Phi3[i,i])/deltaijA*(1/abs(weights[i,j]))
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
      deltaijA <- (tcrossprod(EW) + n*condvar)[j,j]
      thetaijA <- crossprod(EW[j,],(M[,i]-rep(A0[i,], n)))
      omegaijA <- (tcrossprod(EW) + n*condvar)[-j,j]%*%A[i,-j]
      abar <- ((thetaijA-omegaijA))/(deltaijA)
      if(weights[i,j] == 0){
        lambdaA <- (tuningpA*Phi3[i,i])/deltaijA*(1/1e-5)
      }
      else{
        lambdaA <- (tuningpA*Phi3[i,i])/deltaijA*(1/abs(weights[i,j]))
      }
      ACV[i,j] <- lasso(abar, lambdaA)
    }
  }
  
  #################################################################################
  ## back to the while loop to update A0/Gamma/Phi2/Phi3
  #################################################################################
  
  A0CV <- (1/n)*as.matrix(colSums(M - t(ACV%*%EW)))
  
  Phi3CV <- (1/n)*diag(diag(t(M)%*%M - 2*t(M)%*%matrix(A0CV, nrow = n, ncol = p, byrow = T) - 2*t(M)%*%t(EW)%*%t(ACV) + matrix(A0CV, nrow = p, ncol = n, byrow = F)%*%t(matrix(A0CV, nrow = p, ncol = n, byrow = F)) + 
                              2*matrix(A0CV, nrow = p, ncol = n, byrow = F)%*%t(EW)%*%t(ACV) + n*ACV%*%condvar%*%t(ACV) + ACV%*%EW%*%t(EW)%*%t(ACV)))
  
  GammaCV <- as.matrix((n*condvarVW[((s+1):(s+m)),(1:s)] + EW%*%t(EV))%*%solve(n*condvarVW[(1:s),(1:s)] + EV%*%t(EV)))
  
  #Phi3 Update
  if(m == 1){
    Phi2CV <- (1/n)*as.matrix(n*condvarVW[((s+1):(s+m)),((s+1):(s+m))] + EW%*%t(EW) - 2*n*condvarVW[((s+1):(s+m)),(1:s)]%*%t(GammaCV) -
                                2*EW%*%t(EV)%*%t(GammaCV) + n*GammaCV%*%condvarVW[(1:s),(1:s)]%*%t(GammaCV) + GammaCV%*%EV%*%t(EV)%*%t(GammaCV))
  }else{
    Phi2CV <- (1/n)*as.matrix(diag(diag(n*condvarVW[((s+1):(s+m)),((s+1):(s+m))] + EW%*%t(EW) - 2*n*condvarVW[((s+1):(s+m)),(1:s)]%*%t(GammaCV) -
                                          2*EW%*%t(EV)%*%t(GammaCV) + n*GammaCV%*%condvarVW[(1:s),(1:s)]%*%t(GammaCV) + GammaCV%*%EV%*%t(EV)%*%t(GammaCV))))
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