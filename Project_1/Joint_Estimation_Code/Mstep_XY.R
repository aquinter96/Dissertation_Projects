## Mstep_XY.R

Mstep_XY <- function(Data, Old_Par, E_estimates, tuningpB, tuningpA, weightsB, weightsA){
  
  X <- Data$X
  Y <- Data$Y
  A0 <- Old_Par$A0
  A <- Old_Par$A
  B0 <- Old_Par$B0
  B <- Old_Par$B
  Gamma <- Old_Par$Gamma
  Psi <- Old_Par$Psi
  Phi1 <- Old_Par$Phi1
  Phi2 <- Old_Par$Phi2
  Phi3 <- Old_Par$Phi3
  EW <- E_estimates$EW
  EZ <- E_estimates$EZ
  condvarW <- E_estimates$condvarW
  condvarZ <- E_estimates$condvarZ
  condvarWZ <- E_estimates$condvarWZ
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  s <- ncol(Gamma)
  m <- ncol(A)
  
  Aold <- A
  Bold <- B
  
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
        deltaijB <- (tcrossprod(EW) + n*condvarW)[j,j]
        thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0[i], n)))
        omegaijB <- (tcrossprod(EW) + n*condvarW)[-j,j]%*%B[i,-j]
        bbar <- ((thetaijB-omegaijB))/(deltaijB)
        if(weightsB[i,j] == 0){
          lambdaB <- (tuningpB*Phi1[i,i])/deltaijB*(1/1e-5)
        }
        else{
          lambdaB <- (tuningpB*Phi1[i,i])/deltaijB*(1/abs(weightsB[i,j]))
        }
        Bold[i,j] <- lasso(bbar, lambdaB)
      }
    }
  }
  #calculate coordinate descent updates for remaining (q-(s+1))xs submatrix of B
  for(i in (s+1):q){
    for(j in 1:s){
      deltaijB <- 0
      thetaijB <- 0
      omegaijB <- 0
      deltaijB <- (tcrossprod(EW) + n*condvarW)[j,j]
      thetaijB <- crossprod(EW[j,],(X[,i]-rep(B0[i], n)))
      omegaijB <- (tcrossprod(EW) + n*condvarW)[-j,j]%*%B[i,-j]
      bbar <- ((thetaijB-omegaijB))/(deltaijB)
      if(weightsB[i,j] == 0){
        lambdaB <- (tuningpB*Phi1[i,i])/deltaijB*(1/1e-5)
      }
      else{
        lambdaB <- (tuningpB*Phi1[i,i])/deltaijB*(1/abs(weightsB[i,j]))
      }
      Bold[i,j] <- lasso(bbar, lambdaB)
    }
  }
  
  B0old <- (1/n)*as.matrix(colSums(X - t(Bold%*%EW)))
  if(s == 1){
    Psiold <- (1/n)*(n*condvarW + EW%*%t(EW))
  }else{
    Psiold <- (1/n)*diag(diag(n*condvarW + EW%*%t(EW)))
  }
  Phi1old <- (1/n)*diag(diag(t(X)%*%X - 2*t(X)%*%matrix(B0old, nrow = n, ncol = q, byrow = T) - 2*t(X)%*%t(EW)%*%t(Bold) + matrix(B0old, nrow = q, ncol = n, byrow = F)%*%t(matrix(B0old, nrow = q, ncol = n, byrow = F)) + 
                              2*matrix(B0old, nrow = q, ncol = n, byrow = F)%*%t(EW)%*%t(Bold) + n*Bold%*%condvarW%*%t(Bold) + Bold%*%EW%*%t(EW)%*%t(Bold)))
  
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
        deltaijA <- (tcrossprod(EZ) + n*condvarZ)[j,j]
        thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0[i,], n)))
        omegaijA <- (tcrossprod(EZ) + n*condvarZ)[-j,j]%*%A[i,-j]
        abar <- ((thetaijA-omegaijA))/(deltaijA)
        if(weightsA[i,j] == 0){
          lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/1e-5)
        }
        else{
          lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/abs(weightsA[i,j]))
        }
        Aold[i,j] <- lasso(abar, lambdaA)
      }
    }
  }
  for(i in (m+1):p){
    for(j in 1:m){
      deltaijA <- 0
      thetaijA <- 0
      omegaijA <- 0
      deltaijA <- (tcrossprod(EZ) + n*condvarZ)[j,j]
      thetaijA <- crossprod(EZ[j,],(Y[,i]-rep(A0[i,], n)))
      omegaijA <- (tcrossprod(EZ) + n*condvarZ)[-j,j]%*%A[i,-j]
      abar <- ((thetaijA-omegaijA))/(deltaijA)
      if(weightsA[i,j] == 0){
        lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/1e-5)
      }
      else{
        lambdaA <- (tuningpA*Phi2[i,i])/deltaijA*(1/abs(weightsA[i,j]))
      }
      Aold[i,j] <- lasso(abar, lambdaA)
    }
  }
  
  #################################################################################
  ## back to the while loop to update A0/Gamma/Phi2/Phi3
  #################################################################################
  
  A0old <- (1/n)*as.matrix(colSums(Y - t(Aold%*%EZ)))
  
  Phi2old <- (1/n)*diag(diag(t(Y)%*%Y - 2*t(Y)%*%matrix(A0old, nrow = n, ncol = p, byrow = T) - 2*t(Y)%*%t(EZ)%*%t(Aold) + matrix(A0old, nrow = p, ncol = n, byrow = F)%*%t(matrix(A0old, nrow = p, ncol = n, byrow = F)) + 
                              2*matrix(A0old, nrow = p, ncol = n, byrow = F)%*%t(EZ)%*%t(Aold) + n*Aold%*%condvarZ%*%t(Aold) + Aold%*%EZ%*%t(EZ)%*%t(Aold)))
  
  Gammaold <- as.matrix((n*condvarWZ[((s+1):(s+m)),(1:s)] + EZ%*%t(EW))%*%solve(n*condvarWZ[(1:s),(1:s)] + EW%*%t(EW)))
  
  #Phi3 Update
  if(m == 1){
    Phi3old <- (1/n)*as.matrix(n*condvarWZ[((s+1):(s+m)),((s+1):(s+m))] + EZ%*%t(EZ) - 2*n*condvarWZ[((s+1):(s+m)),(1:s)]%*%t(Gammaold) -
                                2*EZ%*%t(EW)%*%t(Gammaold) + n*Gammaold%*%condvarWZ[(1:s),(1:s)]%*%t(Gammaold) + Gammaold%*%EW%*%t(EW)%*%t(Gammaold))
  }else{
    Phi3old <- (1/n)*as.matrix(diag(diag(n*condvarWZ[((s+1):(s+m)),((s+1):(s+m))] + EZ%*%t(EZ) - 2*n*condvarWZ[((s+1):(s+m)),(1:s)]%*%t(Gammaold) -
                                          2*EZ%*%t(EW)%*%t(Gammaold) + n*Gammaold%*%condvarWZ[(1:s),(1:s)]%*%t(Gammaold) + Gammaold%*%EW%*%t(EW)%*%t(Gammaold))))
  }
  
  ests$B0 <- B0old
  ests$B <- Bold
  ests$Psi <- Psiold
  ests$Phi1 <- Phi1old
  ests$A0 <- A0old
  ests$A <- Aold
  ests$Gamma <- Gammaold
  ests$Phi2 <- Phi2old
  ests$Phi3 <- Phi3old
  
  return(ests)
  
}
  
## end of code