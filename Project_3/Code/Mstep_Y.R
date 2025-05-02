## Mstep_M.R

Mstep_Y <- function(Data, Old_Par, E_estimates, tuningpC, tuningpD, weights){
  
  X <- Data$X
  M <- Data$M
  Y <- Data$Y
  C0 <- Old_Par$C0
  C <- Old_Par$C
  D <- Old_Par$D
  Delta <- Old_Par$Delta
  Phi4 <- Old_Par$Phi4
  Phi5 <- Old_Par$Phi5
  EV <- E_estimates$EV
  EW <- E_estimates$EW
  EZ <- E_estimates$EZ
  condvarV <- E_estimates$condvarV
  condvarZ <- E_estimates$condvarZ
  condvarVZ <- E_estimates$condvarVZ
  condvarVWZ <- E_estimates$condvarVWZ
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(M)
  r <- ncol(Y)
  s <- ncol(Old_Par$Gamma)
  m <- ncol(Old_Par$A)
  t <- ncol(C)
  
  Cnew <- C
  Dnew <- D
  
  #ests <- list()
  ests <- Old_Par
  
  ##################################              
  ## for each element in the first part of A: mxm  
  ## S is p*m
  ##################################              
  for(i in 1:r){
    for(j in 1:t){
      ## coordinate descent updates
      ## calculate coordinate descent updates for A1
      if(i > j){
        epsilonijC <- (tcrossprod(EZ) + n*condvarZ)[j,j]
        thetaijC <- crossprod(EZ[j,],Y[,i])
        omegaijC <- (tcrossprod(EZ) + n*condvarZ)[-j,j]%*%as.matrix(C[i,-j])
        tauijC <-  c(t(EZ%*%t(EV) + n*t(condvarVZ))[,j])%*%c(as.matrix(D[i,]))
        cbar <- ((thetaijC-omegaijC-tauijC))/(epsilonijC)
        if(as.matrix(weights$C)[i,j] == 0){
          lambdaC <- (tuningpC*Phi5[i,i])/epsilonijC*(1/1e-5)
        }
        else{
          lambdaC <- (tuningpC*Phi5[i,i])/epsilonijC*(1/abs(as.matrix(weights$C)[i,j]))
        }
        Cnew[i,j] <- lasso(cbar, lambdaC)
      }
    }
  }

  for(i in 1:r){
    for(j in 1:s){
    ## coordinate descent updates
    ## calculate coordinate descent updates for A1
      epsilonijD <- (tcrossprod(EV) + n*condvarV)[j,j]
      thetaijD <- crossprod(EV[j,],Y[,i])
      omegaijD <- (tcrossprod(EV) + n*condvarV)[-j,j]%*%as.matrix(D[i,-j])
      tauijD <-  c(t(EV%*%t(EZ) + n*condvarVZ)[,j])%*%c(as.matrix(C[i,]))
      dbar <- ((thetaijD-omegaijD-tauijD))/(epsilonijD)
      if(as.matrix(weights$D)[i,j] == 0){
        lambdaD <- (tuningpD*Phi5[i,i])/epsilonijD*(1/1e-5)
      }
      else{
        lambdaD <- (tuningpD*Phi5[i,i])/epsilonijD*(1/abs(as.matrix(weights$D)[i,j]))
      }
      Dnew[i,j] <- lasso(dbar, lambdaD)
    }
  }

  #################################################################################
  ## back to the while loop to update A0/Gamma/Phi2/Phi3
  #################################################################################

  C0new <- (1/n)*as.matrix(colSums(Y - t(Cnew%*%EZ) - t(Dnew%*%EV)))

  Phi5new <- (1/n)*diag(diag(t(Y)%*%Y - 2*t(Y)%*%matrix(C0new, nrow = n, ncol = r, byrow = T) - 2*t(Y)%*%t(EZ)%*%t(Cnew) - 2*t(Y)%*%t(EV)%*%t(Dnew) +
                               matrix(C0new, nrow = r, ncol = n, byrow = F)%*%t(matrix(C0new, nrow = r, ncol = n, byrow = F)) +
                               2*matrix(C0new, nrow = r, ncol = n, byrow = F)%*%t(EZ)%*%t(Cnew) +
                               2*matrix(C0new, nrow = r, ncol = n, byrow = F)%*%t(EV)%*%t(Dnew) + n*Cnew%*%condvarZ%*%t(Cnew) +
                               Cnew%*%EZ%*%t(EZ)%*%t(Cnew) + n*Dnew%*%condvarV%*%t(Dnew) + Dnew%*%EV%*%t(EV)%*%t(Dnew) +
                               2*n*Cnew%*%t(condvarVZ)%*%t(Dnew) + 2*Cnew%*%EZ%*%t(Dnew%*%EV)))

  Deltanew <- as.matrix((n*condvarVWZ[s+m+(1:t),s+(1:m)] + EZ%*%t(EW))%*%solve(n*condvarVWZ[s+(1:m),s+(1:m)] + EW%*%t(EW)))

  #Phi3 Update
  if(t == 1){
    Phi4new <- (1/n)*as.matrix(n*condvarVWZ[s+m+(1:t), s+m+(1:t)] + EZ%*%t(EZ) - 2*n*condvarVWZ[s+m+(1:t),s+(1:m)]%*%t(Deltanew) -
                                2*EZ%*%t(EW)%*%t(Deltanew) + n*Deltanew%*%condvarVWZ[s+(1:m),s+(1:m)]%*%t(Deltanew) + Deltanew%*%EW%*%t(EW)%*%t(Deltanew))
  }else{
    Phi4new <- (1/n)*as.matrix(diag(diag(n*condvarVWZ[s+m+(1:t), s+m+(1:t)] + EZ%*%t(EZ) - 2*n*condvarVWZ[s+m+(1:t),s+(1:m)]%*%t(Deltanew) -
                                          2*EZ%*%t(EW)%*%t(Deltanew) + n*Deltanew%*%condvarVWZ[s+(1:m),s+(1:m)]%*%t(Deltanew) + Deltanew%*%EW%*%t(EW)%*%t(Deltanew))))
  }
  
  ests$C0 <- C0new
  ests$C <- Cnew
  ests$D <- Dnew
  ests$Delta <- Deltanew
  ests$Phi4 <- Phi4new
  ests$Phi5 <- Phi5new
  ests$A0 <- Old_Par$A0
  ests$A <- Old_Par$A
  ests$Gamma <- Old_Par$Gamma
  ests$Phi2 <- Old_Par$Phi2
  ests$Phi3 <- Old_Par$Phi3
  ests$B <- Old_Par$B
  ests$B0 <- Old_Par$B0
  ests$Phi1 <- Old_Par$Phi1
  ests$Psi <- Old_Par$Psi
  
  return(ests)
  
}

## end of code