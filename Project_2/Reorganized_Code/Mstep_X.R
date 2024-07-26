## Mstep_X.R

Mstep_X <- function(Data, Old_Par, E_estimates, tuningpA, tuningpB, weights){
  
  ests <- list()
  
  X1 <- Data$X1
  X2 <- Data$X2
  
  PhiR   <-   Old_Par$PhiR
  PhiS1  <-   Old_Par$PhiS1
  PhiS2  <-   Old_Par$PhiS2
  Phi11  <-   Old_Par$Phi11
  Phi12  <-   Old_Par$Phi12
  A1     <-   Old_Par$A1
  A2     <-   Old_Par$A2
  B1     <-   Old_Par$B1
  B2     <-   Old_Par$B2
  
  A1new <- A1
  A2new <- A2
  B1new <- B1
  B2new <- B2
  
  PhiS1S <- cbind(PhiS1, matrix(0, nrow = u1, ncol = u2))
  PhiS2S <- cbind(matrix(0, nrow = u2, ncol = u1), PhiS2)
  PhiS <- diag(c(diag(PhiS1), diag(PhiS2)))
  
  m <- ncol(Old_Par$PhiR)
  u1 <- ncol(Old_Par$PhiS1)
  u2 <-  ncol(Old_Par$PhiS2)   
  q1 <-  ncol(Old_Par$Phi11)
  q2 <-  ncol(Old_Par$Phi12)
  n <- ncol(Data$X1) ## sample size
  
  condvar <- E_estimates$condvar
  
  condvarR <- condvar[(1:m), (1:m)]
  condvarS1 <- condvar[(m+(1:u1)), (m+(1:u1))]
  condvarS2 <- condvar[(m+u1+(1:u2)), (m+u1+(1:u2))]
  condvarRS1 <- condvar[(1:m), (m+(1:u1))]
  condvarRS2 <- condvar[(1:m), (m+u1+(1:u2))]
  
  ERS <- E_estimates$ERS
  
  ER <- ERS[1:m,]
  ES1 <- ERS[m+(1:u1),]
  ES2 <- ERS[m+u1+(1:u2),]
  
  ##################################              
  ## for each element in the first part of A1: mxm  
  ## A1 is q1*m
  ##################################              
  for(i in 1:q1){
    for(j in 1:m){
      ## coordinate descent updates 
      ## calculate coordinate descent updates for A1
      if(i > j){
        epsilonijA1 <- (tcrossprod(ER) + n*condvarR)[j,j]
        thetaijA1 <- crossprod(ER[j,],X1[,i])
        omegaijA1 <- (tcrossprod(ER) + n*condvarR)[-j,j]%*%A1[i,-j]
        tauijA1 <-  c(t(ER%*%t(ES1) + n*condvarRS1)[,j])%*%c(B1[i,])
        a1bar <- ((thetaijA1-omegaijA1-tauijA1))/(epsilonijA1)
        if(weightsA[[1]][i,j] == 0){
          lambdaA <- (tuningpA*Phi11[i,i])/epsilonijA1*(1/1e-5)
        }
        else{
          lambdaA <- (tuningpA*Phi11[i,i])/epsilonijA1*(1/abs(weights$A1[i,j]))
        }
        A1new[i,j] <- lasso(a1bar, lambdaA)
      }
    }
  }
  
  #calculate coordinate descent updates for the q2xm matrix of A2 
  for(i in 1:q2){
    for(j in 1:m){
      epsilonijA2 <- (tcrossprod(ER) + n*condvarR)[j,j]
      thetaijA2 <- crossprod(ER[j,],X2[,i])
      omegaijA2 <- (tcrossprod(ER) + n*condvarR)[-j,j]%*%A2[i,-j]
      tauijA2 <- c(t(ER%*%t(ES2) + n*condvarRS2)[,j])%*%c(B2[i,])
      a2bar <- ((thetaijA2-omegaijA2-tauijA2))/(epsilonijA2)
      if(weightsA[[2]][i,j] == 0){
        lambdaA <- (tuningpA*Phi12[i,i])/epsilonijA2*(1/1e-5)
      }
      else{
        lambdaA <- (tuningpA*Phi12[i,i])/epsilonijA2*(1/abs(weights$A2[i,j]))
      }
      A2new[i,j] <- lasso(a2bar, lambdaA)
    }
  }
  
  # calculate coordinate descent updates for B1 
  # B1 is q1*u1 matrix
  for(i in 1:q1){
    for(j in 1:u1){
      if(i > j){
        epsilonijB1 <- (tcrossprod(ES1) + n*condvarS1)[j,j]
        thetaijB1 <- crossprod(ES1[j,],X1[,i])
        omegaijB1 <- (tcrossprod(ES1) + n*condvarS1)[-j,j]%*%B1[i,-j]
        tauijB1 <- c(t(ES1%*%t(ER) + n*t(condvarRS1))[,j])%*%c(A1[i,])
        b1bar <- ((thetaijB1-omegaijB1-tauijB1))/(epsilonijB1)
        if(weightsB[[1]][i,j] == 0){
          lambdaB <- (tuningpB*Phi11[i,i])/epsilonijB1*(1/1e-5)
        }
        else{
          lambdaB <- (tuningpB*Phi11[i,i])/epsilonijB1*(1/abs(weights$B1[i,j]))
        }
        B1new[i,j] <- lasso(b1bar, lambdaB)
      }
    }
  }

  # calculate coordinate descent updates for B2
  # B1 is q2*u2 matrix
  for(i in 1:q2){
    for(j in 1:u2){
      tauijB2 <- 0
      if(i > j){
        epsilonijB2 <- (tcrossprod(ES2) + n*condvarS2)[j,j]
        thetaijB2 <- crossprod(ES2[j,],X2[,i])
        omegaijB2 <- (tcrossprod(ES2) + n*condvarS2)[-j,j]%*%B2[i,-j]
        tauijB2 <- c(t(ES2%*%t(ER) + n*t(condvarRS2))[,j])%*%c(A2[i,])
        b2bar <- ((thetaijB2-omegaijB2-tauijB2))/(epsilonijB2)
        if(weightsB[[2]][i,j] == 0){
          lambdaB <- (tuningpB*Phi12[i,i])/epsilonijB2*(1/1e-5)
        }
        else{
          lambdaB <- (tuningpB*Phi12[i,i])/epsilonijB2*(1/abs(weights$B2[i,j]))
        }
        B2new[i,j] <- lasso(b2bar, lambdaB)
      }
    }
  }
  
  #################################################################################
  ## update covariance matrices
  #################################################################################
  
  PhiRnew <- diag(diag(condvarR + 1/n*ER%*%t(ER)))
  PhiS1new <- diag(diag(condvarS1 + 1/n*ES1%*%t(ES1)))
  PhiS2new <- diag(diag(condvarS2 + 1/n*ES2%*%t(ES2)))
  Phi11new <- 1/n*diag(diag(t(X1)%*%X1 - 2*t(X1)%*%t(A1new%*%ER) - 2*t(X1)%*%t(B1new%*%ES1) + 2*n*A1new%*%condvarRS1%*%t(B1new)
                          + 2*A1new%*%ER%*%t(B1new%*%ES1) + n*A1new%*%condvarR%*%t(A1new) + A1new%*%ER%*%t(A1new%*%ER) + n*B1new%*%condvarS1%*%t(B1new) + B1new%*%ES1%*%t(B1new%*%ES1)))
  Phi12new <- 1/n*diag(diag(t(X2)%*%X2 - 2*t(X2)%*%t(A2new%*%ER) - 2*t(X2)%*%t(B2new%*%ES2) + 2*n*A2new%*%condvarRS2%*%t(B2new)
                          + 2*A2new%*%ER%*%t(B2new%*%ES2) + n*A2new%*%condvarR%*%t(A2new) + A2new%*%ER%*%t(A2new%*%ER) + n*B2new%*%condvarS2%*%t(B2new) + B2new%*%ES2%*%t(B2new%*%ES2)))

  ests$A1 <- A1new
  ests$A2 <- A2new
  ests$B1 <- B1new
  ests$B2 <- B2new
  
  ests$PhiR <- PhiRnew
  ests$PhiS1 <- PhiS1new
  ests$PhiS2 <- PhiS2new
  ests$Phi11 <- Phi11new
  ests$Phi12 <- Phi12new
  
  return(ests)
}

## end of code