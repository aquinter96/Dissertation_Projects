
Mstep_Y <- function(Data, Old_Par, E_estimates, tuningpG1, tuningpD1, tuningpG2, tuningpD2, weights){
  
  ests <- list()
  
  X1 <- Data$X1
  X2 <- Data$X2
  Y1 <- Data$Y1
  Y2 <- Data$Y2
  
  A1     <-   as.matrix(Old_Par$A1)
  A2     <-   as.matrix(Old_Par$A2)
  B1     <-   as.matrix(Old_Par$B1)
  B2     <-   as.matrix(Old_Par$B2)
  Gamma1 <-  as.matrix(Old_Par$Gamma1)
  Delta1 <-  as.matrix(Old_Par$Delta1)
  Gamma2 <-  as.matrix(Old_Par$Gamma2)
  Delta2 <-  as.matrix(Old_Par$Delta2)
  Theta  <-  as.matrix(Old_Par$Theta)
  Psi    <-  as.matrix(Old_Par$Psi)
  Pi1    <-  as.matrix(Old_Par$Pi1)
  Pi2    <-  as.matrix(Old_Par$Pi2)
  Omega1 <-  as.matrix(Old_Par$Omega1)
  Omega2 <-  as.matrix(Old_Par$Omega2)
  
  Gamma1new <-  Gamma1
  Delta1new <-  Delta1
  Gamma2new <-  Gamma2
  Delta2new <-  Delta2 
  Thetanew  <-  Theta
  Psinew    <-  Psi
  Pi1new    <-  Pi1
  Pi2new    <-  Pi2
  Omega1new <-  Omega1
  Omega2new <-  Omega2
  
  PhiR   <-   as.matrix(Old_Par$PhiR)
  PhiS1  <-   as.matrix(Old_Par$PhiS1)
  PhiS2  <-   as.matrix(Old_Par$PhiS2)
  Phi11  <-   as.matrix(Old_Par$Phi11)
  Phi12  <-   as.matrix(Old_Par$Phi12)
  Phi3   <-  as.matrix(Old_Par$Phi3)
  Phi41  <-  as.matrix(Old_Par$Phi41)
  Phi42  <-  as.matrix(Old_Par$Phi42)
  Phi21  <-  as.matrix(Old_Par$Phi21)
  Phi22  <-  as.matrix(Old_Par$Phi22)
  
  Phi3new   <-  Phi3
  Phi41new  <-  Phi41
  Phi42new  <-  Phi42
  Phi21new  <-  Phi21
  Phi22new  <-  Phi22
  
  m <- ncol(as.matrix(Old_Par$PhiR))
  u1 <- ncol(as.matrix(Old_Par$PhiS1))
  u2 <-  ncol(as.matrix(Old_Par$PhiS2))   
  q1 <-  ncol(as.matrix(Old_Par$Phi11))
  q2 <-  ncol(as.matrix(Old_Par$Phi12))
  p1 <-  ncol(as.matrix(Old_Par$Phi21))
  p2 <- ncol(as.matrix(Old_Par$Phi22))
  j <- ncol(as.matrix(Old_Par$Phi3))
  z1 <- ncol(as.matrix(Old_Par$Phi41))
  z2 <- ncol(as.matrix(Old_Par$Phi42))     
  n <- nrow(Y1)
  
  weights <- lapply(weights, as.matrix)
  
  condvar <- E_estimates$condvar
  
  condvarR <- condvar[(1:m), (1:m)]
  condvarS <- condvar[(m+(1:(u1+u2))), (m+(1:(u1+u2)))]
  condvarV <- condvar[(m+u1+u2+(1:j)), (m+u1+u2+(1:j))]
  condvarW1 <- condvar[(m+u1+u2+j+(1:z1)), (m+u1+u2+j+(1:z1))]
  condvarW2 <- condvar[(m+u1+u2+j+z1+(1:z2)), (m+u1+u2+j+z1+(1:z2))]
  if(m == 1){
    condvarRS <- rbind(condvar[(1:m), (m+(1:(u1+u2)))])
    condvarRW1 <- rbind(condvar[(1:m), (m+u1+u2+j+(1:z1))])
    condvarRW2 <- rbind(condvar[(1:m), (m+u1+u2+j+z1+(1:z2))])
    
  }else{
    condvarRS <- condvar[(1:m), (m+(1:(u1+u2)))]
    condvarRW1 <- condvar[(1:m), (m+u1+u2+j+(1:z1))]
    condvarRW2 <- condvar[(1:m), (m+u1+u2+j+z1+(1:z2))]
  }
  condvarSW1 <- condvar[(m+(1:(u1+u2))), (m+u1+u2+j+(1:z1))]
  condvarSW2 <- condvar[(m+(1:(u1+u2))), (m+u1+u2+j+z1+(1:z2))]
  if(j == 1){
    condvarVR <- rbind(condvar[(m+u1+u2+(1:j)), (1:m)])
    condvarVS <- rbind(condvar[(m+u1+u2+(1:j)), (m+(1:(u1+u2)))])
    condvarVW1 <- rbind(condvar[(m+u1+u2+(1:j)), (m+u1+u2+j+(1:z1))])
    condvarVW2 <- rbind(condvar[(m+u1+u2+(1:j)), (m+u1+u2+j+z1+(1:z2))])
  }else{
    condvarVR <- condvar[(m+u1+u2+(1:j)), (1:m)]
    condvarVS <- condvar[(m+u1+u2+(1:j)), (m+(1:(u1+u2)))]
    condvarVW1 <- condvar[(m+u1+u2+(1:j)), (m+u1+u2+j+(1:z1))]
    condvarVW2 <- condvar[(m+u1+u2+(1:j)), (m+u1+u2+j+z1+(1:z2))]
  }
  
  ERSVW <- E_estimates$ERSVW
  
  if(m == 1){
    ER <- t(as.matrix(ERSVW[1:m,]))
  }else{
    ER <- ERSVW[1:m,]
  }
  
  if(u1 == 1){
    ES1 <- t(as.matrix(ERSVW[m+(1:u1),]))
  }else{
    ES1 <- ERSVW[m+(1:u1),]
  }
  
  if(u2 == 1){
    ES2 <- t(as.matrix(ERSVW[m+u1+(1:u2),]))
  }else{
    ES2 <- ERSVW[m+u1+(1:u2),]
  }
  
  ES <- ERSVW[m+(1:(u1+u2)),]
  
  if(j == 1){
    EV <- t(as.matrix(ERSVW[m+u1+u2+(1:j),]))
  }else{
    EV <- ERSVW[m+u1+u2+(1:j),]
  }
  
  if(z1 == 1){
    EW1 <- t(as.matrix(ERSVW[m+u1+u2+j+(1:z1),]))
  }else{
    EW1 <- ERSVW[m+u1+u2+j+(1:z1),]
  }
  
  if(z2 == 1){
    EW2 <- t(as.matrix(ERSVW[m+u1+u2+j+z1+(1:z2),]))
  }else{
    EW2 <- ERSVW[m+u1+u2+j+z1+(1:z2),]
  }
  
  ##################################              
  ## for each element in the first part of Gamma1: mxm  
  ## Gamma1 is p1*j
  ##################################              
  for(i in 1:p1){
    for(k in 1:j){
      ## coordinate descent updates 
      ## calculate coordinate descent updates for Gamma1
      if(i > k){
        epsilonijG1 <- (tcrossprod(EV) + n*condvarV)[k,k]
        thetaijG1 <- crossprod(EV[k,],Y1[,i])
        omegaijG1 <- (tcrossprod(EV) + n*condvarV)[-k,k]%*%Gamma1[i,-k]
        tauijG1 <-  c(t(EV%*%t(EW1) + n*condvarVW1)[,k])%*%c(Delta1[i,])
        g1bar <- ((thetaijG1-omegaijG1-tauijG1))/(epsilonijG1)
        if(weights$Gamma1[i,k] == 0){
          lambdaG <- (tuningpG1*Phi21[i,i])/epsilonijG1*(1/1e-5)
        }
        else{
          lambdaG <- (tuningpG1*Phi21[i,i])/epsilonijG1*(1/abs(weights$Gamma1[i,k]))
        }
        Gamma1new[i,k] <- lasso(g1bar, lambdaG)
      }
    }
  }
  
  #calculate coordinate descent updates for the q2xm matrix of Gamma2 
  for(i in 1:p2){
    for(k in 1:j){
      epsilonijG2 <- (tcrossprod(EV) + n*condvarV)[k,k]
      thetaijG2 <- crossprod(EV[k,],Y2[,i])
      omegaijG2 <- (tcrossprod(EV) + n*condvarV)[-k,k]%*%Gamma2[i,-k]
      tauijG2 <- c(t(EV%*%t(EW2) + n*condvarVW2)[,k])%*%c(Delta2[i,])
      g2bar <- ((thetaijG2-omegaijG2-tauijG2))/(epsilonijG2)
      if(weights$Gamma2[i,k] == 0){
        lambdaG <- (tuningpG2*Phi22[i,i])/epsilonijG2*(1/1e-5)
      }
      else{
        lambdaG <- (tuningpG2*Phi22[i,i])/epsilonijG2*(1/abs(weights$Gamma2[i,k]))
      }
      Gamma2new[i,k] <- lasso(g2bar, lambdaG)
    }
  }
  
  # calculate coordinate descent updates for Delta1 
  # Delta1 is p1*z1 matrix
  for(i in 1:p1){
    for(k in 1:z1){
      if(i > k){
        epsilonijD1 <- (tcrossprod(EW1) + n*condvarW1)[k,k]
        thetaijD1 <- crossprod(EW1[k,],Y1[,i])
        omegaijD1 <- (tcrossprod(EW1) + n*condvarW1)[-k,k]%*%Delta1[i,-k]
        tauijD1 <- c(t(EW1%*%t(EV) + n*t(condvarVW1))[,k])%*%c(Gamma1[i,])
        d1bar <- ((thetaijD1-omegaijD1-tauijD1))/(epsilonijD1)
        if(weights$Delta1[i,k] == 0){
          lambdaD <- (tuningpD1*Phi21[i,i])/epsilonijD1*(1/1e-5)
        }
        else{
          lambdaD <- (tuningpD1*Phi21[i,i])/epsilonijD1*(1/abs(weights$Delta1[i,k]))
        }
        Delta1new[i,k] <- lasso(d1bar, lambdaD)
      }
    }
  }
  
  # calculate coordinate descent updates for Delta2
  # Delta1 is p2*z2 matrix
  for(i in 1:p2){
    for(k in 1:z2){
      tauijDelta2 <- 0
      if(i > k){
        epsilonijD2 <- (tcrossprod(EW2) + n*condvarW2)[k,k]
        thetaijD2 <- crossprod(EW2[k,],Y2[,i])
        omegaijD2 <- (tcrossprod(EW2) + n*condvarW2)[-k,k]%*%Delta2[i,-k]
        tauijD2 <- c(t(EW2%*%t(EV) + n*t(condvarVW2))[,k])%*%c(Gamma2[i,])
        d2bar <- ((thetaijD2-omegaijD2-tauijD2))/(epsilonijD2)
        if(weights$Delta2[i,k] == 0){
          lambdaD <- (tuningpD2*Phi22[i,i])/epsilonijD2*(1/1e-5)
        }
        else{
          lambdaD <- (tuningpD2*Phi22[i,i])/epsilonijD2*(1/abs(weights$Delta2[i,k]))
        }
        Delta2new[i,k] <- lasso(d2bar, lambdaD)
      }
    }
  }
  
  Thetanew <- (n*condvarVR + EV%*%t(ER) - n*(Old_Par$Psi)%*%t(condvarRS) - (Old_Par$Psi)%*%ES%*%t(ER))%*%solve(n*condvarR + ER%*%t(ER))
  Pi1new <- matrix(0, nrow = z1, ncol = m)
  Pi2new <- (n*t(condvarRW2) + EW2%*%t(ER) - n*Old_Par$Omega2%*%t(condvarRS) - Old_Par$Omega2%*%ES%*%t(ER))%*%solve(n*condvarR + ER%*%t(ER))
  Psinew <- (n*condvarVS + EV%*%t(ES) - n*Thetanew%*%condvarRS - Thetanew%*%ER%*%t(ES))%*%solve(n*condvarS + ES%*%t(ES))
  Omega1new <- matrix(0, nrow = z1, ncol = u1 + u2)
  Omega2new <- (n*t(condvarSW2) + EW2%*%t(ES) - n*Pi2new%*%condvarRS - Pi2new%*%ER%*%t(ES))%*%solve(n*condvarS + ES%*%t(ES))
  
  Phi21new <- 1/n*diag(diag(t(Data$Y1)%*%Data$Y1 - 2*t(Data$Y1)%*%t(Gamma1new%*%EV) - 2*t(Data$Y1)%*%t(Delta1new%*%EW1) + 2*n*Gamma1new%*%condvarVW1%*%t(Delta1new)
                            + 2*Gamma1new%*%EV%*%t(Delta1new%*%EW1) + n*Gamma1new%*%condvarV%*%t(Gamma1new) + Gamma1new%*%EV%*%t(Gamma1new%*%EV) + n*Delta1new%*%condvarW1%*%t(Delta1new) + Delta1new%*%EW1%*%t(Delta1new%*%EW1)))
  Phi22new <- 1/n*diag(diag(t(Data$Y2)%*%Data$Y2 - 2*t(Data$Y2)%*%t(Gamma2new%*%EV) - 2*t(Data$Y2)%*%t(Delta2new%*%EW2) + 2*n*Gamma2new%*%condvarVW2%*%t(Delta2new)
                            + 2*Gamma2new%*%EV%*%t(Delta2new%*%EW2) + n*Gamma2new%*%condvarV%*%t(Gamma2new) + Gamma2new%*%EV%*%t(Gamma2new%*%EV) + n*Delta2new%*%condvarW2%*%t(Delta2new) + Delta2new%*%EW2%*%t(Delta2new%*%EW2)))
  if(j == 1){
    Phi3new <- 1/n*(n*condvarV + EV%*%t(EV) - 2*n*condvarVR%*%t(Thetanew) - 2*EV%*%t(Thetanew%*%ER) - 2*n*condvarVS%*%t(Psinew) - 2*EV%*%t(Psinew%*%ES) + n*Thetanew%*%condvarR%*%t(Thetanew) + Thetanew%*%ER%*%t(Thetanew%*%ER)
                    + 2*n*Thetanew%*%condvarRS%*%t(Psinew) + 2*Thetanew%*%ER%*%t(Psinew%*%ES) + n*Psinew%*%condvarS%*%t(Psinew) + Psinew%*%ES%*%t(Psinew%*%ES))
  }else{
    Phi3new <- 1/n*diag(diag(n*condvarV + EV%*%t(EV) - 2*n*condvarVR%*%t(Thetanew) - 2*EV%*%t(Thetanew%*%ER) - 2*n*condvarVS%*%t(Psinew) - 2*EV%*%t(Psinew%*%ES) + n*Thetanew%*%condvarR%*%t(Thetanew) + Thetanew%*%ER%*%t(Thetanew%*%ER)
                             + 2*n*Thetanew%*%condvarRS%*%t(Psinew) + 2*Thetanew%*%ER%*%t(Psinew%*%ES) + n*Psinew%*%condvarS%*%t(Psinew) + Psinew%*%ES%*%t(Psinew%*%ES)))
  }
  if(z1 == 1){
    Phi41new <- 1/n*(n*condvarW1 + EW1%*%t(EW1) - 2*n*t(Pi1new%*%condvarRW1) -2*EW1%*%t(Pi1new%*%ER) -2*n*t(Omega1new%*%condvarSW1) - 2*EW1%*%t(Omega1new%*%ES) + n*Pi1new%*%condvarR%*%t(Pi1new) + Pi1new%*%ER%*%t(Pi1new%*%ER)
                     + 2*n*Pi1new%*%condvarRS%*%t(Omega1new) + 2*Pi1new%*%ER%*%t(Omega1new%*%ES) + n*Omega1new%*%condvarS%*%t(Omega1new) + Omega1new%*%ES%*%t(Omega1new%*%ES))
  }else{
    Phi41new <- 1/n*diag(diag(n*condvarW1 + EW1%*%t(EW1) - 2*n*t(Pi1new%*%condvarRW1) -2*EW1%*%t(Pi1new%*%ER) -2*n*t(Omega1new%*%condvarSW1) - 2*EW1%*%t(Omega1new%*%ES) + n*Pi1new%*%condvarR%*%t(Pi1new) + Pi1new%*%ER%*%t(Pi1new%*%ER)
                              + 2*n*Pi1new%*%condvarRS%*%t(Omega1new) + 2*Pi1new%*%ER%*%t(Omega1new%*%ES) + n*Omega1new%*%condvarS%*%t(Omega1new) + Omega1new%*%ES%*%t(Omega1new%*%ES)))
  }
  if(z2 == 1){
    Phi42new <- 1/n*(n*condvarW2 + EW2%*%t(EW2) - 2*n*t(Pi2new%*%condvarRW2) -2*EW2%*%t(Pi2new%*%ER) -2*n*t(Omega2new%*%condvarSW2) - 2*EW2%*%t(Omega2new%*%ES) + n*Pi2new%*%condvarR%*%t(Pi2new) + Pi2new%*%ER%*%t(Pi2new%*%ER)
                     + 2*n*Pi2new%*%condvarRS%*%t(Omega2new) + 2*Pi2new%*%ER%*%t(Omega2new%*%ES) + n*Omega2new%*%condvarS%*%t(Omega2new) + Omega2new%*%ES%*%t(Omega2new%*%ES))
  }else{
    Phi42new <- 1/n*diag(diag(n*condvarW2 + EW2%*%t(EW2) - 2*n*t(Pi2new%*%condvarRW2) -2*EW2%*%t(Pi2new%*%ER) -2*n*t(Omega2new%*%condvarSW2) - 2*EW2%*%t(Omega2new%*%ES) + n*Pi2new%*%condvarR%*%t(Pi2new) + Pi2new%*%ER%*%t(Pi2new%*%ER)
                              + 2*n*Pi2new%*%condvarRS%*%t(Omega2new) + 2*Pi2new%*%ER%*%t(Omega2new%*%ES) + n*Omega2new%*%condvarS%*%t(Omega2new) + Omega2new%*%ES%*%t(Omega2new%*%ES)))
  }
  
  
  ests$Gamma1 <- Gamma1new
  ests$Delta1 <- Delta1new
  ests$Gamma2 <- Gamma2new
  ests$Delta2 <- Delta2new
  
  ests$Theta <- Thetanew
  ests$Psi <- Psinew
  ests$Pi1 <- Pi1new
  ests$Omega1 <- Omega1new
  ests$Pi2 <- Pi2new
  ests$Omega2 <-Omega2new
  
  ests$Phi21 <- Phi21new
  ests$Phi22 <- Phi22new
  ests$Phi3 <- Phi3new
  ests$Phi41 <- Phi41new
  ests$Phi42 <- Phi42new
  
  ests$A1 <- Old_Par$A1
  ests$B1 <- Old_Par$B1
  ests$A2 <- Old_Par$A2
  ests$B2 <- Old_Par$B2
  ests$PhiR <- Old_Par$PhiR
  ests$PhiS1 <- Old_Par$PhiS1
  ests$PhiS2 <- Old_Par$PhiS2
  ests$Phi11 <- Old_Par$Phi11
  ests$Phi12 <- Old_Par$Phi12
  
  
  return(ests)
  
}