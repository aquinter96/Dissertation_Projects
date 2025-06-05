## Estep_X.R

Estep_Y <- function(Data, Old_Par){
  
  ests <- list()
  
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
  
  PhiS1S <- cbind(PhiS1, matrix(0, nrow = u1, ncol = u2))
  PhiS2S <- cbind(matrix(0, nrow = u2, ncol = u1), PhiS2)
  PhiS <- diag(c(diag(PhiS1), diag(PhiS2)))

  ests$sigma21 <- rbind(cbind(PhiR%*%t(A1), PhiR%*%t(A2), PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1), PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2)),
                         cbind(PhiS1%*%t(B1), matrix(0, nrow = u1, ncol = q2), t((Gamma1%*%Psi + Delta1%*%Omega1)%*%t(PhiS1S)), t((Gamma2%*%Psi + Delta2%*%Omega2)%*%t(PhiS1S))),
                         cbind(matrix(0, nrow = u2, ncol = q1), PhiS2%*%t(B2), t((Gamma1%*%Psi + Delta1%*%Omega1)%*%t(PhiS2S)), t((Gamma2%*%Psi + Delta2%*%Omega2)%*%t(PhiS2S))),
                         cbind(t(A1%*%PhiR%*%t(Theta) + B1%*%PhiS1S%*%t(Psi)), t(A2%*%PhiR%*%t(Theta) + B2%*%PhiS2S%*%t(Psi)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Theta) + 
                                                                                                                                   (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Psi) + Gamma1%*%Phi3), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Theta) + 
                                                                                                                                                                                                          (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Psi) + Gamma2%*%Phi3)),
                         cbind(t(A1%*%PhiR%*%t(Pi1) + B1%*%PhiS1S%*%t(Omega1)), t(A2%*%PhiR%*%t(Pi1) + B2%*%PhiS2S%*%t(Omega1)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Pi1) + 
                                                                                                                                     (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Omega1) + Delta1%*%Phi41), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Pi1) + (Gamma2%*%Psi + 
                                                                                                                                                                                                                                                                   Delta2%*%Omega2)%*%PhiS%*%t(Omega1))),
                         cbind(t(A1%*%PhiR%*%t(Pi2) + B1%*%PhiS1S%*%t(Omega2)), t(A2%*%PhiR%*%t(Pi2) + B2%*%PhiS2S%*%t(Omega2)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Pi2) + 
                                                                                                                                     (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Omega2)), t((Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Pi2) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Omega2) + 
                                                                                                                                                                                               Delta2%*%Phi42)))
  
  #calculate covariance of observed data, i.e. Cov(Xi1, Xi2, Yi1, Yi2)
  
  ests$sigma11 <- rbind(cbind(A1%*%PhiR%*%t(A1) + B1%*%PhiS1%*%t(B1) + Phi11, A1%*%PhiR%*%t(A2), A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), 
                               A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                         cbind(A2%*%PhiR%*%t(A1), A2%*%PhiR%*%t(A2) + B2%*%PhiS2%*%t(B2) + Phi12, A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1), 
                               A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)),
                         cbind(t(A1%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + B1%*%PhiS1S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), t(A2%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + 
                                                                                                                                      B2%*%PhiS2S%*%t(Gamma1%*%Psi + Delta1%*%Omega1)), (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma1%*%Theta + Delta1%*%Pi1) + 
                                 (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma1%*%Psi + Delta1%*%Omega1) + Gamma1%*%Phi3%*%t(Gamma1) + Delta1%*%Phi41%*%t(Delta1) + Phi21, 
                               (Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + 
                                 Gamma1%*%Phi3%*%t(Gamma2)),
                         cbind(t(A1%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + B1%*%PhiS1S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t(A2%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + 
                                                                                                                                      B2%*%PhiS2S%*%t(Gamma2%*%Psi + Delta2%*%Omega2)), t((Gamma1%*%Theta + Delta1%*%Pi1)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + 
                                                                                                                                                                                            (Gamma1%*%Psi + Delta1%*%Omega1)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + Gamma1%*%Phi3%*%t(Gamma2)), 
                               (Gamma2%*%Theta + Delta2%*%Pi2)%*%PhiR%*%t(Gamma2%*%Theta + Delta2%*%Pi2) + (Gamma2%*%Psi + Delta2%*%Omega2)%*%PhiS%*%t(Gamma2%*%Psi + Delta2%*%Omega2) + 
                                 Gamma2%*%Phi3%*%t(Gamma2) + Delta2%*%Phi42%*%t(Delta2) + Phi22))
  
  #calculate covariance of latent factors, i.e. Cov(Ri, Si1, Si2, Vi, Wi1, Wi2)
  # 
  ests$sigma22 <- rbind(cbind(PhiR, matrix(0, nrow = m, ncol = u1), matrix(0, nrow = m, ncol = u2), PhiR%*%t(Theta), PhiR%*%t(Pi1), PhiR%*%t(Pi2)),
                         cbind(matrix(0, nrow = u1, ncol = m), PhiS1, matrix(0, nrow = u1, ncol = u2), PhiS1S%*%t(Psi), PhiS1S%*%t(Omega1), PhiS1S%*%t(Omega2)),
                         cbind(matrix(0, nrow = u2, ncol = m), matrix(0, nrow = u2, ncol = u1), PhiS2, PhiS2S%*%t(Psi), PhiS2S%*%t(Omega1), PhiS2S%*%t(Omega2)),
                         cbind(Theta%*%PhiR, Psi%*%t(PhiS1S), Psi%*%t(PhiS2S), Theta%*%PhiR%*%t(Theta) + Psi%*%PhiS%*%t(Psi) + Phi3, Theta%*%PhiR%*%t(Pi1) + Psi%*%PhiS%*%t(Omega1), 
                               Theta%*%PhiR%*%t(Pi2) + Psi%*%PhiS%*%t(Omega2)),
                         cbind(Pi1%*%PhiR, Omega1%*%t(PhiS1S), Omega1%*%t(PhiS2S), t(Theta%*%PhiR%*%t(Pi1) + Psi%*%PhiS%*%t(Omega1)), Pi1%*%PhiR%*%t(Pi1) + Omega1%*%PhiS%*%t(Omega1) + 
                                 Phi41, Pi1%*%PhiR%*%t(Pi2) + Omega1%*%PhiS%*%t(Omega2)),
                         cbind(Pi2%*%PhiR, Omega2%*%t(PhiS1S), Omega2%*%t(PhiS2S), t(Theta%*%PhiR%*%t(Pi2) + Psi%*%PhiS%*%t(Omega2)), Pi2%*%PhiR%*%t(Pi1) + Omega2%*%PhiS%*%t(Omega1), 
                               Pi2%*%PhiR%*%t(Pi2) + Omega2%*%PhiS%*%t(Omega2) + Phi42))
  
  #calculate conditional covariance of latent factors given observed data, i.e. Cov(Ei|Fi)
  #and then extract the conditional covariances of each latent factor (Cov(Ri|Fi), Cov(Si1|Fi), etc)
  
  ests$condvar <- ests$sigma22 - ests$sigma21%*%solve(ests$sigma11)%*%t(ests$sigma21)
  
  # calculate conditional expectations of latent factors given observed data
  # 
  # ERSVW <- sigma21%*%solve(sigma11)%*%rbind(t(Xpanel[[1]]), t(Xpanel[[2]]), t(Ypanel[[1]]), t(Ypanel[[2]]))
  # Xpanel = list(t(MyData$X1), t(MyData$X2)), Ypanel = list(t(MyData$Y1), t(MyData$Y2))  
  ests$ERSVW <- ests$sigma21%*%solve(ests$sigma11)%*%rbind( t(Data$X1), t(Data$X2), t(Data$Y1), t(Data$Y2))
  
  return(ests)
}

## end of code