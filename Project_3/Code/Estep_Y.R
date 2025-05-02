## Estep_Y.R

Estep_Y <- function(Data, Old_Par){
  
  ests <- list()
  
  X <- Data$X
  M <- Data$M
  Y <- Data$Y
  A0 <- Old_Par$A0
  B0 <- Old_Par$B0
  C0 <- Old_Par$C0
  C <- Old_Par$C
  D <- Old_Par$D
  Delta <- Old_Par$Delta
  Phi1 <- Old_Par$Phi1
  Phi3 <- Old_Par$Phi3
  Phi4 <- Old_Par$Phi4
  Phi5 <- Old_Par$Phi5
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(M)
  r <- ncol(Y)
  s <- ncol(Old_Par$B)
  m <- ncol(Old_Par$A)
  t <- ncol(Old_Par$C)

  if(s == 1){
    B <- as.matrix(Old_Par$B)
    Psi <- as.matrix(Old_Par$Psi)
  }else{
    B <- Old_Par$B
    Psi <- Old_Par$Psi
  }
  
  if(m == 1){
    A <- as.matrix(Old_Par$A)
    Gamma <- as.matrix(Old_Par$Gamma, nrow = m)
    Phi2 <- as.matrix(Old_Par$Phi2)
  }else{
    A <- Old_Par$A
    Gamma <- Old_Par$Gamma
    Phi2 <- Old_Par$Phi2
  }
  
  sigma21 <- matrix(cbind(rbind(Psi%*%t(B), Gamma%*%Psi%*%t(B), Delta%*%Gamma%*%Psi%*%t(B)), rbind(Psi%*%t(Gamma)%*%t(A),
                   (Phi2 + Gamma%*%Psi%*%t(Gamma))%*%t(A), Delta%*%(Phi2 + Gamma%*%Psi%*%t(Gamma))%*%t(A)), rbind(Psi%*%t(Gamma)%*%t(Delta)%*%t(C) + Psi%*%t(D),
                        (Phi2 + Gamma%*%Psi%*%t(Gamma))%*%t(Delta)%*%t(C) + Gamma%*%Psi%*%t(D), (Delta%*%(Phi2 + Gamma%*%Psi%*%t(Gamma))%*%t(Delta) + Phi4)%*%t(C) + Delta%*%Gamma%*%Psi%*%t(D))),
                            nrow = s+m+t, ncol = q+p+r)

  sigma11 <- matrix(cbind(rbind((Phi1 + B%*%Psi%*%t(B)), A%*%Gamma%*%Psi%*%t(B), (C%*%Delta%*%Gamma + D)%*%Psi%*%t(B)), rbind(B%*%Psi%*%t(A%*%Gamma),
                   (Phi3 + A%*%(Phi2 + Gamma%*%Psi%*%t(Gamma))%*%t(A)), C%*%Delta%*%(Phi2 + Gamma%*%Psi%*%t(Gamma))%*%t(A) + D%*%Psi%*%t(Gamma)%*%t(A)),
                        rbind(B%*%Psi%*%(t(Gamma)%*%t(Delta)%*%t(C) + t(D)), A%*%(Phi2 + Gamma%*%Psi%*%t(Gamma))%*%t(Delta)%*%t(C) + A%*%Gamma%*%Psi%*%t(D),
                              (C%*%(Delta%*%(Gamma%*%Psi%*%t(Gamma) + Phi2)%*%t(Delta) + Phi4)%*%t(C) + D%*%Psi%*%t(D) + 2*C%*%Delta%*%Gamma%*%Psi%*%t(D) + Phi5))), nrow = q+p+r, ncol = q+p+r)
  
  sigma22 <- matrix(cbind(rbind(Psi, Gamma%*%Psi, Delta%*%Gamma%*%Psi), rbind(Psi%*%t(Gamma), (Gamma%*%Psi%*%t(Gamma) + Phi2),
                   Delta%*%(Gamma%*%Psi%*%t(Gamma) + Phi2)), rbind(Psi%*%t(Gamma)%*%t(Delta), (Gamma%*%Psi%*%t(Gamma) + Phi2)%*%t(Delta),
                           (Delta%*%(Gamma%*%Psi%*%t(Gamma) + Phi2)%*%t(Delta) + Phi4))), nrow = s+m+t, ncol = s+m+t)
  
  inv11 <- Matrix::solve(sigma11)
  ests$condvarVWZ <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
  ests$condvarV <- ests$condvarVWZ[(1:s), (1:s)]
  ests$condvarW <- ests$condvarVWZ[s+(1:m), s+(1:m) ]  ## m x m
  ests$condvarZ <- ests$condvarVWZ[s+m+(1:t), s+m+(1:t)] ## t x t
  if(t == 1){
    ests$condvarVZ <- as.matrix(ests$condvarVWZ[(1:s), s+m+(1:t)], nrow = s)
  }else{
    ests$condvarVZ <- ests$condvarVWZ[(1:s), s+m+(1:t)]
  }
  
  EVWZ <- sigma21%*%inv11%*%(rbind(t(X),t(M),t(Y))-matrix(rbind(B0,A0,C0), nrow = q+p+r, ncol = n, byrow = F)) ## (m + s) x n
  
  if(s == 1){
    ests$EV <- t(as.matrix(EVWZ[1:s,]))
  }else{
    ests$EV <- EVWZ[1:s,]
  }
  
  if(m == 1){
    ests$EW <- t(as.matrix(EVWZ[s+(1:m),]))
  }else{
    ests$EW <- EVWZ[s+(1:m),]
  }
  
  if(t == 1){
    ests$EZ <- t(as.matrix(EVWZ[s+m+(1:t),]))
  }else{
    ests$EZ <- EVWZ[s+m+(1:t),]
  }
  
  ests$sigma11 <- sigma11
  ests$inv11 <- inv11
  
  return(ests)
  
}

## end of code