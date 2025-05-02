## Estep_Y.R

Estep_M <- function(Data, Old_Par){
  
  ests <- list()
  
  X <- Data$X
  M <- Data$M
  A0 <- Old_Par$A0
  B0 <- Old_Par$B0
  A <- Old_Par$A
  Gamma <- Old_Par$Gamma
  Phi1 <- Old_Par$Phi1
  Phi2 <- Old_Par$Phi2
  Phi3 <- Old_Par$Phi3

  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(M)
  s <- ncol(Old_Par$B)
  m <- ncol(Old_Par$A)

  if(s == 1){
    B <- as.matrix(Old_Par$B)
    Psi <- as.matrix(Old_Par$Psi)
  }else{
    B <- Old_Par$B
    Psi <- Old_Par$Psi
  }

  sigma21 <-  rbind( matrix( cbind( Psi%*%t(B), Psi%*%t(Gamma)%*% t(A)), nrow = s, ncol = q+p)    , matrix(cbind(Gamma%*%Psi%*%t(B),(Phi2 + Gamma%*%Psi%*%t(Gamma))%*%t(A)), nrow = m, ncol = q+p) )  ## m (m for Z)  s (s for W): Z = Gamma W + E
  
  sigma11 <- matrix(cbind(rbind((Phi1 + B%*%Psi%*%t(B)), A%*%Gamma%*%Psi%*%t(B)), rbind(B%*%Psi%*%t(A%*%Gamma),
                                                                                        (Phi3 + A%*%(Phi2 + Gamma%*%Psi%*%t(Gamma))%*%t(A)))), nrow = q+p, ncol = q+p)
  
  sigma22 <-   rbind( matrix( cbind( Psi, Psi%*%t(Gamma) ), nrow = s, ncol = (m + s)   )  ,matrix( cbind( Gamma%*%Psi  , Phi2 + Gamma%*%Psi%*%t(Gamma) ), nrow = m, ncol = (m + s)   ) ) ## (m +s) x (m +s )
  
  inv11 <- Matrix::solve(sigma11)
  ests$condvarVW <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
  ests$condvar  <- ests$condvarVW[s+(1:m), s+(1:m) ]  ## m x m
 
  EVW <- sigma21%*%inv11%*%(rbind(t(X),t(M))-matrix(rbind(B0,A0), nrow = q+p, ncol = n, byrow = F)) ## (m + s) x n
  ests$EV <-  matrix(EVW[1:s,], nrow = s)   ## s x n
  ests$EW <-  matrix(EVW[s+(1:m), ], nrow = m)  ## m x n
  ests$sigma11 <- sigma11
  ests$inv11 <- inv11

  return(ests)
  
}

## end of code