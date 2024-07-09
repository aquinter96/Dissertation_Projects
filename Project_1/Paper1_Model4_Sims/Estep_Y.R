## Estep_Y.R

Estep_Y <- function(Data, Old_Par){
  
  ests <- list()
  
  X <- Data$X
  Y <- Data$Y
  A0 <- Old_Par$A0
  B0 <- Old_Par$B0
  A <- Old_Par$A
  Gamma <- Old_Par$Gamma
  Phi1 <- Old_Par$Phi1
  Phi2 <- Old_Par$Phi2

  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  s <- ncol(Old_Par$B)
  m <- ncol(A)

  if(s == 1){
    B <- as.matrix(Old_Par$B)
    Phi3 <- as.matrix(Old_Par$Phi3)
  }else{
    B <- Old_Par$B
    Phi3 <- Old_Par$Phi3
  }

  sigma21 <-  rbind( matrix( cbind( t(B), t(Gamma)%*% t(A)), nrow = s, ncol = q+p)    , matrix(cbind(Gamma%*%t(B),(Phi3 + Gamma%*%t(Gamma))%*%t(A)), nrow = m, ncol = q+p) )  ## m (m for Z)  s (s for W): Z = Gamma W + E

  sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), A%*%Gamma%*%t(B)), rbind(B%*%t(A%*%Gamma),
                                                                               (Phi2 + A%*%(Phi3 + tcrossprod(Gamma))%*%t(A)))), nrow = q+p, ncol = q+p)

  sigma22 <-   rbind( matrix( cbind( diag(s), t(Gamma) ), nrow = s, ncol = (m + s)   )  ,matrix( cbind( Gamma  , Phi3 + Gamma%*%t(Gamma) ), nrow = m, ncol = (m + s)   ) ) ## (m +s) x (m +s )

  inv11 <- Matrix::solve(sigma11)
  ests$condvarWZ <- sigma22 - sigma21%*%inv11%*%t(sigma21)  ##  (m +s) x (m +s)
  ests$condvar  <- ests$condvarWZ[s+(1:m), s+(1:m) ]  ## m x m
 
  EWZ <- sigma21%*%inv11%*%(rbind(t(X),t(Y))-matrix(rbind(B0,A0), nrow = q+p, ncol = n, byrow = F)) ## (m + s) x n
  ests$EW <-  matrix(EWZ[1:s,], nrow = s)   ## q x n
  ests$EZ <-  matrix(EWZ[s+(1:m), ], nrow = m)  ## p x n
  ests$sigma11 <- sigma11
  ests$inv11 <- inv11

  return(ests)
  
}

## end of code