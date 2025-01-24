## Estep_X.R

Estep_X <- function(Data, Old_Par){
  
  ests <- list()
  
  B <- Old_Par$B
  B0 <- Old_Par$B0
  Psi <- Old_Par$Psi
  Phi1 <- Old_Par$Phi1
  X <- Data

  n <- nrow(X)
  q <- nrow(B)
  s <- ncol(B)
  
  ests$modcv <- Phi1 + B%*%Psi%*%t(B)
  ests$invmodcv <- Matrix::solve(Phi1 + B%*%Psi%*%t(B))
  ests$condvar <- Psi - Psi%*%t(B)%*%ests$invmodcv%*%B%*%Psi
  ests$EW <-Psi%*%t(B)%*%ests$invmodcv%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))

  return(ests)
}

## end of code