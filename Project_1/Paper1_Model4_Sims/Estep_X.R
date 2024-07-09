## Estep_X.R

Estep_X <- function(Data, Old_Par){
  
  ests <- list()
  
  B <- Old_Par$B
  B0 <- Old_Par$B0
  Phi1 <- Old_Par$Phi1
  X <- Data
  
  n <- nrow(X)
  q <- nrow(B)
  s <- ncol(B)
  
  ests$modcv <- Phi1 + tcrossprod(B)
  ests$invmodcv <- Matrix::solve(Phi1 + tcrossprod(B))
  ests$condvar <- diag(s) - crossprod(B,ests$invmodcv)%*%B
  ests$EW <- crossprod(B,ests$invmodcv)%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))

  return(ests)
}

## end of code