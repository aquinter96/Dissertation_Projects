## Estep_X.R

Estep_X <- function(X, B0, B, Phi1){
  
  n <- nrow(X)
  q <- nrow(B)
  xk_sig <- ncol(B)

  invmodcv <- Matrix::solve(Phi1 + tcrossprod(B))
  condvar <- diag(xk_sig) - crossprod(B,invmodcv)%*%B
  EW <- crossprod(B,invmodcv)%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))

  ests <- list("EW" = EW, "condvar" = condvar)
  
  return(ests)
}

## end of code