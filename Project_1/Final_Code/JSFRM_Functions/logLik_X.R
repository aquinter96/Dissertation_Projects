## logLik_X.R

logLik_X <- function(Data, E_estimates){
  
  modcv <- E_estimates$modcv
  invmodcv <- Matrix::solve(modcv)
  X <- Data
  n <- nrow(X)
  q <- ncol(X)

  log.lik <- -(n/2)*(sum(diag(invmodcv%*%(cov(X)+(colMeans(X)-colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
                       q*log((2*pi))+log(det(modcv)))
  
  return(log.lik)
}

## end of code