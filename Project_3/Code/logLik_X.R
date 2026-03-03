## logLik_X.R

logLik_X <- function(Data, E_estimates){
  
  modcv <- E_estimates$modcv
  invmodcv <- Matrix::solve(modcv)
  X <- Data
  n <- nrow(X)
  q <- ncol(X)
  
  sing.vals <- svd(modcv)$d

  log.lik <- -(1/2)*(sum(diag(invmodcv%*%(cov(X)+(colMeans(X)-
                    colMeans(X))%*%t(colMeans(X)-colMeans(X))))) +
                    n*q*log((2*pi))+n*sum(log(sing.vals)))
  
  return(log.lik)
}

## end of code