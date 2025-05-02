## logLik_M.R

logLik_M <- function(Data, E_estimates){
  
  sigma11 <- E_estimates$sigma11
  inv11 <- E_estimates$inv11
  X <- Data$X
  M <- Data$M
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(M)
  
  log.lik <- -(n/2)*(sum(diag(inv11%*%(cov(cbind(X,M))+(colMeans(cbind(X,M)) - colMeans(cbind(X,M)))%*%t(colMeans(cbind(X,M)) - colMeans(cbind(X,M)))))) + (q+p)*log((2*pi))+log(det(sigma11)))
  
  return(log.lik)
}

## end of code