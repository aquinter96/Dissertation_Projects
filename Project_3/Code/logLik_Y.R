## logLik_Y.R

logLik_Y <- function(Data, E_estimates){
  
  sigma11 <- E_estimates$sigma11
  inv11 <- E_estimates$inv11
  X <- Data$X
  M <- Data$M
  Y <- Data$Y
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(M)
  r <- ncol(Y)
  
  sing.vals <- svd(sigma11)$d
  
  log.lik <- -(1/2)*(sum(diag(inv11%*%(cov(cbind(X,M,Y))+(colMeans(cbind(X,M,Y)) -
                    colMeans(cbind(X,M,Y)))%*%t(colMeans(cbind(X,M,Y)) -
                    colMeans(cbind(X,M,Y)))))) + n*(q+p+r)*log((2*pi))+n*sum(log(sing.vals)))
  
  return(log.lik)
}

## end of code