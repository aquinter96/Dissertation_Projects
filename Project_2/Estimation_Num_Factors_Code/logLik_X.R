## logLik_X.R

logLik_X <- function(Data, E_estimates){
  
  X1 <- Data$X1
  X2 <- Data$X2
  n <- nrow(X1)
  q1 <- ncol(X1)
  q2 <- ncol(X2)

  p1 <- nrow(Data$Y1)
  p2 <- nrow(Data$Y2)
  n <-  ncol(Data$Y2)
  
  sigma11 <- E_estimates$sigma11
  
  log.lik <- -(n/2)*(sum(diag(solve(sigma11)%*%(cov(cbind(X1,X2))+(colMeans(cbind(X1,X2)) - colMeans(cbind(X1,X2)))%*%t(colMeans(cbind(X1,X2)) - colMeans(cbind(X1,X2)))))) + (q1+q2)*log((2*pi))+log(det(sigma11)))
  
  return(log.lik)
}

## end of code