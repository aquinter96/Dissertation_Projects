
logLik_Y <- function(Data, E_estimates){
  ## to caclulate log-likelihood given current parameters   
  
  X1 <- Data$X1
  X2 <- Data$X2
  Y1 <- Data$Y1
  Y2 <- Data$Y2
  
  q1 <- ncol(X1)
  q2 <- ncol(X2)
  p1 <- ncol(Y1)
  p2 <- ncol(Y2)
  n <-  nrow(Y2)
  
  sigma11 <- E_estimates$sigma11
  
  log.lik <- -(n/2)*(sum(diag(solve(sigma11)%*%(cov(cbind(X1, X2, Y1,Y2))+(colMeans(cbind(X1, X2, Y1, Y2)) - 
             colMeans(cbind(X1, X2, Y1, Y2)))%*%t(colMeans(cbind(X1, X2, Y1, Y2)) - 
             colMeans(cbind(X1, X2, Y1, Y2)))))) + (q1+q2+p1+p2)*log((2*pi))+log(det(sigma11)))
  
  return(log.lik)  
}
