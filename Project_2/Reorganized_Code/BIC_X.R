## BIC_X.R

BIC_X <- function(Data, Old_Par, log.lik){
  
  A1 <- Old_Par$A1
  A2 <- Old_Par$A2
  B1 <- Old_Par$B1
  B2 <- Old_Par$B2
  
  n <- nrow(Data$X1)
  m <- ncol(A1)
  u1 <- ncol(B1)
  u2 <- ncol(B2)
  
  BICval <- log(n)*(sum(A1 != 0) + sum(A2 != 0) + sum(B1 != 0) + sum(B2 != 0) - (m+u1+u2)) - 2*log.lik
  
  return(BICval)
  
}

## end of code