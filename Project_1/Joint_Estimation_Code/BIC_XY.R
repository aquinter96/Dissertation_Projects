## BIC_XY.R

BIC_XY <- function(Data, Old_Par, log.lik){
  
  Y <- Data$Y
  A <- Old_Par$A
  B <- Old_Par$B
  
  n <- nrow(Y)
  m <- ncol(A)
  s <- ncol(B)
  
  BICval <- log(n)*(sum(A != 0) - m + sum(B != 0) - s) - 2*log.lik
  
  return(BICval)
  
}

## end of code