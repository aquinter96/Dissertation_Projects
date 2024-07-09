## BIC_Y.R

BIC_Y <- function(Data, Old_Par, log.lik){
  
  Y <- Data$Y
  A <- Old_Par$A
  
  n <- nrow(Y)
  m <- ncol(A)
  
  BICval <- log(n)*(sum(A != 0) - m) - 2*log.lik
  
  return(BICval)
  
}

## end of code