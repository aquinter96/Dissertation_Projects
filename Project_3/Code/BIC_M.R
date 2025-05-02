## BIC_M.R

BIC_M <- function(Data, Old_Par, log.lik){
  
  M <- Data$M
  A <- Old_Par$A
  
  n <- nrow(M)
  m <- ncol(A)
  
  BICval <- log(n)*(sum(A != 0) - m) - 2*log.lik
  
  return(BICval)
  
}

## end of code