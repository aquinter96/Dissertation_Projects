## BIC_Y.R

BIC_Y <- function(Data, Old_Par, log.lik){
  
  Y <- Data$Y
  C <- Old_Par$C
  D <- Old_Par$D
  
  n <- nrow(Y)
  t <- ncol(C)
  s <- ncol(D)
  
  #BICval <- log(n)*(sum(C != 0) - t + sum(D != 0) - s) - 2*log.lik
  BICval <- log(n)*(sum(D != 0) - s) - 2*log.lik
  
  return(BICval)
  
}

## end of code