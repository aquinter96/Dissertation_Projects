## BIC_Y.R

BIC_Y <- function(Data, Old_Par, log.lik){
  
  Gamma1 <- Old_Par$Gamma1
  Delta1 <- Old_Par$Delta1
  Gamma2 <- Old_Par$Gamma2
  Delta2 <- Old_Par$Delta2
  
  n <- nrow(Data$Y1)
  j <- ncol(Gamma1)
  z1 <- ncol(Delta1)
  z2 <- ncol(Delta2)
  
  BICval <- log(n)*(sum(Gamma1 != 0) + sum(Delta1 != 0) + sum(Gamma2 != 0) + sum(Delta2 != 0) - (j+z1+z2)) - 2*log.lik
  
  return(BICval)
  
}

## end of code