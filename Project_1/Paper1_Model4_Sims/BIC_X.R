## BIC_X.R

BIC_X <- function(Data, Old_Par, log.lik){
  
  X <- Data
  B <- Old_Par$B
  
  n <- nrow(X)
  s <- ncol(B)
  
  BICval <- log(n)*(sum(B != 0) - s) - 2*log.lik
 
  return(BICval)
   
}

## end of code