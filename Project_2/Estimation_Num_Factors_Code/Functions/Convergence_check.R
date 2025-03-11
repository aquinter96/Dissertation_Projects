##convergence_check.R

Convergence_check <- function(Old_Par, M_estimates){
  
  myDiffL1 <- NULL  
  
  for (i in 1:length(Old_Par)){
    Param_conv <- c(myDiffL1,   norm(Old_Par[[i]] - M_estimates[[i]], type = "F") ) ## 
  }
  
  return("normDiff" = Param_conv) 
  
}

## end of code