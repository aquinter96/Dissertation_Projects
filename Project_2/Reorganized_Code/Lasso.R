# Lasso.R

lasso <- function(x,y){ #initialize soft-thresholding function
  result <- NULL
  if(abs(x) <= y){
    result <- 0
  } else{
    result <- x - y*sign(x)
  }
  return(result)
}
