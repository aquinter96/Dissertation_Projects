## Data_Generation.R

##### to generate data (X,Y), given model parameters ... 

Data_Generation <- function(n = 1000, Model_Params){
  
  Data <- list()
  
  MyPar <- Model_Params
  
  ak <- cbind(MyPar$A0, MyPar$A)
  bk <- cbind(MyPar$B0, MyPar$B)
  gamma <- MyPar$Gamma

  varvecA <- diag(MyPar$Phi2)
  varvecB <- diag(MyPar$Phi1)
  if((is.atomic(gamma)) & (length(gamma) == 1L)){
    varvecG <- MyPar$Phi3
  }
  else{
    varvecG <- diag(MyPar$Phi3)
  }

  
  
  ## to simulate latent factor and then X/Y
  ## latent factors for X
  
  if((is.atomic(gamma)) & (length(gamma) == 1L)){
    W <- matrix(rnorm(n, mean = 0, sd = 1), nrow = 1, n) ## latent factor for X if s=1
  }
  else{
    W <- matrix(rnorm(n*ncol(gamma), mean = 0, sd = 1), nrow = ncol(gamma), n) ## latent factor for X if s>1
  }
  
  xlin.term  <-  bk %*% rbind(rep(1, ncol(W)),W)
  xlin.term <- t(xlin.term)
  xx <- matrix(0, nrow = n, ncol = nrow(bk))
  for(jj in 1: (ncol(xx))){
    xx[, jj] <- xlin.term[,jj] + cbind(rnorm(n, mean = 0, sd = sqrt(varvecB[jj])))
  }
  colnames(xlin.term) <- paste0('X', 1:(ncol(xx)))
  colnames(xx) <- paste0('X', 1:(ncol(xx)))
  Data$X = xx ## X
  
  ## latent factors for Y
  
  if((is.atomic(gamma)) & (length(gamma) == 1L)){
    Z <- gamma*W + matrix(rnorm(n, mean = 0, sqrt(varvecG)), nrow = 1) ## latent factor for Y if m=1
  }
  else{
    Z <- matrix(0, nrow=nrow(gamma), ncol = n)
    for(i in 1:nrow(gamma)){
      Z[i,] <- gamma[i,]%*%W + matrix(rnorm(n, mean = 0, sd = sqrt(varvecG[i])), nrow = 1) ## latent factor for Y if m>1
    }
  }
  
  ylin.term  <-  ak %*% rbind(rep(1, ncol(Z)),Z)
  ylin.term <- t(ylin.term)
  yy <- matrix(0, nrow = n, ncol = nrow(ak))
  for(jj in 1: (ncol(yy))){
    yy[, jj] <- ylin.term[,jj] + cbind(rnorm(n, mean = 0, sd = sqrt(varvecA[jj])))
  }
  colnames(ylin.term) <- paste0('Y', 1:(ncol(yy)))
  colnames(yy) <- paste0('Y', 1:(ncol(yy)))
  Data$Y = yy ## Y

  return(Data)
  
  ## end of code
  
}

## end of code 
