## Data_Generation.R

##### to generate data (X,Y), given model parameters ... 

Data_Generation <- function(n = 1000, Model_Params){
  
  Data <- list()
  
  MyPar <- Model_Params
  
  ak <- cbind(MyPar$A0, MyPar$A)
  bk <- cbind(MyPar$B0, MyPar$B)
  ck <- cbind(MyPar$C0, MyPar$C)
  dk <- MyPar$D
  
  gamma <- MyPar$Gamma
  delta <- MyPar$Delta

  varvecA <- diag(MyPar$Phi3)
  varvecB <- diag(MyPar$Phi1)
  varvecC <- diag(MyPar$Phi5)

  if((is.atomic(gamma)) & (length(gamma) == 1L)){
    varvecG <- MyPar$Phi2
  }
  else{
    varvecG <- diag(MyPar$Phi2)
  }

  if((is.atomic(delta)) & (length(delta) == 1L)){
    varvecD <- MyPar$Phi4
  }
  else{
    varvecD <- diag(MyPar$Phi4)
  }
  
  
  ## to simulate latent factor and then X/M
  ## latent factors for X
  
  if((is.atomic(gamma)) & (length(gamma) == 1L)){
    V <- matrix(rnorm(n, mean = 0, sd = 1), nrow = 1, n) ## latent factor for X if s=1
  }
  else{
    V <- matrix(rnorm(n*ncol(gamma), mean = 0, sd = 1), nrow = ncol(gamma), n) ## latent factor for X if s>1
  }
  
  xlin.term  <-  bk %*% rbind(rep(1, ncol(V)),V)
  xlin.term <- t(xlin.term)
  xx <- matrix(0, nrow = n, ncol = nrow(bk))
  for(jj in 1: (ncol(xx))){
    xx[, jj] <- xlin.term[,jj] + cbind(rnorm(n, mean = 0, sd = sqrt(varvecB[jj])))
  }
  colnames(xlin.term) <- paste0('X', 1:(ncol(xx)))
  colnames(xx) <- paste0('X', 1:(ncol(xx)))
  Data$X = xx ## X
  
  ## latent factors for M
  
  if((is.atomic(gamma)) & (length(gamma) == 1L)){
    W <- gamma*V + matrix(rnorm(n, mean = 0, sqrt(varvecG)), nrow = 1) ## latent factor for Y if m=1
  }
  else{
    W <- matrix(0, nrow=nrow(gamma), ncol = n)
    for(i in 1:nrow(gamma)){
      W[i,] <- gamma[i,]%*%V + matrix(rnorm(n, mean = 0, sd = sqrt(varvecG[i])), nrow = 1) ## latent factor for Y if m>1
    }
  }
  
  mlin.term  <-  ak %*% rbind(rep(1, ncol(W)),W)
  mlin.term <- t(mlin.term)
  mm <- matrix(0, nrow = n, ncol = nrow(ak))
  for(jj in 1: (ncol(mm))){
    mm[, jj] <- mlin.term[,jj] + cbind(rnorm(n, mean = 0, sd = sqrt(varvecA[jj])))
  }
  colnames(mlin.term) <- paste0('M', 1:(ncol(mm)))
  colnames(mm) <- paste0('M', 1:(ncol(mm)))
  Data$M = mm ## M

  ## latent factors for Y
  
  if((is.atomic(delta)) & (length(delta) == 1L)){
    Z <- delta*W + matrix(rnorm(n, mean = 0, sqrt(varvecD)), nrow = 1) ## latent factor for Y if m=1
  }
  else{
    Z <- matrix(0, nrow=nrow(delta), ncol = n)
    for(i in 1:nrow(delta)){
      Z[i,] <- delta[i,]%*%W + matrix(rnorm(n, mean = 0, sd = sqrt(varvecD[i])), nrow = 1) ## latent factor for Y if m>1
    }
  }
  
  ylin.term  <-  ck %*% rbind(rep(1, ncol(Z)),Z) + dk %*% V
  ylin.term <- t(ylin.term)
  yy <- matrix(0, nrow = n, ncol = nrow(ck))
  for(jj in 1: (ncol(yy))){
    yy[, jj] <- ylin.term[,jj] + cbind(rnorm(n, mean = 0, sd = sqrt(varvecC[jj])))
  }
  colnames(ylin.term) <- paste0('Y', 1:(ncol(yy)))
  colnames(yy) <- paste0('Y', 1:(ncol(yy)))
  Data$Y = yy ## Y
  
  return(Data)
  
  ## end of code
  
}

## end of code 
