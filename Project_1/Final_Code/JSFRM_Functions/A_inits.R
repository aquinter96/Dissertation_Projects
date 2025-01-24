## A_inits.R

A_inits <- function(Data, Best, m){
  
  ############################################################################ 
  ## Obtain estimate of the matrix product A*Gamma using method of moments
  ############################################################################ 
  
  if(is.null(dim(Best$B_final_pars$B))){
    B <- as.matrix(Best$B_final_pars$B)
  }else{
    B <- Best$B_final_pars$B
  }
  Psi <- Best$B_final_pars$Psi
  X <- Data$X
  Y <- Data$Y
  
  n <- nrow(X)
  p <- ncol(Y)
  s <- ncol(B)
  
  inits <- list()
  
  A0init <- as.matrix(colMeans(Y))

  AGamma <- t(solve(t(B)%*%B)%*%t(B)%*%cov(X, Y))%*%Matrix::solve(Psi)

  ############################################################################ 
  ## Use sample covariance matrix and estimate of A*Gamma to get the initial
  ## estimates of A/Phi2 using the nlminb() function
  ############################################################################ 
  
  Aobj <-function(a){
    amat <- matrix(0, nrow = p, ncol = m)
    amat[(1:m),(1:m)][lower.tri(amat[(1:m),(1:m)], diag = T)] <- a[1:(m*(m+1)/2)]
    amat[(m+1):p,] <- a[(1+m*(m+1)/2):(p*m-m*(m-1)/2)]
    phi2mat <- diag(a[(p*m-m*(m-1)/2)+(1:p)])
    norm((var(Y) - AGamma%*%Psi%*%t(AGamma) - amat%*%t(amat) - phi2mat), type = "F")
  }
  
  Ainit <- matrix(0, nrow = p, ncol = m)
  if(m==1){
    Ainit[1,1] <- 1
  }
  else{
    diag(Ainit[(1:m),(1:m)]) <- 1
  }
  
  
  Phi2init <- diag(p)
  Ainitvec <- c(Ainit[(1:m),(1:m)][lower.tri(Ainit[(1:m),(1:m)], diag = T)], Ainit[(m+1):p,])
  Phi2initvec <- diag(Phi2init)
  Anewvec <- nlminb(c(Ainitvec, Phi2initvec), Aobj)
  Ainit <- matrix(0, nrow = p, ncol = m)
  Ainit[(1:m),(1:m)][lower.tri(Ainit[(1:m),(1:m)], diag = T)] <- Anewvec$par[1:(m*(m+1)/2)]
  Ainit[(m+1):p,] <- Anewvec$par[(1+m*(m+1)/2):(p*m-m*(m-1)/2)]
  
  Phi2init <- diag(Anewvec$par[(p*m-m*(m-1)/2)+(1:p)])
  A3 <- Ainit%*%(eigen(1/n*t(Ainit)%*%Ainit)$vectors)
  A5 <- A3%*%qr.Q(qr(t(A3[(1:m),(1:m)])))
  
  ############################################################################ 
  ## Use estimates of A/A*Gamma/Phi2 to get initial estimates of Gamma and Phi3
  ############################################################################ 
  
  if(m==1){
    Phi3init <- A5[1,1]^2
    Ainit <- A5 %*% solve(A5[(1:m),(1:m)])
    Ainit[1,1] <- 1
  }
  else{
    Phi3init <- diag(diag(A5[(1:m),(1:m)]))%*%diag(diag(A5[(1:m),(1:m)]))
    Ainit <- A5 %*% solve(diag(diag(A5[(1:m),(1:m)])))
    Ainit[(1:m),(1:(m))][upper.tri(Ainit[(1:m),(1:(m))], diag = F)] <- 0
    diag(Ainit[(1:m),(1:(m))]) <- 1
  }
  AAinv <- Matrix::solve(t(Ainit)%*%Ainit)
  Gammainit <- AAinv%*%t(Ainit)%*%AGamma%*%Matrix::solve(Psi)
  
  inits$A0 <- A0init
  inits$A <- Ainit
  inits$Gamma <- Gammainit
  inits$Phi2 <- Phi2init
  inits$Phi3 <- Phi3init

  return(inits)
}

## end of code